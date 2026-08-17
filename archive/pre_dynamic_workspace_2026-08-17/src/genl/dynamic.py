from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .form_factors import ELEMENT_SYMBOLS, form_factors, read_form_factor_coefficients
from .kinematic import Control, Instrument, Layer, generate_strain, matlab_round
from .poscar import PoscarStructure, read_poscar

R0_ANGSTROM = 2.814042735053330e-5

try:
    from numba import njit, prange
except ImportError:  # pragma: no cover - depends on optional dependency
    njit = None
    prange = range


@dataclass(frozen=True)
class DynamicResult:
    refl: np.ndarray
    z: np.ndarray
    rho_e: np.ndarray
    amplitude_s: np.ndarray | None = None
    amplitude_p: np.ndarray | None = None


def calc_dynamic_density(
    q: np.ndarray,
    wavelength: float,
    stack: list[Layer],
    control: Control | None = None,
    instrument: Instrument | None = None,
    poscar_dir: str | Path | None = None,
    form_factor_dir: str | Path | None = None,
    vacuum_thick: float = 5.0,
    slices: int = 50,
    max_q0: float = 30.0,
    step_q0: float = 0.1,
    propagation_backend: str | None = None,
) -> DynamicResult:
    """Density-based dynamic scattering calculation ported from GenL MATLAB.

    This first Python port covers substrate-only and substrate + film paths
    without roughness averaging. It mirrors `parratt_matrix_repeated_rhobuiltin`
    enough to support dynamic fitting while keeping the implementation testable.
    The backend selects both electron-density generation and matrix propagation.
    """

    if len(stack) < 1:
        raise ValueError("dynamic density calculation expects at least one layer")
    if any(layer.roughness for layer in stack):
        raise NotImplementedError("dynamic roughness averaging is not ported yet")

    q = np.asarray(q, dtype=float)
    control = control or Control(pol=2, model="density")
    instrument = instrument or Instrument(theta_m=2.0)
    backend = (propagation_backend or os.environ.get("GENL_DYNAMIC_BACKEND", "auto")).lower()
    if backend not in {"auto", "fused", "legacy"}:
        raise ValueError("propagation_backend must be 'auto', 'fused', or 'legacy'")
    add_atom_density = (
        _add_atom_density_numba
        if backend != "legacy" and _add_atom_density_numba is not None
        else _add_atom_density_numpy
    )

    q0 = np.arange(-max_q0, max_q0 + step_q0 * 0.5, step_q0)
    prepared = [
        _prepare_layer(layer, wavelength, q0, poscar_dir, form_factor_dir)
        for layer in stack
    ]

    strain_positions: dict[int, np.ndarray] = {}
    sorted_indices: dict[int, np.ndarray] = {}
    last_atom_z = prepared[0]["lat_par"] * 2.0 + float(np.max(prepared[0]["z_s"]))
    for idx in range(1, len(stack)):
        layer = stack[idx]
        startz = last_atom_z + layer.dinterface
        positions, sorted_idx, _ = generate_strain(
            prepared[idx]["z_s"],
            prepared[idx]["lat_par"],
            matlab_round(layer.n),
            idx,
            startz,
            stack,
        )
        strain_positions[idx] = positions
        sorted_indices[idx] = sorted_idx
        last_atom_z = float(positions[-1])

    if len(stack) == 1:
        totthick = float(last_atom_z)
    else:
        totthick = float(max(strain_positions[len(stack) - 1]))
    substrate = prepared[0]
    dz = substrate["lat_par"] / slices
    vacuum_slices = matlab_round(vacuum_thick / dz)
    vacuum_thick_exact = dz * vacuum_slices
    z = np.arange(-vacuum_thick_exact, totthick + vacuum_thick_exact + dz * 0.5, dz)
    rho_e = np.zeros_like(z, dtype=complex)

    for repeat in range(3):
        for atom_idx, z_atom in enumerate(substrate["z_s"]):
            pos_z = z_atom + substrate["lat_par"] * repeat
            add_atom_density(
                rho_e,
                z,
                pos_z,
                substrate["lat_par"],
                substrate["area"],
                q0,
                substrate["ff"][:, atom_idx],
            )

    start_idx = slices
    rho_0 = rho_e[: vacuum_slices + start_idx]
    rho_1 = rho_e[vacuum_slices + start_idx : vacuum_slices + start_idx * 2]
    substrate_end = vacuum_slices + start_idx * 2

    for idx in range(1, len(stack)):
        layer_data = prepared[idx]
        positions = strain_positions[idx]
        sorted_idx = sorted_indices[idx]
        for atom_counter, pos_z in enumerate(positions):
            form_factor_idx = sorted_idx[atom_counter % len(sorted_idx)]
            add_atom_density(
                rho_e,
                z,
                pos_z,
                layer_data["lat_par"],
                layer_data["area"],
                q0,
                layer_data["ff"][:, form_factor_idx],
            )

    rho_rest = rho_e[substrate_end:]
    rho_0r = rho_0[::-1]
    rho_1r = rho_1[::-1]
    rho_restr = rho_rest[::-1]
    substrate_repeats = matlab_round(stack[0].n)

    if control.pol in (0, 1):
        amplitude = propagate_amplitude_vectorized(
            q, wavelength, rho_0r, rho_1r, substrate_repeats, rho_restr, dz, control.pol, backend
        )
        refl = np.abs(amplitude) ** 2
        amplitude_s = amplitude if control.pol == 0 else None
        amplitude_p = amplitude if control.pol == 1 else None
    elif control.pol == 2:
        amplitude_s = propagate_amplitude_vectorized(
            q, wavelength, rho_0r, rho_1r, substrate_repeats, rho_restr, dz, 0, backend
        )
        amplitude_p = propagate_amplitude_vectorized(
            q, wavelength, rho_0r, rho_1r, substrate_repeats, rho_restr, dz, 1, backend
        )
        mono = np.cos(np.deg2rad(instrument.theta_m * 2.0)) ** 2
        refl = (np.abs(amplitude_s) ** 2 + mono * np.abs(amplitude_p) ** 2) / (1.0 + mono)
    else:
        raise ValueError("control.pol must be 0, 1, or 2")

    return DynamicResult(
        refl=np.real_if_close(refl).real,
        z=z,
        rho_e=rho_e,
        amplitude_s=amplitude_s,
        amplitude_p=amplitude_p,
    )


def propagate_vectorized(
    q: np.ndarray,
    wavelength: float,
    rho_0r: np.ndarray,
    rho_1r: np.ndarray,
    substrate_repeats: int,
    rho_restr: np.ndarray,
    dz: float,
    pol: int,
    backend: str = "auto",
) -> np.ndarray:
    amplitude = propagate_amplitude_vectorized(
        q, wavelength, rho_0r, rho_1r, substrate_repeats, rho_restr, dz, pol, backend
    )
    return np.abs(amplitude) ** 2


def propagate_amplitude_vectorized(
    q: np.ndarray,
    wavelength: float,
    rho_0r: np.ndarray,
    rho_1r: np.ndarray,
    substrate_repeats: int,
    rho_restr: np.ndarray,
    dz: float,
    pol: int,
    backend: str = "auto",
) -> np.ndarray:
    if backend in {"auto", "fused"} and _propagate_vectorized_fused_numba is not None:
        return _propagate_vectorized_fused_numba(
            q, wavelength, rho_0r, rho_1r, max(0, substrate_repeats - 2), rho_restr, dz, pol
        )
    if backend == "fused":
        raise RuntimeError("fused dynamic propagation requires numba")

    a1, a2, a3, a4 = make_a_matrix(q, wavelength, rho_1r, dz, pol)
    unit_cell_matrix = do_matrix_propagation(a1, a2, a3, a4)
    exponent = max(0, substrate_repeats - 2)

    a1, a2, a3, a4 = make_a_matrix(q, wavelength, rho_0r, dz, pol)
    cap_matrix = do_matrix_propagation(a1, a2, a3, a4)

    a1, a2, a3, a4 = make_a_matrix(q, wavelength, rho_restr, dz, pol)
    film_matrix = do_matrix_propagation(a1, a2, a3, a4)

    if _combine_transfer_matrices_numba is not None:
        return _combine_transfer_matrices_numba(unit_cell_matrix, cap_matrix, film_matrix, exponent)

    substrate_matrix = np.empty_like(unit_cell_matrix)
    for idx in range(len(q)):
        substrate_matrix[:, :, idx] = fast_matrix_power(unit_cell_matrix[:, :, idx], exponent)
        substrate_matrix[:, :, idx] = substrate_matrix[:, :, idx] @ cap_matrix[:, :, idx]
        substrate_matrix[:, :, idx] = film_matrix[:, :, idx] @ substrate_matrix[:, :, idx]

    return substrate_matrix[1, 0, :] / substrate_matrix[0, 0, :]


if njit is not None:

    @njit(cache=True)
    def _add_atom_density_numba(
        rho_e: np.ndarray,
        z: np.ndarray,
        pos_z: float,
        lat_par: float,
        area: float,
        q0: np.ndarray,
        f_q0: np.ndarray,
    ) -> None:
        if len(z) < 2 or len(q0) < 2:
            return
        dz = z[1] - z[0]
        lower = pos_z - lat_par
        upper = pos_z + lat_par
        start = max(0, int(np.ceil((lower - z[0]) / dz)))
        stop = min(len(z) - 1, int(np.floor((upper - z[0]) / dz)))

        while start > 0 and z[start - 1] >= lower:
            start -= 1
        while start < len(z) and z[start] < lower:
            start += 1
        while stop + 1 < len(z) and z[stop + 1] <= upper:
            stop += 1
        while stop >= 0 and z[stop] > upper:
            stop -= 1
        if stop < start:
            return

        delta = z[start] - pos_z
        for q_idx in range(len(q0)):
            if q_idx == 0:
                weight = 0.5 * (q0[1] - q0[0])
            elif q_idx == len(q0) - 1:
                weight = 0.5 * (q0[-1] - q0[-2])
            else:
                weight = 0.5 * (q0[q_idx + 1] - q0[q_idx - 1])

            value = (
                weight
                * f_q0[q_idx]
                * np.exp(1j * q0[q_idx] * delta)
                / area
            )
            phase_step = np.exp(1j * q0[q_idx] * dz)
            for z_idx in range(start, stop + 1):
                rho_e[z_idx] += value
                value *= phase_step

    @njit(parallel=True, cache=True)
    def _combine_transfer_matrices_numba(
        unit_cell_matrix: np.ndarray,
        cap_matrix: np.ndarray,
        film_matrix: np.ndarray,
        exponent: int,
    ) -> np.ndarray:
        q_count = unit_cell_matrix.shape[2]
        amplitude = np.empty(q_count, dtype=np.complex128)
        for idx in prange(q_count):
            r00 = 1.0 + 0.0j
            r01 = 0.0 + 0.0j
            r10 = 0.0 + 0.0j
            r11 = 1.0 + 0.0j

            b00 = unit_cell_matrix[0, 0, idx]
            b01 = unit_cell_matrix[0, 1, idx]
            b10 = unit_cell_matrix[1, 0, idx]
            b11 = unit_cell_matrix[1, 1, idx]
            k = exponent
            while k > 0:
                if k % 2 == 1:
                    n00 = r00 * b00 + r01 * b10
                    n01 = r00 * b01 + r01 * b11
                    n10 = r10 * b00 + r11 * b10
                    n11 = r10 * b01 + r11 * b11
                    r00 = n00
                    r01 = n01
                    r10 = n10
                    r11 = n11

                n00 = b00 * b00 + b01 * b10
                n01 = b00 * b01 + b01 * b11
                n10 = b10 * b00 + b11 * b10
                n11 = b10 * b01 + b11 * b11
                b00 = n00
                b01 = n01
                b10 = n10
                b11 = n11

                scale = max(max(b00.real, b01.real), max(b10.real, b11.real))
                if scale != 0.0:
                    b00 = b00 / scale
                    b01 = b01 / scale
                    b10 = b10 / scale
                    b11 = b11 / scale
                k = k // 2

            c00 = cap_matrix[0, 0, idx]
            c01 = cap_matrix[0, 1, idx]
            c10 = cap_matrix[1, 0, idx]
            c11 = cap_matrix[1, 1, idx]
            n00 = r00 * c00 + r01 * c10
            n01 = r00 * c01 + r01 * c11
            n10 = r10 * c00 + r11 * c10
            n11 = r10 * c01 + r11 * c11
            r00 = n00
            r01 = n01
            r10 = n10
            r11 = n11

            f00 = film_matrix[0, 0, idx]
            f01 = film_matrix[0, 1, idx]
            f10 = film_matrix[1, 0, idx]
            f11 = film_matrix[1, 1, idx]
            n00 = f00 * r00 + f01 * r10
            n10 = f10 * r00 + f11 * r10
            amplitude[idx] = n10 / n00
        return amplitude

    @njit(cache=True)
    def _a_values_at_layer(
        q_value: float,
        wavelength: float,
        rho_current: complex,
        rho_previous: complex,
        dz: float,
        pol: int,
    ) -> tuple[complex, complex, complex, complex]:
        k0 = 2.0 * np.pi / wavelength
        factor = 8.0 * k0**2 * wavelength**2 / (2.0 * np.pi)
        sld_current = rho_current * R0_ANGSTROM / (2.0 * np.pi)
        sld_previous = rho_previous * R0_ANGSTROM / (2.0 * np.pi)

        kz_current = np.sqrt(
            q_value**2
            - factor * sld_current.real
            + 1j * factor * sld_current.imag
        )
        kz_previous = np.sqrt(
            q_value**2
            - factor * sld_previous.real
            + 1j * factor * sld_previous.imag
        )
        exp_pos = np.exp(kz_current * (1j * dz / 2.0))

        if pol == 0:
            denom = 1.0 / (kz_previous + kz_current)
            r = (kz_previous - kz_current) * denom
            t = 2.0 * kz_previous * denom
            invtp = exp_pos / t
            invtm = 1.0 / (t * exp_pos)
            return invtp, invtm * r, invtp * r, invtm

        deltan_current = wavelength**2 * sld_current.real / (2.0 * np.pi)
        betan_current = 1j * wavelength**2 * sld_current.imag / (2.0 * np.pi)
        deltan_previous = wavelength**2 * sld_previous.real / (2.0 * np.pi)
        betan_previous = 1j * wavelength**2 * sld_previous.imag / (2.0 * np.pi)
        n_current = 1.0 - deltan_current - betan_current
        n_previous = 1.0 - deltan_previous - betan_previous
        fact1 = n_previous**2 * kz_current
        fact2 = n_current**2 * kz_previous
        denom = fact1 + fact2
        r_pi = (fact1 - fact2) / denom
        t_pi = 2.0 * fact1 / denom
        invtp_pi = exp_pos / t_pi
        invtm_pi = 1.0 / (t_pi * exp_pos)
        return invtp_pi, r_pi * invtm_pi, r_pi * invtp_pi, invtm_pi

    @njit(cache=True)
    def _profile_matrix_at_q(
        q_value: float,
        wavelength: float,
        rho_e: np.ndarray,
        dz: float,
        pol: int,
    ) -> tuple[complex, complex, complex, complex]:
        b1 = 1.0 + 0.0j
        b2 = 0.0 + 0.0j
        b3 = 0.0 + 0.0j
        b4 = 1.0 + 0.0j

        for layer in range(len(rho_e) - 1, 0, -1):
            a1, a2, a3, a4 = _a_values_at_layer(
                q_value, wavelength, rho_e[layer], rho_e[layer - 1], dz, pol
            )
            wa1 = a3 - a1
            wa2 = a3 + a4
            wa3 = wa1 + a4
            wa4 = a2 - wa3

            t = a1 * b1
            b1_minus = b2 - b4
            u = wa1 * b1_minus
            v = wa2 * (b2 - b1)
            w = t + wa3 * (b1 - b1_minus)
            w1 = w + u
            b3_old = b3
            b3 = w1 + a4 * (b3_old - b1 + b1_minus)
            b1 = t + a2 * b3_old
            b2 = w + v + wa4 * b4
            b4 = w1 + v

        return b1, b2, b3, b4

    @njit(parallel=True, cache=True)
    def _propagate_vectorized_fused_numba(
        q: np.ndarray,
        wavelength: float,
        rho_0r: np.ndarray,
        rho_1r: np.ndarray,
        exponent: int,
        rho_restr: np.ndarray,
        dz: float,
        pol: int,
    ) -> np.ndarray:
        amplitude = np.empty(len(q), dtype=np.complex128)
        for idx in prange(len(q)):
            q_value = q[idx]
            r00 = 1.0 + 0.0j
            r01 = 0.0 + 0.0j
            r10 = 0.0 + 0.0j
            r11 = 1.0 + 0.0j

            b00, b01, b10, b11 = _profile_matrix_at_q(q_value, wavelength, rho_1r, dz, pol)
            k = exponent
            while k > 0:
                if k % 2 == 1:
                    n00 = r00 * b00 + r01 * b10
                    n01 = r00 * b01 + r01 * b11
                    n10 = r10 * b00 + r11 * b10
                    n11 = r10 * b01 + r11 * b11
                    r00 = n00
                    r01 = n01
                    r10 = n10
                    r11 = n11

                n00 = b00 * b00 + b01 * b10
                n01 = b00 * b01 + b01 * b11
                n10 = b10 * b00 + b11 * b10
                n11 = b10 * b01 + b11 * b11
                b00 = n00
                b01 = n01
                b10 = n10
                b11 = n11

                scale = max(max(b00.real, b01.real), max(b10.real, b11.real))
                if scale != 0.0:
                    b00 = b00 / scale
                    b01 = b01 / scale
                    b10 = b10 / scale
                    b11 = b11 / scale
                k = k // 2

            c00, c01, c10, c11 = _profile_matrix_at_q(q_value, wavelength, rho_0r, dz, pol)
            n00 = r00 * c00 + r01 * c10
            n01 = r00 * c01 + r01 * c11
            n10 = r10 * c00 + r11 * c10
            n11 = r10 * c01 + r11 * c11
            r00 = n00
            r01 = n01
            r10 = n10
            r11 = n11

            f00, f01, f10, f11 = _profile_matrix_at_q(q_value, wavelength, rho_restr, dz, pol)
            n00 = f00 * r00 + f01 * r10
            n10 = f10 * r00 + f11 * r10
            amplitude[idx] = n10 / n00
        return amplitude

else:
    _add_atom_density_numba = None
    _combine_transfer_matrices_numba = None
    _propagate_vectorized_fused_numba = None


def make_a_matrix(
    q: np.ndarray,
    wavelength: float,
    rho_e: np.ndarray,
    dz: float,
    pol: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    k0 = 2.0 * np.pi / wavelength
    factor = 8.0 * k0**2 * wavelength**2 / (2.0 * np.pi)
    sld = rho_e * R0_ANGSTROM / (2.0 * np.pi)
    delta = factor * np.real(sld)
    beta = 1j * factor * np.imag(sld)

    kz_all = np.sqrt(q[np.newaxis, :] ** 2 - delta[:, np.newaxis] + beta[:, np.newaxis])
    kz_shifted = np.roll(kz_all, 1, axis=0)
    exp_pos = np.exp(kz_all * (1j * dz / 2.0))

    if pol == 0:
        denom = 1.0 / (kz_shifted + kz_all)
        r = (kz_shifted - kz_all) * denom
        t = 2.0 * kz_shifted * denom
        invtp = exp_pos / t
        invtm = 1.0 / (t * exp_pos)
        return invtp, invtm * r, invtp * r, invtm

    if pol == 1:
        deltan = wavelength**2 * np.real(sld) / (2.0 * np.pi)
        betan = 1j * wavelength**2 * np.imag(sld) / (2.0 * np.pi)
        n = 1.0 - deltan[:, np.newaxis] - betan[:, np.newaxis]
        n_shifted = np.roll(n, 1, axis=0)
        fact1 = n_shifted**2 * kz_all
        fact2 = n**2 * kz_shifted
        denom = fact1 + fact2
        r_pi = (fact1 - fact2) / denom
        t_pi = 2.0 * fact1 / denom
        invtp_pi = exp_pos / t_pi
        invtm_pi = 1.0 / (t_pi * exp_pos)
        return invtp_pi, r_pi * invtm_pi, r_pi * invtp_pi, invtm_pi

    raise ValueError("pol must be 0 or 1")


def do_matrix_propagation(
    a1: np.ndarray,
    a2: np.ndarray,
    a3: np.ndarray,
    a4: np.ndarray,
) -> np.ndarray:
    layers, q_count = a1.shape
    b1 = np.ones(q_count, dtype=complex)
    b2 = np.zeros(q_count, dtype=complex)
    b3 = np.zeros(q_count, dtype=complex)
    b4 = np.ones(q_count, dtype=complex)

    wa1 = a3 - a1
    wa2 = a3 + a4
    wa3 = wa1 + a4
    wa4 = a2 - wa3

    for layer in range(layers - 1, 0, -1):
        t = a1[layer, :] * b1
        b1_minus = b2 - b4
        u = wa1[layer, :] * b1_minus
        v = wa2[layer, :] * (b2 - b1)
        w = t + wa3[layer, :] * (b1 - b1_minus)
        w1 = w + u
        b3_old = b3
        b3 = w1 + a4[layer, :] * (b3_old - b1 + b1_minus)
        b1 = t + a2[layer, :] * b3_old
        b2 = w + v + wa4[layer, :] * b4
        b4 = w1 + v

    result = np.zeros((2, 2, q_count), dtype=complex)
    result[0, 0, :] = b1
    result[0, 1, :] = b2
    result[1, 0, :] = b3
    result[1, 1, :] = b4
    return result


def fast_matrix_power(matrix: np.ndarray, exponent: int) -> np.ndarray:
    if exponent == 0:
        return np.eye(2, dtype=complex)
    result = np.eye(2, dtype=complex)
    base = matrix.copy()
    k = int(exponent)
    while k > 0:
        if k % 2 == 1:
            result = result @ base
        base = base @ base
        scale = np.max(np.real(base))
        if scale != 0:
            base = base / scale
        k //= 2
    return result


def _prepare_layer(
    layer: Layer,
    wavelength: float,
    q0: np.ndarray,
    poscar_dir: str | Path | None,
    form_factor_dir: str | Path | None,
) -> dict[str, np.ndarray | float | PoscarStructure]:
    structure = read_poscar(layer.filename, poscar_dir)
    scaling, area, _ = _layer_geometry(structure, layer)
    z_s = structure.positions @ scaling
    ff = np.zeros((len(q0), len(z_s)), dtype=complex)

    cursor = 0
    for element, count in zip(structure.types, structure.type_counts):
        atomic_number = ELEMENT_SYMBOLS.index(element) + 1
        f_q, _, _ = form_factors(
            q0,
            read_form_factor_coefficients(atomic_number, wavelength, form_factor_dir),
            1.0,
        )
        for _ in range(int(count)):
            ff[:, cursor] = f_q
            cursor += 1

    return {
        "structure": structure,
        "scaling": scaling,
        "area": area,
        "lat_par": float(np.linalg.norm(scaling)),
        "z_s": z_s,
        "ff": ff,
    }


def _layer_geometry(structure: PoscarStructure, layer: Layer) -> tuple[np.ndarray, float, float]:
    area_scale = float(np.sqrt(layer.area_scale))
    axes = [structure.a1, structure.a2, structure.a3]
    direction_index = layer.direction - 1
    scaling = axes[direction_index] * layer.scale
    transverse = [axis * area_scale for idx, axis in enumerate(axes) if idx != direction_index]
    area = float(abs(np.linalg.norm(np.cross(transverse[0], transverse[1]))))
    unit_cell_volume = float(abs(np.dot(scaling, np.cross(transverse[0], transverse[1]))))
    return scaling, area, unit_cell_volume


def _add_atom_density_numpy(
    rho_e: np.ndarray,
    z: np.ndarray,
    pos_z: float,
    lat_par: float,
    area: float,
    q0: np.ndarray,
    f_q0: np.ndarray,
) -> None:
    idx = (z >= pos_z - lat_par) & (z <= pos_z + lat_par)
    if not np.any(idx):
        return
    phase = np.exp(1j * q0[:, np.newaxis] * (z[idx][np.newaxis, :] - pos_z))
    rho_e[idx] += np.trapezoid(f_q0[:, np.newaxis] * phase, q0, axis=0) / area
