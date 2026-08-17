from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .debye import debye_waller_prefactor
from .form_factors import ELEMENT_SYMBOLS, form_factors, read_form_factor_coefficients
from .poscar import PoscarStructure, read_poscar

R0_ANGSTROM = 2.814042735053330e-5


def matlab_round(value: float) -> int:
    return int(np.floor(value + 0.5)) if value >= 0 else int(np.ceil(value - 0.5))


@dataclass
class Layer:
    direction: int
    n: float
    filename: str | Path
    dinterface: float = 0.0
    scale: float = 1.0
    area_scale: float = 1.0
    roughness: bool = False
    sigma: float = 0.0
    bottom_strain_amplitude: float = 0.0
    bottom_strain_end: float = 0.0
    top_strain_amplitude: float = 0.0
    top_strain_end: float = 0.0
    pre_calc_f: dict[str, Any] | None = field(default=None, repr=False)


@dataclass
class Control:
    pol: int = 0
    model: str = "kinematic"


@dataclass
class Instrument:
    theta_m_path: int = 1
    theta_m: float = 2.0
    theta: np.ndarray | None = None


@dataclass(frozen=True)
class KinematicResult:
    refl: np.ndarray


def calc_kinematic(
    q: np.ndarray,
    wavelength: float,
    stack: list[Layer],
    control: Control | None = None,
    instrument: Instrument | None = None,
    poscar_dir: str | Path | None = None,
    form_factor_dir: str | Path | None = None,
) -> KinematicResult:
    """Calculate the GenL kinematic approximation for a stack of layers."""

    q = np.asarray(q, dtype=float)
    control = control or Control()
    instrument = instrument or Instrument()
    total_amplitude = np.zeros_like(q, dtype=complex)
    strain: dict[int, np.ndarray] = {}

    for layer_index, layer in enumerate(stack):
        pre_calc = layer.pre_calc_f
        if pre_calc is None:
            structure = read_poscar(layer.filename, poscar_dir)
            pre_calc = {"structure": structure, "elements": {}}
        else:
            structure = pre_calc["structure"]

        scaling, area, unit_cell_volume = _layer_geometry(structure, layer)
        lattice_parameter = float(np.linalg.norm(scaling))
        atom_z = structure.positions @ scaling
        f_by_atom, substrate_f, substrate_f0 = _atom_form_factors(
            q,
            wavelength,
            structure,
            atom_z,
            scaling,
            pre_calc,
            layer_index == 0,
            form_factor_dir,
        )
        n_cells = matlab_round(layer.n)

        if layer_index == 0:
            g = 4.0 * np.pi * substrate_f / unit_cell_volume * lattice_parameter * R0_ANGSTROM
            with np.errstate(divide="ignore", invalid="ignore"):
                g0 = (
                    4.0
                    * np.pi
                    * substrate_f0
                    / unit_cell_volume
                    * lattice_parameter
                    * R0_ANGSTROM
                    / q
                )
            layer_amplitude = _substrate_like_amplitude(
                q, g, g0, lattice_parameter, n_cells, layer
            )
            total_amplitude = total_amplitude + layer_amplitude
            last_atom_z = n_cells * lattice_parameter + float(np.max(atom_z))
        else:
            if layer_index > 1:
                last_atom_z = float(strain[layer_index - 1][-1])
            startz = last_atom_z + layer.dinterface
            layer_amplitude, strained = _layer_amplitude(
                q,
                atom_z,
                lattice_parameter,
                n_cells,
                layer_index,
                startz,
                stack,
                f_by_atom,
                unit_cell_volume,
            )
            strain[layer_index] = strained
            total_amplitude = total_amplitude + np.exp(1j * q * layer.dinterface) * layer_amplitude

        layer.pre_calc_f = pre_calc

    intensity = np.abs(total_amplitude) ** 2
    theta = instrument.theta
    if theta is None:
        theta = np.rad2deg(np.arcsin(np.clip(q * wavelength / (4.0 * np.pi), -1.0, 1.0)))

    if instrument.theta_m_path == 1:
        if control.pol == 1:
            intensity = intensity * np.cos(np.deg2rad(theta * 2.0)) ** 2
        elif control.pol == 2:
            numerator = 1.0 + (
                np.cos(np.deg2rad(instrument.theta_m * 2.0)) ** 2
                * np.cos(np.deg2rad(theta * 2.0)) ** 2
            )
            denominator = 1.0 + np.cos(np.deg2rad(instrument.theta_m * 2.0)) ** 2
            intensity = intensity * numerator / denominator

    return KinematicResult(refl=intensity)


def generate_strain(
    atom_z: np.ndarray,
    lattice_parameter: float,
    n_cells: int,
    layer_index: int,
    startz: float,
    stack: list[Layer],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    pos_vector = np.zeros(n_cells * len(atom_z), dtype=float)
    f_index = np.zeros(n_cells * len(atom_z), dtype=int)
    sorted_idx = np.argsort(atom_z)
    sorted_z = atom_z[sorted_idx]

    cursor = 0
    for cell_index in range(n_cells):
        for atom_index, z_value in zip(sorted_idx, sorted_z):
            pos_vector[cursor] = z_value + lattice_parameter * cell_index
            f_index[cursor] = atom_index
            cursor += 1

    pos_vector = pos_vector - pos_vector[0] + startz
    strained_position = np.zeros_like(pos_vector)
    displacement = np.zeros_like(pos_vector)
    strained_position[0] = pos_vector[0]

    layer = stack[layer_index]
    bottom_boundary_matlab = matlab_round(layer.bottom_strain_end)
    if bottom_boundary_matlab > len(pos_vector):
        bottom_boundary_matlab = len(pos_vector)
    bottom_boundary = max(bottom_boundary_matlab - 1, 0)

    top_boundary_matlab = len(pos_vector) - matlab_round(layer.top_strain_end)
    if top_boundary_matlab <= 0:
        top_boundary_matlab = len(pos_vector)
    top_boundary = max(top_boundary_matlab - 1, 0)

    for idx in range(1, len(pos_vector)):
        delta = pos_vector[idx] - pos_vector[idx - 1]
        if idx + 1 > layer.bottom_strain_end:
            displacement[idx] = delta
        elif delta == 0:
            displacement[idx] = 0.0
        else:
            displacement[idx] = delta + layer.bottom_strain_amplitude * (
                pos_vector[bottom_boundary] - pos_vector[idx]
            )

        if idx + 1 > top_boundary_matlab:
            if delta == 0:
                displacement[idx] = 0.0
            else:
                displacement[idx] = displacement[idx] + layer.top_strain_amplitude * (
                    pos_vector[idx] - pos_vector[top_boundary]
                )

        strained_position[idx] = strained_position[idx - 1] + displacement[idx]

    return strained_position, sorted_idx, f_index


def _layer_geometry(structure: PoscarStructure, layer: Layer) -> tuple[np.ndarray, float, float]:
    area_scale = float(np.sqrt(layer.area_scale))
    axes = [structure.a1, structure.a2, structure.a3]
    direction_index = layer.direction - 1
    if direction_index not in (0, 1, 2):
        raise ValueError("Layer direction must be 1, 2, or 3")

    scaling = axes[direction_index] * layer.scale
    transverse = [axis * area_scale for idx, axis in enumerate(axes) if idx != direction_index]
    area = float(abs(np.linalg.norm(np.cross(transverse[0], transverse[1]))))
    unit_cell_volume = float(abs(np.dot(scaling, np.cross(transverse[0], transverse[1]))))
    return scaling, area, unit_cell_volume


def _atom_form_factors(
    q: np.ndarray,
    wavelength: float,
    structure: PoscarStructure,
    atom_z: np.ndarray,
    scaling: np.ndarray,
    pre_calc: dict[str, Any],
    calculate_substrate_terms: bool,
    form_factor_dir: str | Path | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    f_by_atom = np.zeros((len(q), len(atom_z)), dtype=complex)
    substrate_f = np.zeros_like(q, dtype=complex)
    substrate_f0 = np.zeros_like(q, dtype=complex)

    atom_cursor = 0
    for element, count in zip(structure.types, structure.type_counts):
        atomic_number = ELEMENT_SYMBOLS.index(element) + 1
        element_cache = pre_calc["elements"].get(element)
        if element_cache is None:
            f_q, _, _ = form_factors(
                q,
                read_form_factor_coefficients(atomic_number, wavelength, form_factor_dir),
                1.0,
            )
            b_factor = debye_waller_prefactor(atomic_number)
            with np.errstate(divide="ignore", invalid="ignore"):
                f_q = f_q * np.exp(-b_factor * (q / (4.0 * np.pi)) ** 2) / q
            f0, _, _ = form_factors(
                0.0,
                read_form_factor_coefficients(atomic_number, wavelength, form_factor_dir),
                1.0,
            )
            element_cache = {"f": f_q, "f0": f0[0]}
            pre_calc["elements"][element] = element_cache

        for _ in range(int(count)):
            f_by_atom[:, atom_cursor] = element_cache["f"]
            if calculate_substrate_terms:
                substrate_f = substrate_f + element_cache["f"] * np.exp(
                    1j * q * (structure.positions[atom_cursor] @ scaling)
                )
                substrate_f0 = substrate_f0 + element_cache["f0"]
            atom_cursor += 1

    return f_by_atom, substrate_f, substrate_f0


def _substrate_like_amplitude(
    q: np.ndarray,
    g: np.ndarray,
    g0: np.ndarray,
    lattice_parameter: float,
    n_cells: int,
    layer: Layer,
) -> np.ndarray:
    if layer.roughness:
        nvector, weights = _roughness_distribution(n_cells, layer.sigma)
        amplitudes = np.zeros((len(q), len(nvector)), dtype=complex)
        tally = 0.0
        for idx, (n_value, weight) in enumerate(zip(nvector, weights)):
            amplitudes[:, idx] = (
                np.exp(1j * q * layer.dinterface)
                * -1j
                * g
                * (1.0 - np.exp(1j * (q * lattice_parameter - 2.0 * g0) * n_value))
                / (1.0 - np.exp(1j * (q * lattice_parameter - 2.0 * g0)))
            )
            amplitudes[:, idx] *= np.sqrt(weight)
            tally += np.sqrt(weight)
        return np.sum(amplitudes / tally, axis=1)

    return (
        -1j
        * g
        * (1.0 - np.exp(1j * q * lattice_parameter * n_cells))
        / (1.0 - np.exp(1j * q * lattice_parameter))
    )


def _layer_amplitude(
    q: np.ndarray,
    atom_z: np.ndarray,
    lattice_parameter: float,
    n_cells: int,
    layer_index: int,
    startz: float,
    stack: list[Layer],
    f_by_atom: np.ndarray,
    unit_cell_volume: float,
) -> tuple[np.ndarray, np.ndarray]:
    layer = stack[layer_index]
    if layer.roughness:
        nvector, weights = _roughness_distribution(n_cells, layer.sigma)
        amplitudes = np.zeros((len(q), len(nvector)), dtype=complex)
        tally = 0.0
        last_strain = np.array([], dtype=float)
        for idx, (n_value, weight) in enumerate(zip(nvector, weights)):
            layer_total, last_strain = _calc_f(
                q,
                atom_z,
                lattice_parameter,
                int(n_value),
                layer_index,
                startz,
                stack,
                f_by_atom,
                unit_cell_volume,
            )
            amplitudes[:, idx] = np.exp(1j * q * layer.dinterface) * layer_total
            amplitudes[:, idx] *= np.sqrt(weight)
            tally += np.sqrt(weight)
        return np.sum(amplitudes / tally, axis=1), last_strain

    return _calc_f(
        q,
        atom_z,
        lattice_parameter,
        n_cells,
        layer_index,
        startz,
        stack,
        f_by_atom,
        unit_cell_volume,
    )


def _calc_f(
    q: np.ndarray,
    atom_z: np.ndarray,
    lattice_parameter: float,
    n_cells: int,
    layer_index: int,
    startz: float,
    stack: list[Layer],
    f_by_atom: np.ndarray,
    unit_cell_volume: float,
) -> tuple[np.ndarray, np.ndarray]:
    position_vector, _, f_index = generate_strain(
        atom_z, lattice_parameter, n_cells, layer_index, startz, stack
    )
    pre_factor = -1j * 4.0 * np.pi * R0_ANGSTROM * lattice_parameter / unit_cell_volume
    weighted_f = f_by_atom * pre_factor
    total = np.sum(
        weighted_f[:, f_index] * np.exp(1j * q[:, np.newaxis] * position_vector[np.newaxis, :]),
        axis=1,
    )
    return total, position_vector


def _roughness_distribution(n_cells: int, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    nvector = np.arange(
        max(1, matlab_round(n_cells - 3.0 * sigma)),
        max(1, matlab_round(n_cells + 3.0 * sigma)) + 1,
    )
    if sigma == 0:
        weights = np.ones_like(nvector, dtype=float)
    else:
        a = 1.0 / (np.sqrt(2.0 * np.pi) * sigma)
        weights = a * np.exp(-((nvector - n_cells) ** 2) / (2.0 * sigma**2))
    return nvector, weights
