from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genl import Control, Instrument, Layer, calc_kinematic  # noqa: E402
from genl.debye import debye_waller_prefactor  # noqa: E402
from genl.form_factors import form_factors, read_form_factor_coefficients  # noqa: E402
from genl.paths import FORM_FACTOR_DIR, STRUCTURE_DIR  # noqa: E402
from genl.poscar import read_poscar  # noqa: E402

R0_ANGSTROM = 2.814042735053330e-5


def brute_force_fe_kinematic(
    theta: np.ndarray,
    wavelength: float,
    n_cells: int,
) -> np.ndarray:
    """Independent finite sum for the single-layer Fe example.

    This intentionally does not call `calc_kinematic`; it expands every unit
    cell and atom explicitly and sums the scattering amplitude atom by atom.
    """

    q = 4.0 * np.pi / wavelength * np.sin(np.deg2rad(theta))
    structure = read_poscar(
        "Fe_fractional.vasp", STRUCTURE_DIR
    )

    scaling = structure.a3
    lattice_parameter = np.linalg.norm(scaling)
    unit_cell_volume = abs(np.dot(scaling, np.cross(structure.a1, structure.a2)))
    atom_z = structure.positions @ scaling

    coeffs = read_form_factor_coefficients(
        26,
        wavelength,
        FORM_FACTOR_DIR,
    )
    f_q, _, _ = form_factors(q, coeffs, 1.0)
    b_factor = debye_waller_prefactor(26)
    f_q = f_q * np.exp(-b_factor * (q / (4.0 * np.pi)) ** 2) / q

    positions = np.concatenate(
        [atom_z + cell_index * lattice_parameter for cell_index in range(n_cells)]
    )
    amplitude = (
        -1j
        * 4.0
        * np.pi
        * R0_ANGSTROM
        * lattice_parameter
        / unit_cell_volume
        * np.sum(f_q[:, np.newaxis] * np.exp(1j * q[:, np.newaxis] * positions), axis=1)
    )
    return np.abs(amplitude) ** 2


def translated_fe_kinematic(
    theta: np.ndarray,
    wavelength: float,
    n_cells: int,
) -> np.ndarray:
    q = 4.0 * np.pi / wavelength * np.sin(np.deg2rad(theta))
    layer = Layer(direction=3, n=n_cells, filename="Fe_fractional.vasp")
    result = calc_kinematic(
        q,
        wavelength,
        [layer],
        Control(pol=0),
        Instrument(theta=theta),
        poscar_dir=STRUCTURE_DIR,
        form_factor_dir=FORM_FACTOR_DIR,
    )
    return result.refl


def validate_fe_kinematic() -> dict[str, float]:
    wavelength = 1.54056
    n_cells = 40
    theta = np.arange(20.0, 90.0, 0.02)

    translated = translated_fe_kinematic(theta, wavelength, n_cells)
    brute_force = brute_force_fe_kinematic(theta, wavelength, n_cells)

    abs_error = np.abs(translated - brute_force)
    rel_error = abs_error / np.maximum(np.abs(brute_force), np.finfo(float).tiny)
    worst_rel_idx = int(np.argmax(rel_error))

    return {
        "points": float(theta.size),
        "max_abs_error": float(np.max(abs_error)),
        "max_rel_error": float(np.max(rel_error)),
        "rms_rel_error": float(np.sqrt(np.mean(rel_error**2))),
        "theta_at_max_rel_error": float(theta[worst_rel_idx]),
        "intensity_at_max_rel_error": float(brute_force[worst_rel_idx]),
        "translated_min": float(np.min(translated)),
        "translated_max": float(np.max(translated)),
        "brute_force_min": float(np.min(brute_force)),
        "brute_force_max": float(np.max(brute_force)),
    }


def main() -> int:
    metrics = validate_fe_kinematic()
    print("Fe kinematic validation: translated closed form vs brute-force finite sum")
    print(f"points: {int(metrics['points'])}")
    print(f"max_abs_error: {metrics['max_abs_error']:.6e}")
    print(f"max_rel_error: {metrics['max_rel_error']:.6e}")
    print(f"rms_rel_error: {metrics['rms_rel_error']:.6e}")
    print(
        "worst relative point: "
        f"theta={metrics['theta_at_max_rel_error']:.4f} deg, "
        f"I={metrics['intensity_at_max_rel_error']:.6e}"
    )
    print(
        "translated range: "
        f"{metrics['translated_min']:.6e} to {metrics['translated_max']:.6e}"
    )
    print(
        "brute-force range: "
        f"{metrics['brute_force_min']:.6e} to {metrics['brute_force_max']:.6e}"
    )

    if metrics["max_abs_error"] > 1e-14 or metrics["max_rel_error"] > 1e-9:
        print("FAILED: error exceeded validation tolerance")
        return 1
    print("PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
