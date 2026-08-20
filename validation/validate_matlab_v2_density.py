from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genl import Control, DynamicWorkspace, Instrument, Layer, calc_dynamic_density  # noqa: E402
from genl.paths import FORM_FACTOR_DIR, STRUCTURE_DIR  # noqa: E402

REFERENCE_DIR = ROOT / "validation" / "matlab_v2"
WAVELENGTH = 1.54056


def layer(
    direction: int,
    count: float,
    filename: str,
    interface: float,
    scale: float,
    area_scale: float,
) -> Layer:
    return Layer(
        direction=direction,
        n=count,
        filename=filename,
        dinterface=interface,
        scale=scale,
        area_scale=area_scale,
    )


def cases() -> dict[str, tuple[np.ndarray, list[Layer]]]:
    fe_v = [layer(1, 1e6, "MgO_001_fractional.vasp", 0.279226, 1, 1)]
    for _ in range(11):
        fe_v.extend(
            [
                layer(1, 13, "V_fractional.vasp", 1.033, 1.02552, 1),
                layer(1, 2, "Fe_fractional.vasp", 1.64691, 0.97458, 1),
            ]
        )
    fe_v.append(layer(1, 14, "V_fractional.vasp", 0.997391, 1.02715, 1))
    return {
        "gaas": (
            np.arange(20, 35.0001, 0.02),
            [layer(3, 1e6, "GaAs_alt_fractional.vasp", 0, 1.001, 1.001)],
        ),
        "fe_film": (
            np.arange(58.92, 68.0001, 0.02),
            [
                layer(1, 1e6, "MgO_001_fractional.vasp", 0, 1, 1),
                layer(1, 28.5, "Fe_fractional.vasp", 1.4, 1.04, 1.1927),
            ],
        ),
        "w_film": (
            np.arange(81, 93.0001, 0.02),
            [
                layer(1, 1e6, "Al2O3_11-20_fractional.vasp", 0, 1, 1),
                layer(3, 45, "W_110_fractional.vasp", 1.4, 1, 1),
            ],
        ),
        "fe_v_superlattice": (np.arange(50, 75.0001, 0.02), fe_v),
    }


def compare_case(name: str, twotheta: np.ndarray, stack: list[Layer]) -> bool:
    reflectivity_path = REFERENCE_DIR / f"{name}_reflectivity.csv"
    density_path = REFERENCE_DIR / f"{name}_density.csv"
    if not reflectivity_path.exists() or not density_path.exists():
        print(f"{name}: missing MATLAB reference")
        return False

    reference = np.loadtxt(reflectivity_path, delimiter=",")
    density = np.loadtxt(density_path, delimiter=",")
    q = 4.0 * np.pi / WAVELENGTH * np.sin(np.deg2rad(twotheta / 2.0))
    result = calc_dynamic_density(
        q,
        WAVELENGTH,
        stack,
        Control(pol=2, model="density"),
        Instrument(theta_m=2),
        poscar_dir=STRUCTURE_DIR,
        form_factor_dir=FORM_FACTOR_DIR,
        vacuum_thick=5,
        slices=100,
        max_q0=30,
        propagation_backend="reflection",
        workspace=DynamicWorkspace(),
        density_method="analytic",
    )

    reference_refl = reference[:, 2]
    scale = max(float(np.max(np.abs(reference_refl))), float(np.max(result.refl)), 1e-18)
    reflectivity_relative = float(np.max(np.abs(result.refl - reference_refl)) / scale)
    floor = scale * 1e-14
    mask = np.maximum(result.refl, reference_refl) > floor
    reflectivity_log = float(
        np.max(
            np.abs(
                np.log10(np.maximum(result.refl[mask], floor))
                - np.log10(np.maximum(reference_refl[mask], floor))
            )
        )
    )

    reference_rho = density[:, 1] + 1j * density[:, 2]
    python_rho = np.interp(density[:, 0], result.z, result.rho_e.real) + 1j * np.interp(
        density[:, 0], result.z, result.rho_e.imag
    )
    density_scale = max(float(np.max(np.abs(reference_rho))), 1e-18)
    density_relative = float(np.max(np.abs(python_rho - reference_rho)) / density_scale)
    print(
        f"{name}: reflectivity={reflectivity_relative:.3e}, "
        f"log={reflectivity_log:.3e} decades, density={density_relative:.3e}"
    )
    return (
        reflectivity_relative <= 1e-3
        and reflectivity_log <= 1e-2
        and density_relative <= 1e-6
    )


def main() -> int:
    available = list(REFERENCE_DIR.glob("*_reflectivity.csv"))
    if not available:
        print(
            "MATLAB references are missing. Run "
            "matlab/kinematic_and_dynamic/export_subroutines_v2_reference.m first."
        )
        return 2
    results = [compare_case(name, *case) for name, case in cases().items()]
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
