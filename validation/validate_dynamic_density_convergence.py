from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genl import DynamicWorkspace, Layer, calc_dynamic_density  # noqa: E402
from genl.paths import EXAMPLE_DATA_DIR, FORM_FACTOR_DIR, STRUCTURE_DIR  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Check analytic density slice convergence")
    parser.add_argument("--slices", type=int, default=200)
    args = parser.parse_args()

    data = np.loadtxt(EXAMPLE_DATA_DIR / "Example_data_10nmFe.txt")
    twotheta = data[(data[:, 0] >= 58.92) & (data[:, 0] <= 68.0), 0]
    wavelength = 1.54056
    q = 4.0 * np.pi / wavelength * np.sin(np.deg2rad(twotheta / 2.0))
    stack = [
        Layer(direction=1, n=1e6, filename="MgO_001_fractional.vasp"),
        Layer(
            direction=1,
            n=28.5,
            filename="Fe_fractional.vasp",
            dinterface=1.4,
            scale=1.04,
            area_scale=1.1927,
        ),
    ]
    common = {
        "poscar_dir": STRUCTURE_DIR,
        "form_factor_dir": FORM_FACTOR_DIR,
        "max_q0": 30.0,
        "density_method": "analytic",
        "propagation_backend": "reflection",
        "workspace": DynamicWorkspace(),
    }
    coarse = calc_dynamic_density(q, wavelength, stack, slices=args.slices, **common)
    fine = calc_dynamic_density(q, wavelength, stack, slices=2 * args.slices, **common)
    scale = max(float(np.max(coarse.refl)), float(np.max(fine.refl)), 1e-18)
    relative_peak = float(np.max(np.abs(coarse.refl - fine.refl)) / scale)
    floor = scale * 1e-14
    mask = np.maximum(coarse.refl, fine.refl) > floor
    log_change = float(
        np.max(
            np.abs(
                np.log10(np.maximum(coarse.refl[mask], floor))
                - np.log10(np.maximum(fine.refl[mask], floor))
            )
        )
    )
    print(
        f"{args.slices} -> {2 * args.slices} slices/cell: "
        f"relative peak change={relative_peak:.6g}, log change={log_change:.6g} decades"
    )
    return 0 if relative_peak <= 1e-3 and log_change <= 1e-2 else 1


if __name__ == "__main__":
    raise SystemExit(main())
