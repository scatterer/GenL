from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

cache_dir = ROOT / ".matplotlib-cache"
cache_dir.mkdir(exist_ok=True)
(cache_dir / "xdg").mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))

from genl import Control, Instrument, Layer, calc_dynamic_density  # noqa: E402
from genl.paths import FORM_FACTOR_DIR, STRUCTURE_DIR  # noqa: E402


def sind(x: np.ndarray | float) -> np.ndarray | float:
    return np.sin(np.deg2rad(x))


def run_gaas_substrate(
    theta_step: float = 0.01,
    slices: int = 400,
    max_q0: float = 75.0,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    wavelength = 1.54056
    theta = np.arange(0.0, 90.0 + theta_step * 0.5, theta_step)
    q = 4.0 * np.pi / wavelength * sind(theta)
    q_calc = np.maximum(q, 1e-12)
    stack = [
        Layer(
            direction=3,
            n=1e8,
            filename="GaAs_alt_fractional.vasp",
            dinterface=0.0,
            scale=1.001,
            area_scale=1.001,
            roughness=False,
            sigma=0.0,
        )
    ]
    common = dict(
        q=q_calc,
        wavelength=wavelength,
        stack=stack,
        instrument=Instrument(theta_m=2.0),
        poscar_dir=STRUCTURE_DIR,
        form_factor_dir=FORM_FACTOR_DIR,
        vacuum_thick=20.0,
        slices=slices,
        max_q0=max_q0,
    )
    reflectivity = {
        "sigma": calc_dynamic_density(control=Control(pol=0, model="density"), **common).refl,
        "pi": calc_dynamic_density(control=Control(pol=1, model="density"), **common).refl,
        "unpolarized": calc_dynamic_density(control=Control(pol=2, model="density"), **common).refl,
    }
    return theta, q, reflectivity


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run the GaAs substrate dynamic scattering example from "
            "matlab/kinematic_and_dynamic/run_gaas_substrate.m."
        )
    )
    parser.add_argument("--theta-step", type=float, default=0.01)
    parser.add_argument("--slices", type=int, default=400)
    parser.add_argument("--max-q0", type=float, default=75.0)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    theta, q, reflectivity = run_gaas_substrate(
        theta_step=args.theta_step,
        slices=args.slices,
        max_q0=args.max_q0,
    )

    output_csv = ROOT / "validation" / "gaas_substrate_dynamic.csv"
    np.savetxt(
        output_csv,
        np.column_stack(
            [
                theta,
                q,
                reflectivity["sigma"],
                reflectivity["pi"],
                reflectivity["unpolarized"],
            ]
        ),
        delimiter=",",
        header="theta_deg,Q_inv_A,refl_sigma,refl_pi,refl_unpolarized",
        comments="",
    )

    if not args.no_plot:
        import matplotlib.pyplot as plt

        output_png = ROOT / "validation" / "gaas_substrate_dynamic.png"
        fig, axis = plt.subplots(figsize=(8, 6), constrained_layout=True)
        axis.plot(theta, reflectivity["sigma"], "-k", linewidth=2, label="density sigma")
        axis.plot(theta, reflectivity["pi"], "--k", linewidth=2, label="density pi")
        axis.plot(
            theta,
            reflectivity["unpolarized"],
            "-.k",
            linewidth=2,
            label="density partial pol",
        )
        axis.set_yscale("log")
        axis.set_xlabel("alpha_i (deg)")
        axis.set_ylabel("I/I_i")
        axis.legend(loc="best")
        axis.grid(True, alpha=0.25)
        fig.savefig(output_png, dpi=180)
        print(f"PNG: {output_png}")

    print(f"CSV: {output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
