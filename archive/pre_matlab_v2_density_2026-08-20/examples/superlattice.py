from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

CACHE_DIR = ROOT / ".matplotlib-cache"
CACHE_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(CACHE_DIR))

from genl import StackDefinition, StackModel  # noqa: E402
from genl.paths import STACK_DIR  # noqa: E402


def read_data(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    try:
        values = np.loadtxt(path)
    except ValueError:
        values = np.loadtxt(path, delimiter=",")
    if values.ndim != 2 or values.shape[1] < 2:
        raise ValueError("Data file must contain 2theta and intensity columns")
    return values[:, 0], values[:, 1]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Simulate or fit a GenL superlattice definition."
    )
    parser.add_argument(
        "--superlattice",
        type=Path,
        default=STACK_DIR / "fe_v_4_28_x11.json",
        help="versioned GenL superlattice JSON file",
    )
    parser.add_argument(
        "--stack",
        dest="superlattice",
        type=Path,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--data", type=Path, help="optional 2theta/intensity data file")
    parser.add_argument(
        "--fit", action="store_true", help="fit parameters listed in the superlattice file"
    )
    parser.add_argument("--model", choices=("kinematic", "dynamic"))
    parser.add_argument("--maxiter", type=int, default=20)
    parser.add_argument("--popsize", type=int, default=8)
    parser.add_argument("--seed", type=int, default=3)
    parser.add_argument("--output-prefix", type=Path)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    definition = StackDefinition.load(args.superlattice)
    calculation = definition.calculation()
    observed = None
    if args.data:
        twotheta, observed = read_data(args.data)
        mask = (twotheta >= float(calculation["twotheta_min"])) & (
            twotheta <= float(calculation["twotheta_max"])
        )
        twotheta, observed = twotheta[mask], observed[mask]
        if not len(twotheta):
            raise ValueError("No experimental points are inside the superlattice 2theta range")
    else:
        step = float(calculation["twotheta_step"])
        twotheta = np.arange(
            float(calculation["twotheta_min"]),
            float(calculation["twotheta_max"]) + step * 0.5,
            step,
        )
    if args.fit and observed is None:
        parser.error("--fit requires --data")
    if args.fit and not len(definition.parameter_names):
        parser.error("the superlattice file does not select any fit parameters")

    model = StackModel(definition, twotheta, observed, args.model)
    parameters = definition.start
    if args.fit:
        iteration = 0

        def progress(_values: np.ndarray, convergence: float) -> bool:
            nonlocal iteration
            iteration += 1
            print(f"DE iteration {iteration}: convergence={convergence:.4g}")
            return False

        result = differential_evolution(
            model.objective,
            definition.bounds,
            x0=definition.start,
            seed=args.seed,
            maxiter=args.maxiter,
            popsize=args.popsize,
            callback=progress,
            polish=True,
            updating="immediate",
            workers=1,
        )
        parameters = result.x
        print(f"Fit cost: {result.fun:.8g}")

    calculated = model.predict(parameters)
    prefix = args.output_prefix or ROOT / "validation" / (
        f"{args.superlattice.stem}_{model.model}"
    )
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    pattern_path = Path(f"{prefix}_pattern.csv")
    resolved_path = Path(f"{prefix}_resolved_superlattice.json")
    np.savetxt(
        pattern_path,
        np.column_stack(
            [
                twotheta,
                model.q,
                np.full_like(twotheta, np.nan) if observed is None else observed,
                calculated,
            ]
        ),
        delimiter=",",
        header="twotheta_deg,Q_inv_A,observed,calculated",
        comments="",
    )
    resolved = definition.resolved_document(parameters)
    resolved["model"] = model.model
    with resolved_path.open("w", encoding="utf-8") as handle:
        json.dump(resolved, handle, indent=2)
        handle.write("\n")

    if model.last_dynamic_result is not None:
        density_path = Path(f"{prefix}_density.csv")
        density = model.last_dynamic_result
        np.savetxt(
            density_path,
            np.column_stack([density.z, density.rho_e.real, density.rho_e.imag]),
            delimiter=",",
            header="z_A,rho_real,rho_imag",
            comments="",
        )
        print(f"Density: {density_path}")

    if not args.no_plot:
        import matplotlib.pyplot as plt

        figure, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
        if observed is not None:
            axis.plot(twotheta, observed, ".", color="black", markersize=3, label="data")
        axis.plot(twotheta, calculated, color="#c62828", linewidth=1.5, label="GenL")
        axis.set_yscale("log")
        axis.set_xlabel("2theta (deg)")
        axis.set_ylabel("intensity")
        axis.grid(True, alpha=0.25)
        axis.legend()

        def theta_to_q(value):
            return 4.0 * np.pi / definition.wavelength * np.sin(np.deg2rad(value / 2.0))

        def q_to_theta(value):
            return 2.0 * np.rad2deg(
                np.arcsin(np.clip(value * definition.wavelength / (4.0 * np.pi), -1.0, 1.0))
            )

        top = axis.secondary_xaxis("top", functions=(theta_to_q, q_to_theta))
        top.set_xlabel("q (1/A)")
        plot_path = Path(f"{prefix}_pattern.png")
        figure.savefig(plot_path, dpi=180)
        plt.close(figure)
        print(f"Plot: {plot_path}")

    print(f"Pattern: {pattern_path}")
    print(f"Resolved superlattice: {resolved_path}")
    for target, value in definition.overrides(parameters).items():
        print(f"{target} = {value:.8g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
