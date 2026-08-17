from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, least_squares, minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genl.convolution import gauss_conv  # noqa: E402
from genl.form_factors import form_factors, read_form_factor_coefficients  # noqa: E402
from genl.paths import EXAMPLE_DATA_DIR, FORM_FACTOR_DIR  # noqa: E402

DEFAULT_PROGRESS_PLOT = ROOT / "validation" / "fe_100a_fit_progress.png"


def sind(x: np.ndarray | float) -> np.ndarray | float:
    return np.sin(np.deg2rad(x))


def fe_film_shape(
    q: np.ndarray,
    wavelength: float,
    f: np.ndarray,
    mu: np.ndarray,
    d_spacing: float,
    n_planes: float,
    resolution: float,
    debye_waller_coeff: float = -0.3328,
    theta_m: float = 2.0,
) -> np.ndarray:
    n = int(round(n_planes))
    positions = np.arange(n + 1, dtype=float) * d_spacing
    debye_waller = np.exp(debye_waller_coeff * (q / (4.0 * np.pi)) ** 2)
    structure_factor = np.sum(np.exp(1j * q[:, np.newaxis] * positions), axis=1)
    intensity = np.abs(debye_waller * f * structure_factor) ** 2

    theta = np.rad2deg(np.arcsin(np.clip(q * wavelength / (4.0 * np.pi), -1.0, 1.0)))
    tau = positions[-1] * 1e-10
    absorption = 1.0 - np.exp(-2.0 * np.abs(mu) * tau / sind(theta))
    polarization = (
        1.0 + np.cos(np.deg2rad(2.0 * theta_m)) ** 2 * np.cos(np.deg2rad(2.0 * theta)) ** 2
    ) / (1.0 + np.cos(np.deg2rad(2.0 * theta_m)) ** 2)
    lorentz = 1.0 / sind(2.0 * theta)

    broadened = gauss_conv(q, intensity, resolution)
    return broadened * absorption * polarization * lorentz


def log_cost(predicted: np.ndarray, observed: np.ndarray) -> float:
    floor = max(np.min(observed[observed > 0]) * 0.1, 1e-12)
    predicted = np.maximum(predicted, floor)
    return float(np.mean(np.abs(np.log10(predicted) - np.log10(observed))))


def predict(
    params: np.ndarray,
    q: np.ndarray,
    wavelength: float,
    f: np.ndarray,
    mu: np.ndarray,
) -> np.ndarray:
    d_spacing, n_planes, resolution, amplitude, bkg_a, bkg_b = params
    shape = fe_film_shape(q, wavelength, f, mu, d_spacing, n_planes, resolution)
    return amplitude * shape + bkg_a * q + bkg_b


def residual_vector(
    params: np.ndarray,
    q: np.ndarray,
    observed: np.ndarray,
    wavelength: float,
    f: np.ndarray,
    mu: np.ndarray,
) -> np.ndarray:
    floor = max(np.min(observed[observed > 0]) * 0.1, 1e-12)
    predicted = np.maximum(predict(params, q, wavelength, f, mu), floor)
    return np.log10(predicted) - np.log10(observed)


def objective(
    params: np.ndarray,
    q: np.ndarray,
    observed: np.ndarray,
    wavelength: float,
    f: np.ndarray,
    mu: np.ndarray,
) -> float:
    predicted = predict(params, q, wavelength, f, mu)
    if not np.all(np.isfinite(predicted)):
        return np.inf
    return log_cost(predicted, observed)


class FitProgressPlotter:
    def __init__(
        self,
        twotheta: np.ndarray,
        q: np.ndarray,
        observed: np.ndarray,
        wavelength: float,
        f: np.ndarray,
        mu: np.ndarray,
        output_path: Path,
        interval: int,
    ) -> None:
        cache_dir = ROOT / ".matplotlib-cache"
        cache_dir.mkdir(exist_ok=True)
        (cache_dir / "xdg").mkdir(exist_ok=True)
        os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
        os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))
        try:
            import matplotlib.pyplot as plt
            import matplotlib
        except ImportError as exc:
            raise RuntimeError(
                "Progress plotting requires matplotlib. Install it with "
                "`.venv/bin/python -m pip install matplotlib`."
            ) from exc

        self.plt = plt
        self.interactive_backend = "agg" not in matplotlib.get_backend().lower()
        self.twotheta = twotheta
        self.q = q
        self.observed = observed
        self.wavelength = wavelength
        self.f = f
        self.mu = mu
        self.output_path = output_path
        self.interval = max(1, interval)
        self.count = 0
        self.evaluations: list[int] = []
        self.costs: list[float] = []

        self.plt.ion()
        self.figure, (self.loss_axis, self.fit_axis) = self.plt.subplots(
            2, 1, figsize=(9, 8), constrained_layout=True
        )

    def update(self, phase: str, params: np.ndarray, force: bool = False) -> None:
        self.count += 1
        if not force and self.count % self.interval != 0:
            return

        predicted = predict(params, self.q, self.wavelength, self.f, self.mu)
        cost = objective(params, self.q, self.observed, self.wavelength, self.f, self.mu)
        self.evaluations.append(self.count)
        self.costs.append(cost)

        self.loss_axis.clear()
        self.loss_axis.plot(self.evaluations, self.costs, color="tab:blue", linewidth=1.5)
        self.loss_axis.set_xlabel("progress callback")
        self.loss_axis.set_ylabel("mean abs log10 error")
        self.loss_axis.set_title(f"{phase}: cost={cost:.5g}")
        self.loss_axis.grid(True, alpha=0.25)

        self.fit_axis.clear()
        self.fit_axis.plot(self.twotheta, self.observed, ".", color="black", markersize=3, label="data")
        self.fit_axis.plot(self.twotheta, predicted, color="tab:red", linewidth=1.4, label="current fit")
        self.fit_axis.set_yscale("log")
        self.fit_axis.set_xlabel("2theta (deg)")
        self.fit_axis.set_ylabel("intensity (cps)")
        self.fit_axis.legend(loc="best")
        self.fit_axis.grid(True, alpha=0.25)

        self.figure.canvas.draw_idle()
        if self.interactive_backend:
            self.plt.pause(0.001)

    def save(self) -> None:
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self.figure.savefig(self.output_path, dpi=180)


def fit_fe_100a(
    seed: int = 20260601,
    plot_progress: bool = False,
    progress_interval: int = 1,
    progress_plot_path: Path = DEFAULT_PROGRESS_PLOT,
) -> dict[str, float]:
    wavelength = 1.5406
    peak = 64.92
    fit_range = 6.0
    thickness = 100.0
    density = 2.0 / (2.866**3)

    data = np.loadtxt(EXAMPLE_DATA_DIR / "Example_data_10nmFe.txt")
    mask = (data[:, 0] >= peak - fit_range) & (data[:, 0] <= peak + fit_range)
    twotheta = data[mask, 0]
    observed = data[mask, 1]
    q = 4.0 * np.pi / wavelength * sind(twotheta / 2.0)

    coefficients = read_form_factor_coefficients(
        26,
        wavelength,
        FORM_FACTOR_DIR,
    )
    f, _, _ = form_factors(q, coefficients, 1.0)
    classical_electron_radius = 2.81794e-15 * 1e10
    mu = 2.0 * classical_electron_radius * density * wavelength * np.imag(f) * 1e10

    d0 = wavelength / (2.0 * sind(peak / 2.0))
    bounds_array = np.array(
        [
            [
                wavelength / (2.0 * sind((peak * 1.0035) / 2.0)),
                wavelength / (2.0 * sind((peak * 0.9965) / 2.0)),
            ],
            [0.80 * thickness / d0, thickness / d0],
            [0.005 / 20.0, 0.005 * 5.0],
            [1e-3, 0.05],
            [0.0, 1.0],
            [0.0, 3.0],
        ],
        dtype=float,
    )
    bounds = [tuple(row) for row in bounds_array]

    start = np.array([d0, thickness / d0, 0.005, 0.0069, 0.0, 0.1], dtype=float)
    plotter = (
        FitProgressPlotter(
            twotheta,
            q,
            observed,
            wavelength,
            f,
            mu,
            progress_plot_path,
            progress_interval,
        )
        if plot_progress
        else None
    )

    if plotter is not None:
        plotter.update("initial", start, force=True)

    def differential_evolution_callback(xk: np.ndarray, convergence: float | None = None) -> bool:
        if plotter is not None:
            plotter.update("differential evolution", np.asarray(xk))
        return False

    global_result = differential_evolution(
        objective,
        bounds,
        args=(q, observed, wavelength, f, mu),
        seed=seed,
        maxiter=120,
        popsize=14,
        tol=1e-7,
        polish=False,
        updating="immediate",
        workers=1,
        callback=differential_evolution_callback,
    )
    local_start = global_result.x if global_result.fun < objective(start, q, observed, wavelength, f, mu) else start
    local_result = least_squares(
        residual_vector,
        local_start,
        bounds=(bounds_array[:, 0], bounds_array[:, 1]),
        args=(q, observed, wavelength, f, mu),
        xtol=1e-11,
        ftol=1e-11,
        gtol=1e-11,
        max_nfev=3000,
    )

    if plotter is not None:
        plotter.update("least squares", local_result.x, force=True)

    def polish_callback(xk: np.ndarray) -> None:
        if plotter is not None:
            plotter.update("Powell polish", np.asarray(xk))

    polish_result = minimize(
        objective,
        global_result.x,
        args=(q, observed, wavelength, f, mu),
        method="Powell",
        bounds=bounds,
        callback=polish_callback,
        options={"maxiter": 1200, "xtol": 1e-9, "ftol": 1e-9},
    )

    candidates = [global_result.x, local_result.x, polish_result.x, start]
    best_params = min(
        candidates,
        key=lambda params: objective(params, q, observed, wavelength, f, mu),
    )
    best_cost = objective(best_params, q, observed, wavelength, f, mu)
    best_predicted = predict(best_params, q, wavelength, f, mu)

    if plotter is not None:
        plotter.update("best fit", best_params, force=True)
        plotter.save()

    residual = best_predicted - observed
    output = np.column_stack([twotheta, q, observed, best_predicted, residual])
    output_path = ROOT / "validation" / "fe_100a_fit.csv"
    np.savetxt(
        output_path,
        output,
        delimiter=",",
        header="twotheta_deg,Q_inv_A,observed_cps,fitted_cps,residual_cps",
        comments="",
    )

    return {
        "points": float(len(observed)),
        "cost_log10_mae": float(best_cost),
        "d_spacing_A": float(best_params[0]),
        "n_planes": float(round(best_params[1])),
        "coherent_thickness_A": float(round(best_params[1]) * best_params[0]),
        "resolution_fwhm_Q": float(best_params[2]),
        "scale": float(best_params[3]),
        "background_slope": float(best_params[4]),
        "background_intercept": float(best_params[5]),
        "rmse_cps": float(np.sqrt(np.mean(residual**2))),
        "mean_abs_log10_error": float(best_cost),
        "global_optimizer_cost": float(global_result.fun),
        "local_optimizer_cost": float(objective(local_result.x, q, observed, wavelength, f, mu)),
        "polish_optimizer_cost": float(objective(polish_result.x, q, observed, wavelength, f, mu)),
        "output_path": str(output_path),
        "progress_plot_path": str(progress_plot_path) if plot_progress else "",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Fit the bundled 100 A Fe film data.")
    parser.add_argument("--seed", type=int, default=20260601)
    parser.add_argument(
        "--plot-progress",
        action="store_true",
        help="show and save an updating progress plot during the fit",
    )
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=1,
        help="plot every N optimizer callbacks",
    )
    parser.add_argument(
        "--progress-plot-path",
        type=Path,
        default=DEFAULT_PROGRESS_PLOT,
        help="where to save the final progress plot PNG",
    )
    args = parser.parse_args()

    metrics = fit_fe_100a(
        seed=args.seed,
        plot_progress=args.plot_progress,
        progress_interval=args.progress_interval,
        progress_plot_path=args.progress_plot_path,
    )
    print("100 A Fe film fit attempt: kinematic-only single-layer model")
    print(f"points fitted: {int(metrics['points'])}")
    print(f"d spacing: {metrics['d_spacing_A']:.6f} A")
    print(f"coherent planes: {metrics['n_planes']:.0f}")
    print(f"coherent thickness: {metrics['coherent_thickness_A']:.3f} A")
    print(f"resolution FWHM in Q: {metrics['resolution_fwhm_Q']:.6e} 1/A")
    print(f"scale: {metrics['scale']:.6e}")
    print(
        "linear background: "
        f"{metrics['background_slope']:.6e} * Q + {metrics['background_intercept']:.6e}"
    )
    print(f"mean abs log10 error: {metrics['mean_abs_log10_error']:.6e}")
    print(f"RMSE: {metrics['rmse_cps']:.6e} cps")
    print(f"global optimizer cost: {metrics['global_optimizer_cost']:.6e}")
    print(f"local optimizer cost: {metrics['local_optimizer_cost']:.6e}")
    print(f"polish optimizer cost: {metrics['polish_optimizer_cost']:.6e}")
    print(f"fit curve CSV: {metrics['output_path']}")
    if metrics["progress_plot_path"]:
        print(f"progress plot PNG: {metrics['progress_plot_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
