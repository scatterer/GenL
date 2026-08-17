from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, least_squares, minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genl import Control, Instrument, Layer, calc_dynamic_density  # noqa: E402
from genl.convolution import gauss_conv  # noqa: E402
from genl.paths import EXAMPLE_DATA_DIR, FORM_FACTOR_DIR, STRUCTURE_DIR  # noqa: E402

DEFAULT_PROGRESS_PLOT = ROOT / "validation" / "fe_100a_dynamic_fit_progress.png"


def sind(x: np.ndarray | float) -> np.ndarray | float:
    return np.sin(np.deg2rad(x))


class DynamicModel:
    def __init__(
        self,
        twotheta: np.ndarray,
        observed: np.ndarray,
        include_strain: bool = False,
        include_roughness: bool = False,
        film_filename: str = "Fe_fractional.vasp",
        film_direction: int = 1,
        substrate_filename: str = "MgO_001_fractional.vasp",
        substrate_direction: int = 1,
        substrate_n: float = 1e6,
        substrate_dinterface: float = 0.0,
        substrate_scale: float = 1.0,
        substrate_area_scale: float = 1.0,
    ) -> None:
        self.twotheta = twotheta
        self.observed = observed
        self.include_strain = include_strain
        self.include_roughness = include_roughness
        self.film_filename = film_filename
        self.film_direction = film_direction
        self.substrate_filename = substrate_filename
        self.substrate_direction = substrate_direction
        self.substrate_n = substrate_n
        self.substrate_dinterface = substrate_dinterface
        self.substrate_scale = substrate_scale
        self.substrate_area_scale = substrate_area_scale
        self.wavelength = 1.54056
        self.q = 4.0 * np.pi / self.wavelength * sind(twotheta / 2.0)
        self.poscar_dir = STRUCTURE_DIR
        self.form_factor_dir = FORM_FACTOR_DIR
        self.last_z: np.ndarray | None = None
        self.last_rho_e: np.ndarray | None = None

    def _parse_optional_params(self, params: np.ndarray) -> tuple[float, float, float, float, float, float]:
        offset = 8
        if self.include_strain:
            bottom_amp, bottom_end, top_amp, top_end = params[offset : offset + 4]
            offset += 4
        else:
            bottom_amp = bottom_end = top_amp = top_end = 0.0
        if self.include_roughness:
            film_sigma, substrate_sigma = params[offset : offset + 2]
        else:
            film_sigma = substrate_sigma = 0.0
        return bottom_amp, bottom_end, top_amp, top_end, film_sigma, substrate_sigma

    def _single_reflectivity(
        self,
        params: np.ndarray,
        n_planes: float,
        bottom_amp: float,
        bottom_end: float,
        top_amp: float,
        top_end: float,
    ) -> np.ndarray:
        _, scale, area_scale, dinterface = params[:4]
        stack = [
            Layer(
                direction=self.substrate_direction,
                n=self.substrate_n,
                filename=self.substrate_filename,
                dinterface=self.substrate_dinterface,
                scale=self.substrate_scale,
                area_scale=self.substrate_area_scale,
            ),
            Layer(
                direction=self.film_direction,
                n=n_planes,
                filename=self.film_filename,
                dinterface=dinterface,
                scale=scale,
                area_scale=area_scale,
                bottom_strain_amplitude=bottom_amp,
                bottom_strain_end=bottom_end,
                top_strain_amplitude=top_amp,
                top_strain_end=top_end,
            ),
        ]
        result = calc_dynamic_density(
            self.q,
            self.wavelength,
            stack,
            Control(pol=2, model="density"),
            Instrument(theta_m=2.0),
            poscar_dir=self.poscar_dir,
            form_factor_dir=self.form_factor_dir,
            vacuum_thick=5.0,
            slices=50,
            max_q0=30.0,
            step_q0=0.1,
        )
        self.last_z = result.z
        self.last_rho_e = result.rho_e
        return result.refl

    def reflectivity(self, params: np.ndarray) -> np.ndarray:
        n_planes = params[0]
        bottom_amp, bottom_end, top_amp, top_end, film_sigma, substrate_sigma = (
            self._parse_optional_params(params)
        )
        if film_sigma > 0:
            n_values, weights = roughness_distribution(n_planes, film_sigma)
            reflectivity = np.zeros_like(self.q)
            for n_value, weight in zip(n_values, weights):
                reflectivity += weight * self._single_reflectivity(
                    params, n_value, bottom_amp, bottom_end, top_amp, top_end
                )
        else:
            reflectivity = self._single_reflectivity(
                params, n_planes, bottom_amp, bottom_end, top_amp, top_end
            )

        if substrate_sigma > 0:
            reflectivity = reflectivity * np.exp(-((self.q * substrate_sigma) ** 2))
        return reflectivity

    def predict(self, params: np.ndarray) -> np.ndarray:
        resolution, amplitude, bkg_a, bkg_b = params[4:8]
        reflectivity = self.reflectivity(params)
        broadened = gauss_conv(self.q, reflectivity, resolution)
        return amplitude * broadened + bkg_a * self.q + bkg_b

    def residual_vector(self, params: np.ndarray) -> np.ndarray:
        floor = max(np.min(self.observed[self.observed > 0]) * 0.1, 1e-12)
        predicted = np.maximum(self.predict(params), floor)
        return np.log10(predicted) - np.log10(self.observed)

    def objective(self, params: np.ndarray) -> float:
        predicted = self.predict(params)
        if not np.all(np.isfinite(predicted)):
            return np.inf
        floor = max(np.min(self.observed[self.observed > 0]) * 0.1, 1e-12)
        predicted = np.maximum(predicted, floor)
        return float(np.mean(np.abs(np.log10(predicted) - np.log10(self.observed))))


def roughness_distribution(n_planes: float, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    if sigma <= 0:
        return np.array([int(round(n_planes))]), np.array([1.0])
    n_values = np.arange(
        max(1, int(round(n_planes - 3.0 * sigma))),
        max(1, int(round(n_planes + 3.0 * sigma))) + 1,
    )
    weights = np.exp(-((n_values - n_planes) ** 2) / (2.0 * sigma**2))
    weights = weights / np.sum(weights)
    return n_values, weights


class DynamicProgressPlotter:
    def __init__(
        self,
        model: DynamicModel,
        output_path: Path,
        interval: int,
    ) -> None:
        cache_dir = ROOT / ".matplotlib-cache"
        cache_dir.mkdir(exist_ok=True)
        (cache_dir / "xdg").mkdir(exist_ok=True)
        os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
        os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))
        try:
            import matplotlib
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise RuntimeError(
                "Progress plotting requires matplotlib. Install it with "
                "`.venv/bin/python -m pip install matplotlib`."
            ) from exc

        self.model = model
        self.output_path = output_path
        self.interval = max(1, interval)
        self.count = 0
        self.evaluations: list[int] = []
        self.costs: list[float] = []
        self.plt = plt
        self.interactive_backend = "agg" not in matplotlib.get_backend().lower()
        self.plt.ion()
        self.figure, (self.loss_axis, self.fit_axis) = self.plt.subplots(
            2, 1, figsize=(9, 8), constrained_layout=True
        )

    def update(self, phase: str, params: np.ndarray, force: bool = False) -> None:
        self.count += 1
        if not force and self.count % self.interval != 0:
            return

        predicted = self.model.predict(params)
        cost = self.model.objective(params)
        self.evaluations.append(self.count)
        self.costs.append(cost)

        self.loss_axis.clear()
        self.loss_axis.plot(self.evaluations, self.costs, color="tab:blue", linewidth=1.5)
        self.loss_axis.set_xlabel("progress callback")
        self.loss_axis.set_ylabel("mean abs log10 error")
        self.loss_axis.set_title(f"{phase}: cost={cost:.5g}")
        self.loss_axis.grid(True, alpha=0.25)

        self.fit_axis.clear()
        self.fit_axis.plot(
            self.model.twotheta,
            self.model.observed,
            ".",
            color="black",
            markersize=3,
            label="data",
        )
        self.fit_axis.plot(
            self.model.twotheta,
            predicted,
            color="tab:red",
            linewidth=1.4,
            label="current dynamic fit",
        )
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


def fit_fe_100a_dynamic(
    seed: int = 20260601,
    plot_progress: bool = False,
    progress_interval: int = 1,
    progress_plot_path: Path = DEFAULT_PROGRESS_PLOT,
    maxiter: int = 18,
    popsize: int = 6,
    local_max_nfev: int = 180,
    polish_maxiter: int = 350,
    twotheta_min: float = 58.92,
    twotheta_max: float = 68.0,
) -> dict[str, float | str]:
    data = np.loadtxt(EXAMPLE_DATA_DIR / "Example_data_10nmFe.txt")
    mask = (data[:, 0] >= twotheta_min) & (data[:, 0] <= twotheta_max)
    if not np.any(mask):
        raise ValueError(
            f"No data points found in fit window {twotheta_min} <= 2theta <= {twotheta_max}"
        )
    model = DynamicModel(data[mask, 0], data[mask, 1])

    # [N, Fe scale, Fe area_scale, interface, resolution, I0, bkg_a, bkg_b]
    bounds_array = np.array(
        [
            [10.0, 40.0],
            [0.90, 1.06],
            [0.50, 1.50],
            [1.30, 1.43],
            [0.0001, 0.006],
            [50.0, 1_000_000.0],
            [0.0, 1.0],
            [0.0, 0.1],
        ],
        dtype=float,
    )
    bounds = [tuple(row) for row in bounds_array]
    start = np.array([28.5, 1.04, 1.1927, 1.40, 0.0054205, 5000.8437, 5.1469e-4, 1.2366e-7])

    plotter = (
        DynamicProgressPlotter(model, progress_plot_path, progress_interval)
        if plot_progress
        else None
    )
    if plotter is not None:
        plotter.update("initial", start, force=True)

    def de_callback(xk: np.ndarray, convergence: float | None = None) -> bool:
        if plotter is not None:
            plotter.update("differential evolution", np.asarray(xk))
        return False

    global_result = differential_evolution(
        model.objective,
        bounds,
        seed=seed,
        maxiter=maxiter,
        popsize=popsize,
        tol=2e-4,
        polish=False,
        updating="immediate",
        workers=1,
        callback=de_callback,
    )

    local_start = global_result.x if global_result.fun < model.objective(start) else start
    local_result = least_squares(
        model.residual_vector,
        local_start,
        bounds=(bounds_array[:, 0], bounds_array[:, 1]),
        xtol=1e-8,
        ftol=1e-8,
        gtol=1e-8,
        max_nfev=local_max_nfev,
    )
    if plotter is not None:
        plotter.update("least squares", local_result.x, force=True)

    def polish_callback(xk: np.ndarray) -> None:
        if plotter is not None:
            plotter.update("Powell polish", np.asarray(xk))

    polish_result = minimize(
        model.objective,
        global_result.x,
        method="Powell",
        bounds=bounds,
        callback=polish_callback,
        options={"maxiter": polish_maxiter, "xtol": 1e-6, "ftol": 1e-6},
    )

    candidates = [start, global_result.x, local_result.x, polish_result.x]
    best_params = min(candidates, key=model.objective)
    best_cost = model.objective(best_params)
    best_predicted = model.predict(best_params)
    residual = best_predicted - model.observed

    if plotter is not None:
        plotter.update("best dynamic fit", best_params, force=True)
        plotter.save()

    output_path = ROOT / "validation" / "fe_100a_dynamic_fit.csv"
    np.savetxt(
        output_path,
        np.column_stack([model.twotheta, model.q, model.observed, best_predicted, residual]),
        delimiter=",",
        header="twotheta_deg,Q_inv_A,observed_cps,fitted_cps,residual_cps",
        comments="",
    )

    return {
        "points": float(len(model.observed)),
        "twotheta_min": float(np.min(model.twotheta)),
        "twotheta_max": float(np.max(model.twotheta)),
        "n_planes": float(round(best_params[0])),
        "fe_scale": float(best_params[1]),
        "fe_area_scale": float(best_params[2]),
        "dinterface_A": float(best_params[3]),
        "resolution_fwhm_Q": float(best_params[4]),
        "scale": float(best_params[5]),
        "background_slope": float(best_params[6]),
        "background_intercept": float(best_params[7]),
        "coherent_thickness_A": float(round(best_params[0]) * 2.866 * best_params[1]),
        "mean_abs_log10_error": float(best_cost),
        "rmse_cps": float(np.sqrt(np.mean(residual**2))),
        "start_cost": float(model.objective(start)),
        "global_optimizer_cost": float(global_result.fun),
        "local_optimizer_cost": float(model.objective(local_result.x)),
        "polish_optimizer_cost": float(model.objective(polish_result.x)),
        "output_path": str(output_path),
        "progress_plot_path": str(progress_plot_path) if plot_progress else "",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Fit the 100 A Fe data with the dynamic density model.")
    parser.add_argument("--seed", type=int, default=20260601)
    parser.add_argument("--plot-progress", action="store_true")
    parser.add_argument("--progress-interval", type=int, default=1)
    parser.add_argument("--progress-plot-path", type=Path, default=DEFAULT_PROGRESS_PLOT)
    parser.add_argument("--maxiter", type=int, default=18, help="differential-evolution iterations")
    parser.add_argument("--popsize", type=int, default=6, help="differential-evolution population factor")
    parser.add_argument("--local-max-nfev", type=int, default=180, help="least-squares evaluation budget")
    parser.add_argument("--polish-maxiter", type=int, default=350, help="Powell polishing iteration budget")
    parser.add_argument("--twotheta-min", type=float, default=58.92, help="lower 2theta fit limit")
    parser.add_argument("--twotheta-max", type=float, default=68.0, help="upper 2theta fit limit")
    args = parser.parse_args()

    metrics = fit_fe_100a_dynamic(
        seed=args.seed,
        plot_progress=args.plot_progress,
        progress_interval=args.progress_interval,
        progress_plot_path=args.progress_plot_path,
        maxiter=args.maxiter,
        popsize=args.popsize,
        local_max_nfev=args.local_max_nfev,
        polish_maxiter=args.polish_maxiter,
        twotheta_min=args.twotheta_min,
        twotheta_max=args.twotheta_max,
    )
    print("100 A Fe film fit attempt: dynamic density model")
    print(f"points fitted: {int(metrics['points'])}")
    print(
        "2theta fit window: "
        f"{metrics['twotheta_min']:.3f} to {metrics['twotheta_max']:.3f} deg"
    )
    print(f"Fe coherent planes: {metrics['n_planes']:.0f}")
    print(f"Fe scale: {metrics['fe_scale']:.6f}")
    print(f"Fe area scale: {metrics['fe_area_scale']:.6f}")
    print(f"interface distance: {metrics['dinterface_A']:.6f} A")
    print(f"coherent thickness: {metrics['coherent_thickness_A']:.3f} A")
    print(f"resolution FWHM in Q: {metrics['resolution_fwhm_Q']:.6e} 1/A")
    print(f"scale: {metrics['scale']:.6e}")
    print(
        "linear background: "
        f"{metrics['background_slope']:.6e} * Q + {metrics['background_intercept']:.6e}"
    )
    print(f"mean abs log10 error: {metrics['mean_abs_log10_error']:.6e}")
    print(f"RMSE: {metrics['rmse_cps']:.6e} cps")
    print(f"start cost: {metrics['start_cost']:.6e}")
    print(f"global optimizer cost: {metrics['global_optimizer_cost']:.6e}")
    print(f"local optimizer cost: {metrics['local_optimizer_cost']:.6e}")
    print(f"polish optimizer cost: {metrics['polish_optimizer_cost']:.6e}")
    print(f"fit curve CSV: {metrics['output_path']}")
    if metrics["progress_plot_path"]:
        print(f"progress plot PNG: {metrics['progress_plot_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
