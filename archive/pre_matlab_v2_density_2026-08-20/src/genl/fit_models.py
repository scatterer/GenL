from __future__ import annotations

import numpy as np

from . import Control, DynamicResult, DynamicWorkspace, Instrument, Layer, calc_dynamic_density
from .convolution import gauss_conv
from .kinematic import matlab_round
from .paths import FORM_FACTOR_DIR, STRUCTURE_DIR


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
        wavelength: float = 1.5406,
        propagation_backend: str = "auto",
        density_slices: int = 100,
        density_max_q0: float = 30.0,
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
        self.wavelength = wavelength
        self.propagation_backend = propagation_backend
        self.density_slices = density_slices
        self.density_max_q0 = density_max_q0
        self.q = 4.0 * np.pi / self.wavelength * sind(twotheta / 2.0)
        self.poscar_dir = STRUCTURE_DIR
        self.form_factor_dir = FORM_FACTOR_DIR
        self.last_z: np.ndarray | None = None
        self.last_rho_e: np.ndarray | None = None
        self.workspace = DynamicWorkspace()

    def _parse_optional_params(self, params: np.ndarray) -> tuple[float, float, float, float, float, float]:
        offset = 9
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

    def _single_result(
        self,
        params: np.ndarray,
        n_planes: float,
        bottom_amp: float,
        bottom_end: float,
        top_amp: float,
        top_end: float,
    ) -> DynamicResult:
        _, scale, area_scale, dinterface = params[:4]
        stack = [
            Layer(
                direction=self.substrate_direction,
                n=self.substrate_n,
                filename=self.substrate_filename,
                dinterface=self.substrate_dinterface,
                scale=params[8],
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
        return calc_dynamic_density(
            self.q,
            self.wavelength,
            stack,
            Control(pol=2, model="density"),
            Instrument(theta_m=2.0),
            poscar_dir=self.poscar_dir,
            form_factor_dir=self.form_factor_dir,
            vacuum_thick=5.0,
            slices=self.density_slices,
            max_q0=self.density_max_q0,
            step_q0=0.1,
            propagation_backend=self.propagation_backend,
            workspace=self.workspace,
        )

    def reflectivity(self, params: np.ndarray) -> np.ndarray:
        n_planes = params[0]
        bottom_amp, bottom_end, top_amp, top_end, film_sigma, substrate_sigma = (
            self._parse_optional_params(params)
        )
        if film_sigma > 0:
            n_values, weights = roughness_distribution(n_planes, film_sigma)
            results: list[DynamicResult] = []
            for n_value, weight in zip(n_values, weights):
                results.append(
                    self._single_result(
                        params, n_value, bottom_amp, bottom_end, top_amp, top_end
                    )
                )

            amplitude_s = sum(
                weight * result.amplitude_s
                for weight, result in zip(weights, results)
            )
            amplitude_p = sum(
                weight * result.amplitude_p
                for weight, result in zip(weights, results)
            )
            mono = np.cos(np.deg2rad(4.0)) ** 2
            reflectivity = (
                np.abs(amplitude_s) ** 2 + mono * np.abs(amplitude_p) ** 2
            ) / (1.0 + mono)

            longest = max(results, key=lambda result: len(result.z))
            self.last_z = longest.z
            self.last_rho_e = np.zeros_like(longest.rho_e)
            for weight, result in zip(weights, results):
                self.last_rho_e[: len(result.rho_e)] += weight * result.rho_e
        else:
            result = self._single_result(
                params, n_planes, bottom_amp, bottom_end, top_amp, top_end
            )
            reflectivity = result.refl
            self.last_z = result.z
            self.last_rho_e = result.rho_e

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
    center = max(1, matlab_round(n_planes))
    if sigma <= 0:
        return np.array([center]), np.array([1.0])
    n_values = np.arange(
        max(1, matlab_round(center - 3.0 * sigma)),
        max(1, matlab_round(center + 3.0 * sigma)) + 1,
    )
    weights = np.exp(-((n_values - center) ** 2) / (2.0 * sigma**2))
    weights = weights / np.sum(weights)
    return n_values, weights
