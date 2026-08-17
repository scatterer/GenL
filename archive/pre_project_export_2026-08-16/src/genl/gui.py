from __future__ import annotations

import os
import queue
import threading
import tkinter as tk
import webbrowser
from dataclasses import dataclass
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CU_K_ALPHA_WAVELENGTH = 1.5406
STRAIN_KEYS = ("bottom_amplitude", "bottom_extent", "top_amplitude", "top_extent")
GENL_DOI_URL = "https://doi.org/10.1107/S1600576726002566"
GENL_CITATION = (
    "GenL Python version by Vassilios Kapaklis and Gunnar K. Palsson. "
    f"Please cite: J. Appl. Cryst. 59, 968-977 (2026), {GENL_DOI_URL}"
)
UI_COLORS = {
    "window": "#f3f5f7",
    "panel": "#ffffff",
    "setup": "#7a8793",
    "kinematic": "#4f8cc9",
    "film": "#4f9f68",
    "fit": "#d5972c",
    "substrate": "#8a7bb8",
    "strain": "#3f8f8a",
    "roughness": "#c9677d",
    "simulate": "#2f6fad",
    "run": "#2e8540",
    "stop": "#b83232",
    "entry": "#ffffff",
    "entry_disabled": "#eceff1",
    "watermark": "#6b7280",
}

cache_dir = ROOT / ".matplotlib-cache"
cache_dir.mkdir(exist_ok=True)
(cache_dir / "xdg").mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from scipy.optimize import differential_evolution, least_squares, minimize

from .convolution import gauss_conv  # noqa: E402
from .fit_models import DynamicModel, roughness_distribution  # noqa: E402
from .form_factors import form_factors, read_form_factor_coefficients  # noqa: E402
from .kinematic import matlab_round  # noqa: E402
from .poscar import read_poscar  # noqa: E402


def sind(x: np.ndarray | float) -> np.ndarray | float:
    return np.sin(np.deg2rad(x))


def q_from_twotheta(twotheta: np.ndarray | float, wavelength: float) -> np.ndarray | float:
    return 4.0 * np.pi / wavelength * sind(np.asarray(twotheta) / 2.0)


def twotheta_from_q(q: np.ndarray | float, wavelength: float) -> np.ndarray | float:
    q_array = np.asarray(q, dtype=float)
    argument = q_array * wavelength / (4.0 * np.pi)
    if np.any(argument > 1.0):
        raise ValueError("q limit is too large for the selected wavelength")
    return 2.0 * np.rad2deg(np.arcsin(np.clip(argument, -1.0, 1.0)))


@dataclass(frozen=True)
class SampleConfig:
    name: str
    data_file: str
    film_filename: str
    atomic_number: int
    d_spacing: float
    thickness: float
    debye_waller_coeff: float
    dynamic_direction: int
    twotheta_min: float
    twotheta_max: float
    substrate_filename: str = "MgO_001_fractional.vasp"
    substrate_direction: int = 1
    substrate_n: float = 1e6
    substrate_dinterface: float = 0.0
    substrate_scale: float = 1.0
    substrate_area_scale: float = 1.0

    @property
    def peak_twotheta(self) -> float:
        return float(
            2.0
            * np.rad2deg(
                np.arcsin(CU_K_ALPHA_WAVELENGTH / (2.0 * self.d_spacing))
            )
        )


SAMPLES: dict[str, SampleConfig] = {
    "Fe 10 nm": SampleConfig(
        name="Fe 10 nm",
        data_file="Example_data_10nmFe.txt",
        film_filename="Fe_fractional.vasp",
        atomic_number=26,
        d_spacing=1.433,
        thickness=100.0,
        debye_waller_coeff=-0.34,
        dynamic_direction=1,
        twotheta_min=58.92,
        twotheta_max=68.0,
    ),
    "W 10 nm": SampleConfig(
        name="W 10 nm",
        data_file="Example_data_10nmW.txt",
        film_filename="W_110_fractional.vasp",
        atomic_number=74,
        d_spacing=1.119,
        thickness=100.0,
        debye_waller_coeff=-0.161,
        dynamic_direction=3,
        twotheta_min=81.0,
        twotheta_max=93.0,
        substrate_filename="Al2O3_11-20_fractional.vasp",
        substrate_direction=1,
    ),
    "V 10 nm": SampleConfig(
        name="V 10 nm",
        data_file="Example_data_10nmV.txt",
        film_filename="V_fractional.vasp",
        atomic_number=23,
        d_spacing=1.515,
        thickness=105.0,
        debye_waller_coeff=-0.58,
        dynamic_direction=1,
        twotheta_min=55.0,
        twotheta_max=65.0,
    ),
}


POSCAR_FILES = tuple(
    sorted(path.name for path in (ROOT / "kinematic_and_dynamic" / "POSCAR").glob("*.vasp"))
)


def default_data_path(sample: SampleConfig) -> Path:
    return ROOT / "kinematic_and_dynamic" / "examples" / sample.data_file


def sample_from_data_path(path: str | Path, fallback: SampleConfig) -> SampleConfig:
    data_name = Path(path).name.lower()
    for sample in SAMPLES.values():
        if data_name == sample.data_file.lower():
            return sample
    return fallback


def resolve_data_path(path_text: str) -> Path:
    path = Path(path_text).expanduser()
    if not path.is_absolute():
        path = ROOT / path
    return path


def unit_cell_density(sample: SampleConfig) -> float:
    structure = read_poscar(
        sample.film_filename, ROOT / "kinematic_and_dynamic" / "POSCAR"
    )
    volume = abs(np.dot(structure.a1, np.cross(structure.a2, structure.a3)))
    return float(np.sum(structure.type_counts) / volume)


def propagation_period(sample: SampleConfig) -> float:
    return propagation_period_for(sample.film_filename, sample.dynamic_direction)


def propagation_period_for(filename: str, direction: int) -> float:
    structure = read_poscar(filename, ROOT / "kinematic_and_dynamic" / "POSCAR")
    if direction not in (1, 2, 3):
        raise ValueError("Layer direction must be 1, 2, or 3")
    axis = [structure.a1, structure.a2, structure.a3][direction - 1]
    return float(np.linalg.norm(axis))


@dataclass(frozen=True)
class FitUpdate:
    phase: str
    cost: float
    twotheta: np.ndarray
    q: np.ndarray
    observed: np.ndarray
    predicted: np.ndarray
    params: np.ndarray
    density_z: np.ndarray | None = None
    density_rho_e: np.ndarray | None = None


@dataclass
class RangeIndicator:
    name: str
    canvas: tk.Canvas
    start_var: tk.StringVar
    min_var: tk.StringVar
    max_var: tk.StringVar
    fit_var: tk.StringVar
    fit_enabled_var: tk.BooleanVar | None = None


class FitCancelled(Exception):
    pass


def strained_plane_positions(
    d_spacing: float,
    n_planes: float,
    bottom_amplitude: float,
    bottom_end: float,
    top_amplitude: float,
    top_end: float,
) -> np.ndarray:
    positions = np.arange(matlab_round(n_planes) + 1, dtype=float) * d_spacing
    if len(positions) <= 1:
        return positions

    strained = np.zeros_like(positions)
    displacement = np.zeros_like(positions)
    bottom_boundary = min(max(matlab_round(bottom_end), 0), len(positions) - 1)
    top_boundary = len(positions) - min(max(matlab_round(top_end), 0), len(positions))
    if top_boundary <= 0:
        top_boundary = len(positions) - 1

    for idx in range(1, len(positions)):
        delta = positions[idx] - positions[idx - 1]
        if idx + 1 > bottom_end or bottom_end <= 0:
            displacement[idx] = delta
        else:
            displacement[idx] = delta + bottom_amplitude * (
                positions[bottom_boundary] - positions[idx]
            )

        if top_end > 0 and idx + 1 > top_boundary:
            displacement[idx] += top_amplitude * (positions[idx] - positions[top_boundary])

        strained[idx] = strained[idx - 1] + displacement[idx]
    return strained


class KinematicModel:
    def __init__(
        self,
        twotheta: np.ndarray,
        observed: np.ndarray,
        sample: SampleConfig,
        debye_waller_coeff: float | None = None,
        include_strain: bool = False,
        include_roughness: bool = False,
        wavelength: float = CU_K_ALPHA_WAVELENGTH,
    ) -> None:
        self.twotheta = twotheta
        self.observed = observed
        self.sample = sample
        self.debye_waller_coeff = (
            sample.debye_waller_coeff
            if debye_waller_coeff is None
            else debye_waller_coeff
        )
        self.include_strain = include_strain
        self.include_roughness = include_roughness
        self.wavelength = wavelength
        self.thickness = sample.thickness
        self.density = unit_cell_density(sample)
        self.q = q_from_twotheta(twotheta, self.wavelength)
        coefficients = read_form_factor_coefficients(
            sample.atomic_number,
            self.wavelength,
            ROOT / "kinematic_and_dynamic" / "Form_Factor_and_Elemental_data",
        )
        self.f, _, _ = form_factors(self.q, coefficients, 1.0)
        classical_electron_radius = 2.81794e-15 * 1e10
        self.mu = (
            2.0
            * classical_electron_radius
            * self.density
            * self.wavelength
            * np.imag(self.f)
            * 1e10
        )

    def _parse_optional_params(self, params: np.ndarray) -> tuple[float, float, float, float, float, float]:
        offset = 6
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

    def _single_shape(
        self,
        params: np.ndarray,
        n_planes: float,
        bottom_amp: float,
        bottom_end: float,
        top_amp: float,
        top_end: float,
    ) -> np.ndarray:
        d_spacing, _, resolution, _, _, _ = params[:6]
        positions = strained_plane_positions(
            d_spacing, n_planes, bottom_amp, bottom_end, top_amp, top_end
        )
        debye_waller = np.exp(self.debye_waller_coeff * (self.q / (4.0 * np.pi)) ** 2)
        structure_factor = np.sum(
            np.exp(1j * self.q[:, np.newaxis] * positions[np.newaxis, :]), axis=1
        )
        intensity = np.abs(debye_waller * self.f * structure_factor) ** 2
        theta = np.rad2deg(
            np.arcsin(np.clip(self.q * self.wavelength / (4.0 * np.pi), -1.0, 1.0))
        )
        tau = positions[-1] * 1e-10
        absorption = 1.0 - np.exp(-2.0 * np.abs(self.mu) * tau / sind(theta))
        polarization = (
            1.0
            + np.cos(np.deg2rad(4.0)) ** 2 * np.cos(np.deg2rad(2.0 * theta)) ** 2
        ) / (1.0 + np.cos(np.deg2rad(4.0)) ** 2)
        lorentz = 1.0 / sind(2.0 * theta)
        return gauss_conv(self.q, intensity, resolution) * absorption * polarization * lorentz

    def predict(self, params: np.ndarray) -> np.ndarray:
        d_spacing, n_planes, resolution, amplitude, bkg_a, bkg_b = params[:6]
        bottom_amp, bottom_end, top_amp, top_end, film_sigma, substrate_sigma = (
            self._parse_optional_params(params)
        )
        if film_sigma > 0:
            n_values, weights = roughness_distribution(n_planes, film_sigma)
            shape = np.zeros_like(self.q)
            for n_value, weight in zip(n_values, weights):
                shape += weight * self._single_shape(
                    params, n_value, bottom_amp, bottom_end, top_amp, top_end
                )
        else:
            shape = self._single_shape(params, n_planes, bottom_amp, bottom_end, top_amp, top_end)

        if substrate_sigma > 0:
            shape = shape * np.exp(-((self.q * substrate_sigma) ** 2))
        return amplitude * shape + bkg_a * self.q + bkg_b

    def _residual_for_prediction(self, predicted: np.ndarray) -> np.ndarray:
        floor = max(np.min(self.observed[self.observed > 0]) * 0.1, 1e-12)
        predicted = np.maximum(predicted, floor)
        return np.log10(predicted) - np.log10(self.observed)

    def _objective_for_prediction(self, predicted: np.ndarray) -> float:
        if not np.all(np.isfinite(predicted)):
            return np.inf
        return float(np.mean(np.abs(self._residual_for_prediction(predicted))))

    def residual_vector(self, params: np.ndarray) -> np.ndarray:
        return self._residual_for_prediction(self.predict(params))

    def objective(self, params: np.ndarray) -> float:
        return self._objective_for_prediction(self.predict(params))


def read_experimental_data(data_path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    path = resolve_data_path(str(data_path))
    if not path.is_file():
        raise ValueError(f"Data file not found: {path}")
    try:
        data = np.loadtxt(path)
    except ValueError:
        data = np.loadtxt(path, delimiter=",")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("Data file must contain at least two columns: 2theta and intensity")
    return data[:, 0], data[:, 1]


def load_sample_data(
    twotheta_min: float,
    twotheta_max: float,
    data_path: str | Path,
) -> tuple[np.ndarray, np.ndarray]:
    twotheta, observed = read_experimental_data(data_path)
    mask = (twotheta >= twotheta_min) & (twotheta <= twotheta_max)
    if not np.any(mask):
        raise ValueError("No data points found inside the selected 2theta window")
    return twotheta[mask], observed[mask]


def kinematic_bounds_and_start(
    sample: SampleConfig,
    include_strain: bool,
    include_roughness: bool,
    roughness_settings: dict[str, tuple[float, float, float]] | None = None,
    kinematic_settings: dict[str, tuple[float, float, float]] | None = None,
    strain_settings: dict[str, tuple[float, float, float]] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    if kinematic_settings is None:
        wavelength = CU_K_ALPHA_WAVELENGTH
        d0 = sample.d_spacing
        peak = float(
            2.0
            * np.rad2deg(
                np.arcsin(np.clip(wavelength / (2.0 * d0), -1.0, 1.0))
            )
        )
        d_spacing = (
            d0,
            wavelength / (2.0 * sind((peak * 1.0035) / 2.0)),
            wavelength / (2.0 * sind((peak * 0.9965) / 2.0)),
        )
        planes = (sample.thickness / d0, 0.60 * sample.thickness / d0, 1.20 * sample.thickness / d0)
        resolution = (0.005, 0.005 / 20.0, 0.005 * 5.0)
        scale = (0.0069, 1e-3, 0.05)
        bkg_a = (0.0, 0.0, 1.0)
        bkg_b = (0.1, 0.0, 3.0)
    else:
        d_spacing = kinematic_settings["d_spacing"]
        planes = kinematic_settings["planes"]
        resolution = kinematic_settings["resolution"]
        scale = kinematic_settings["scale"]
        bkg_a = kinematic_settings["bkg_a"]
        bkg_b = kinematic_settings["bkg_b"]
    bounds_array = np.array(
        [
            [d_spacing[1], d_spacing[2]],
            [planes[1], planes[2]],
            [resolution[1], resolution[2]],
            [scale[1], scale[2]],
            [bkg_a[1], bkg_a[2]],
            [bkg_b[1], bkg_b[2]],
        ],
        dtype=float,
    )
    start = np.array(
        [d_spacing[0], planes[0], resolution[0], scale[0], bkg_a[0], bkg_b[0]],
        dtype=float,
    )
    if include_strain:
        settings = strain_settings or {
            "bottom_amplitude": (0.0, -0.2, 0.2),
            "bottom_extent": (0.0, 0.0, 20.0),
            "top_amplitude": (0.0, -0.2, 0.2),
            "top_extent": (0.0, 0.0, 20.0),
        }
        strain_values = [settings[key] for key in STRAIN_KEYS]
        strain_bounds = np.array([[value[1], value[2]] for value in strain_values], dtype=float)
        bounds_array = np.vstack([bounds_array, strain_bounds])
        start = np.concatenate([start, [value[0] for value in strain_values]])
    if include_roughness:
        settings = roughness_settings or {
            "film": (0.0, 0.0, 5.0),
            "substrate": (0.0, 0.0, 5.0),
        }
        roughness_bounds = np.array(
            [[settings["film"][1], settings["film"][2]], [settings["substrate"][1], settings["substrate"][2]]],
            dtype=float,
        )
        bounds_array = np.vstack([bounds_array, roughness_bounds])
        start = np.concatenate([start, [settings["film"][0], settings["substrate"][0]]])
    return bounds_array, start


def dynamic_bounds_and_start(
    sample: SampleConfig,
    include_strain: bool,
    include_roughness: bool,
    roughness_settings: dict[str, tuple[float, float, float]] | None = None,
    film_settings: dict[str, object] | None = None,
    dynamic_fit_settings: dict[str, tuple[float, float, float]] | None = None,
    strain_settings: dict[str, tuple[float, float, float]] | None = None,
    substrate_settings: dict[str, object] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    if film_settings is None:
        n0 = sample.thickness / propagation_period(sample)
        film_n = (n0, max(1.0, 0.45 * n0), 1.35 * n0)
        film_scale = (1.0, 0.90, 1.10)
        film_area_scale = (1.0, 0.50, 1.50)
        film_dinterface = (1.35, 1.30, 1.43)
    else:
        film_n = film_settings["n"]
        film_scale = film_settings["scale"]
        film_area_scale = film_settings["area_scale"]
        film_dinterface = film_settings["dinterface"]
    if dynamic_fit_settings is None:
        resolution = (0.0054205, 0.0001, 0.006)
        intensity_scale = (5000.8437, 50.0, 1_000_000.0)
        bkg_a = (5.1469e-4, 0.0, 1.0)
        bkg_b = (1.2366e-7, 0.0, 0.1)
    else:
        resolution = dynamic_fit_settings["resolution"]
        intensity_scale = dynamic_fit_settings["intensity_scale"]
        bkg_a = dynamic_fit_settings["bkg_a"]
        bkg_b = dynamic_fit_settings["bkg_b"]
    substrate_scale = (
        (sample.substrate_scale, sample.substrate_scale * 0.995, sample.substrate_scale * 1.005)
        if substrate_settings is None
        else substrate_settings["scale"]
    )
    bounds_array = np.array(
        [
            [film_n[1], film_n[2]],
            [film_scale[1], film_scale[2]],
            [film_area_scale[1], film_area_scale[2]],
            [film_dinterface[1], film_dinterface[2]],
            [resolution[1], resolution[2]],
            [intensity_scale[1], intensity_scale[2]],
            [bkg_a[1], bkg_a[2]],
            [bkg_b[1], bkg_b[2]],
            [substrate_scale[1], substrate_scale[2]],
        ],
        dtype=float,
    )
    start = np.array(
        [
            film_n[0],
            film_scale[0],
            film_area_scale[0],
            film_dinterface[0],
            resolution[0],
            intensity_scale[0],
            bkg_a[0],
            bkg_b[0],
            substrate_scale[0],
        ],
        dtype=float,
    )
    if include_strain:
        settings = strain_settings or {
            "bottom_amplitude": (0.0, -0.4, 0.4),
            "bottom_extent": (0.0, 0.0, 20.0),
            "top_amplitude": (0.0, -0.4, 0.4),
            "top_extent": (0.0, 0.0, 20.0),
        }
        strain_values = [settings[key] for key in STRAIN_KEYS]
        strain_bounds = np.array([[value[1], value[2]] for value in strain_values], dtype=float)
        bounds_array = np.vstack([bounds_array, strain_bounds])
        start = np.concatenate([start, [value[0] for value in strain_values]])
    if include_roughness:
        settings = roughness_settings or {
            "film": (0.0, 0.0, 5.0),
            "substrate": (0.0, 0.0, 5.0),
        }
        roughness_bounds = np.array(
            [[settings["film"][1], settings["film"][2]], [settings["substrate"][1], settings["substrate"][2]]],
            dtype=float,
        )
        bounds_array = np.vstack([bounds_array, roughness_bounds])
        start = np.concatenate([start, [settings["film"][0], settings["substrate"][0]]])
    return bounds_array, start


def validate_start_min_max(name: str, start: float, lower: float, upper: float) -> None:
    if lower > upper:
        raise ValueError(f"{name}: min must be <= max")
    if not lower <= start <= upper:
        raise ValueError(f"{name}: start must be between min and max")


def validate_direction(name: str, direction: int) -> None:
    if direction not in (1, 2, 3):
        raise ValueError(f"{name}: direction must be 1, 2, or 3")


def validate_poscar(name: str, filename: str) -> None:
    if filename not in POSCAR_FILES:
        raise ValueError(f"{name}: choose an available POSCAR file")


def optional_param_summary(
    params: np.ndarray,
    base_len: int,
    include_strain: bool,
    include_roughness: bool,
) -> str:
    offset = base_len
    text = ""
    if include_strain:
        text += (
            f"bottom strain: amp={params[offset]:.6g}, end={params[offset + 1]:.3f}\n"
            f"top strain: amp={params[offset + 2]:.6g}, end={params[offset + 3]:.3f}\n"
        )
        offset += 4
    if include_roughness:
        film_unit = "planes" if base_len == 6 else "layers"
        text += (
            f"film roughness sigma: {params[offset]:.6g} {film_unit}\n"
            f"substrate/interface roughness sigma: {params[offset + 1]:.6g} A\n"
        )
    return text


def summarize_fit(
    sample: SampleConfig,
    model_name: str,
    params: np.ndarray,
    cost: float,
    rmse: float,
    include_strain: bool,
    include_roughness: bool,
    film_settings: dict[str, object] | None = None,
    substrate_settings: dict[str, object] | None = None,
) -> str:
    if model_name == "Kinematic":
        return (
            f"{sample.name} kinematic fit\n"
            f"d spacing: {params[0]:.6f} A\n"
            f"coherent planes: {matlab_round(params[1])}\n"
            f"coherent thickness: {matlab_round(params[1]) * params[0]:.3f} A\n"
            f"resolution FWHM Q: {params[2]:.6e} 1/A\n"
            f"scale: {params[3]:.6e}\n"
            f"background: {params[4]:.6e} * Q + {params[5]:.6e}\n"
            + optional_param_summary(params, 6, include_strain, include_roughness)
            + f"mean abs log10 error: {cost:.6e}\n"
            f"RMSE: {rmse:.6e} cps"
        )

    dynamic_period = (
        propagation_period_for(str(film_settings["filename"]), int(film_settings["direction"]))
        if film_settings is not None
        else propagation_period(sample)
    )
    layer_text = ""
    if film_settings is not None:
        layer_text += (
            f"film POSCAR: {film_settings['filename']} direction {film_settings['direction']}\n"
        )
    if substrate_settings is not None:
        layer_text += (
            f"substrate POSCAR: {substrate_settings['filename']} "
            f"direction {substrate_settings['direction']}\n"
        )

    return (
        f"{sample.name} dynamic fit\n"
        + layer_text
        + f"coherent planes: {matlab_round(params[0])}\n"
        f"film scale: {params[1]:.6f}\n"
        f"film area scale: {params[2]:.6f}\n"
        f"interface distance: {params[3]:.6f} A\n"
        f"coherent thickness: {matlab_round(params[0]) * dynamic_period * params[1]:.3f} A\n"
        f"resolution FWHM Q: {params[4]:.6e} 1/A\n"
        f"scale: {params[5]:.6e}\n"
        f"background: {params[6]:.6e} * Q + {params[7]:.6e}\n"
        + f"substrate lattice scale: {params[8]:.8f}\n"
        + optional_param_summary(params, 9, include_strain, include_roughness)
        + f"mean abs log10 error: {cost:.6e}\n"
        f"RMSE: {rmse:.6e} cps"
    )


class FitApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("GenL Film Fitting")
        self.queue: queue.Queue[tuple[str, object]] = queue.Queue()
        self.running = False
        self.stop_event = threading.Event()
        self.history_x: list[int] = []
        self.history_y: list[float] = []
        self.last_update: FitUpdate | None = None
        self.active_fit_config: dict[str, object] | None = None
        self.preview_after_id: str | None = None
        self.preview_data_path: Path | None = None
        self.updating_twotheta_window = False

        self.sample_var = tk.StringVar(value="Fe 10 nm")
        self.model_var = tk.StringVar(value="Dynamic")
        self.data_path_var = tk.StringVar(value=str(default_data_path(SAMPLES["Fe 10 nm"])))
        self.wavelength_var = tk.StringVar(value=f"{CU_K_ALPHA_WAVELENGTH:.5g}")
        self.axis_var = tk.StringVar(value="2\u03b8")
        self.axis_mode = "twotheta"
        self.min_label_var = tk.StringVar(value="2\u03b8 min (deg)")
        self.max_label_var = tk.StringVar(value="2\u03b8 max (deg)")
        self.min_var = tk.StringVar(value="")
        self.max_var = tk.StringVar(value="")
        self.seed_var = tk.StringVar(value="20260601")
        self.maxiter_var = tk.StringVar(value="18")
        self.popsize_var = tk.StringVar(value="6")
        self.local_var = tk.StringVar(value="180")
        self.polish_var = tk.StringVar(value="350")
        self.interval_var = tk.StringVar(value="1")
        self.strain_var = tk.BooleanVar(value=False)
        self.roughness_var = tk.BooleanVar(value=False)
        self.kin_d_start_var = tk.StringVar(value="")
        self.kin_d_min_var = tk.StringVar(value="")
        self.kin_d_max_var = tk.StringVar(value="")
        self.kin_planes_start_var = tk.StringVar(value="")
        self.kin_planes_min_var = tk.StringVar(value="")
        self.kin_planes_max_var = tk.StringVar(value="")
        self.kin_resolution_start_var = tk.StringVar(value="0.005")
        self.kin_resolution_min_var = tk.StringVar(value="0.00025")
        self.kin_resolution_max_var = tk.StringVar(value="0.025")
        self.kin_scale_start_var = tk.StringVar(value="0.0069")
        self.kin_scale_min_var = tk.StringVar(value="0.001")
        self.kin_scale_max_var = tk.StringVar(value="0.05")
        self.kin_bkg_a_start_var = tk.StringVar(value="0.0")
        self.kin_bkg_a_min_var = tk.StringVar(value="0.0")
        self.kin_bkg_a_max_var = tk.StringVar(value="1.0")
        self.kin_bkg_b_start_var = tk.StringVar(value="0.1")
        self.kin_bkg_b_min_var = tk.StringVar(value="0.0")
        self.kin_bkg_b_max_var = tk.StringVar(value="3.0")
        self.kin_debye_var = tk.StringVar(value=f"{SAMPLES['Fe 10 nm'].debye_waller_coeff:.6g}")
        self.kin_d_fit_var = tk.StringVar(value="")
        self.kin_planes_fit_var = tk.StringVar(value="")
        self.kin_resolution_fit_var = tk.StringVar(value="")
        self.kin_scale_fit_var = tk.StringVar(value="")
        self.kin_bkg_a_fit_var = tk.StringVar(value="")
        self.kin_bkg_b_fit_var = tk.StringVar(value="")
        self.kin_debye_fit_var = tk.StringVar(value="")
        self.kin_d_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_planes_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_resolution_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_scale_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_bkg_a_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_bkg_b_fit_enabled_var = tk.BooleanVar(value=True)
        self.film_filename_var = tk.StringVar(value=SAMPLES["Fe 10 nm"].film_filename)
        self.film_direction_var = tk.StringVar(value=str(SAMPLES["Fe 10 nm"].dynamic_direction))
        self.film_n_start_var = tk.StringVar(value="")
        self.film_n_min_var = tk.StringVar(value="")
        self.film_n_max_var = tk.StringVar(value="")
        self.film_scale_start_var = tk.StringVar(value="1.0")
        self.film_scale_min_var = tk.StringVar(value="0.9")
        self.film_scale_max_var = tk.StringVar(value="1.1")
        self.film_area_start_var = tk.StringVar(value="1.0")
        self.film_area_min_var = tk.StringVar(value="0.5")
        self.film_area_max_var = tk.StringVar(value="1.5")
        self.film_interface_start_var = tk.StringVar(value="1.35")
        self.film_interface_min_var = tk.StringVar(value="1.3")
        self.film_interface_max_var = tk.StringVar(value="1.43")
        self.film_n_fit_var = tk.StringVar(value="")
        self.film_scale_fit_var = tk.StringVar(value="")
        self.film_area_fit_var = tk.StringVar(value="")
        self.film_interface_fit_var = tk.StringVar(value="")
        self.film_n_fit_enabled_var = tk.BooleanVar(value=True)
        self.film_scale_fit_enabled_var = tk.BooleanVar(value=True)
        self.film_area_fit_enabled_var = tk.BooleanVar(value=True)
        self.film_interface_fit_enabled_var = tk.BooleanVar(value=True)
        self.dynamic_resolution_start_var = tk.StringVar(value="0.0054205")
        self.dynamic_resolution_min_var = tk.StringVar(value="0.0001")
        self.dynamic_resolution_max_var = tk.StringVar(value="0.006")
        self.dynamic_resolution_fit_var = tk.StringVar(value="")
        self.dynamic_intensity_start_var = tk.StringVar(value="5000.8437")
        self.dynamic_intensity_min_var = tk.StringVar(value="50.0")
        self.dynamic_intensity_max_var = tk.StringVar(value="1000000.0")
        self.dynamic_intensity_fit_var = tk.StringVar(value="")
        self.dynamic_bkg_a_start_var = tk.StringVar(value="5.1469e-4")
        self.dynamic_bkg_a_min_var = tk.StringVar(value="0.0")
        self.dynamic_bkg_a_max_var = tk.StringVar(value="1.0")
        self.dynamic_bkg_a_fit_var = tk.StringVar(value="")
        self.dynamic_bkg_b_start_var = tk.StringVar(value="1.2366e-7")
        self.dynamic_bkg_b_min_var = tk.StringVar(value="0.0")
        self.dynamic_bkg_b_max_var = tk.StringVar(value="0.1")
        self.dynamic_bkg_b_fit_var = tk.StringVar(value="")
        self.dynamic_backend_var = tk.StringVar(value="auto")
        self.dynamic_resolution_fit_enabled_var = tk.BooleanVar(value=True)
        self.dynamic_intensity_fit_enabled_var = tk.BooleanVar(value=True)
        self.dynamic_bkg_a_fit_enabled_var = tk.BooleanVar(value=True)
        self.dynamic_bkg_b_fit_enabled_var = tk.BooleanVar(value=True)
        self.substrate_filename_var = tk.StringVar(value=SAMPLES["Fe 10 nm"].substrate_filename)
        self.substrate_direction_var = tk.StringVar(value=str(SAMPLES["Fe 10 nm"].substrate_direction))
        self.substrate_n_var = tk.StringVar(value=f"{SAMPLES['Fe 10 nm'].substrate_n:.6g}")
        self.substrate_interface_var = tk.StringVar(value=f"{SAMPLES['Fe 10 nm'].substrate_dinterface:.6g}")
        self.substrate_scale_var = tk.StringVar(value=f"{SAMPLES['Fe 10 nm'].substrate_scale:.6g}")
        self.substrate_scale_min_var = tk.StringVar(value="0.995")
        self.substrate_scale_max_var = tk.StringVar(value="1.005")
        self.substrate_scale_fit_var = tk.StringVar(value="")
        self.substrate_scale_fit_enabled_var = tk.BooleanVar(value=False)
        self.substrate_area_var = tk.StringVar(value=f"{SAMPLES['Fe 10 nm'].substrate_area_scale:.6g}")
        self.film_rough_start_var = tk.StringVar(value="0.0")
        self.film_rough_min_var = tk.StringVar(value="0.0")
        self.film_rough_max_var = tk.StringVar(value="5.0")
        self.substrate_rough_start_var = tk.StringVar(value="0.0")
        self.substrate_rough_min_var = tk.StringVar(value="0.0")
        self.substrate_rough_max_var = tk.StringVar(value="5.0")
        self.film_rough_fit_var = tk.StringVar(value="")
        self.substrate_rough_fit_var = tk.StringVar(value="")
        self.film_rough_fit_enabled_var = tk.BooleanVar(value=True)
        self.substrate_rough_fit_enabled_var = tk.BooleanVar(value=True)
        self.strain_values_by_model = {
            "Kinematic": [
                ("0.0", "-0.2", "0.2"),
                ("0.0", "0.0", "20.0"),
                ("0.0", "-0.2", "0.2"),
                ("0.0", "0.0", "20.0"),
            ],
            "Dynamic": [
                ("0.0", "-0.4", "0.4"),
                ("0.0", "0.0", "20.0"),
                ("0.0", "-0.4", "0.4"),
                ("0.0", "0.0", "20.0"),
            ],
        }
        self.active_strain_model = self.model_var.get()
        strain_values = self.strain_values_by_model[self.active_strain_model]
        self.strain_start_vars = [tk.StringVar(value=value[0]) for value in strain_values]
        self.strain_min_vars = [tk.StringVar(value=value[1]) for value in strain_values]
        self.strain_max_vars = [tk.StringVar(value=value[2]) for value in strain_values]
        self.strain_fit_vars = [tk.StringVar(value="") for _ in STRAIN_KEYS]
        self.strain_fit_enabled_vars = [tk.BooleanVar(value=True) for _ in STRAIN_KEYS]
        self.bottom_extent_label_var = tk.StringVar(value="Bottom extent (atomic positions)")
        self.top_extent_label_var = tk.StringVar(value="Top extent (atomic positions)")
        self.status_var = tk.StringVar(value="Ready")
        self.kinematic_widgets: list[tk.Widget] = []
        self.dynamic_widgets: list[tk.Widget] = []
        self.strain_widgets: list[tk.Widget] = []
        self.range_indicators: list[RangeIndicator] = []

        self._set_sample_layer_defaults(SAMPLES[self.sample_var.get()])
        self._set_sample_kinematic_defaults(SAMPLES[self.sample_var.get()])
        self._configure_styles()
        self._build_layout()
        self.status_var.trace_add("write", self._sync_status_style)
        self._sync_status_style()
        self.min_var.trace_add("write", self._schedule_data_preview)
        self.max_var.trace_add("write", self._schedule_data_preview)
        self.wavelength_var.trace_add("write", self._on_wavelength_changed)
        self.axis_var.trace_add("write", self._on_axis_mode_changed)
        self.data_path_var.trace_add("write", self._on_data_path_changed)
        self.model_var.trace_add("write", self._on_model_changed)
        self.strain_var.trace_add("write", self._on_strain_changed)
        self.root.after(0, self._draw_experimental_preview)
        self.root.after(150, self._process_queue)

    def _configure_styles(self) -> None:
        style = ttk.Style(self.root)
        if "clam" in style.theme_names():
            style.theme_use("clam")
        self.root.configure(background=UI_COLORS["window"])
        style.configure(".", background=UI_COLORS["window"])
        style.configure("TFrame", background=UI_COLORS["window"])
        style.configure("Panel.TFrame", background=UI_COLORS["panel"])
        style.configure("TLabelframe", background=UI_COLORS["panel"])
        style.configure("TLabelframe.Label", background=UI_COLORS["window"])
        style.configure("TLabel", background=UI_COLORS["panel"])
        style.configure("Status.TLabel", background=UI_COLORS["panel"], foreground="#374151")
        style.configure("StatusActive.TLabel", background=UI_COLORS["panel"], foreground=UI_COLORS["simulate"])
        style.configure("StatusSuccess.TLabel", background=UI_COLORS["panel"], foreground=UI_COLORS["run"])
        style.configure("StatusWarning.TLabel", background=UI_COLORS["panel"], foreground=UI_COLORS["fit"])
        style.configure("StatusError.TLabel", background=UI_COLORS["panel"], foreground=UI_COLORS["stop"])
        style.configure(
            "Watermark.TLabel",
            background=UI_COLORS["window"],
            foreground=UI_COLORS["watermark"],
            font=("TkDefaultFont", 9),
        )
        style.configure("TCheckbutton", background=UI_COLORS["panel"])
        style.configure("TEntry", fieldbackground=UI_COLORS["entry"])
        style.map(
            "TEntry",
            fieldbackground=[("disabled", UI_COLORS["entry_disabled"])],
            foreground=[("disabled", "#6b7280")],
        )
        for name, color in (
            ("Simulate", UI_COLORS["simulate"]),
            ("Run", UI_COLORS["run"]),
            ("Stop", UI_COLORS["stop"]),
        ):
            style.configure(f"{name}.TButton", background=color, foreground="#ffffff")
            style.map(
                f"{name}.TButton",
                background=[("disabled", "#c8cdd2"), ("active", color)],
                foreground=[("disabled", "#6b7280"), ("active", "#ffffff")],
            )

    def _add_accent_strip(self, parent: tk.Widget, color: str) -> None:
        strip = tk.Frame(parent, height=4, bg=color, highlightthickness=0)
        strip.place(relx=0.0, rely=0.0, relwidth=1.0)

    def _sync_status_style(self, *_args: object) -> None:
        if not hasattr(self, "status_label"):
            return
        text = self.status_var.get().lower()
        if any(marker in text for marker in ("failed", "unavailable", "invalid", "stopped", "cancelled")):
            style = "StatusError.TLabel"
        elif "warning" in text:
            style = "StatusWarning.TLabel"
        elif any(marker in text for marker in ("running", "stopping", "cost=")):
            style = "StatusActive.TLabel"
        elif "complete" in text:
            style = "StatusSuccess.TLabel"
        else:
            style = "Status.TLabel"
        self.status_label.configure(style=style)

    def _build_layout(self) -> None:
        self.root.geometry("1400x850")
        self.root.minsize(1180, 720)

        main_pane = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_pane.pack(fill=tk.BOTH, expand=True)

        sidebar = ttk.Frame(main_pane, padding=10)
        main_pane.add(sidebar, weight=0)

        run_frame = ttk.LabelFrame(sidebar, text="Simulation and fit setup")
        run_frame.pack(fill=tk.X)
        self._add_accent_strip(run_frame, UI_COLORS["setup"])

        ttk.Label(run_frame, text="Scattering model").grid(row=0, column=0, sticky="w")
        self.model_combo = ttk.Combobox(
            run_frame,
            textvariable=self.model_var,
            values=("Dynamic", "Kinematic"),
            state="readonly",
            width=14,
        )
        self.model_combo.grid(row=0, column=1, sticky="ew", padx=(4, 8), pady=2)

        ttk.Label(run_frame, text="Experimental data file").grid(row=1, column=0, sticky="w")
        ttk.Entry(run_frame, textvariable=self.data_path_var, width=42).grid(
            row=1, column=1, columnspan=2, sticky="ew", padx=(4, 8), pady=2
        )
        ttk.Button(run_frame, text="Browse...", command=self._browse_data_file).grid(
            row=1, column=3, sticky="ew", pady=2
        )

        ttk.Label(run_frame, text="X-ray wavelength (\u00c5)").grid(row=2, column=0, sticky="w")
        ttk.Entry(run_frame, textvariable=self.wavelength_var, width=10).grid(
            row=2, column=1, sticky="ew", padx=(4, 8), pady=2
        )
        ttk.Label(run_frame, text="Horizontal axis").grid(row=2, column=2, sticky="w")
        ttk.Combobox(
            run_frame,
            textvariable=self.axis_var,
            values=("2\u03b8", "q"),
            state="readonly",
            width=10,
        ).grid(row=2, column=3, sticky="ew", padx=(4, 8), pady=2)

        run_controls = [
            (self.min_label_var, self.min_var, 3, 0),
            (self.max_label_var, self.max_var, 3, 2),
            ("Seed", self.seed_var, 4, 0),
            ("Progress update interval", self.interval_var, 4, 2),
            ("DE max iterations", self.maxiter_var, 5, 0),
            ("DE population size", self.popsize_var, 5, 2),
            ("Local max evaluations", self.local_var, 6, 0),
            ("Polish iterations", self.polish_var, 6, 2),
        ]
        for label, var, row, column in run_controls:
            if isinstance(label, tk.StringVar):
                ttk.Label(run_frame, textvariable=label).grid(row=row, column=column, sticky="w")
            else:
                ttk.Label(run_frame, text=label).grid(row=row, column=column, sticky="w")
            ttk.Entry(run_frame, textvariable=var, width=10).grid(
                row=row, column=column + 1, sticky="ew", padx=(4, 8), pady=2
            )

        ttk.Checkbutton(run_frame, text="Include strain", variable=self.strain_var).grid(
            row=7, column=0, columnspan=2, sticky="w", pady=(4, 2)
        )
        ttk.Checkbutton(run_frame, text="Include roughness", variable=self.roughness_var).grid(
            row=7, column=2, columnspan=2, sticky="w", pady=(4, 2)
        )

        self.simulate_button = ttk.Button(
            run_frame,
            text="Simulate",
            command=self.simulate_pattern,
            style="Simulate.TButton",
        )
        self.simulate_button.grid(row=8, column=0, sticky="ew", padx=(0, 4), pady=(6, 2))
        self.run_button = ttk.Button(
            run_frame,
            text="Run fit",
            command=self.run_fit,
            style="Run.TButton",
        )
        self.run_button.grid(row=8, column=1, columnspan=2, sticky="ew", padx=(0, 4), pady=(6, 2))
        self.stop_button = ttk.Button(
            run_frame,
            text="Stop fit",
            command=self.stop_fit,
            state=tk.DISABLED,
            style="Stop.TButton",
        )
        self.stop_button.grid(row=8, column=3, sticky="ew", padx=(4, 0), pady=(6, 2))

        self.status_label = ttk.Label(
            run_frame,
            textvariable=self.status_var,
            wraplength=390,
            style="Status.TLabel",
        )
        self.status_label.grid(row=9, column=0, columnspan=4, sticky="ew", pady=(4, 2))
        for column in (1, 3):
            run_frame.columnconfigure(column, weight=1)

        parameter_tabs = ttk.Notebook(sidebar)
        parameter_tabs.pack(fill=tk.BOTH, expand=True, pady=(8, 8))

        self.kinematic_frame = ttk.Frame(parameter_tabs, padding=8, style="Panel.TFrame")
        parameter_tabs.add(self.kinematic_frame, text="Kinematic parameters")
        self._add_accent_strip(self.kinematic_frame, UI_COLORS["kinematic"])
        ttk.Label(self.kinematic_frame, text="Fit").grid(row=0, column=1, sticky="ew")
        ttk.Label(self.kinematic_frame, text="Value").grid(row=0, column=2, sticky="ew")
        ttk.Label(self.kinematic_frame, text="Min").grid(row=0, column=3, sticky="ew")
        ttk.Label(self.kinematic_frame, text="Max").grid(row=0, column=4, sticky="ew")
        ttk.Label(self.kinematic_frame, text="Range").grid(row=0, column=5, sticky="ew")
        kinematic_controls = [
            ("Plane spacing d (\u00c5)", self.kin_d_fit_enabled_var, self.kin_d_start_var, self.kin_d_min_var, self.kin_d_max_var, self.kin_d_fit_var),
            ("Number of planes", self.kin_planes_fit_enabled_var, self.kin_planes_start_var, self.kin_planes_min_var, self.kin_planes_max_var, self.kin_planes_fit_var),
            (
                "Resolution (deg)",
                self.kin_resolution_fit_enabled_var,
                self.kin_resolution_start_var,
                self.kin_resolution_min_var,
                self.kin_resolution_max_var,
                self.kin_resolution_fit_var,
            ),
            ("Intensity scale", self.kin_scale_fit_enabled_var, self.kin_scale_start_var, self.kin_scale_min_var, self.kin_scale_max_var, self.kin_scale_fit_var),
            ("Background offset", self.kin_bkg_a_fit_enabled_var, self.kin_bkg_a_start_var, self.kin_bkg_a_min_var, self.kin_bkg_a_max_var, self.kin_bkg_a_fit_var),
            ("Background slope", self.kin_bkg_b_fit_enabled_var, self.kin_bkg_b_start_var, self.kin_bkg_b_min_var, self.kin_bkg_b_max_var, self.kin_bkg_b_fit_var),
        ]
        for row, (label, fit_enabled_var, start_var, min_var, max_var, fit_var) in enumerate(kinematic_controls, start=1):
            ttk.Label(self.kinematic_frame, text=label).grid(row=row, column=0, sticky="w")
            ttk.Checkbutton(self.kinematic_frame, variable=fit_enabled_var).grid(row=row, column=1, sticky="ew")
            ttk.Entry(self.kinematic_frame, textvariable=start_var, width=7).grid(row=row, column=2, sticky="ew")
            ttk.Entry(self.kinematic_frame, textvariable=min_var, width=7).grid(row=row, column=3, sticky="ew")
            ttk.Entry(self.kinematic_frame, textvariable=max_var, width=7).grid(row=row, column=4, sticky="ew")
            self._add_range_indicator(self.kinematic_frame, row, label, start_var, min_var, max_var, fit_var, fit_enabled_var)
        ttk.Label(self.kinematic_frame, text="Debye-Waller coefficient").grid(row=7, column=0, sticky="w")
        ttk.Entry(self.kinematic_frame, textvariable=self.kin_debye_var, width=7).grid(
            row=7, column=2, columnspan=3, sticky="ew"
        )
        self.kinematic_frame.columnconfigure(2, weight=1)
        self.kinematic_frame.columnconfigure(3, weight=1)
        self.kinematic_frame.columnconfigure(4, weight=1)
        self.kinematic_frame.columnconfigure(5, weight=1)

        self.dynamic_film_frame = ttk.Frame(parameter_tabs, padding=8, style="Panel.TFrame")
        parameter_tabs.add(self.dynamic_film_frame, text="Dynamic film parameters")
        self._add_accent_strip(self.dynamic_film_frame, UI_COLORS["film"])
        ttk.Label(self.dynamic_film_frame, text="Structure file").grid(row=0, column=0, sticky="w")
        ttk.Combobox(
            self.dynamic_film_frame,
            textvariable=self.film_filename_var,
            values=POSCAR_FILES,
            state="readonly",
            width=18,
        ).grid(row=0, column=2, columnspan=3, sticky="ew")
        ttk.Label(self.dynamic_film_frame, text="Layer direction").grid(row=1, column=0, sticky="w")
        ttk.Combobox(
            self.dynamic_film_frame,
            textvariable=self.film_direction_var,
            values=("1", "2", "3"),
            state="readonly",
            width=6,
        ).grid(row=1, column=2, sticky="w")
        ttk.Label(self.dynamic_film_frame, text="Fit").grid(row=2, column=1, sticky="ew")
        ttk.Label(self.dynamic_film_frame, text="Value").grid(row=2, column=2, sticky="ew")
        ttk.Label(self.dynamic_film_frame, text="Min").grid(row=2, column=3, sticky="ew")
        ttk.Label(self.dynamic_film_frame, text="Max").grid(row=2, column=4, sticky="ew")
        ttk.Label(self.dynamic_film_frame, text="Range").grid(row=2, column=5, sticky="ew")
        film_controls = [
            ("Number of layers", self.film_n_fit_enabled_var, self.film_n_start_var, self.film_n_min_var, self.film_n_max_var, self.film_n_fit_var),
            ("Lattice scale", self.film_scale_fit_enabled_var, self.film_scale_start_var, self.film_scale_min_var, self.film_scale_max_var, self.film_scale_fit_var),
            ("Area scale", self.film_area_fit_enabled_var, self.film_area_start_var, self.film_area_min_var, self.film_area_max_var, self.film_area_fit_var),
            (
                "Interface spacing (\u00c5)",
                self.film_interface_fit_enabled_var,
                self.film_interface_start_var,
                self.film_interface_min_var,
                self.film_interface_max_var,
                self.film_interface_fit_var,
            ),
        ]
        for row, (label, fit_enabled_var, start_var, min_var, max_var, fit_var) in enumerate(film_controls, start=3):
            ttk.Label(self.dynamic_film_frame, text=label).grid(row=row, column=0, sticky="w")
            ttk.Checkbutton(self.dynamic_film_frame, variable=fit_enabled_var).grid(row=row, column=1, sticky="ew")
            ttk.Entry(self.dynamic_film_frame, textvariable=start_var, width=7).grid(row=row, column=2, sticky="ew")
            ttk.Entry(self.dynamic_film_frame, textvariable=min_var, width=7).grid(row=row, column=3, sticky="ew")
            ttk.Entry(self.dynamic_film_frame, textvariable=max_var, width=7).grid(row=row, column=4, sticky="ew")
            self._add_range_indicator(self.dynamic_film_frame, row, label, start_var, min_var, max_var, fit_var, fit_enabled_var)
        self.dynamic_film_frame.columnconfigure(2, weight=1)
        self.dynamic_film_frame.columnconfigure(3, weight=1)
        self.dynamic_film_frame.columnconfigure(4, weight=1)
        self.dynamic_film_frame.columnconfigure(5, weight=1)

        self.dynamic_fit_frame = ttk.Frame(parameter_tabs, padding=8, style="Panel.TFrame")
        parameter_tabs.add(self.dynamic_fit_frame, text="Dynamic fit parameters")
        self._add_accent_strip(self.dynamic_fit_frame, UI_COLORS["fit"])
        ttk.Label(self.dynamic_fit_frame, text="Dynamic backend").grid(row=0, column=0, sticky="w")
        ttk.Combobox(
            self.dynamic_fit_frame,
            textvariable=self.dynamic_backend_var,
            values=("auto", "fused", "legacy"),
            state="readonly",
            width=10,
        ).grid(row=0, column=2, sticky="w")
        ttk.Label(self.dynamic_fit_frame, text="Fit").grid(row=1, column=1, sticky="ew")
        ttk.Label(self.dynamic_fit_frame, text="Value").grid(row=1, column=2, sticky="ew")
        ttk.Label(self.dynamic_fit_frame, text="Min").grid(row=1, column=3, sticky="ew")
        ttk.Label(self.dynamic_fit_frame, text="Max").grid(row=1, column=4, sticky="ew")
        ttk.Label(self.dynamic_fit_frame, text="Range").grid(row=1, column=5, sticky="ew")
        dynamic_fit_controls = [
            (
                "Resolution (deg)",
                self.dynamic_resolution_fit_enabled_var,
                self.dynamic_resolution_start_var,
                self.dynamic_resolution_min_var,
                self.dynamic_resolution_max_var,
                self.dynamic_resolution_fit_var,
            ),
            (
                "Intensity scale",
                self.dynamic_intensity_fit_enabled_var,
                self.dynamic_intensity_start_var,
                self.dynamic_intensity_min_var,
                self.dynamic_intensity_max_var,
                self.dynamic_intensity_fit_var,
            ),
            (
                "Background offset",
                self.dynamic_bkg_a_fit_enabled_var,
                self.dynamic_bkg_a_start_var,
                self.dynamic_bkg_a_min_var,
                self.dynamic_bkg_a_max_var,
                self.dynamic_bkg_a_fit_var,
            ),
            (
                "Background slope",
                self.dynamic_bkg_b_fit_enabled_var,
                self.dynamic_bkg_b_start_var,
                self.dynamic_bkg_b_min_var,
                self.dynamic_bkg_b_max_var,
                self.dynamic_bkg_b_fit_var,
            ),
        ]
        for row, (label, fit_enabled_var, start_var, min_var, max_var, fit_var) in enumerate(dynamic_fit_controls, start=2):
            ttk.Label(self.dynamic_fit_frame, text=label).grid(row=row, column=0, sticky="w")
            ttk.Checkbutton(self.dynamic_fit_frame, variable=fit_enabled_var).grid(row=row, column=1, sticky="ew")
            ttk.Entry(self.dynamic_fit_frame, textvariable=start_var, width=7).grid(row=row, column=2, sticky="ew")
            ttk.Entry(self.dynamic_fit_frame, textvariable=min_var, width=7).grid(row=row, column=3, sticky="ew")
            ttk.Entry(self.dynamic_fit_frame, textvariable=max_var, width=7).grid(row=row, column=4, sticky="ew")
            self._add_range_indicator(self.dynamic_fit_frame, row, label, start_var, min_var, max_var, fit_var, fit_enabled_var)
        self.dynamic_fit_frame.columnconfigure(2, weight=1)
        self.dynamic_fit_frame.columnconfigure(3, weight=1)
        self.dynamic_fit_frame.columnconfigure(4, weight=1)
        self.dynamic_fit_frame.columnconfigure(5, weight=1)

        self.dynamic_substrate_frame = ttk.Frame(parameter_tabs, padding=8, style="Panel.TFrame")
        parameter_tabs.add(self.dynamic_substrate_frame, text="Dynamic substrate setup")
        self._add_accent_strip(self.dynamic_substrate_frame, UI_COLORS["substrate"])
        ttk.Label(self.dynamic_substrate_frame, text="Structure file").grid(row=0, column=0, sticky="w")
        ttk.Combobox(
            self.dynamic_substrate_frame,
            textvariable=self.substrate_filename_var,
            values=POSCAR_FILES,
            state="readonly",
            width=18,
        ).grid(row=0, column=1, sticky="ew")
        substrate_controls = [
            ("Layer direction", self.substrate_direction_var, ("1", "2", "3")),
            ("Number of layers", self.substrate_n_var, None),
            ("Interface spacing (\u00c5)", self.substrate_interface_var, None),
            ("Area scale", self.substrate_area_var, None),
        ]
        for row, (label, var, values) in enumerate(substrate_controls, start=1):
            ttk.Label(self.dynamic_substrate_frame, text=label).grid(row=row, column=0, sticky="w")
            if values is None:
                ttk.Entry(self.dynamic_substrate_frame, textvariable=var, width=12).grid(row=row, column=1, sticky="ew")
            else:
                ttk.Combobox(
                    self.dynamic_substrate_frame,
                    textvariable=var,
                    values=values,
                    state="readonly",
                    width=8,
                ).grid(row=row, column=1, sticky="w")
        ttk.Label(self.dynamic_substrate_frame, text="Fit").grid(row=5, column=1, sticky="ew")
        ttk.Label(self.dynamic_substrate_frame, text="Value").grid(row=5, column=2, sticky="ew")
        ttk.Label(self.dynamic_substrate_frame, text="Min").grid(row=5, column=3, sticky="ew")
        ttk.Label(self.dynamic_substrate_frame, text="Max").grid(row=5, column=4, sticky="ew")
        ttk.Label(self.dynamic_substrate_frame, text="Range").grid(row=5, column=5, sticky="ew")
        ttk.Label(self.dynamic_substrate_frame, text="Lattice scale").grid(row=6, column=0, sticky="w")
        ttk.Checkbutton(
            self.dynamic_substrate_frame,
            variable=self.substrate_scale_fit_enabled_var,
        ).grid(row=6, column=1, sticky="ew")
        ttk.Entry(self.dynamic_substrate_frame, textvariable=self.substrate_scale_var, width=7).grid(
            row=6, column=2, sticky="ew"
        )
        ttk.Entry(self.dynamic_substrate_frame, textvariable=self.substrate_scale_min_var, width=7).grid(
            row=6, column=3, sticky="ew"
        )
        ttk.Entry(self.dynamic_substrate_frame, textvariable=self.substrate_scale_max_var, width=7).grid(
            row=6, column=4, sticky="ew"
        )
        self._add_range_indicator(
            self.dynamic_substrate_frame,
            6,
            "substrate lattice scale",
            self.substrate_scale_var,
            self.substrate_scale_min_var,
            self.substrate_scale_max_var,
            self.substrate_scale_fit_var,
            self.substrate_scale_fit_enabled_var,
        )
        for column in range(1, 6):
            self.dynamic_substrate_frame.columnconfigure(column, weight=1)

        self.strain_frame = ttk.Frame(parameter_tabs, padding=8, style="Panel.TFrame")
        parameter_tabs.add(self.strain_frame, text="Strain parameters")
        self._add_accent_strip(self.strain_frame, UI_COLORS["strain"])
        ttk.Label(self.strain_frame, text="").grid(row=0, column=0, sticky="w")
        ttk.Label(self.strain_frame, text="Fit").grid(row=0, column=1, sticky="ew")
        ttk.Label(self.strain_frame, text="Value").grid(row=0, column=2, sticky="ew")
        ttk.Label(self.strain_frame, text="Min").grid(row=0, column=3, sticky="ew")
        ttk.Label(self.strain_frame, text="Max").grid(row=0, column=4, sticky="ew")
        ttk.Label(self.strain_frame, text="Range").grid(row=0, column=5, sticky="ew")
        strain_controls = [
            "Bottom strain amplitude",
            self.bottom_extent_label_var,
            "Top strain amplitude",
            self.top_extent_label_var,
        ]
        for row, (label, fit_enabled_var, start_var, min_var, max_var, fit_var) in enumerate(
            zip(
                strain_controls,
                self.strain_fit_enabled_vars,
                self.strain_start_vars,
                self.strain_min_vars,
                self.strain_max_vars,
                self.strain_fit_vars,
            ),
            start=1,
        ):
            name = label.get() if isinstance(label, tk.StringVar) else label
            if isinstance(label, tk.StringVar):
                ttk.Label(self.strain_frame, textvariable=label).grid(row=row, column=0, sticky="w")
            else:
                ttk.Label(self.strain_frame, text=label).grid(row=row, column=0, sticky="w")
            ttk.Checkbutton(self.strain_frame, variable=fit_enabled_var).grid(row=row, column=1, sticky="ew")
            ttk.Entry(self.strain_frame, textvariable=start_var, width=7).grid(row=row, column=2, sticky="ew")
            ttk.Entry(self.strain_frame, textvariable=min_var, width=7).grid(row=row, column=3, sticky="ew")
            ttk.Entry(self.strain_frame, textvariable=max_var, width=7).grid(row=row, column=4, sticky="ew")
            self._add_range_indicator(
                self.strain_frame,
                row,
                name,
                start_var,
                min_var,
                max_var,
                fit_var,
                fit_enabled_var,
            )
        for column in range(2, 6):
            self.strain_frame.columnconfigure(column, weight=1)

        rough_frame = ttk.Frame(parameter_tabs, padding=8, style="Panel.TFrame")
        parameter_tabs.add(rough_frame, text="Roughness parameters")
        self._add_accent_strip(rough_frame, UI_COLORS["roughness"])
        ttk.Label(rough_frame, text="").grid(row=0, column=0, sticky="w")
        ttk.Label(rough_frame, text="Fit").grid(row=0, column=1, sticky="ew")
        ttk.Label(rough_frame, text="Value").grid(row=0, column=2, sticky="ew")
        ttk.Label(rough_frame, text="Min").grid(row=0, column=3, sticky="ew")
        ttk.Label(rough_frame, text="Max").grid(row=0, column=4, sticky="ew")
        ttk.Label(rough_frame, text="Range").grid(row=0, column=5, sticky="ew")
        rough_controls = [
            (
                "Film thickness \u03c3 (planes/layers)",
                self.film_rough_fit_enabled_var,
                self.film_rough_start_var,
                self.film_rough_min_var,
                self.film_rough_max_var,
                self.film_rough_fit_var,
            ),
            (
                "Substrate/interface roughness \u03c3 (\u00c5)",
                self.substrate_rough_fit_enabled_var,
                self.substrate_rough_start_var,
                self.substrate_rough_min_var,
                self.substrate_rough_max_var,
                self.substrate_rough_fit_var,
            ),
        ]
        for row, (label, fit_enabled_var, start_var, min_var, max_var, fit_var) in enumerate(rough_controls, start=1):
            ttk.Label(rough_frame, text=label).grid(row=row, column=0, sticky="w")
            ttk.Checkbutton(rough_frame, variable=fit_enabled_var).grid(row=row, column=1, sticky="ew")
            ttk.Entry(rough_frame, textvariable=start_var, width=7).grid(row=row, column=2, sticky="ew")
            ttk.Entry(rough_frame, textvariable=min_var, width=7).grid(row=row, column=3, sticky="ew")
            ttk.Entry(rough_frame, textvariable=max_var, width=7).grid(row=row, column=4, sticky="ew")
            self._add_range_indicator(rough_frame, row, label, start_var, min_var, max_var, fit_var, fit_enabled_var)
        rough_frame.columnconfigure(2, weight=1)
        rough_frame.columnconfigure(3, weight=1)
        rough_frame.columnconfigure(4, weight=1)
        rough_frame.columnconfigure(5, weight=1)

        summary_frame = ttk.LabelFrame(sidebar, text="Fit results")
        summary_frame.pack(fill=tk.X)
        self._add_accent_strip(summary_frame, UI_COLORS["setup"])
        self.summary_text = tk.Text(
            summary_frame,
            width=48,
            height=9,
            background=UI_COLORS["panel"],
            foreground="#111827",
            insertbackground="#111827",
            relief=tk.FLAT,
            borderwidth=1,
            highlightthickness=1,
            highlightbackground="#d1d5db",
            highlightcolor="#9ca3af",
        )
        summary_scrollbar = ttk.Scrollbar(
            summary_frame,
            orient=tk.VERTICAL,
            command=self.summary_text.yview,
        )
        self.summary_text.configure(yscrollcommand=summary_scrollbar.set)
        self.summary_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        summary_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        plot_frame = ttk.Frame(main_pane, padding=(0, 10, 10, 10))
        main_pane.add(plot_frame, weight=1)
        self.figure = Figure(figsize=(9.5, 8), constrained_layout=True)
        self.loss_axis = self.figure.add_subplot(311)
        self.fit_axis = self.figure.add_subplot(312)
        self.density_axis = self.figure.add_subplot(313)
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        watermark = ttk.Label(
            plot_frame,
            text=GENL_CITATION,
            style="Watermark.TLabel",
            cursor="hand2",
            anchor="center",
            justify=tk.CENTER,
            wraplength=900,
        )
        watermark.pack(fill=tk.X, pady=(4, 0))
        watermark.bind("<Button-1>", lambda _event: webbrowser.open_new(GENL_DOI_URL))
        self._draw_empty_plot()
        self.kinematic_widgets = self._collect_children(self.kinematic_frame)
        self.dynamic_widgets = (
            self._collect_children(self.dynamic_film_frame)
            + self._collect_children(self.dynamic_fit_frame)
            + self._collect_children(self.dynamic_substrate_frame)
        )
        self.strain_widgets = self._collect_children(self.strain_frame)
        self._sync_model_controls()
        self._sync_optional_controls()
        self._redraw_range_indicators()

    def _collect_children(self, parent: tk.Widget) -> list[tk.Widget]:
        widgets = list(parent.winfo_children())
        for child in list(widgets):
            widgets.extend(self._collect_children(child))
        return widgets

    def _add_range_indicator(
        self,
        parent: tk.Widget,
        row: int,
        name: str,
        start_var: tk.StringVar,
        min_var: tk.StringVar,
        max_var: tk.StringVar,
        fit_var: tk.StringVar,
        fit_enabled_var: tk.BooleanVar | None = None,
    ) -> None:
        canvas = tk.Canvas(parent, width=96, height=20, highlightthickness=0, bg="#f5f5f5")
        canvas.grid(row=row, column=5, sticky="ew", padx=(4, 0))
        indicator = RangeIndicator(name, canvas, start_var, min_var, max_var, fit_var, fit_enabled_var)
        self.range_indicators.append(indicator)

        def redraw(*_args: object) -> None:
            self._draw_range_indicator(indicator)

        for variable in (start_var, min_var, max_var, fit_var, fit_enabled_var):
            if variable is None:
                continue
            variable.trace_add("write", redraw)
        canvas.bind("<Configure>", lambda _event: self._draw_range_indicator(indicator))
        self._draw_range_indicator(indicator)

    def _parse_float_var(self, variable: tk.StringVar) -> float | None:
        try:
            return float(variable.get())
        except ValueError:
            return None

    def _draw_range_indicator(self, indicator: RangeIndicator) -> None:
        canvas = indicator.canvas
        canvas.delete("all")
        width = max(canvas.winfo_width(), 96)
        height = max(canvas.winfo_height(), 20)
        pad = 8
        y = height // 2
        lower = self._parse_float_var(indicator.min_var)
        upper = self._parse_float_var(indicator.max_var)
        value = self._parse_float_var(indicator.start_var)
        fitted_value = self._parse_float_var(indicator.fit_var)
        if indicator.fit_enabled_var is not None and not indicator.fit_enabled_var.get():
            canvas.create_text(width // 2, y, text="fixed", fill="#777777", font=("TkDefaultFont", 8))
            return
        if lower is None or upper is None or upper <= lower:
            canvas.create_text(width // 2, y, text="limits", fill="#777777", font=("TkDefaultFont", 8))
            return

        canvas.create_line(pad, y, width - pad, y, fill="#b7b7b7", width=2)
        canvas.create_line(pad, y - 4, pad, y + 4, fill="#777777")
        canvas.create_line(width - pad, y - 4, width - pad, y + 4, fill="#777777")

        def x_for(value: float) -> float:
            ratio = (value - lower) / (upper - lower)
            ratio = min(max(ratio, 0.0), 1.0)
            return pad + ratio * (width - 2 * pad)

        if value is not None:
            x = x_for(value)
            ratio = (value - lower) / (upper - lower)
            if ratio <= 0.02 or ratio >= 0.98:
                color = "#d62728"
            elif ratio <= 0.05 or ratio >= 0.95:
                color = "#ffbf00"
            else:
                color = "#2ca02c" if fitted_value is not None else "#1f77b4"
            canvas.create_oval(x - 4, y - 4, x + 4, y + 4, fill=color, outline=color)
        elif indicator.fit_var.get() == "off":
            canvas.create_text(width // 2, y, text="off", fill="#777777", font=("TkDefaultFont", 8))

    def _redraw_range_indicators(self) -> None:
        for indicator in self.range_indicators:
            self._draw_range_indicator(indicator)

    def _set_widgets_enabled(self, widgets: list[tk.Widget], enabled: bool) -> None:
        for widget in widgets:
            try:
                if isinstance(widget, ttk.Combobox):
                    widget.configure(state="readonly" if enabled else tk.DISABLED)
                else:
                    widget.configure(state=tk.NORMAL if enabled else tk.DISABLED)
            except tk.TclError:
                pass

    def _sync_model_controls(self) -> None:
        is_dynamic = self.model_var.get() == "Dynamic"
        self._set_widgets_enabled(self.dynamic_widgets, is_dynamic)
        self._set_widgets_enabled(self.kinematic_widgets, not is_dynamic)

    def _sync_optional_controls(self) -> None:
        self._set_widgets_enabled(self.strain_widgets, bool(self.strain_var.get()))

    def _store_strain_values(self) -> None:
        self.strain_values_by_model[self.active_strain_model] = list(
            zip(
                (variable.get() for variable in self.strain_start_vars),
                (variable.get() for variable in self.strain_min_vars),
                (variable.get() for variable in self.strain_max_vars),
            )
        )

    def _load_strain_values(self, model_name: str) -> None:
        for values, start_var, min_var, max_var in zip(
            self.strain_values_by_model[model_name],
            self.strain_start_vars,
            self.strain_min_vars,
            self.strain_max_vars,
        ):
            start_var.set(values[0])
            min_var.set(values[1])
            max_var.set(values[2])
        extent_unit = "atomic positions" if model_name == "Dynamic" else "planes"
        self.bottom_extent_label_var.set(f"Bottom extent ({extent_unit})")
        self.top_extent_label_var.set(f"Top extent ({extent_unit})")
        self.active_strain_model = model_name

    def _on_model_changed(self, *_args: object) -> None:
        self._store_strain_values()
        self._load_strain_values(self.model_var.get())
        self._sync_model_controls()
        self._clear_fit_result_values()

    def _on_strain_changed(self, *_args: object) -> None:
        self._sync_optional_controls()
        self._clear_fit_result_values()

    def _browse_data_file(self) -> None:
        initial_path = resolve_data_path(self.data_path_var.get())
        selected = filedialog.askopenfilename(
            title="Choose experimental data file",
            initialdir=str(initial_path.parent if initial_path.parent.exists() else ROOT),
            filetypes=[
                ("Text data", "*.txt *.dat *.csv"),
                ("All files", "*"),
            ],
        )
        if selected:
            self.data_path_var.set(selected)

    def _axis_is_q(self) -> bool:
        return self.axis_var.get() == "q"

    def _axis_label(self) -> str:
        return "q (1/\u00c5)" if self._axis_is_q() else "2\u03b8 (deg)"

    def _axis_limit_labels(self) -> tuple[str, str]:
        return ("q min (1/\u00c5)", "q max (1/\u00c5)") if self._axis_is_q() else ("2\u03b8 min (deg)", "2\u03b8 max (deg)")

    def _parse_wavelength(self) -> float:
        wavelength = float(self.wavelength_var.get())
        if wavelength <= 0:
            raise ValueError("X-ray wavelength must be positive")
        return wavelength

    def _x_values(self, twotheta: np.ndarray, q: np.ndarray | None, wavelength: float) -> np.ndarray:
        if self._axis_is_q():
            return q if q is not None else np.asarray(q_from_twotheta(twotheta, wavelength), dtype=float)
        return twotheta

    def _window_twotheta_limits(self) -> tuple[float, float]:
        lower = float(self.min_var.get())
        upper = float(self.max_var.get())
        if lower > upper:
            raise ValueError(f"{self._axis_label()} min must be <= max")
        if self._axis_is_q():
            if lower < 0 or upper < 0:
                raise ValueError("q limits must be non-negative")
            wavelength = self._parse_wavelength()
            return float(twotheta_from_q(lower, wavelength)), float(twotheta_from_q(upper, wavelength))
        return lower, upper

    def _on_axis_mode_changed(self, *_args: object) -> None:
        old_mode = self.axis_mode
        new_mode = "q" if self._axis_is_q() else "twotheta"
        if old_mode == new_mode:
            return

        self.updating_twotheta_window = True
        try:
            wavelength = self._parse_wavelength()
            lower_text = self.min_var.get().strip()
            upper_text = self.max_var.get().strip()
            if lower_text and upper_text:
                lower = float(lower_text)
                upper = float(upper_text)
                if old_mode == "twotheta" and new_mode == "q":
                    self.min_var.set(self._format_window_limit(float(q_from_twotheta(lower, wavelength))))
                    self.max_var.set(self._format_window_limit(float(q_from_twotheta(upper, wavelength))))
                elif old_mode == "q" and new_mode == "twotheta":
                    self.min_var.set(self._format_window_limit(float(twotheta_from_q(lower, wavelength))))
                    self.max_var.set(self._format_window_limit(float(twotheta_from_q(upper, wavelength))))
        except ValueError:
            pass
        finally:
            self.axis_mode = new_mode
            min_label, max_label = self._axis_limit_labels()
            self.min_label_var.set(min_label)
            self.max_label_var.set(max_label)
            self.updating_twotheta_window = False
        self._schedule_data_preview()

    def _on_wavelength_changed(self, *_args: object) -> None:
        self._clear_fit_result_values()
        self._schedule_data_preview()

    def _set_sample_kinematic_defaults(self, sample: SampleConfig) -> None:
        try:
            wavelength = self._parse_wavelength()
        except ValueError:
            wavelength = CU_K_ALPHA_WAVELENGTH
        d0 = sample.d_spacing
        peak = float(
            2.0
            * np.rad2deg(
                np.arcsin(np.clip(wavelength / (2.0 * d0), -1.0, 1.0))
            )
        )
        n0 = sample.thickness / d0
        self.kin_d_start_var.set(f"{d0:.6g}")
        self.kin_d_min_var.set(f"{wavelength / (2.0 * sind((peak * 1.0035) / 2.0)):.6g}")
        self.kin_d_max_var.set(f"{wavelength / (2.0 * sind((peak * 0.9965) / 2.0)):.6g}")
        self.kin_planes_start_var.set(f"{n0:.6g}")
        self.kin_planes_min_var.set(f"{0.60 * n0:.6g}")
        self.kin_planes_max_var.set(f"{1.20 * n0:.6g}")
        self.kin_resolution_start_var.set("0.005")
        self.kin_resolution_min_var.set("0.00025")
        self.kin_resolution_max_var.set("0.025")
        self.kin_scale_start_var.set("0.0069")
        self.kin_scale_min_var.set("0.001")
        self.kin_scale_max_var.set("0.05")
        self.kin_bkg_a_start_var.set("0.0")
        self.kin_bkg_a_min_var.set("0.0")
        self.kin_bkg_a_max_var.set("1.0")
        self.kin_bkg_b_start_var.set("0.1")
        self.kin_bkg_b_min_var.set("0.0")
        self.kin_bkg_b_max_var.set("3.0")
        self.kin_debye_var.set(f"{sample.debye_waller_coeff:.6g}")

    def _set_sample_layer_defaults(self, sample: SampleConfig) -> None:
        self.film_filename_var.set(sample.film_filename)
        self.film_direction_var.set(str(sample.dynamic_direction))
        n0 = sample.thickness / propagation_period(sample)
        self.film_n_start_var.set(f"{n0:.6g}")
        self.film_n_min_var.set(f"{max(1.0, 0.45 * n0):.6g}")
        self.film_n_max_var.set(f"{1.35 * n0:.6g}")
        self.film_scale_start_var.set("1.0")
        self.film_scale_min_var.set("0.9")
        self.film_scale_max_var.set("1.1")
        self.film_area_start_var.set("1.0")
        self.film_area_min_var.set("0.5")
        self.film_area_max_var.set("1.5")
        self.film_interface_start_var.set("1.35")
        self.film_interface_min_var.set("1.3")
        self.film_interface_max_var.set("1.43")

        self.substrate_filename_var.set(sample.substrate_filename)
        self.substrate_direction_var.set(str(sample.substrate_direction))
        self.substrate_n_var.set(f"{sample.substrate_n:.6g}")
        self.substrate_interface_var.set(f"{sample.substrate_dinterface:.6g}")
        self.substrate_scale_var.set(f"{sample.substrate_scale:.6g}")
        self.substrate_scale_min_var.set(f"{sample.substrate_scale * 0.995:.8g}")
        self.substrate_scale_max_var.set(f"{sample.substrate_scale * 1.005:.8g}")
        self.substrate_scale_fit_var.set("")
        self.substrate_area_var.set(f"{sample.substrate_area_scale:.6g}")

    def _format_window_limit(self, value: float) -> str:
        return f"{value:.6f}".rstrip("0").rstrip(".")

    def _set_twotheta_window(self, twotheta: np.ndarray) -> None:
        self.updating_twotheta_window = True
        try:
            if self._axis_is_q():
                wavelength = self._parse_wavelength()
                values = q_from_twotheta(twotheta, wavelength)
            else:
                values = twotheta
            self.min_var.set(self._format_window_limit(float(np.min(values))))
            self.max_var.set(self._format_window_limit(float(np.max(values))))
        finally:
            self.updating_twotheta_window = False

    def _apply_sample_profile(self, sample: SampleConfig) -> None:
        self.sample_var.set(sample.name)
        self._set_sample_layer_defaults(sample)
        self._set_sample_kinematic_defaults(sample)

    def _on_data_path_changed(self, *_args: object) -> None:
        self.preview_data_path = None
        self._clear_fit_result_values()
        self._schedule_data_preview()

    def _schedule_data_preview(self, *_args: object) -> None:
        if self.updating_twotheta_window:
            return
        if self.running:
            return
        if self.preview_after_id is not None:
            self.root.after_cancel(self.preview_after_id)
        self.preview_after_id = self.root.after(250, self._draw_experimental_preview)

    def _draw_empty_plot(self) -> None:
        self.loss_axis.clear()
        self.loss_axis.set_xlabel("progress callback")
        self.loss_axis.set_ylabel("mean abs log10 error")
        self.loss_axis.grid(True, alpha=0.25)
        self.fit_axis.clear()
        self.fit_axis.set_xlabel(self._axis_label())
        self.fit_axis.set_ylabel("intensity (cps)")
        self.fit_axis.set_yscale("log")
        self.fit_axis.grid(True, alpha=0.25)
        self.density_axis.clear()
        self.density_axis.set_xlabel("z (A)")
        self.density_axis.set_ylabel("Re density")
        self.density_axis.grid(True, alpha=0.25)
        self.canvas.draw_idle()

    def _draw_experimental_preview(self) -> None:
        self.preview_after_id = None
        if self.running:
            return

        self._draw_empty_plot()
        try:
            wavelength = self._parse_wavelength()
            data_path = resolve_data_path(self.data_path_var.get())
            all_twotheta, all_observed = read_experimental_data(data_path)
            if data_path != self.preview_data_path:
                sample = sample_from_data_path(data_path, SAMPLES[self.sample_var.get()])
                self._apply_sample_profile(sample)
                self._set_twotheta_window(all_twotheta)
                self.preview_data_path = data_path
            else:
                sample = SAMPLES[self.sample_var.get()]

            twotheta_min_text = self.min_var.get().strip()
            twotheta_max_text = self.max_var.get().strip()
            if not twotheta_min_text or not twotheta_max_text:
                self._set_twotheta_window(all_twotheta)
                twotheta_min_text = self.min_var.get().strip()
                twotheta_max_text = self.max_var.get().strip()

            twotheta_min, twotheta_max = self._window_twotheta_limits()
            mask = (all_twotheta >= twotheta_min) & (all_twotheta <= twotheta_max)
            if not np.any(mask):
                raise ValueError("No data points found inside the selected 2theta window")
            twotheta = all_twotheta[mask]
            observed = all_observed[mask]
            x_values = self._x_values(twotheta, None, wavelength)
        except ValueError as exc:
            self.status_var.set(f"Data preview unavailable: {exc}")
            return

        self.fit_axis.clear()
        self.fit_axis.plot(
            x_values,
            observed,
            ".",
            color="black",
            markersize=3,
            label="data",
        )
        self.fit_axis.set_yscale("log")
        self.fit_axis.set_xlabel(self._axis_label())
        self.fit_axis.set_ylabel("intensity (cps)")
        self.fit_axis.set_title(f"{data_path.name}: experimental data")
        self.fit_axis.legend(loc="best")
        self.fit_axis.grid(True, alpha=0.25)

        self.density_axis.clear()
        self.density_axis.text(
            0.5,
            0.5,
            "electron density profile appears during dynamic fits",
            transform=self.density_axis.transAxes,
            ha="center",
            va="center",
        )
        self.density_axis.set_xlabel("z (A)")
        self.density_axis.set_ylabel("Re density")
        self.density_axis.grid(True, alpha=0.25)

        self.status_var.set(
            f"Loaded {data_path.name}: showing {len(twotheta)} of {len(all_twotheta)} data points from "
            f"{x_values[0]:.3f} to {x_values[-1]:.3f} {self._axis_label()}"
        )
        self.canvas.draw_idle()

    def _format_fit_result(self, value: float) -> str:
        return f"{value:.6g}"

    def _clear_fit_result_values(self) -> None:
        for variable in (
            self.kin_d_fit_var,
            self.kin_planes_fit_var,
            self.kin_resolution_fit_var,
            self.kin_scale_fit_var,
            self.kin_bkg_a_fit_var,
            self.kin_bkg_b_fit_var,
            self.kin_debye_fit_var,
            self.film_n_fit_var,
            self.film_scale_fit_var,
            self.film_area_fit_var,
            self.film_interface_fit_var,
            self.dynamic_resolution_fit_var,
            self.dynamic_intensity_fit_var,
            self.dynamic_bkg_a_fit_var,
            self.dynamic_bkg_b_fit_var,
            self.substrate_scale_fit_var,
            *self.strain_fit_vars,
            self.film_rough_fit_var,
            self.substrate_rough_fit_var,
        ):
            variable.set("")

    def _set_fit_result_values(self, config: dict[str, object], params: np.ndarray) -> None:
        self._clear_fit_result_values()
        if config["model"] == "Kinematic":
            self.kin_d_start_var.set(self._format_fit_result(float(params[0])))
            self.kin_planes_start_var.set(self._format_fit_result(float(params[1])))
            self.kin_resolution_start_var.set(self._format_fit_result(float(params[2])))
            self.kin_scale_start_var.set(self._format_fit_result(float(params[3])))
            self.kin_bkg_a_start_var.set(self._format_fit_result(float(params[4])))
            self.kin_bkg_b_start_var.set(self._format_fit_result(float(params[5])))
            self.kin_d_fit_var.set(self.kin_d_start_var.get())
            self.kin_planes_fit_var.set(self.kin_planes_start_var.get())
            self.kin_resolution_fit_var.set(self.kin_resolution_start_var.get())
            self.kin_scale_fit_var.set(self.kin_scale_start_var.get())
            self.kin_bkg_a_fit_var.set(self.kin_bkg_a_start_var.get())
            self.kin_bkg_b_fit_var.set(self.kin_bkg_b_start_var.get())
            self.kin_debye_fit_var.set(self.kin_debye_var.get())
            offset = 6
        else:
            self.film_n_start_var.set(self._format_fit_result(float(params[0])))
            self.film_scale_start_var.set(self._format_fit_result(float(params[1])))
            self.film_area_start_var.set(self._format_fit_result(float(params[2])))
            self.film_interface_start_var.set(self._format_fit_result(float(params[3])))
            self.dynamic_resolution_start_var.set(self._format_fit_result(float(params[4])))
            self.dynamic_intensity_start_var.set(self._format_fit_result(float(params[5])))
            self.dynamic_bkg_a_start_var.set(self._format_fit_result(float(params[6])))
            self.dynamic_bkg_b_start_var.set(self._format_fit_result(float(params[7])))
            self.film_n_fit_var.set(self.film_n_start_var.get())
            self.film_scale_fit_var.set(self.film_scale_start_var.get())
            self.film_area_fit_var.set(self.film_area_start_var.get())
            self.film_interface_fit_var.set(self.film_interface_start_var.get())
            self.dynamic_resolution_fit_var.set(self.dynamic_resolution_start_var.get())
            self.dynamic_intensity_fit_var.set(self.dynamic_intensity_start_var.get())
            self.dynamic_bkg_a_fit_var.set(self.dynamic_bkg_a_start_var.get())
            self.dynamic_bkg_b_fit_var.set(self.dynamic_bkg_b_start_var.get())
            self.substrate_scale_var.set(self._format_fit_result(float(params[8])))
            self.substrate_scale_fit_var.set(self.substrate_scale_var.get())
            offset = 9

        if bool(config["include_strain"]):
            for index, (start_var, fit_var) in enumerate(
                zip(self.strain_start_vars, self.strain_fit_vars)
            ):
                start_var.set(self._format_fit_result(float(params[offset + index])))
                fit_var.set(start_var.get())
            offset += 4
        else:
            for fit_var in self.strain_fit_vars:
                fit_var.set("off")
        if bool(config["include_roughness"]):
            self.film_rough_start_var.set(self._format_fit_result(float(params[offset])))
            self.substrate_rough_start_var.set(self._format_fit_result(float(params[offset + 1])))
            self.film_rough_fit_var.set(self.film_rough_start_var.get())
            self.substrate_rough_fit_var.set(self.substrate_rough_start_var.get())
        else:
            self.film_rough_fit_var.set("off")
            self.substrate_rough_fit_var.set("off")
        self._redraw_range_indicators()

    def _fit_boundary_warnings(self, config: dict[str, object], params: np.ndarray) -> list[str]:
        rows: list[tuple[str, float, float, float]] = []
        if config["model"] == "Kinematic":
            settings = config["kinematic_settings"]
            base = (
                ("d spacing", params[0], settings["d_spacing"]),
                ("planes", params[1], settings["planes"]),
                ("resolution", params[2], settings["resolution"]),
                ("scale", params[3], settings["scale"]),
                ("bkg a", params[4], settings["bkg_a"]),
                ("bkg b", params[5], settings["bkg_b"]),
            )
            rows.extend(
                (name, float(value), float(bounds[1]), float(bounds[2]))
                for (name, value, bounds), enabled in zip(base, config["kinematic_fit_flags"])
                if enabled
            )
            offset = 6
        else:
            settings = config["film_settings"]
            base = (
                ("film N", params[0], settings["n"]),
                ("film scale", params[1], settings["scale"]),
                ("film area_scale", params[2], settings["area_scale"]),
                ("film dinterface", params[3], settings["dinterface"]),
            )
            rows.extend(
                (name, float(value), float(bounds[1]), float(bounds[2]))
                for (name, value, bounds), enabled in zip(base, config["film_fit_flags"])
                if enabled
            )
            settings = config["dynamic_fit_settings"]
            base = (
                ("dynamic resolution", params[4], settings["resolution"]),
                ("dynamic intensity scale", params[5], settings["intensity_scale"]),
                ("dynamic bkg a", params[6], settings["bkg_a"]),
                ("dynamic bkg b", params[7], settings["bkg_b"]),
            )
            rows.extend(
                (name, float(value), float(bounds[1]), float(bounds[2]))
                for (name, value, bounds), enabled in zip(base, config["dynamic_fit_flags"])
                if enabled
            )
            if bool(config["substrate_scale_fit_flag"]):
                bounds = config["substrate_settings"]["scale"]
                rows.append(
                    (
                        "substrate lattice scale",
                        float(params[8]),
                        float(bounds[1]),
                        float(bounds[2]),
                    )
                )
            offset = 9

        if bool(config["include_strain"]):
            settings = config["strain_settings"]
            strain_rows = tuple(
                (
                    name.replace("_", " "),
                    float(params[offset + index]),
                    float(settings[name][1]),
                    float(settings[name][2]),
                )
                for index, name in enumerate(STRAIN_KEYS)
            )
            rows.extend(
                row
                for row, enabled in zip(strain_rows, config["strain_fit_flags"])
                if enabled
            )
            offset += 4
        if bool(config["include_roughness"]):
            settings = config["roughness_settings"]
            rough_rows = (
                ("film roughness", float(params[offset]), float(settings["film"][1]), float(settings["film"][2])),
                (
                    "substrate/interface roughness",
                    float(params[offset + 1]),
                    float(settings["substrate"][1]),
                    float(settings["substrate"][2]),
                ),
            )
            rows.extend(row for row, enabled in zip(rough_rows, config["roughness_fit_flags"]) if enabled)

        warnings: list[str] = []
        for name, value, lower, upper in rows:
            if upper <= lower:
                continue
            ratio = (value - lower) / (upper - lower)
            if ratio <= 0.02:
                warnings.append(f"{name} reached the lower bound region: {value:.6g} [{lower:.6g}, {upper:.6g}]")
            elif ratio >= 0.98:
                warnings.append(f"{name} reached the upper bound region: {value:.6g} [{lower:.6g}, {upper:.6g}]")
            elif ratio <= 0.05:
                warnings.append(f"{name} is close to the lower bound: {value:.6g} [{lower:.6g}, {upper:.6g}]")
            elif ratio >= 0.95:
                warnings.append(f"{name} is close to the upper bound: {value:.6g} [{lower:.6g}, {upper:.6g}]")
        return warnings

    def _build_run_config(self) -> dict[str, object]:
        data_path = resolve_data_path(self.data_path_var.get())
        wavelength = self._parse_wavelength()
        sample = sample_from_data_path(data_path, SAMPLES[self.sample_var.get()])
        if sample.name != self.sample_var.get():
            self._apply_sample_profile(sample)
        if not self.min_var.get().strip() or not self.max_var.get().strip():
            all_twotheta, _ = read_experimental_data(data_path)
            self._set_twotheta_window(all_twotheta)

        kinematic_fit_flags = (
            bool(self.kin_d_fit_enabled_var.get()),
            bool(self.kin_planes_fit_enabled_var.get()),
            bool(self.kin_resolution_fit_enabled_var.get()),
            bool(self.kin_scale_fit_enabled_var.get()),
            bool(self.kin_bkg_a_fit_enabled_var.get()),
            bool(self.kin_bkg_b_fit_enabled_var.get()),
        )
        film_fit_flags = (
            bool(self.film_n_fit_enabled_var.get()),
            bool(self.film_scale_fit_enabled_var.get()),
            bool(self.film_area_fit_enabled_var.get()),
            bool(self.film_interface_fit_enabled_var.get()),
        )
        dynamic_fit_flags = (
            bool(self.dynamic_resolution_fit_enabled_var.get()),
            bool(self.dynamic_intensity_fit_enabled_var.get()),
            bool(self.dynamic_bkg_a_fit_enabled_var.get()),
            bool(self.dynamic_bkg_b_fit_enabled_var.get()),
        )
        strain_fit_flags = tuple(
            bool(variable.get()) for variable in self.strain_fit_enabled_vars
        )
        substrate_scale_fit_flag = bool(self.substrate_scale_fit_enabled_var.get())
        roughness_fit_flags = (
            bool(self.film_rough_fit_enabled_var.get()),
            bool(self.substrate_rough_fit_enabled_var.get()),
        )

        kinematic_d = (
            float(self.kin_d_start_var.get()),
            float(self.kin_d_min_var.get()),
            float(self.kin_d_max_var.get()),
        )
        kinematic_planes = (
            float(self.kin_planes_start_var.get()),
            float(self.kin_planes_min_var.get()),
            float(self.kin_planes_max_var.get()),
        )
        kinematic_resolution = (
            float(self.kin_resolution_start_var.get()),
            float(self.kin_resolution_min_var.get()),
            float(self.kin_resolution_max_var.get()),
        )
        kinematic_scale = (
            float(self.kin_scale_start_var.get()),
            float(self.kin_scale_min_var.get()),
            float(self.kin_scale_max_var.get()),
        )
        kinematic_bkg_a = (
            float(self.kin_bkg_a_start_var.get()),
            float(self.kin_bkg_a_min_var.get()),
            float(self.kin_bkg_a_max_var.get()),
        )
        kinematic_bkg_b = (
            float(self.kin_bkg_b_start_var.get()),
            float(self.kin_bkg_b_min_var.get()),
            float(self.kin_bkg_b_max_var.get()),
        )
        kinematic_debye = float(self.kin_debye_var.get())
        for enabled, name, values in zip(
            kinematic_fit_flags,
            (
                "kinematic d spacing",
                "kinematic planes",
                "kinematic resolution",
                "kinematic scale",
                "kinematic bkg a",
                "kinematic bkg b",
            ),
            (
                kinematic_d,
                kinematic_planes,
                kinematic_resolution,
                kinematic_scale,
                kinematic_bkg_a,
                kinematic_bkg_b,
            ),
        ):
            if enabled:
                validate_start_min_max(name, *values)
        if kinematic_d[0] <= 0 or (kinematic_fit_flags[0] and kinematic_d[1] <= 0):
            raise ValueError("kinematic d spacing must be positive")
        if kinematic_planes[0] <= 0 or (kinematic_fit_flags[1] and kinematic_planes[1] <= 0):
            raise ValueError("kinematic planes must be positive")
        if kinematic_resolution[0] <= 0 or (kinematic_fit_flags[2] and kinematic_resolution[1] <= 0):
            raise ValueError("kinematic resolution must be positive")
        if kinematic_scale[0] <= 0 or (kinematic_fit_flags[3] and kinematic_scale[1] <= 0):
            raise ValueError("kinematic scale must be positive")

        film_filename = self.film_filename_var.get()
        substrate_filename = self.substrate_filename_var.get()
        validate_poscar("film POSCAR", film_filename)
        validate_poscar("substrate POSCAR", substrate_filename)

        film_direction = int(self.film_direction_var.get())
        substrate_direction = int(self.substrate_direction_var.get())
        validate_direction("film", film_direction)
        validate_direction("substrate", substrate_direction)

        film_n = (
            float(self.film_n_start_var.get()),
            float(self.film_n_min_var.get()),
            float(self.film_n_max_var.get()),
        )
        film_scale = (
            float(self.film_scale_start_var.get()),
            float(self.film_scale_min_var.get()),
            float(self.film_scale_max_var.get()),
        )
        film_area_scale = (
            float(self.film_area_start_var.get()),
            float(self.film_area_min_var.get()),
            float(self.film_area_max_var.get()),
        )
        film_dinterface = (
            float(self.film_interface_start_var.get()),
            float(self.film_interface_min_var.get()),
            float(self.film_interface_max_var.get()),
        )
        for enabled, name, values in zip(
            film_fit_flags,
            ("film N", "film scale", "film area_scale", "film dinterface"),
            (film_n, film_scale, film_area_scale, film_dinterface),
        ):
            if enabled:
                validate_start_min_max(name, *values)
        if film_n[0] <= 0 or (film_fit_flags[0] and film_n[1] <= 0):
            raise ValueError("film N must be positive")
        if film_scale[0] <= 0 or (film_fit_flags[1] and film_scale[1] <= 0):
            raise ValueError("film scale must be positive")
        if film_area_scale[0] <= 0 or (film_fit_flags[2] and film_area_scale[1] <= 0):
            raise ValueError("film area_scale must be positive")

        dynamic_resolution = (
            float(self.dynamic_resolution_start_var.get()),
            float(self.dynamic_resolution_min_var.get()),
            float(self.dynamic_resolution_max_var.get()),
        )
        dynamic_intensity_scale = (
            float(self.dynamic_intensity_start_var.get()),
            float(self.dynamic_intensity_min_var.get()),
            float(self.dynamic_intensity_max_var.get()),
        )
        dynamic_bkg_a = (
            float(self.dynamic_bkg_a_start_var.get()),
            float(self.dynamic_bkg_a_min_var.get()),
            float(self.dynamic_bkg_a_max_var.get()),
        )
        dynamic_bkg_b = (
            float(self.dynamic_bkg_b_start_var.get()),
            float(self.dynamic_bkg_b_min_var.get()),
            float(self.dynamic_bkg_b_max_var.get()),
        )
        for enabled, name, values in zip(
            dynamic_fit_flags,
            (
                "dynamic resolution",
                "dynamic intensity scale",
                "dynamic bkg a",
                "dynamic bkg b",
            ),
            (
                dynamic_resolution,
                dynamic_intensity_scale,
                dynamic_bkg_a,
                dynamic_bkg_b,
            ),
        ):
            if enabled:
                validate_start_min_max(name, *values)
        if dynamic_resolution[0] <= 0 or (dynamic_fit_flags[0] and dynamic_resolution[1] <= 0):
            raise ValueError("dynamic resolution must be positive")
        if dynamic_intensity_scale[0] <= 0 or (dynamic_fit_flags[1] and dynamic_intensity_scale[1] <= 0):
            raise ValueError("dynamic intensity scale must be positive")

        substrate_n = float(self.substrate_n_var.get())
        substrate_dinterface = float(self.substrate_interface_var.get())
        substrate_scale = (
            float(self.substrate_scale_var.get()),
            float(self.substrate_scale_min_var.get()),
            float(self.substrate_scale_max_var.get()),
        )
        substrate_area_scale = float(self.substrate_area_var.get())
        if substrate_n <= 0:
            raise ValueError("substrate N must be positive")
        if substrate_scale_fit_flag:
            validate_start_min_max("substrate lattice scale", *substrate_scale)
        if substrate_scale[0] <= 0 or (substrate_scale_fit_flag and substrate_scale[1] <= 0):
            raise ValueError("substrate scale must be positive")
        if substrate_area_scale <= 0:
            raise ValueError("substrate area_scale must be positive")

        include_strain = bool(self.strain_var.get())
        if include_strain:
            strain_settings = {
                key: (
                    float(self.strain_start_vars[index].get()),
                    float(self.strain_min_vars[index].get()),
                    float(self.strain_max_vars[index].get()),
                )
                for index, key in enumerate(STRAIN_KEYS)
            }
            for key, enabled in zip(STRAIN_KEYS, strain_fit_flags):
                if enabled:
                    validate_start_min_max(key.replace("_", " "), *strain_settings[key])
            for key in ("bottom_extent", "top_extent"):
                value, lower, _ = strain_settings[key]
                index = STRAIN_KEYS.index(key)
                if value < 0 or (strain_fit_flags[index] and lower < 0):
                    raise ValueError(f"{key.replace('_', ' ')} must be non-negative")
        else:
            strain_settings = {
                key: (0.0, 0.0, 0.0) for key in STRAIN_KEYS
            }

        film_roughness = (
            float(self.film_rough_start_var.get()),
            float(self.film_rough_min_var.get()),
            float(self.film_rough_max_var.get()),
        )
        substrate_roughness = (
            float(self.substrate_rough_start_var.get()),
            float(self.substrate_rough_min_var.get()),
            float(self.substrate_rough_max_var.get()),
        )
        if roughness_fit_flags[0]:
            validate_start_min_max("film roughness", *film_roughness)
        if roughness_fit_flags[1]:
            validate_start_min_max("substrate/interface roughness", *substrate_roughness)
        twotheta_min, twotheta_max = self._window_twotheta_limits()
        dynamic_backend = self.dynamic_backend_var.get().strip().lower()
        if dynamic_backend not in {"auto", "fused", "legacy"}:
            raise ValueError("dynamic backend must be auto, fused, or legacy")
        return {
            "sample_profile": self.sample_var.get(),
            "data_path": str(data_path),
            "model": self.model_var.get(),
            "wavelength": wavelength,
            "dynamic_backend": dynamic_backend,
            "twotheta_min": twotheta_min,
            "twotheta_max": twotheta_max,
            "seed": int(self.seed_var.get()),
            "maxiter": int(self.maxiter_var.get()),
            "popsize": int(self.popsize_var.get()),
            "local_nfev": int(self.local_var.get()),
            "polish_iter": int(self.polish_var.get()),
            "interval": max(1, int(self.interval_var.get())),
            "include_strain": include_strain,
            "include_roughness": bool(self.roughness_var.get()),
            "kinematic_fit_flags": kinematic_fit_flags,
            "film_fit_flags": film_fit_flags,
            "dynamic_fit_flags": dynamic_fit_flags,
            "strain_fit_flags": strain_fit_flags,
            "substrate_scale_fit_flag": substrate_scale_fit_flag,
            "roughness_fit_flags": roughness_fit_flags,
            "kinematic_settings": {
                "d_spacing": kinematic_d,
                "planes": kinematic_planes,
                "resolution": kinematic_resolution,
                "scale": kinematic_scale,
                "bkg_a": kinematic_bkg_a,
                "bkg_b": kinematic_bkg_b,
                "debye_waller_coeff": kinematic_debye,
            },
            "film_settings": {
                "filename": film_filename,
                "direction": film_direction,
                "n": film_n,
                "scale": film_scale,
                "area_scale": film_area_scale,
                "dinterface": film_dinterface,
            },
            "dynamic_fit_settings": {
                "resolution": dynamic_resolution,
                "intensity_scale": dynamic_intensity_scale,
                "bkg_a": dynamic_bkg_a,
                "bkg_b": dynamic_bkg_b,
            },
            "substrate_settings": {
                "filename": substrate_filename,
                "direction": substrate_direction,
                "n": substrate_n,
                "dinterface": substrate_dinterface,
                "scale": substrate_scale,
                "area_scale": substrate_area_scale,
            },
            "strain_settings": strain_settings,
            "roughness_settings": {
                "film": film_roughness,
                "substrate": substrate_roughness,
            },
        }

    def _make_model_and_start(
        self, config: dict[str, object]
    ) -> tuple[SampleConfig, KinematicModel | DynamicModel, np.ndarray, np.ndarray, dict[str, object], dict[str, object]]:
        sample = SAMPLES[str(config["sample_profile"])]
        kinematic_settings = config["kinematic_settings"]
        film_settings = config["film_settings"]
        dynamic_fit_settings = config["dynamic_fit_settings"]
        substrate_settings = config["substrate_settings"]
        twotheta, observed = load_sample_data(
            float(config["twotheta_min"]), float(config["twotheta_max"]),
            str(config["data_path"]),
        )
        if config["model"] == "Kinematic":
            model = KinematicModel(
                twotheta,
                observed,
                sample,
                debye_waller_coeff=float(kinematic_settings["debye_waller_coeff"]),
                include_strain=bool(config["include_strain"]),
                include_roughness=bool(config["include_roughness"]),
                wavelength=float(config["wavelength"]),
            )
            bounds_array, start = kinematic_bounds_and_start(
                sample,
                bool(config["include_strain"]),
                bool(config["include_roughness"]),
                config["roughness_settings"],
                kinematic_settings,
                config["strain_settings"],
            )
        else:
            model = DynamicModel(
                twotheta,
                observed,
                include_strain=bool(config["include_strain"]),
                include_roughness=bool(config["include_roughness"]),
                film_filename=str(film_settings["filename"]),
                film_direction=int(film_settings["direction"]),
                substrate_filename=str(substrate_settings["filename"]),
                substrate_direction=int(substrate_settings["direction"]),
                substrate_n=float(substrate_settings["n"]),
                substrate_dinterface=float(substrate_settings["dinterface"]),
                substrate_scale=float(substrate_settings["scale"][0]),
                substrate_area_scale=float(substrate_settings["area_scale"]),
                wavelength=float(config["wavelength"]),
                propagation_backend=str(config["dynamic_backend"]),
            )
            bounds_array, start = dynamic_bounds_and_start(
                sample,
                bool(config["include_strain"]),
                bool(config["include_roughness"]),
                config["roughness_settings"],
                film_settings,
                dynamic_fit_settings,
                config["strain_settings"],
                substrate_settings,
            )
        return sample, model, bounds_array, start, film_settings, substrate_settings

    def _fit_mask_for_config(self, config: dict[str, object], n_params: int) -> np.ndarray:
        if config["model"] == "Kinematic":
            flags = list(config["kinematic_fit_flags"])
        else:
            flags = list(config["film_fit_flags"]) + list(config["dynamic_fit_flags"])
            flags.append(bool(config["substrate_scale_fit_flag"]))

        if bool(config["include_strain"]):
            flags.extend(list(config["strain_fit_flags"]))
        if bool(config["include_roughness"]):
            flags.extend(list(config["roughness_fit_flags"]))

        mask = np.asarray(flags, dtype=bool)
        if len(mask) != n_params:
            raise ValueError("Internal fit-mask size does not match the parameter vector")
        if not np.any(mask):
            raise ValueError("Select at least one parameter to fit")
        return mask

    def _has_selected_fit_parameter(self, config: dict[str, object]) -> bool:
        if config["model"] == "Kinematic":
            flags = list(config["kinematic_fit_flags"])
        else:
            flags = list(config["film_fit_flags"]) + list(config["dynamic_fit_flags"])
            flags.append(bool(config["substrate_scale_fit_flag"]))
        if bool(config["include_strain"]):
            flags.extend(list(config["strain_fit_flags"]))
        if bool(config["include_roughness"]):
            flags.extend(list(config["roughness_fit_flags"]))
        return any(flags)

    def simulate_pattern(self) -> None:
        if self.running:
            messagebox.showinfo("Fit running", "Stop the fit before simulating a new pattern.")
            return

        try:
            config = self._build_run_config()
            sample, model, _bounds_array, start, film_settings, substrate_settings = self._make_model_and_start(config)
            predicted = model.predict(start)
            cost = model.objective(start)
            residual = predicted - model.observed
            rmse = float(np.sqrt(np.mean(residual**2)))
            density_z = getattr(model, "last_z", None)
            density_rho_e = getattr(model, "last_rho_e", None)
            self.history_x.clear()
            self.history_y.clear()
            self._clear_fit_result_values()
            self.summary_text.delete("1.0", tk.END)
            self._apply_update(
                FitUpdate(
                    phase="simulation",
                    cost=cost,
                    twotheta=model.twotheta,
                    q=model.q,
                    observed=model.observed,
                    predicted=predicted,
                    params=start,
                    density_z=None if density_z is None else np.asarray(density_z),
                    density_rho_e=None if density_rho_e is None else np.asarray(density_rho_e),
                )
            )
            summary = summarize_fit(
                sample,
                str(config["model"]),
                start,
                cost,
                rmse,
                bool(config["include_strain"]),
                bool(config["include_roughness"]),
                film_settings if config["model"] == "Dynamic" else None,
                substrate_settings if config["model"] == "Dynamic" else None,
            )
            backend_text = (
                f"Dynamic backend: {config['dynamic_backend']}\n"
                if config["model"] == "Dynamic"
                else ""
            )
            summary = (
                f"Simulation using current parameter values\n"
                f"Data file: {resolve_data_path(str(config['data_path']))}\n"
                f"X-ray wavelength: {float(config['wavelength']):.6g} A\n"
                f"{backend_text}"
                + summary
            )
            self.summary_text.insert(tk.END, summary)
            self.status_var.set(f"Simulation complete: cost={cost:.5g}")
        except ValueError as exc:
            messagebox.showerror("Invalid input", str(exc))
        except Exception as exc:  # noqa: BLE001 - report simulation exceptions to GUI
            messagebox.showerror("Simulation failed", str(exc))

    def run_fit(self) -> None:
        if self.running:
            messagebox.showinfo("Fit running", "A fit is already running.")
            return

        try:
            config = self._build_run_config()
        except ValueError as exc:
            messagebox.showerror("Invalid input", str(exc))
            return
        if not self._has_selected_fit_parameter(config):
            messagebox.showerror("Invalid input", "Select at least one parameter to fit.")
            return

        self.running = True
        self.active_fit_config = config
        self.stop_event.clear()
        self._clear_fit_result_values()
        self.simulate_button.configure(state=tk.DISABLED)
        self.run_button.configure(state=tk.DISABLED)
        self.stop_button.configure(state=tk.NORMAL)
        self.history_x.clear()
        self.history_y.clear()
        self.summary_text.delete("1.0", tk.END)
        self.status_var.set("Running fit...")
        self._draw_empty_plot()

        thread = threading.Thread(target=self._fit_worker, args=(config,), daemon=True)
        thread.start()

    def stop_fit(self) -> None:
        if self.running:
            self.stop_event.set()
            self.status_var.set("Stopping after current optimizer evaluation...")

    def _fit_worker(self, config: dict[str, object]) -> None:
        try:
            sample, model, bounds_array, start, film_settings, substrate_settings = self._make_model_and_start(config)
            fit_mask = self._fit_mask_for_config(config, len(start))
            fitted_start = start[fit_mask]
            fitted_bounds_array = bounds_array[fit_mask]
            bounds = [tuple(row) for row in fitted_bounds_array]
            counter = {"value": 0}

            def full_params(fitted_params: np.ndarray) -> np.ndarray:
                params = start.copy()
                params[fit_mask] = fitted_params
                return params

            def check_stop() -> None:
                if self.stop_event.is_set():
                    raise FitCancelled("Fit stopped by user")

            def objective(fitted_params: np.ndarray) -> float:
                check_stop()
                return model.objective(full_params(fitted_params))

            def residual_vector(fitted_params: np.ndarray) -> np.ndarray:
                check_stop()
                return model.residual_vector(full_params(fitted_params))

            def send_update(phase: str, params: np.ndarray, force: bool = False) -> None:
                check_stop()
                counter["value"] += 1
                if not force and counter["value"] % int(config["interval"]) != 0:
                    return
                predicted = model.predict(params)
                cost = model.objective(params)
                density_z = getattr(model, "last_z", None)
                density_rho_e = getattr(model, "last_rho_e", None)
                self.queue.put(
                    (
                        "update",
                        FitUpdate(
                            phase=phase,
                            cost=cost,
                            twotheta=model.twotheta,
                            q=model.q,
                            observed=model.observed,
                            predicted=predicted,
                            params=np.asarray(params, dtype=float),
                            density_z=None if density_z is None else np.asarray(density_z),
                            density_rho_e=None if density_rho_e is None else np.asarray(density_rho_e),
                        ),
                    )
                )

            send_update("initial", start, force=True)

            def de_callback(xk: np.ndarray, convergence: float | None = None) -> bool:
                send_update("differential evolution", full_params(np.asarray(xk)))
                return self.stop_event.is_set()

            global_result = differential_evolution(
                objective,
                bounds,
                seed=int(config["seed"]),
                maxiter=int(config["maxiter"]),
                popsize=int(config["popsize"]),
                tol=2e-4,
                polish=False,
                updating="immediate",
                workers=1,
                callback=de_callback,
            )
            check_stop()

            local_start = (
                global_result.x if global_result.fun < objective(fitted_start) else fitted_start
            )
            local_result = least_squares(
                residual_vector,
                local_start,
                bounds=(fitted_bounds_array[:, 0], fitted_bounds_array[:, 1]),
                xtol=1e-8,
                ftol=1e-8,
                gtol=1e-8,
                max_nfev=int(config["local_nfev"]),
            )
            check_stop()
            send_update("least squares", full_params(local_result.x), force=True)

            def polish_callback(xk: np.ndarray) -> None:
                send_update("Powell polish", full_params(np.asarray(xk)))

            polish_result = minimize(
                objective,
                global_result.x,
                method="Powell",
                bounds=bounds,
                callback=polish_callback,
                options={
                    "maxiter": int(config["polish_iter"]),
                    "xtol": 1e-6,
                    "ftol": 1e-6,
                },
            )
            check_stop()

            candidates = [
                start,
                full_params(global_result.x),
                full_params(local_result.x),
                full_params(polish_result.x),
            ]
            best_params = min(candidates, key=model.objective)
            predicted = model.predict(best_params)
            residual = predicted - model.observed
            cost = model.objective(best_params)
            rmse = float(np.sqrt(np.mean(residual**2)))
            send_update("best fit", best_params, force=True)

            model_slug = str(config["model"]).lower()
            sample_slug = resolve_data_path(str(config["data_path"])).stem.lower().replace(" ", "_").replace("/", "_")
            output_path = ROOT / "validation" / f"{sample_slug}_gui_{model_slug}_fit.csv"
            np.savetxt(
                output_path,
                np.column_stack([model.twotheta, model.q, model.observed, predicted, residual]),
                delimiter=",",
                header="twotheta_deg,Q_inv_A,observed_cps,fitted_cps,residual_cps",
                comments="",
            )
            figure_path = ROOT / "validation" / f"{sample_slug}_gui_{model_slug}_fit_progress.png"
            summary = summarize_fit(
                sample,
                str(config["model"]),
                best_params,
                cost,
                rmse,
                bool(config["include_strain"]),
                bool(config["include_roughness"]),
                film_settings if config["model"] == "Dynamic" else None,
                substrate_settings if config["model"] == "Dynamic" else None,
            )
            backend_text = (
                f"Dynamic backend: {config['dynamic_backend']}\n"
                if config["model"] == "Dynamic"
                else ""
            )
            summary = (
                f"Data file: {resolve_data_path(str(config['data_path']))}\n"
                f"X-ray wavelength: {float(config['wavelength']):.6g} A\n"
                f"{backend_text}"
                + summary
            )
            summary += f"\nCSV: {output_path}\nPNG: {figure_path}"
            self.queue.put(("done", (summary, figure_path, config, best_params)))
        except FitCancelled as exc:
            self.queue.put(("cancelled", str(exc)))
        except Exception as exc:  # noqa: BLE001 - report worker exceptions to GUI
            self.queue.put(("error", exc))

    def _process_queue(self) -> None:
        while True:
            try:
                message_type, payload = self.queue.get_nowait()
            except queue.Empty:
                break

            if message_type == "update":
                self._apply_update(payload)
            elif message_type == "done":
                summary, figure_path, config, best_params = payload
                self.figure.savefig(figure_path, dpi=180)
                best_params = np.asarray(best_params, dtype=float)
                warnings = self._fit_boundary_warnings(config, best_params)
                if warnings:
                    summary += "\n\nBoundary warnings:\n" + "\n".join(f"- {warning}" for warning in warnings)
                self.summary_text.delete("1.0", tk.END)
                self.summary_text.insert(tk.END, summary)
                self._set_fit_result_values(config, best_params)
                if warnings:
                    self.status_var.set(f"Fit complete with {len(warnings)} boundary warning(s)")
                else:
                    self.status_var.set("Fit complete")
                self.running = False
                self.active_fit_config = None
                self.simulate_button.configure(state=tk.NORMAL)
                self.run_button.configure(state=tk.NORMAL)
                self.stop_button.configure(state=tk.DISABLED)
            elif message_type == "cancelled":
                self.status_var.set(str(payload))
                self.running = False
                self.active_fit_config = None
                self.simulate_button.configure(state=tk.NORMAL)
                self.run_button.configure(state=tk.NORMAL)
                self.stop_button.configure(state=tk.DISABLED)
            elif message_type == "error":
                self.status_var.set("Fit failed")
                self.running = False
                self.active_fit_config = None
                self.simulate_button.configure(state=tk.NORMAL)
                self.run_button.configure(state=tk.NORMAL)
                self.stop_button.configure(state=tk.DISABLED)
                messagebox.showerror("Fit failed", str(payload))

        self.root.after(150, self._process_queue)

    def _apply_update(self, update: FitUpdate) -> None:
        self.last_update = update
        if self.active_fit_config is not None:
            self._set_fit_result_values(self.active_fit_config, update.params)
        self.history_x.append(len(self.history_x) + 1)
        self.history_y.append(update.cost)
        self.status_var.set(f"{update.phase}: cost={update.cost:.5g}")

        self.loss_axis.clear()
        self.loss_axis.plot(self.history_x, self.history_y, color="tab:blue", linewidth=1.5)
        self.loss_axis.set_xlabel("progress callback")
        self.loss_axis.set_ylabel("mean abs log10 error")
        self.loss_axis.set_title(f"{update.phase}: cost={update.cost:.5g}")
        self.loss_axis.grid(True, alpha=0.25)

        self.fit_axis.clear()
        x_values = self._x_values(update.twotheta, update.q, CU_K_ALPHA_WAVELENGTH)
        self.fit_axis.plot(x_values, update.observed, ".", color="black", markersize=3, label="data")
        self.fit_axis.plot(
            x_values,
            update.predicted,
            color="tab:red",
            linewidth=1.4,
            label="simulation" if update.phase == "simulation" else "current fit",
        )
        self.fit_axis.set_yscale("log")
        self.fit_axis.set_xlabel(self._axis_label())
        self.fit_axis.set_ylabel("intensity (cps)")
        self.fit_axis.legend(loc="best")
        self.fit_axis.grid(True, alpha=0.25)

        self.density_axis.clear()
        if update.density_z is not None and update.density_rho_e is not None:
            self.density_axis.plot(
                update.density_z,
                np.real(update.density_rho_e),
                color="tab:purple",
                linewidth=1.2,
                label="Re(rho_e)",
            )
            self.density_axis.axhline(0.0, color="black", linewidth=0.7, alpha=0.35)
            self.density_axis.legend(loc="best")
            self.density_axis.set_title("electron density profile")
        else:
            self.density_axis.text(
                0.5,
                0.5,
                "electron density profile is available for dynamic fits",
                transform=self.density_axis.transAxes,
                ha="center",
                va="center",
            )
        self.density_axis.set_xlabel("z (A)")
        self.density_axis.set_ylabel("Re density")
        self.density_axis.grid(True, alpha=0.25)
        self.canvas.draw_idle()


def main() -> int:
    root = tk.Tk()
    FitApp(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
