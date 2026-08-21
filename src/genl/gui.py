from __future__ import annotations

import copy
import json
import os
import queue
import secrets
import tempfile
import threading
import tkinter as tk
import webbrowser
from dataclasses import dataclass
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

import numpy as np

from .paths import EXAMPLE_DATA_DIR, FORM_FACTOR_DIR, REPOSITORY_ROOT, STACK_DIR, STRUCTURE_DIR

ROOT = REPOSITORY_ROOT
CU_K_ALPHA_WAVELENGTH = 1.5406
STRAIN_KEYS = ("bottom_amplitude", "bottom_extent", "top_amplitude", "top_extent")
STACK_STRAIN_FIELDS = (
    ("bottom_strain_amplitude", "Bottom amplitude", -0.4, 0.4),
    ("bottom_strain_end", "Bottom extent (atomic positions)", 0.0, 20.0),
    ("top_strain_amplitude", "Top amplitude", -0.4, 0.4),
    ("top_strain_end", "Top extent (atomic positions)", 0.0, 20.0),
)
GENL_DOI_URL = "https://doi.org/10.1107/S1600576726002566"
PROJECT_FORMAT = "GenL GUI project"
PROJECT_VERSION = 1
GENL_CITATION = (
    "GenL Python version by Vassilios Kapaklis and Gunnar K. Palsson. "
    "Please cite: J. Appl. Cryst. 59, 968-977 (2026)\n"
    f"{GENL_DOI_URL}"
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
    "stack": "#6f843f",
    "simulate": "#2f6fad",
    "run": "#2e8540",
    "stop": "#b83232",
    "entry": "#ffffff",
    "entry_disabled": "#eceff1",
    "watermark": "#6b7280",
}
TOOLTIPS = {
    "scattering_model": (
        "Kinematic uses the single-scattering approximation. Dynamic includes "
        "multiple scattering through the layered structure."
    ),
    "sample_configuration": (
        "Choose between a single film on a substrate and a repeated superlattice structure."
    ),
    "wavelength": (
        "Incident X-ray wavelength in angstroms. The default is Cu K\u03b1 radiation, 1.5406 \u00c5."
    ),
    "horizontal_axis": (
        "Display and define the calculation range using either 2\u03b8 in degrees or "
        "scattering vector q in \u00c5\u207b\u00b9."
    ),
    "include_strain": "Apply the strain profile defined in the Strain tab.",
    "include_roughness": (
        "Apply the film-thickness and substrate/interface roughness defined in the Roughness tab."
    ),
    "fit": "Allow this parameter to vary during fitting.",
    "plane_spacing": "Spacing between consecutive diffracting planes in angstroms.",
    "number_of_planes": "Number of coherently scattering atomic planes in the film.",
    "resolution": "Gaussian instrumental broadening expressed as FWHM in degrees.",
    "intensity_scale": (
        "Multiplicative scale factor applied to the calculated diffraction intensity."
    ),
    "background_a": (
        "Linear coefficient of the background about the center of the fitted range."
    ),
    "background_b": "Constant background intensity.",
    "background_c": "Quadratic coefficient of the centered polynomial background.",
    "debye_waller": (
        "Coefficient controlling attenuation of scattering with increasing q due to atomic displacement."
    ),
    "include_substrate_peak": (
        "Add a Lorentzian substrate reflection to the kinematic calculation."
    ),
    "substrate_intensity": "Integrated intensity of the kinematic substrate reflection.",
    "substrate_fwhm": "Full width at half maximum of the substrate reflection in degrees.",
    "substrate_d_spacing": (
        "Substrate diffracting-plane spacing used to determine the peak position."
    ),
    "structure_file": "VASP structure file defining the unit cell and atomic positions.",
    "layer_direction": (
        "POSCAR lattice-vector direction used as the surface-normal propagation direction: 1, 2, or 3."
    ),
    "number_of_layers": (
        "Number of unit cells propagated along the selected layer direction."
    ),
    "lattice_scale": (
        "Multiplier applied to the lattice spacing along the propagation direction."
    ),
    "area_scale": (
        "Multiplier applied to the in-plane unit-cell area and therefore the projected density."
    ),
    "interface_spacing": (
        "Additional separation between this structure and the preceding layer, in angstroms."
    ),
    "density_slices": (
        "Number of depth slices used to discretize the electron-density profile of each unit cell."
    ),
    "density_q_max": (
        "Maximum Fourier-space sampling limit used to construct the electron-density profile."
    ),
    "substrate_layers": (
        "Effective number of substrate unit cells used in the propagation calculation."
    ),
    "bottom_strain_amplitude": "Strain amplitude applied at the bottom of the layer.",
    "bottom_strain_extent": (
        "Distance over which the strain extends from the bottom into the layer."
    ),
    "top_strain_amplitude": "Strain amplitude applied at the top of the layer.",
    "top_strain_extent": "Distance over which the strain extends from the top into the layer.",
    "film_roughness": (
        "Standard deviation of the film-thickness distribution, expressed in planes or layers."
    ),
    "substrate_roughness": (
        "Gaussian width of the substrate or interface electron-density transition in angstroms."
    ),
    "seed": (
        "Random-number seed used by differential evolution. Reusing a seed makes an optimization reproducible."
    ),
    "new_seed": "Generate a new random seed for the next optimization.",
    "progress_interval": "Number of objective evaluations between GUI and plot updates.",
    "de_iterations": "Maximum number of differential-evolution generations.",
    "de_population": (
        "Differential-evolution population multiplier. Larger values explore more candidates "
        "but require more calculations."
    ),
    "local_evaluations": "Maximum objective evaluations allowed during local refinement.",
    "polish_iterations": "Maximum iterations used for the final polishing step.",
    "repetitions": "Number of times the repeated layer sequence is propagated.",
    "capping_layer": "Add one non-repeated layer after the repeated sequence.",
    "layer_structure": (
        "VASP structure file defining this substrate, repeated layer, or capping layer."
    ),
    "unit_cells": "Number of unit cells of this material included in each repetition.",
    "points_per_unit": (
        "Number of calculated sampling points per selected horizontal-axis unit."
    ),
}

cache_dir = ROOT / ".matplotlib-cache"
cache_dir.mkdir(exist_ok=True)
(cache_dir / "xdg").mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from scipy.optimize import differential_evolution, least_squares, minimize

from .background import centered_polynomial_background  # noqa: E402
from .convolution import gauss_conv  # noqa: E402
from .dynamic import validate_density_sampling  # noqa: E402
from .fit_models import DynamicModel, roughness_distribution  # noqa: E402
from .form_factors import form_factors, read_form_factor_coefficients  # noqa: E402
from .kinematic import matlab_round  # noqa: E402
from .poscar import read_poscar  # noqa: E402
from .stack import StackDefinition, StackModel  # noqa: E402


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


def stack_simulation_grid(
    lower: float,
    upper: float,
    points_per_unit: float,
    q_axis: bool,
    wavelength: float,
) -> tuple[np.ndarray, np.ndarray]:
    if upper <= lower:
        raise ValueError("Superlattice simulation maximum must be greater than minimum")
    if points_per_unit <= 0:
        raise ValueError("Superlattice simulation points per axis unit must be positive")
    point_count = int(np.ceil((upper - lower) * points_per_unit)) + 1
    if point_count > 1_000_000:
        raise ValueError("Superlattice simulation grid exceeds 1,000,000 points")
    axis_values = np.linspace(lower, upper, max(2, point_count))
    if q_axis:
        return np.asarray(twotheta_from_q(axis_values, wavelength)), axis_values
    return axis_values, np.asarray(q_from_twotheta(axis_values, wavelength))


def kinematic_substrate_peak(
    twotheta: np.ndarray,
    integrated_intensity: float,
    fwhm: float,
    d_spacing: float,
    wavelength: float,
) -> np.ndarray:
    center = 2.0 * np.rad2deg(np.arcsin(wavelength / (2.0 * d_spacing)))
    half_width = 0.5 * fwhm
    return integrated_intensity / np.pi * half_width / (
        (np.asarray(twotheta, dtype=float) - center) ** 2 + half_width**2
    )


def _least_squares_residual(residual: np.ndarray) -> np.ndarray:
    """Make least-squares minimize the displayed mean absolute error."""
    residual = np.asarray(residual, dtype=float)
    return np.sign(residual) * np.sqrt(np.abs(residual) / residual.size)


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
    sorted(path.name for path in STRUCTURE_DIR.glob("*.vasp"))
)


def default_data_path(sample: SampleConfig) -> Path:
    return EXAMPLE_DATA_DIR / sample.data_file


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
    structure = read_poscar(sample.film_filename, STRUCTURE_DIR)
    volume = abs(np.dot(structure.a1, np.cross(structure.a2, structure.a3)))
    return float(np.sum(structure.type_counts) / volume)


def propagation_period(sample: SampleConfig) -> float:
    return propagation_period_for(sample.film_filename, sample.dynamic_direction)


def propagation_period_for(filename: str, direction: int) -> float:
    structure = read_poscar(filename, STRUCTURE_DIR)
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
    show_observed: bool = True


def export_result_data(update: FitUpdate, path: Path) -> list[Path]:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    if update.show_observed:
        columns = [
            update.twotheta,
            update.q,
            update.observed,
            update.predicted,
            update.predicted - update.observed,
        ]
        headings = ["twotheta_deg", "q_inv_A", "observed_cps", "fitted_cps", "residual_cps"]
    else:
        columns = [update.twotheta, update.q, update.predicted]
        headings = ["twotheta_deg", "q_inv_A", "simulated_cps"]
    np.savetxt(
        path,
        np.column_stack(columns),
        delimiter=delimiter,
        header=delimiter.join(headings),
        comments="",
    )
    written = [path]
    if update.density_z is not None and update.density_rho_e is not None:
        density_path = path.with_name(f"{path.stem}_density{path.suffix}")
        np.savetxt(
            density_path,
            np.column_stack(
                [update.density_z, np.real(update.density_rho_e), np.imag(update.density_rho_e)]
            ),
            delimiter=delimiter,
            header=delimiter.join(["z_A", "density_real", "density_imag"]),
            comments="",
        )
        written.append(density_path)
    return written


def save_result_plots(update: FitUpdate, wavelength: float, path: Path) -> list[Path]:
    def style_axis(axis, *, top_ticks: bool = True) -> None:
        axis.tick_params(which="both", direction="in", top=top_ticks, right=True)
        axis.minorticks_on()
        for spine in axis.spines.values():
            spine.set_color("black")
            spine.set_linewidth(0.8)

    suffix = path.suffix
    diffraction_path = path.with_name(f"{path.stem}_diffraction{suffix}")
    diffraction_figure = Figure(figsize=(5.2, 4.2), constrained_layout=True)
    diffraction_axis = diffraction_figure.add_subplot(111)
    curve_label = "GenL simulation" if update.phase == "simulation" else "GenL fit"
    diffraction_axis.plot(
        update.twotheta,
        update.predicted,
        color="red",
        linewidth=1.5,
        label=curve_label,
    )
    if update.show_observed:
        diffraction_axis.plot(
            update.twotheta,
            update.observed,
            linestyle="none",
            marker="o",
            markersize=3.5,
            markerfacecolor="#aeb5dc",
            markeredgecolor="#6873a8",
            markeredgewidth=0.5,
            alpha=0.75,
            label="Data",
        )
    diffraction_axis.set_yscale("log")
    diffraction_axis.set_xlim(float(np.min(update.twotheta)), float(np.max(update.twotheta)))
    diffraction_axis.set_xlabel(r"$2\theta$ [degrees]")
    diffraction_axis.set_ylabel("Intensity [cps]")
    diffraction_axis.legend(loc="best", frameon=True, edgecolor="0.75")
    style_axis(diffraction_axis, top_ticks=False)

    def theta_to_q(values: np.ndarray) -> np.ndarray:
        return np.asarray(q_from_twotheta(values, wavelength), dtype=float)

    def q_to_theta(values: np.ndarray) -> np.ndarray:
        argument = np.asarray(values, dtype=float) * wavelength / (4.0 * np.pi)
        return 2.0 * np.rad2deg(np.arcsin(np.clip(argument, -1.0, 1.0)))

    q_axis = diffraction_axis.secondary_xaxis("top", functions=(theta_to_q, q_to_theta))
    q_axis.set_xlabel(r"$Q$ [$\mathrm{\AA}^{-1}$]")
    q_axis.tick_params(which="both", direction="in")
    q_axis.minorticks_on()
    diffraction_figure.savefig(diffraction_path, dpi=300, facecolor="white")
    written = [diffraction_path]

    if update.density_z is not None and update.density_rho_e is not None:
        density_path = path.with_name(f"{path.stem}_density{suffix}")
        density_figure = Figure(figsize=(5.2, 4.2), constrained_layout=True)
        density_axis = density_figure.add_subplot(111)
        density_axis.plot(
            update.density_z,
            np.real(update.density_rho_e),
            color="red",
            linewidth=1.5,
        )
        density_axis.set_xlim(float(np.min(update.density_z)), float(np.max(update.density_z)))
        density_axis.set_xlabel(r"Depth $z$ [$\mathrm{\AA}$]")
        density_axis.set_ylabel(r"Electron density [$e\,\mathrm{\AA}^{-3}$]")
        style_axis(density_axis)
        density_figure.savefig(density_path, dpi=300, facecolor="white")
        written.append(density_path)

    return written


@dataclass(frozen=True)
class FitCheckpoint:
    signature: tuple[object, ...]
    fit_mask: np.ndarray
    params: np.ndarray
    population: np.ndarray | None
    rng_state: dict[str, object]
    phase: str
    cost: float


def fit_resume_signature(
    config: dict[str, object],
    start: np.ndarray,
    fit_mask: np.ndarray,
) -> tuple[object, ...]:
    if bool(config.get("stack_enabled", False)):
        document = copy.deepcopy(config["stack_document"])
        targets = [
            str(parameter["target"])
            for parameter in document.get("fit_parameters", [])
        ]
        document["fit_parameters"] = targets
        document.pop("calculation_ranges", None)
        layers = {
            "substrate": document["substrate"],
            **{
                str(layer["name"]): layer
                for layer in document["sequence"]["layers"]
            },
            "calculation": document["calculation"],
        }
        if document.get("capping_layer") is not None:
            layers["capping"] = document["capping_layer"]
        for target in targets:
            prefix, field = target.split(".", 1)
            layers[prefix][field] = "<fitted>"
        return (
            "superlattice",
            str(resolve_data_path(str(config["data_path"]))),
            float(config["twotheta_min"]),
            float(config["twotheta_max"]),
            int(config["seed"]),
            int(config["popsize"]),
            json.dumps(document, sort_keys=True),
            tuple(bool(value) for value in fit_mask),
            tuple(float(value) for value in start[~fit_mask]),
        )
    film = config["film_settings"]
    substrate = config["substrate_settings"]
    kinematic = config["kinematic_settings"]
    return (
        config["sample_profile"],
        str(resolve_data_path(str(config["data_path"]))),
        config["model"],
        float(config["wavelength"]),
        int(config.get("density_slices", 100)),
        float(config.get("density_max_q0", 30.0)),
        float(config["twotheta_min"]),
        float(config["twotheta_max"]),
        bool(config["include_strain"]),
        bool(config["include_roughness"]),
        bool(config.get("include_kinematic_substrate", False)),
        int(config["seed"]),
        int(config["popsize"]),
        tuple(bool(value) for value in fit_mask),
        tuple(float(value) for value in start[~fit_mask]),
        str(film["filename"]),
        int(film["direction"]),
        str(substrate["filename"]),
        int(substrate["direction"]),
        float(substrate["n"]),
        float(substrate["dinterface"]),
        float(substrate["area_scale"]),
        float(kinematic["debye_waller_coeff"]),
    )


def clipped_checkpoint_population(
    population: np.ndarray | None,
    fitted_bounds: np.ndarray,
) -> np.ndarray | None:
    if population is None:
        return None
    population = np.asarray(population, dtype=float)
    if population.ndim != 2 or population.shape[1] != len(fitted_bounds):
        raise ValueError("Saved optimizer population does not match the selected fit parameters")
    return np.clip(population, fitted_bounds[:, 0], fitted_bounds[:, 1])


def fit_update_to_dict(update: FitUpdate) -> dict[str, object]:
    density = None if update.density_rho_e is None else np.asarray(update.density_rho_e)
    return {
        "phase": update.phase,
        "cost": update.cost,
        "parameters": update.params.tolist(),
        "twotheta_deg": update.twotheta.tolist(),
        "q_inv_A": update.q.tolist(),
        "show_observed": update.show_observed,
        "observed_cps": update.observed.tolist() if update.show_observed else None,
        "predicted_cps": update.predicted.tolist(),
        "residual_cps": (
            (update.predicted - update.observed).tolist() if update.show_observed else None
        ),
        "density_z_A": None if update.density_z is None else update.density_z.tolist(),
        "density_real": None if density is None else np.real(density).tolist(),
        "density_imag": None if density is None else np.imag(density).tolist(),
    }


def write_json_document(path: Path, document: dict[str, object]) -> None:
    payload = json.dumps(document, indent=2, ensure_ascii=False) + "\n"
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            handle.write(payload)
        os.replace(temporary_path, path)
    except Exception:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
        raise


def fit_update_from_dict(data: dict[str, object]) -> FitUpdate:
    twotheta = np.asarray(data["twotheta_deg"], dtype=float)
    q = np.asarray(data["q_inv_A"], dtype=float)
    show_observed = bool(data.get("show_observed", data.get("observed_cps") is not None))
    observed_value = data.get("observed_cps")
    observed = (
        np.asarray(observed_value, dtype=float)
        if observed_value is not None
        else np.zeros_like(twotheta)
    )
    predicted = np.asarray(data["predicted_cps"], dtype=float)
    if not (len(twotheta) == len(q) == len(observed) == len(predicted)):
        raise ValueError("Saved result arrays have inconsistent lengths")

    density_z_value = data.get("density_z_A")
    density_real_value = data.get("density_real")
    density_z = None if density_z_value is None else np.asarray(density_z_value, dtype=float)
    density = None
    if density_real_value is not None:
        density_real = np.asarray(density_real_value, dtype=float)
        density_imag = np.asarray(data.get("density_imag", np.zeros_like(density_real)), dtype=float)
        if density_z is None or len(density_z) != len(density_real) or len(density_real) != len(density_imag):
            raise ValueError("Saved density arrays have inconsistent lengths")
        density = density_real + 1j * density_imag

    return FitUpdate(
        phase=str(data["phase"]),
        cost=float(data["cost"]),
        twotheta=twotheta,
        q=q,
        observed=observed,
        predicted=predicted,
        params=np.asarray(data["parameters"], dtype=float),
        density_z=density_z,
        density_rho_e=density,
        show_observed=show_observed,
    )


class ToolTip:
    def __init__(self, widget: tk.Widget, text: str, delay_ms: int = 500) -> None:
        self.widget = widget
        self.text = text
        self.delay_ms = delay_ms
        self.after_id: str | None = None
        self.window: tk.Toplevel | None = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")
        widget.bind("<Destroy>", self._hide, add="+")

    def _schedule(self, _event: tk.Event | None = None) -> None:
        self._cancel()
        self.after_id = self.widget.after(self.delay_ms, self._show)

    def _cancel(self) -> None:
        if self.after_id is not None:
            self.widget.after_cancel(self.after_id)
            self.after_id = None

    def _show(self) -> None:
        self.after_id = None
        if self.window is not None or not self.widget.winfo_exists():
            return
        self.window = tk.Toplevel(self.widget)
        self.window.overrideredirect(True)
        label = tk.Label(
            self.window,
            text=self.text,
            justify=tk.LEFT,
            wraplength=340,
            background="#fffbe8",
            foreground="#111827",
            relief=tk.SOLID,
            borderwidth=1,
            padx=7,
            pady=5,
        )
        label.pack()
        self.window.update_idletasks()
        x = min(
            self.widget.winfo_pointerx() + 12,
            self.widget.winfo_screenwidth() - self.window.winfo_reqwidth() - 4,
        )
        y = min(
            self.widget.winfo_pointery() + 18,
            self.widget.winfo_screenheight() - self.window.winfo_reqheight() - 4,
        )
        self.window.geometry(f"+{max(0, x)}+{max(0, y)}")

    def _hide(self, _event: tk.Event | None = None) -> None:
        self._cancel()
        if self.window is not None:
            self.window.destroy()
            self.window = None


def add_tooltip(widget: tk.Widget, key: str) -> tk.Widget:
    ToolTip(widget, TOOLTIPS[key])
    return widget


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
        include_substrate: bool = False,
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
        self.include_substrate = include_substrate
        self.include_strain = include_strain
        self.include_roughness = include_roughness
        self.wavelength = wavelength
        self.thickness = sample.thickness
        self.density = unit_cell_density(sample)
        self.q = q_from_twotheta(twotheta, self.wavelength)
        coefficients = read_form_factor_coefficients(
            sample.atomic_number,
            self.wavelength,
            FORM_FACTOR_DIR,
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

    def _parse_optional_params(
        self, params: np.ndarray
    ) -> tuple[float, float, float, float, float, float, float, float, float]:
        offset = 7
        if self.include_substrate:
            substrate_intensity, substrate_width, substrate_d = params[offset : offset + 3]
            offset += 3
        else:
            substrate_intensity = substrate_width = substrate_d = 0.0
        if self.include_strain:
            bottom_amp, bottom_end, top_amp, top_end = params[offset : offset + 4]
            offset += 4
        else:
            bottom_amp = bottom_end = top_amp = top_end = 0.0
        if self.include_roughness:
            film_sigma, substrate_sigma = params[offset : offset + 2]
        else:
            film_sigma = substrate_sigma = 0.0
        return (
            substrate_intensity,
            substrate_width,
            substrate_d,
            bottom_amp,
            bottom_end,
            top_amp,
            top_end,
            film_sigma,
            substrate_sigma,
        )

    def _single_shape(
        self,
        params: np.ndarray,
        n_planes: float,
        bottom_amp: float,
        bottom_end: float,
        top_amp: float,
        top_end: float,
    ) -> np.ndarray:
        d_spacing, _, resolution, _, _, _, _ = params[:7]
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
        d_spacing, n_planes, resolution, amplitude, bkg_a, bkg_b, bkg_c = params[:7]
        (
            substrate_intensity,
            substrate_width,
            substrate_d,
            bottom_amp,
            bottom_end,
            top_amp,
            top_end,
            film_sigma,
            substrate_sigma,
        ) = self._parse_optional_params(params)
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
        predicted = amplitude * shape + centered_polynomial_background(
            self.q, bkg_a, bkg_b, bkg_c
        )
        if self.include_substrate:
            predicted = predicted + kinematic_substrate_peak(
                self.twotheta,
                substrate_intensity,
                substrate_width,
                substrate_d,
                self.wavelength,
            )
        return predicted

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
        raise ValueError(
            "Data file must contain 2theta/intensity or 2theta/q/intensity columns"
        )
    intensity_column = 2 if data.shape[1] >= 3 else 1
    return data[:, 0], data[:, intensity_column]


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
    include_substrate: bool = False,
    substrate_peak_settings: dict[str, tuple[float, float, float]] | None = None,
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
        bkg_a = (0.0, -3.0, 3.0)
        bkg_b = (0.1, 0.0, 3.0)
        bkg_c = (0.0, -3.0, 3.0)
    else:
        d_spacing = kinematic_settings["d_spacing"]
        planes = kinematic_settings["planes"]
        resolution = kinematic_settings["resolution"]
        scale = kinematic_settings["scale"]
        bkg_a = kinematic_settings["bkg_a"]
        bkg_b = kinematic_settings["bkg_b"]
        bkg_c = kinematic_settings.get("bkg_c", (0.0, -3.0, 3.0))
    bounds_array = np.array(
        [
            [d_spacing[1], d_spacing[2]],
            [planes[1], planes[2]],
            [resolution[1], resolution[2]],
            [scale[1], scale[2]],
            [bkg_a[1], bkg_a[2]],
            [bkg_b[1], bkg_b[2]],
            [bkg_c[1], bkg_c[2]],
        ],
        dtype=float,
    )
    start = np.array(
        [
            d_spacing[0],
            planes[0],
            resolution[0],
            scale[0],
            bkg_a[0],
            bkg_b[0],
            bkg_c[0],
        ],
        dtype=float,
    )
    if include_substrate:
        settings = substrate_peak_settings or {
            "intensity": (51.0, 0.0, 1e6),
            "width": (0.004, 0.0001, 0.05),
            "d_spacing": (sample.d_spacing, sample.d_spacing * 0.995, sample.d_spacing * 1.005),
        }
        substrate_values = [settings[key] for key in ("intensity", "width", "d_spacing")]
        bounds_array = np.vstack(
            [bounds_array, [[value[1], value[2]] for value in substrate_values]]
        )
        start = np.concatenate([start, [value[0] for value in substrate_values]])
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
        bkg_a = (5.1469e-4, -0.1, 0.1)
        bkg_b = (1.2366e-7, 0.0, 0.1)
        bkg_c = (0.0, -0.1, 0.1)
    else:
        resolution = dynamic_fit_settings["resolution"]
        intensity_scale = dynamic_fit_settings["intensity_scale"]
        bkg_a = dynamic_fit_settings["bkg_a"]
        bkg_b = dynamic_fit_settings["bkg_b"]
        bkg_c = dynamic_fit_settings.get("bkg_c", (0.0, -0.1, 0.1))
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
            [bkg_c[1], bkg_c[2]],
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
            bkg_c[0],
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


def validate_start_min_max(
    name: str,
    start: float,
    lower: float,
    upper: float,
    allow_outside_start: bool = False,
) -> None:
    if lower > upper:
        raise ValueError(f"{name}: min must be <= max")
    if not allow_outside_start and not lower <= start <= upper:
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
    film_unit: str,
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
    include_kinematic_substrate: bool = False,
    film_settings: dict[str, object] | None = None,
    substrate_settings: dict[str, object] | None = None,
) -> str:
    if model_name == "Kinematic":
        offset = 7
        substrate_text = ""
        if include_kinematic_substrate:
            substrate_text = (
                f"substrate integrated intensity: {params[offset]:.6e}\n"
                f"substrate FWHM: {params[offset + 1]:.6e} deg\n"
                f"substrate d spacing: {params[offset + 2]:.6f} A\n"
            )
            offset += 3
        return (
            f"{sample.name} kinematic fit\n"
            f"d spacing: {params[0]:.6f} A\n"
            f"coherent planes: {matlab_round(params[1])}\n"
            f"coherent thickness: {matlab_round(params[1]) * params[0]:.3f} A\n"
            f"resolution FWHM Q: {params[2]:.6e} 1/A\n"
            f"scale: {params[3]:.6e}\n"
            f"centered background: {params[5]:.6e} + {params[4]:.6e} x "
            f"+ {params[6]:.6e} x^2\n"
            + substrate_text
            + optional_param_summary(params, offset, include_strain, include_roughness, "planes")
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
        f"centered background: {params[7]:.6e} + {params[6]:.6e} x "
        f"+ {params[8]:.6e} x^2\n"
        + f"substrate lattice scale: {params[9]:.8f}\n"
        + optional_param_summary(params, 10, include_strain, include_roughness, "layers")
        + f"mean abs log10 error: {cost:.6e}\n"
        f"RMSE: {rmse:.6e} cps"
    )


class FitApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("GenL: Fitting Laue oscillation patterns")
        self.queue: queue.Queue[tuple[str, object]] = queue.Queue()
        self.running = False
        self.stop_event = threading.Event()
        self.pause_event = threading.Event()
        self.paused = False
        self.fit_checkpoint: FitCheckpoint | None = None
        self.pending_resume_config: dict[str, object] | None = None
        self.history_x: list[int] = []
        self.history_y: list[float] = []
        self.history_phase: list[str] = []
        self.last_update: FitUpdate | None = None
        self.active_fit_config: dict[str, object] | None = None
        self.preview_after_id: str | None = None
        self.preview_data_path: Path | None = None
        self.superlattice_data_preview = False
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
        self.seed_var = tk.StringVar(value=str(secrets.randbits(32)))
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
        self.kin_bkg_a_min_var = tk.StringVar(value="-3.0")
        self.kin_bkg_a_max_var = tk.StringVar(value="3.0")
        self.kin_bkg_b_start_var = tk.StringVar(value="0.1")
        self.kin_bkg_b_min_var = tk.StringVar(value="0.0")
        self.kin_bkg_b_max_var = tk.StringVar(value="3.0")
        self.kin_bkg_c_start_var = tk.StringVar(value="0.0")
        self.kin_bkg_c_min_var = tk.StringVar(value="-3.0")
        self.kin_bkg_c_max_var = tk.StringVar(value="3.0")
        self.kin_debye_var = tk.StringVar(value=f"{SAMPLES['Fe 10 nm'].debye_waller_coeff:.6g}")
        self.kin_d_fit_var = tk.StringVar(value="")
        self.kin_planes_fit_var = tk.StringVar(value="")
        self.kin_resolution_fit_var = tk.StringVar(value="")
        self.kin_scale_fit_var = tk.StringVar(value="")
        self.kin_bkg_a_fit_var = tk.StringVar(value="")
        self.kin_bkg_b_fit_var = tk.StringVar(value="")
        self.kin_bkg_c_fit_var = tk.StringVar(value="")
        self.kin_debye_fit_var = tk.StringVar(value="")
        self.kin_d_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_planes_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_resolution_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_scale_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_bkg_a_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_bkg_b_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_bkg_c_fit_enabled_var = tk.BooleanVar(value=False)
        self.kin_substrate_var = tk.BooleanVar(value=False)
        self.kin_substrate_intensity_start_var = tk.StringVar(value="51.0")
        self.kin_substrate_intensity_min_var = tk.StringVar(value="0.0")
        self.kin_substrate_intensity_max_var = tk.StringVar(value="1000000.0")
        self.kin_substrate_width_start_var = tk.StringVar(value="0.004")
        self.kin_substrate_width_min_var = tk.StringVar(value="0.0001")
        self.kin_substrate_width_max_var = tk.StringVar(value="0.05")
        self.kin_substrate_d_start_var = tk.StringVar(value="")
        self.kin_substrate_d_min_var = tk.StringVar(value="")
        self.kin_substrate_d_max_var = tk.StringVar(value="")
        self.kin_substrate_intensity_fit_var = tk.StringVar(value="")
        self.kin_substrate_width_fit_var = tk.StringVar(value="")
        self.kin_substrate_d_fit_var = tk.StringVar(value="")
        self.kin_substrate_intensity_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_substrate_width_fit_enabled_var = tk.BooleanVar(value=True)
        self.kin_substrate_d_fit_enabled_var = tk.BooleanVar(value=True)
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
        self.dynamic_bkg_a_min_var = tk.StringVar(value="-0.1")
        self.dynamic_bkg_a_max_var = tk.StringVar(value="0.1")
        self.dynamic_bkg_a_fit_var = tk.StringVar(value="")
        self.dynamic_bkg_b_start_var = tk.StringVar(value="1.2366e-7")
        self.dynamic_bkg_b_min_var = tk.StringVar(value="0.0")
        self.dynamic_bkg_b_max_var = tk.StringVar(value="0.1")
        self.dynamic_bkg_b_fit_var = tk.StringVar(value="")
        self.dynamic_bkg_c_start_var = tk.StringVar(value="0.0")
        self.dynamic_bkg_c_min_var = tk.StringVar(value="-0.1")
        self.dynamic_bkg_c_max_var = tk.StringVar(value="0.1")
        self.dynamic_bkg_c_fit_var = tk.StringVar(value="")
        self.density_slices_var = tk.StringVar(value="100")
        self.density_max_q0_var = tk.StringVar(value="30.0")
        self.dynamic_resolution_fit_enabled_var = tk.BooleanVar(value=True)
        self.dynamic_intensity_fit_enabled_var = tk.BooleanVar(value=True)
        self.dynamic_bkg_a_fit_enabled_var = tk.BooleanVar(value=True)
        self.dynamic_bkg_b_fit_enabled_var = tk.BooleanVar(value=True)
        self.dynamic_bkg_c_fit_enabled_var = tk.BooleanVar(value=False)
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
        self.kinematic_substrate_widgets: list[tk.Widget] = []
        self.dynamic_widgets: list[tk.Widget] = []
        self.strain_widgets: list[tk.Widget] = []
        self.range_indicators: list[RangeIndicator] = []
        default_stack_path = STACK_DIR / "fe_v_4_28_x11.json"
        self.stack_enabled_var = tk.BooleanVar(value=False)
        self.stack_path_var = tk.StringVar(value=str(default_stack_path))
        self.stack_name_var = tk.StringVar()
        self.stack_repetitions_var = tk.StringVar()
        self.stack_capping_enabled_var = tk.BooleanVar(value=False)
        self.stack_points_per_unit_var = tk.StringVar(value="50")
        self.stack_sampling_label_var = tk.StringVar(value="Points per degree")
        self.stack_calculation_fit_enabled_vars = {
            key: tk.BooleanVar(value=False)
            for key in (
                "resolution",
                "intensity_scale",
                "background_a",
                "background_b",
                "background_c",
            )
        }
        self.stack_calculation_fit_vars = {
            key: tk.StringVar(value="")
            for key in self.stack_calculation_fit_enabled_vars
        }
        self.stack_document = copy.deepcopy(StackDefinition.load(default_stack_path).document)
        self.stack_row_vars: list[dict[str, tk.Variable]] = []
        self.stack_row_specs: list[dict[str, object]] = []
        self.stack_row_roles: list[str] = []
        self.stack_structure_indicators: list[RangeIndicator] = []

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
        self.kin_substrate_var.trace_add("write", self._on_kinematic_substrate_changed)
        self.stack_enabled_var.trace_add("write", self._on_stack_enabled_changed)
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
        style.configure("TRadiobutton", background=UI_COLORS["panel"])
        style.configure("TEntry", fieldbackground=UI_COLORS["entry"])
        style.map(
            "TEntry",
            fieldbackground=[("disabled", UI_COLORS["entry_disabled"])],
            foreground=[("disabled", "#6b7280")],
        )
        for name, color in (
            ("Simulate", UI_COLORS["simulate"]),
            ("Run", UI_COLORS["run"]),
            ("Pause", UI_COLORS["fit"]),
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
        elif "warning" in text or "paused" in text or "pause requested" in text:
            style = "StatusWarning.TLabel"
        elif any(marker in text for marker in ("running", "stopping", "cost=")):
            style = "StatusActive.TLabel"
        elif "complete" in text:
            style = "StatusSuccess.TLabel"
        else:
            style = "Status.TLabel"
        self.status_label.configure(style=style)

    def _build_layout(self) -> None:
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        self.root.geometry(f"{screen_width}x{screen_height}+0+0")
        self.root.minsize(1180, 720)

        main_pane = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_pane.pack(fill=tk.BOTH, expand=True)

        sidebar = ttk.Frame(main_pane, padding=10)
        main_pane.add(sidebar, weight=0)

        run_frame = ttk.LabelFrame(sidebar, text="Simulation and fit setup")
        run_frame.pack(fill=tk.X)
        self._add_accent_strip(run_frame, UI_COLORS["setup"])

        add_tooltip(ttk.Label(run_frame, text="Scattering model"), "scattering_model").grid(
            row=0, column=0, sticky="w"
        )
        self.model_combo = ttk.Combobox(
            run_frame,
            textvariable=self.model_var,
            values=("Dynamic", "Kinematic"),
            state="readonly",
            width=14,
        )
        self.model_combo.grid(row=0, column=1, sticky="ew", padx=(4, 8), pady=2)
        add_tooltip(self.model_combo, "scattering_model")

        add_tooltip(
            ttk.Label(run_frame, text="Sample configuration"), "sample_configuration"
        ).grid(row=0, column=2, sticky="w")
        sample_configuration = ttk.Frame(run_frame, style="Panel.TFrame")
        sample_configuration.grid(row=0, column=3, sticky="ew", pady=2)
        for text, value in (("Film", False), ("Superlattice", True)):
            add_tooltip(
                ttk.Radiobutton(
                    sample_configuration,
                    text=text,
                    variable=self.stack_enabled_var,
                    value=value,
                ),
                "sample_configuration",
            ).pack(side=tk.LEFT, padx=(8, 0) if value else 0)

        add_tooltip(ttk.Label(run_frame, text="X-ray wavelength (\u00c5)"), "wavelength").grid(
            row=1, column=0, sticky="w"
        )
        add_tooltip(
            ttk.Entry(run_frame, textvariable=self.wavelength_var, width=10), "wavelength"
        ).grid(
            row=1, column=1, sticky="ew", padx=(4, 8), pady=2
        )
        add_tooltip(ttk.Label(run_frame, text="Horizontal axis"), "horizontal_axis").grid(
            row=1, column=2, sticky="w"
        )
        add_tooltip(
            ttk.Combobox(
                run_frame,
                textvariable=self.axis_var,
                values=("2\u03b8", "q"),
                state="readonly",
                width=10,
            ),
            "horizontal_axis",
        ).grid(row=1, column=3, sticky="ew", padx=(4, 8), pady=2)

        run_controls = [
            (self.min_label_var, self.min_var, 2, 0),
            (self.max_label_var, self.max_var, 2, 2),
        ]
        for label, var, row, column in run_controls:
            if isinstance(label, tk.StringVar):
                ttk.Label(run_frame, textvariable=label).grid(row=row, column=column, sticky="w")
            else:
                ttk.Label(run_frame, text=label).grid(row=row, column=column, sticky="w")
            ttk.Entry(run_frame, textvariable=var, width=10).grid(
                row=row, column=column + 1, sticky="ew", padx=(4, 8), pady=2
            )

        self.simulate_button = ttk.Button(
            run_frame,
            text="Simulate",
            command=self.simulate_pattern,
            style="Simulate.TButton",
        )
        self.simulate_button.grid(row=3, column=0, sticky="ew", padx=(0, 4), pady=(6, 2))
        self.run_button = ttk.Button(
            run_frame,
            text="Run fit",
            command=self.run_fit,
            style="Run.TButton",
        )
        self.run_button.grid(row=3, column=1, sticky="ew", padx=(0, 4), pady=(6, 2))
        self.pause_button = ttk.Button(
            run_frame,
            text="Pause fit",
            command=self.pause_or_resume_fit,
            state=tk.DISABLED,
            style="Pause.TButton",
        )
        self.pause_button.grid(row=3, column=2, sticky="ew", padx=(0, 4), pady=(6, 2))
        self.stop_button = ttk.Button(
            run_frame,
            text="Stop fit",
            command=self.stop_fit,
            state=tk.DISABLED,
            style="Stop.TButton",
        )
        self.stop_button.grid(row=3, column=3, sticky="ew", padx=(4, 0), pady=(6, 2))

        ttk.Button(run_frame, text="Save setup/results...", command=self._save_project).grid(
            row=4, column=0, columnspan=2, sticky="ew", padx=(0, 4), pady=2
        )
        ttk.Button(run_frame, text="Load setup/results...", command=self._load_project).grid(
            row=4, column=2, columnspan=2, sticky="ew", padx=(4, 0), pady=2
        )
        ttk.Button(run_frame, text="Save plots...", command=self._save_plot_image).grid(
            row=5, column=0, columnspan=2, sticky="ew", padx=(0, 4), pady=2
        )
        ttk.Button(run_frame, text="Export graph data...", command=self._export_graph_data).grid(
            row=5, column=2, columnspan=2, sticky="ew", padx=(4, 0), pady=2
        )

        self.status_label = ttk.Label(
            run_frame,
            textvariable=self.status_var,
            wraplength=390,
            style="Status.TLabel",
        )
        self.status_label.grid(row=6, column=0, columnspan=4, sticky="ew", pady=(4, 2))
        for column in (1, 3):
            run_frame.columnconfigure(column, weight=1)

        self.workspace_frame = ttk.Frame(sidebar)
        self.workspace_frame.pack(fill=tk.BOTH, expand=True, pady=(8, 8))
        self.workspace_frame.rowconfigure(0, weight=1)
        self.workspace_frame.columnconfigure(0, weight=1)

        self.film_container = ttk.LabelFrame(
            self.workspace_frame, text="Film simulation and fitting"
        )
        self.film_container.grid(row=0, column=0, sticky="nsew")
        film_setup_frame = ttk.Frame(
            self.film_container, padding=(8, 6), style="Panel.TFrame"
        )
        film_setup_frame.pack(fill=tk.X, padx=4, pady=(4, 0))
        ttk.Label(film_setup_frame, text="Film data").grid(row=0, column=0, sticky="w")
        ttk.Entry(film_setup_frame, textvariable=self.data_path_var, width=42).grid(
            row=0, column=1, columnspan=2, sticky="ew", padx=(4, 8), pady=2
        )
        ttk.Button(film_setup_frame, text="Browse...", command=self._browse_data_file).grid(
            row=0, column=3, sticky="ew", pady=2
        )
        self.strain_checkbutton = ttk.Checkbutton(
            film_setup_frame, text="Include strain", variable=self.strain_var
        )
        self.strain_checkbutton.grid(
            row=1, column=0, columnspan=2, sticky="w", pady=(4, 0)
        )
        add_tooltip(self.strain_checkbutton, "include_strain")
        self.roughness_checkbutton = ttk.Checkbutton(
            film_setup_frame, text="Include roughness", variable=self.roughness_var
        )
        self.roughness_checkbutton.grid(
            row=1, column=2, columnspan=2, sticky="w", pady=(4, 0)
        )
        add_tooltip(self.roughness_checkbutton, "include_roughness")
        film_setup_frame.columnconfigure(1, weight=1)
        parameter_tabs = self.film_parameter_tabs = ttk.Notebook(self.film_container)
        parameter_tabs.pack(fill=tk.BOTH, expand=True, padx=4, pady=(4, 4))

        kinematic_panel = ttk.Frame(parameter_tabs, style="Panel.TFrame")
        parameter_tabs.add(kinematic_panel, text="Kinematic")
        kinematic_tabs = ttk.Notebook(kinematic_panel)
        kinematic_tabs.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.kinematic_frame = ttk.Frame(kinematic_tabs, padding=8, style="Panel.TFrame")
        kinematic_tabs.add(self.kinematic_frame, text="Film and fit")
        self._add_accent_strip(self.kinematic_frame, UI_COLORS["kinematic"])
        ttk.Label(self.kinematic_frame, text="Fit").grid(row=0, column=1, sticky="ew")
        ttk.Label(self.kinematic_frame, text="Value").grid(row=0, column=2, sticky="ew")
        ttk.Label(self.kinematic_frame, text="Min").grid(row=0, column=3, sticky="ew")
        ttk.Label(self.kinematic_frame, text="Max").grid(row=0, column=4, sticky="ew")
        ttk.Label(self.kinematic_frame, text="Range").grid(row=0, column=5, sticky="ew")
        kinematic_controls = [
            ("Plane spacing d (\u00c5)", self.kin_d_fit_enabled_var, self.kin_d_start_var, self.kin_d_min_var, self.kin_d_max_var, self.kin_d_fit_var, "plane_spacing"),
            ("Number of planes", self.kin_planes_fit_enabled_var, self.kin_planes_start_var, self.kin_planes_min_var, self.kin_planes_max_var, self.kin_planes_fit_var, "number_of_planes"),
            (
                "Resolution (deg)",
                self.kin_resolution_fit_enabled_var,
                self.kin_resolution_start_var,
                self.kin_resolution_min_var,
                self.kin_resolution_max_var,
                self.kin_resolution_fit_var,
                "resolution",
            ),
            ("Intensity scale", self.kin_scale_fit_enabled_var, self.kin_scale_start_var, self.kin_scale_min_var, self.kin_scale_max_var, self.kin_scale_fit_var, "intensity_scale"),
            ("Background slope (a)", self.kin_bkg_a_fit_enabled_var, self.kin_bkg_a_start_var, self.kin_bkg_a_min_var, self.kin_bkg_a_max_var, self.kin_bkg_a_fit_var, "background_a"),
            ("Background offset (b)", self.kin_bkg_b_fit_enabled_var, self.kin_bkg_b_start_var, self.kin_bkg_b_min_var, self.kin_bkg_b_max_var, self.kin_bkg_b_fit_var, "background_b"),
            ("Background curvature (c)", self.kin_bkg_c_fit_enabled_var, self.kin_bkg_c_start_var, self.kin_bkg_c_min_var, self.kin_bkg_c_max_var, self.kin_bkg_c_fit_var, "background_c"),
        ]
        for row, (label, fit_enabled_var, start_var, min_var, max_var, fit_var, tooltip_key) in enumerate(kinematic_controls, start=1):
            add_tooltip(ttk.Label(self.kinematic_frame, text=label), tooltip_key).grid(row=row, column=0, sticky="w")
            add_tooltip(ttk.Checkbutton(self.kinematic_frame, variable=fit_enabled_var), "fit").grid(row=row, column=1, sticky="ew")
            ttk.Entry(self.kinematic_frame, textvariable=start_var, width=7).grid(row=row, column=2, sticky="ew")
            ttk.Entry(self.kinematic_frame, textvariable=min_var, width=7).grid(row=row, column=3, sticky="ew")
            ttk.Entry(self.kinematic_frame, textvariable=max_var, width=7).grid(row=row, column=4, sticky="ew")
            self._add_range_indicator(self.kinematic_frame, row, label, start_var, min_var, max_var, fit_var, fit_enabled_var)
        add_tooltip(
            ttk.Label(self.kinematic_frame, text="Debye-Waller coefficient"), "debye_waller"
        ).grid(row=7, column=0, sticky="w")
        ttk.Entry(self.kinematic_frame, textvariable=self.kin_debye_var, width=7).grid(
            row=7, column=2, columnspan=3, sticky="ew"
        )
        self.kinematic_frame.columnconfigure(2, weight=1)
        self.kinematic_frame.columnconfigure(3, weight=1)
        self.kinematic_frame.columnconfigure(4, weight=1)
        self.kinematic_frame.columnconfigure(5, weight=1)

        self.kinematic_substrate_frame = ttk.Frame(
            kinematic_tabs, padding=8, style="Panel.TFrame"
        )
        kinematic_tabs.add(self.kinematic_substrate_frame, text="Substrate")
        self._add_accent_strip(self.kinematic_substrate_frame, UI_COLORS["substrate"])
        add_tooltip(
            ttk.Checkbutton(
                self.kinematic_substrate_frame,
                text="Include substrate peak",
                variable=self.kin_substrate_var,
            ),
            "include_substrate_peak",
        ).grid(row=0, column=0, columnspan=5, sticky="w", pady=(0, 6))
        self.kinematic_substrate_controls_frame = ttk.Frame(
            self.kinematic_substrate_frame, style="Panel.TFrame"
        )
        self.kinematic_substrate_controls_frame.grid(
            row=1, column=0, columnspan=6, sticky="nsew"
        )
        ttk.Label(self.kinematic_substrate_controls_frame, text="Fit").grid(
            row=0, column=1, sticky="ew"
        )
        ttk.Label(self.kinematic_substrate_controls_frame, text="Value").grid(
            row=0, column=2, sticky="ew"
        )
        ttk.Label(self.kinematic_substrate_controls_frame, text="Min").grid(
            row=0, column=3, sticky="ew"
        )
        ttk.Label(self.kinematic_substrate_controls_frame, text="Max").grid(
            row=0, column=4, sticky="ew"
        )
        ttk.Label(self.kinematic_substrate_controls_frame, text="Range").grid(
            row=0, column=5, sticky="ew"
        )
        substrate_peak_controls = [
            (
                "Integrated intensity",
                self.kin_substrate_intensity_fit_enabled_var,
                self.kin_substrate_intensity_start_var,
                self.kin_substrate_intensity_min_var,
                self.kin_substrate_intensity_max_var,
                self.kin_substrate_intensity_fit_var,
                "substrate_intensity",
            ),
            (
                "FWHM (deg)",
                self.kin_substrate_width_fit_enabled_var,
                self.kin_substrate_width_start_var,
                self.kin_substrate_width_min_var,
                self.kin_substrate_width_max_var,
                self.kin_substrate_width_fit_var,
                "substrate_fwhm",
            ),
            (
                "d spacing (Å)",
                self.kin_substrate_d_fit_enabled_var,
                self.kin_substrate_d_start_var,
                self.kin_substrate_d_min_var,
                self.kin_substrate_d_max_var,
                self.kin_substrate_d_fit_var,
                "substrate_d_spacing",
            ),
        ]
        for row, (label, fit_enabled_var, start_var, min_var, max_var, fit_var, tooltip_key) in enumerate(
            substrate_peak_controls, start=1
        ):
            add_tooltip(ttk.Label(self.kinematic_substrate_controls_frame, text=label), tooltip_key).grid(
                row=row, column=0, sticky="w"
            )
            add_tooltip(
                ttk.Checkbutton(self.kinematic_substrate_controls_frame, variable=fit_enabled_var),
                "fit",
            ).grid(row=row, column=1, sticky="ew")
            ttk.Entry(
                self.kinematic_substrate_controls_frame, textvariable=start_var, width=7
            ).grid(row=row, column=2, sticky="ew")
            ttk.Entry(
                self.kinematic_substrate_controls_frame, textvariable=min_var, width=7
            ).grid(row=row, column=3, sticky="ew")
            ttk.Entry(
                self.kinematic_substrate_controls_frame, textvariable=max_var, width=7
            ).grid(row=row, column=4, sticky="ew")
            self._add_range_indicator(
                self.kinematic_substrate_controls_frame,
                row,
                label,
                start_var,
                min_var,
                max_var,
                fit_var,
                fit_enabled_var,
            )
        for column in range(2, 6):
            self.kinematic_substrate_controls_frame.columnconfigure(column, weight=1)
        self.kinematic_substrate_frame.columnconfigure(0, weight=1)

        dynamic_panel = ttk.Frame(parameter_tabs, style="Panel.TFrame")
        parameter_tabs.add(dynamic_panel, text="Dynamic")
        dynamic_tabs = ttk.Notebook(dynamic_panel)
        dynamic_tabs.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.dynamic_film_frame = ttk.Frame(dynamic_tabs, padding=8, style="Panel.TFrame")
        dynamic_tabs.add(self.dynamic_film_frame, text="Film")
        self._add_accent_strip(self.dynamic_film_frame, UI_COLORS["film"])
        add_tooltip(
            ttk.Label(self.dynamic_film_frame, text="Structure file"), "structure_file"
        ).grid(row=0, column=0, sticky="w")
        add_tooltip(
            ttk.Combobox(
                self.dynamic_film_frame,
                textvariable=self.film_filename_var,
                values=POSCAR_FILES,
                state="readonly",
                width=18,
            ),
            "structure_file",
        ).grid(row=0, column=2, columnspan=3, sticky="ew")
        add_tooltip(
            ttk.Label(self.dynamic_film_frame, text="Layer direction"), "layer_direction"
        ).grid(row=1, column=0, sticky="w")
        add_tooltip(
            ttk.Combobox(
                self.dynamic_film_frame,
                textvariable=self.film_direction_var,
                values=("1", "2", "3"),
                state="readonly",
                width=6,
            ),
            "layer_direction",
        ).grid(row=1, column=2, sticky="w")
        ttk.Label(self.dynamic_film_frame, text="Fit").grid(row=2, column=1, sticky="ew")
        ttk.Label(self.dynamic_film_frame, text="Value").grid(row=2, column=2, sticky="ew")
        ttk.Label(self.dynamic_film_frame, text="Min").grid(row=2, column=3, sticky="ew")
        ttk.Label(self.dynamic_film_frame, text="Max").grid(row=2, column=4, sticky="ew")
        ttk.Label(self.dynamic_film_frame, text="Range").grid(row=2, column=5, sticky="ew")
        film_controls = [
            ("Number of layers", self.film_n_fit_enabled_var, self.film_n_start_var, self.film_n_min_var, self.film_n_max_var, self.film_n_fit_var, "number_of_layers"),
            ("Lattice scale", self.film_scale_fit_enabled_var, self.film_scale_start_var, self.film_scale_min_var, self.film_scale_max_var, self.film_scale_fit_var, "lattice_scale"),
            ("Area scale", self.film_area_fit_enabled_var, self.film_area_start_var, self.film_area_min_var, self.film_area_max_var, self.film_area_fit_var, "area_scale"),
            (
                "Interface spacing (\u00c5)",
                self.film_interface_fit_enabled_var,
                self.film_interface_start_var,
                self.film_interface_min_var,
                self.film_interface_max_var,
                self.film_interface_fit_var,
                "interface_spacing",
            ),
        ]
        for row, (label, fit_enabled_var, start_var, min_var, max_var, fit_var, tooltip_key) in enumerate(film_controls, start=3):
            add_tooltip(ttk.Label(self.dynamic_film_frame, text=label), tooltip_key).grid(row=row, column=0, sticky="w")
            add_tooltip(ttk.Checkbutton(self.dynamic_film_frame, variable=fit_enabled_var), "fit").grid(row=row, column=1, sticky="ew")
            ttk.Entry(self.dynamic_film_frame, textvariable=start_var, width=7).grid(row=row, column=2, sticky="ew")
            ttk.Entry(self.dynamic_film_frame, textvariable=min_var, width=7).grid(row=row, column=3, sticky="ew")
            ttk.Entry(self.dynamic_film_frame, textvariable=max_var, width=7).grid(row=row, column=4, sticky="ew")
            self._add_range_indicator(self.dynamic_film_frame, row, label, start_var, min_var, max_var, fit_var, fit_enabled_var)
        self.dynamic_film_frame.columnconfigure(2, weight=1)
        self.dynamic_film_frame.columnconfigure(3, weight=1)
        self.dynamic_film_frame.columnconfigure(4, weight=1)
        self.dynamic_film_frame.columnconfigure(5, weight=1)

        self.dynamic_fit_frame = ttk.Frame(dynamic_tabs, padding=8, style="Panel.TFrame")
        dynamic_tabs.add(self.dynamic_fit_frame, text="Calculation and fit")
        self._add_accent_strip(self.dynamic_fit_frame, UI_COLORS["fit"])
        add_tooltip(
            ttk.Label(self.dynamic_fit_frame, text="Density slices per cell"), "density_slices"
        ).grid(
            row=0, column=0, sticky="w"
        )
        add_tooltip(
            ttk.Entry(self.dynamic_fit_frame, textvariable=self.density_slices_var, width=10),
            "density_slices",
        ).grid(row=0, column=2, sticky="w")
        add_tooltip(
            ttk.Label(self.dynamic_fit_frame, text="Density Q max (1/Å)"), "density_q_max"
        ).grid(
            row=1, column=0, sticky="w"
        )
        add_tooltip(
            ttk.Entry(self.dynamic_fit_frame, textvariable=self.density_max_q0_var, width=10),
            "density_q_max",
        ).grid(row=1, column=2, sticky="w")
        ttk.Label(self.dynamic_fit_frame, text="Fit").grid(row=2, column=1, sticky="ew")
        ttk.Label(self.dynamic_fit_frame, text="Value").grid(row=2, column=2, sticky="ew")
        ttk.Label(self.dynamic_fit_frame, text="Min").grid(row=2, column=3, sticky="ew")
        ttk.Label(self.dynamic_fit_frame, text="Max").grid(row=2, column=4, sticky="ew")
        ttk.Label(self.dynamic_fit_frame, text="Range").grid(row=2, column=5, sticky="ew")
        dynamic_fit_controls = [
            (
                "Resolution (deg)",
                self.dynamic_resolution_fit_enabled_var,
                self.dynamic_resolution_start_var,
                self.dynamic_resolution_min_var,
                self.dynamic_resolution_max_var,
                self.dynamic_resolution_fit_var,
                "resolution",
            ),
            (
                "Intensity scale",
                self.dynamic_intensity_fit_enabled_var,
                self.dynamic_intensity_start_var,
                self.dynamic_intensity_min_var,
                self.dynamic_intensity_max_var,
                self.dynamic_intensity_fit_var,
                "intensity_scale",
            ),
            (
                "Background slope (a)",
                self.dynamic_bkg_a_fit_enabled_var,
                self.dynamic_bkg_a_start_var,
                self.dynamic_bkg_a_min_var,
                self.dynamic_bkg_a_max_var,
                self.dynamic_bkg_a_fit_var,
                "background_a",
            ),
            (
                "Background offset (b)",
                self.dynamic_bkg_b_fit_enabled_var,
                self.dynamic_bkg_b_start_var,
                self.dynamic_bkg_b_min_var,
                self.dynamic_bkg_b_max_var,
                self.dynamic_bkg_b_fit_var,
                "background_b",
            ),
            (
                "Background curvature (c)",
                self.dynamic_bkg_c_fit_enabled_var,
                self.dynamic_bkg_c_start_var,
                self.dynamic_bkg_c_min_var,
                self.dynamic_bkg_c_max_var,
                self.dynamic_bkg_c_fit_var,
                "background_c",
            ),
        ]
        for row, (label, fit_enabled_var, start_var, min_var, max_var, fit_var, tooltip_key) in enumerate(dynamic_fit_controls, start=3):
            add_tooltip(ttk.Label(self.dynamic_fit_frame, text=label), tooltip_key).grid(row=row, column=0, sticky="w")
            add_tooltip(ttk.Checkbutton(self.dynamic_fit_frame, variable=fit_enabled_var), "fit").grid(row=row, column=1, sticky="ew")
            ttk.Entry(self.dynamic_fit_frame, textvariable=start_var, width=7).grid(row=row, column=2, sticky="ew")
            ttk.Entry(self.dynamic_fit_frame, textvariable=min_var, width=7).grid(row=row, column=3, sticky="ew")
            ttk.Entry(self.dynamic_fit_frame, textvariable=max_var, width=7).grid(row=row, column=4, sticky="ew")
            self._add_range_indicator(self.dynamic_fit_frame, row, label, start_var, min_var, max_var, fit_var, fit_enabled_var)
        self.dynamic_fit_frame.columnconfigure(2, weight=1)
        self.dynamic_fit_frame.columnconfigure(3, weight=1)
        self.dynamic_fit_frame.columnconfigure(4, weight=1)
        self.dynamic_fit_frame.columnconfigure(5, weight=1)

        self.dynamic_substrate_frame = ttk.Frame(dynamic_tabs, padding=8, style="Panel.TFrame")
        dynamic_tabs.add(self.dynamic_substrate_frame, text="Substrate")
        self._add_accent_strip(self.dynamic_substrate_frame, UI_COLORS["substrate"])
        add_tooltip(
            ttk.Label(self.dynamic_substrate_frame, text="Structure file"), "structure_file"
        ).grid(row=0, column=0, sticky="w")
        add_tooltip(
            ttk.Combobox(
                self.dynamic_substrate_frame,
                textvariable=self.substrate_filename_var,
                values=POSCAR_FILES,
                state="readonly",
                width=18,
            ),
            "structure_file",
        ).grid(row=0, column=1, sticky="ew")
        substrate_controls = [
            ("Layer direction", self.substrate_direction_var, ("1", "2", "3"), "layer_direction"),
            ("Number of layers", self.substrate_n_var, None, "substrate_layers"),
            ("Interface spacing (\u00c5)", self.substrate_interface_var, None, "interface_spacing"),
            ("Area scale", self.substrate_area_var, None, "area_scale"),
        ]
        for row, (label, var, values, tooltip_key) in enumerate(substrate_controls, start=1):
            add_tooltip(ttk.Label(self.dynamic_substrate_frame, text=label), tooltip_key).grid(row=row, column=0, sticky="w")
            if values is None:
                add_tooltip(
                    ttk.Entry(self.dynamic_substrate_frame, textvariable=var, width=12),
                    tooltip_key,
                ).grid(row=row, column=1, sticky="ew")
            else:
                add_tooltip(
                    ttk.Combobox(
                        self.dynamic_substrate_frame,
                        textvariable=var,
                        values=values,
                        state="readonly",
                        width=8,
                    ),
                    tooltip_key,
                ).grid(row=row, column=1, sticky="w")
        ttk.Label(self.dynamic_substrate_frame, text="Fit").grid(row=5, column=1, sticky="ew")
        ttk.Label(self.dynamic_substrate_frame, text="Value").grid(row=5, column=2, sticky="ew")
        ttk.Label(self.dynamic_substrate_frame, text="Min").grid(row=5, column=3, sticky="ew")
        ttk.Label(self.dynamic_substrate_frame, text="Max").grid(row=5, column=4, sticky="ew")
        ttk.Label(self.dynamic_substrate_frame, text="Range").grid(row=5, column=5, sticky="ew")
        add_tooltip(
            ttk.Label(self.dynamic_substrate_frame, text="Lattice scale"), "lattice_scale"
        ).grid(row=6, column=0, sticky="w")
        add_tooltip(
            ttk.Checkbutton(
                self.dynamic_substrate_frame,
                variable=self.substrate_scale_fit_enabled_var,
            ),
            "fit",
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

        optional_panel = ttk.Frame(parameter_tabs, style="Panel.TFrame")
        parameter_tabs.add(optional_panel, text="Strain / roughness")
        optional_tabs = ttk.Notebook(optional_panel)
        optional_tabs.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.strain_frame = ttk.Frame(optional_tabs, padding=8, style="Panel.TFrame")
        optional_tabs.add(self.strain_frame, text="Strain")
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
        strain_tooltips = (
            "bottom_strain_amplitude",
            "bottom_strain_extent",
            "top_strain_amplitude",
            "top_strain_extent",
        )
        for row, (label, fit_enabled_var, start_var, min_var, max_var, fit_var, tooltip_key) in enumerate(
            zip(
                strain_controls,
                self.strain_fit_enabled_vars,
                self.strain_start_vars,
                self.strain_min_vars,
                self.strain_max_vars,
                self.strain_fit_vars,
                strain_tooltips,
            ),
            start=1,
        ):
            name = label.get() if isinstance(label, tk.StringVar) else label
            if isinstance(label, tk.StringVar):
                label_widget = ttk.Label(self.strain_frame, textvariable=label)
            else:
                label_widget = ttk.Label(self.strain_frame, text=label)
            add_tooltip(label_widget, tooltip_key).grid(row=row, column=0, sticky="w")
            add_tooltip(ttk.Checkbutton(self.strain_frame, variable=fit_enabled_var), "fit").grid(row=row, column=1, sticky="ew")
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

        rough_frame = ttk.Frame(optional_tabs, padding=8, style="Panel.TFrame")
        optional_tabs.add(rough_frame, text="Roughness")
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
                "film_roughness",
            ),
            (
                "Substrate/interface roughness \u03c3 (\u00c5)",
                self.substrate_rough_fit_enabled_var,
                self.substrate_rough_start_var,
                self.substrate_rough_min_var,
                self.substrate_rough_max_var,
                self.substrate_rough_fit_var,
                "substrate_roughness",
            ),
        ]
        for row, (label, fit_enabled_var, start_var, min_var, max_var, fit_var, tooltip_key) in enumerate(rough_controls, start=1):
            add_tooltip(ttk.Label(rough_frame, text=label), tooltip_key).grid(row=row, column=0, sticky="w")
            add_tooltip(ttk.Checkbutton(rough_frame, variable=fit_enabled_var), "fit").grid(row=row, column=1, sticky="ew")
            ttk.Entry(rough_frame, textvariable=start_var, width=7).grid(row=row, column=2, sticky="ew")
            ttk.Entry(rough_frame, textvariable=min_var, width=7).grid(row=row, column=3, sticky="ew")
            ttk.Entry(rough_frame, textvariable=max_var, width=7).grid(row=row, column=4, sticky="ew")
            self._add_range_indicator(rough_frame, row, label, start_var, min_var, max_var, fit_var, fit_enabled_var)
        rough_frame.columnconfigure(2, weight=1)
        rough_frame.columnconfigure(3, weight=1)
        rough_frame.columnconfigure(4, weight=1)
        rough_frame.columnconfigure(5, weight=1)

        optimization_frame = ttk.Frame(parameter_tabs, padding=8, style="Panel.TFrame")
        parameter_tabs.add(optimization_frame, text="Optimization settings")
        self._add_accent_strip(optimization_frame, UI_COLORS["fit"])
        optimization_controls = [
            ("Seed", self.seed_var, 0, 0, "seed"),
            ("Progress update interval", self.interval_var, 0, 2, "progress_interval"),
            ("DE max iterations", self.maxiter_var, 1, 0, "de_iterations"),
            ("DE population size", self.popsize_var, 1, 2, "de_population"),
            ("Local max evaluations", self.local_var, 2, 0, "local_evaluations"),
            ("Polish iterations", self.polish_var, 2, 2, "polish_iterations"),
        ]

        def add_optimization_control(
            parent: ttk.Frame,
            label: str,
            variable: tk.StringVar,
            row: int,
            column: int,
            tooltip_key: str,
        ) -> None:
            add_tooltip(ttk.Label(parent, text=label), tooltip_key).grid(
                row=row, column=column, sticky="w"
            )
            if variable is self.seed_var:
                editor = ttk.Frame(parent, style="Panel.TFrame")
                editor.grid(
                    row=row,
                    column=column + 1,
                    sticky="ew",
                    padx=(4, 8),
                    pady=2,
                )
                add_tooltip(
                    ttk.Entry(editor, textvariable=variable, width=10), tooltip_key
                ).pack(
                    side=tk.LEFT, fill=tk.X, expand=True
                )
                add_tooltip(
                    ttk.Button(
                        editor,
                        text="New seed",
                        command=self._generate_optimizer_seed,
                    ),
                    "new_seed",
                ).pack(side=tk.LEFT, padx=(4, 0))
            else:
                add_tooltip(
                    ttk.Entry(parent, textvariable=variable, width=10), tooltip_key
                ).grid(
                    row=row,
                    column=column + 1,
                    sticky="ew",
                    padx=(4, 8),
                    pady=2,
                )

        for label, var, row, column, tooltip_key in optimization_controls:
            add_optimization_control(
                optimization_frame, label, var, row, column, tooltip_key
            )
        optimization_frame.columnconfigure(1, weight=1)
        optimization_frame.columnconfigure(3, weight=1)

        self.stack_container = ttk.LabelFrame(
            self.workspace_frame, text="Superlattice simulation and fitting"
        )
        self.stack_container.grid(row=0, column=0, sticky="nsew")
        self._add_accent_strip(self.stack_container, UI_COLORS["stack"])
        stack_tabs = ttk.Notebook(self.stack_container)
        stack_tabs.pack(fill=tk.BOTH, expand=True, padx=4, pady=(4, 4))

        stack_structure_panel = ttk.Frame(stack_tabs, padding=8, style="Panel.TFrame")
        stack_tabs.add(stack_structure_panel, text="Structure")
        ttk.Label(stack_structure_panel, text="Superlattice file").grid(
            row=0, column=0, sticky="w"
        )
        ttk.Entry(stack_structure_panel, textvariable=self.stack_path_var).grid(
            row=0, column=1, columnspan=3, sticky="ew", padx=(4, 0)
        )
        ttk.Button(stack_structure_panel, text="Browse...", command=self._browse_stack_file).grid(
            row=1, column=1, sticky="ew", padx=(4, 2), pady=3
        )
        ttk.Button(stack_structure_panel, text="Reload", command=self._load_stack_file).grid(
            row=1, column=2, sticky="ew", padx=2, pady=3
        )
        ttk.Button(stack_structure_panel, text="Save as...", command=self._save_stack_file).grid(
            row=1, column=3, sticky="ew", padx=(2, 0), pady=3
        )
        ttk.Label(stack_structure_panel, text="Superlattice data").grid(
            row=2, column=0, sticky="w"
        )
        ttk.Entry(stack_structure_panel, textvariable=self.data_path_var).grid(
            row=2, column=1, columnspan=3, sticky="ew", padx=(4, 4), pady=2
        )
        ttk.Button(
            stack_structure_panel,
            text="Load data...",
            command=self._browse_superlattice_data_file,
        ).grid(row=2, column=4, sticky="ew", pady=2)
        ttk.Label(stack_structure_panel, text="Superlattice name").grid(
            row=3, column=0, sticky="w"
        )
        ttk.Entry(stack_structure_panel, textvariable=self.stack_name_var).grid(
            row=3, column=1, columnspan=2, sticky="ew", padx=(4, 8), pady=2
        )
        add_tooltip(
            ttk.Label(stack_structure_panel, text="Repetitions"), "repetitions"
        ).grid(row=3, column=3, sticky="e")
        add_tooltip(
            ttk.Spinbox(
                stack_structure_panel,
                textvariable=self.stack_repetitions_var,
                from_=1,
                to=10000,
                width=6,
            ),
            "repetitions",
        ).grid(row=3, column=4, sticky="ew", padx=(4, 0), pady=2)
        add_tooltip(
            ttk.Checkbutton(
                stack_structure_panel,
                text="Include capping layer",
                variable=self.stack_capping_enabled_var,
                command=self._on_stack_capping_changed,
            ),
            "capping_layer",
        ).grid(row=4, column=1, columnspan=3, sticky="w", padx=(4, 0), pady=2)
        stack_layers_scroller = ttk.Frame(stack_structure_panel, style="Panel.TFrame")
        stack_layers_scroller.grid(
            row=5, column=0, columnspan=5, sticky="nsew", pady=(6, 4)
        )
        self.stack_layers_canvas = tk.Canvas(
            stack_layers_scroller,
            height=260,
            background=UI_COLORS["panel"],
            highlightthickness=0,
        )
        stack_layers_scrollbar = ttk.Scrollbar(
            stack_layers_scroller,
            orient=tk.VERTICAL,
            command=self.stack_layers_canvas.yview,
        )
        self.stack_layers_canvas.configure(yscrollcommand=stack_layers_scrollbar.set)
        self.stack_layers_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        stack_layers_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.stack_layers_frame = ttk.Frame(
            self.stack_layers_canvas, style="Panel.TFrame"
        )
        stack_layers_window = self.stack_layers_canvas.create_window(
            (0, 0), window=self.stack_layers_frame, anchor="nw"
        )
        self.stack_layers_frame.bind(
            "<Configure>",
            lambda _event: self.stack_layers_canvas.configure(
                scrollregion=self.stack_layers_canvas.bbox("all")
            ),
        )
        self.stack_layers_canvas.bind(
            "<Configure>",
            lambda event: self.stack_layers_canvas.itemconfigure(
                stack_layers_window, width=event.width
            ),
        )
        ttk.Button(
            stack_structure_panel,
            text="Add repeated layer",
            command=self._add_stack_layer,
        ).grid(
            row=6, column=1, sticky="ew", padx=(4, 2), pady=2
        )
        ttk.Button(
            stack_structure_panel,
            text="Remove last repeated layer",
            command=self._remove_stack_layer,
        ).grid(row=6, column=2, columnspan=2, sticky="ew", padx=(2, 0), pady=2)
        for column in (1, 2, 3):
            stack_structure_panel.columnconfigure(column, weight=1)
        stack_structure_panel.rowconfigure(5, weight=1)

        stack_strain_panel = ttk.Frame(stack_tabs, padding=8, style="Panel.TFrame")
        stack_tabs.add(stack_strain_panel, text="Strain")
        stack_strain_scroller = ttk.Frame(stack_strain_panel, style="Panel.TFrame")
        stack_strain_scroller.pack(fill=tk.BOTH, expand=True)
        self.stack_strain_canvas = tk.Canvas(
            stack_strain_scroller,
            height=260,
            background=UI_COLORS["panel"],
            highlightthickness=0,
        )
        stack_strain_scrollbar = ttk.Scrollbar(
            stack_strain_scroller,
            orient=tk.VERTICAL,
            command=self.stack_strain_canvas.yview,
        )
        self.stack_strain_canvas.configure(yscrollcommand=stack_strain_scrollbar.set)
        self.stack_strain_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        stack_strain_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.stack_strain_frame = ttk.Frame(
            self.stack_strain_canvas, style="Panel.TFrame"
        )
        stack_strain_window = self.stack_strain_canvas.create_window(
            (0, 0), window=self.stack_strain_frame, anchor="nw"
        )
        self.stack_strain_frame.bind(
            "<Configure>",
            lambda _event: self.stack_strain_canvas.configure(
                scrollregion=self.stack_strain_canvas.bbox("all")
            ),
        )
        self.stack_strain_canvas.bind(
            "<Configure>",
            lambda event: self.stack_strain_canvas.itemconfigure(
                stack_strain_window, width=event.width
            ),
        )
        self._populate_stack_rows(self.stack_document)

        stack_calculation_panel = ttk.Frame(stack_tabs, padding=8, style="Panel.TFrame")
        stack_tabs.add(stack_calculation_panel, text="Calculation")
        for column, label in ((1, "Fit"), (2, "Value"), (3, "Min"), (4, "Max"), (5, "Range")):
            ttk.Label(stack_calculation_panel, text=label).grid(
                row=0, column=column, sticky="ew"
            )

        def add_stack_calculation_row(
            row: int,
            key: str,
            label: str | tk.StringVar,
            value_var: tk.StringVar,
            min_var: tk.StringVar,
            max_var: tk.StringVar,
            tooltip_key: str,
        ) -> tuple[list[tk.Widget], tuple[ttk.Entry, ttk.Entry, ttk.Entry], RangeIndicator]:
            label_widget = ttk.Label(
                stack_calculation_panel,
                textvariable=label if isinstance(label, tk.StringVar) else None,
                text="" if isinstance(label, tk.StringVar) else label,
            )
            label_widget.grid(row=row, column=0, sticky="w", pady=2)
            add_tooltip(label_widget, tooltip_key)
            fit_checkbutton = ttk.Checkbutton(
                stack_calculation_panel,
                variable=self.stack_calculation_fit_enabled_vars[key],
            )
            fit_checkbutton.grid(row=row, column=1, sticky="ew", pady=2)
            add_tooltip(fit_checkbutton, "fit")
            entries = tuple(
                ttk.Entry(stack_calculation_panel, textvariable=variable, width=9)
                for variable in (value_var, min_var, max_var)
            )
            for column, entry in zip((2, 3, 4), entries):
                entry.grid(row=row, column=column, sticky="ew", padx=(4, 0), pady=2)
            self._add_range_indicator(
                stack_calculation_panel,
                row,
                str(label),
                value_var,
                min_var,
                max_var,
                self.stack_calculation_fit_vars[key],
                self.stack_calculation_fit_enabled_vars[key],
            )
            indicator = self.range_indicators[-1]
            return [label_widget, fit_checkbutton, *entries, indicator.canvas], entries, indicator

        add_tooltip(
            ttk.Label(stack_calculation_panel, textvariable=self.stack_sampling_label_var),
            "points_per_unit",
        ).grid(row=1, column=0, sticky="w", pady=2)
        add_tooltip(
            ttk.Entry(
                stack_calculation_panel,
                textvariable=self.stack_points_per_unit_var,
                width=9,
            ),
            "points_per_unit",
        ).grid(row=1, column=2, sticky="ew", padx=(4, 0), pady=2)

        self.stack_calculation_entries: list[ttk.Entry] = []
        self.stack_calculation_min_entries: list[ttk.Entry] = []
        self.stack_calculation_max_entries: list[ttk.Entry] = []
        self.stack_model_calculation_indicators: list[RangeIndicator] = []
        dynamic_rows = (
            ("resolution", "Resolution (deg)", self.dynamic_resolution_start_var, self.dynamic_resolution_min_var, self.dynamic_resolution_max_var, "resolution"),
            ("intensity_scale", "Intensity scale", self.dynamic_intensity_start_var, self.dynamic_intensity_min_var, self.dynamic_intensity_max_var, "intensity_scale"),
            ("background_a", "Background slope (a)", self.dynamic_bkg_a_start_var, self.dynamic_bkg_a_min_var, self.dynamic_bkg_a_max_var, "background_a"),
            ("background_b", "Background offset (b)", self.dynamic_bkg_b_start_var, self.dynamic_bkg_b_min_var, self.dynamic_bkg_b_max_var, "background_b"),
            ("background_c", "Background curvature (c)", self.dynamic_bkg_c_start_var, self.dynamic_bkg_c_min_var, self.dynamic_bkg_c_max_var, "background_c"),
        )
        for row, values in enumerate(dynamic_rows, start=2):
            _widgets, entries, indicator = add_stack_calculation_row(row, *values)
            self.stack_calculation_entries.append(entries[0])
            self.stack_calculation_min_entries.append(entries[1])
            self.stack_calculation_max_entries.append(entries[2])
            self.stack_model_calculation_indicators.append(indicator)

        self.stack_dynamic_calculation_widgets: list[tk.Widget] = []
        for row, label, variable in (
            (7, "Density slices per cell", self.density_slices_var),
            (8, "Density Q max (1/Å)", self.density_max_q0_var),
        ):
            label_widget = ttk.Label(stack_calculation_panel, text=label)
            label_widget.grid(row=row, column=0, sticky="w", pady=2)
            entry = ttk.Entry(stack_calculation_panel, textvariable=variable, width=9)
            entry.grid(row=row, column=2, sticky="ew", padx=(4, 0), pady=2)
            tooltip_key = "density_slices" if variable is self.density_slices_var else "density_q_max"
            add_tooltip(label_widget, tooltip_key)
            add_tooltip(entry, tooltip_key)
            self.stack_dynamic_calculation_widgets.extend((label_widget, entry))
        for column in (2, 3, 4, 5):
            stack_calculation_panel.columnconfigure(column, weight=1)

        stack_optimization_frame = ttk.Frame(
            stack_tabs, padding=8, style="Panel.TFrame"
        )
        stack_tabs.add(stack_optimization_frame, text="Optimization settings")
        self._add_accent_strip(stack_optimization_frame, UI_COLORS["fit"])
        for label, variable, row, column, tooltip_key in optimization_controls:
            add_optimization_control(
                stack_optimization_frame,
                label,
                variable,
                row,
                column,
                tooltip_key,
            )
        stack_optimization_frame.columnconfigure(1, weight=1)
        stack_optimization_frame.columnconfigure(3, weight=1)

        kinematic_trace_vars = (
            self.kin_resolution_start_var,
            self.kin_resolution_min_var,
            self.kin_resolution_max_var,
            self.kin_resolution_fit_var,
            self.kin_scale_start_var,
            self.kin_scale_min_var,
            self.kin_scale_max_var,
            self.kin_scale_fit_var,
            self.kin_bkg_a_start_var,
            self.kin_bkg_a_min_var,
            self.kin_bkg_a_max_var,
            self.kin_bkg_a_fit_var,
            self.kin_bkg_b_start_var,
            self.kin_bkg_b_min_var,
            self.kin_bkg_b_max_var,
            self.kin_bkg_b_fit_var,
            self.kin_bkg_c_start_var,
            self.kin_bkg_c_min_var,
            self.kin_bkg_c_max_var,
            self.kin_bkg_c_fit_var,
        )
        for variable in kinematic_trace_vars:
            variable.trace_add("write", self._redraw_stack_model_calculation_indicators)

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
        self.root.after(
            100, lambda: main_pane.sashpos(0, round(main_pane.winfo_width() * 0.4))
        )
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
        )
        watermark.pack(fill=tk.X, pady=(4, 0))
        watermark.bind("<Button-1>", lambda _event: webbrowser.open_new(GENL_DOI_URL))
        self._draw_empty_plot()
        self.kinematic_widgets = (
            self._collect_children(self.kinematic_frame)
            + self._collect_children(self.kinematic_substrate_frame)
        )
        self.kinematic_substrate_widgets = self._collect_children(
            self.kinematic_substrate_controls_frame
        )
        self.dynamic_widgets = (
            self._collect_children(self.dynamic_film_frame)
            + self._collect_children(self.dynamic_fit_frame)
            + self._collect_children(self.dynamic_substrate_frame)
        )
        self.strain_widgets = self._collect_children(self.strain_frame)
        self._sync_workspace_panels()
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
        column: int = 5,
    ) -> None:
        canvas = tk.Canvas(parent, width=96, height=20, highlightthickness=0, bg="#f5f5f5")
        canvas.grid(row=row, column=column, sticky="ew", padx=(4, 0))
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

    def _generate_optimizer_seed(self) -> None:
        seed = secrets.randbits(32)
        self.seed_var.set(str(seed))
        self.status_var.set(f"Generated optimizer seed: {seed}")

    def _populate_stack_rows(self, document: dict[str, object]) -> None:
        for indicator in self.stack_structure_indicators:
            if indicator in self.range_indicators:
                self.range_indicators.remove(indicator)
        self.stack_structure_indicators.clear()
        for child in self.stack_layers_frame.winfo_children():
            child.destroy()
        for child in self.stack_strain_frame.winfo_children():
            child.destroy()
        self.stack_document = copy.deepcopy(document)
        self.stack_name_var.set(str(document["name"]))
        sequence = document["sequence"]
        self.stack_repetitions_var.set(str(sequence["repetitions"]))
        capping_layer = document.get("capping_layer")
        self.stack_capping_enabled_var.set(capping_layer is not None)
        self.stack_row_specs = [
            copy.deepcopy(document["substrate"]),
            *(copy.deepcopy(layer) for layer in sequence["layers"]),
            *([] if capping_layer is None else [copy.deepcopy(capping_layer)]),
        ]
        self.stack_row_roles = [
            "substrate",
            *("repeat" for _layer in sequence["layers"]),
            *([] if capping_layer is None else ["capping"]),
        ]
        self.stack_row_vars.clear()
        fit_parameters = {
            str(parameter["target"]): parameter
            for parameter in document.get("fit_parameters", [])
        }
        for column, (label, width) in enumerate(
            (
                ("Role", 8),
                ("Name / parameter", 13),
                ("Structure / fit", 14),
                ("Dir / value", 8),
                ("Cells / min", 8),
                ("Max", 7),
                ("Range", 12),
            )
        ):
            ttk.Label(self.stack_layers_frame, text=label, width=width).grid(
                row=0, column=column, sticky="ew", padx=1
            )

        def default_bounds(key: str, value: float) -> tuple[float, float]:
            if key == "dinterface":
                return max(0.0, value - 0.5), value + 0.5
            if key == "scale":
                return 0.9, 1.1
            for strain_key, _label, lower, upper in STACK_STRAIN_FIELDS:
                if key == strain_key:
                    return lower, upper
            return 0.5, 1.5

        for index, (role, spec) in enumerate(
            zip(self.stack_row_roles, self.stack_row_specs)
        ):
            row = 1 + index * 4
            is_substrate = role == "substrate"
            prefix = role if role in {"substrate", "capping"} else str(spec["name"])
            variables = {
                "name": tk.StringVar(value="substrate" if is_substrate else str(spec["name"])),
                "filename": tk.StringVar(value=str(spec["filename"])),
                "direction": tk.StringVar(value=str(spec["direction"])),
                "unit_cells": tk.StringVar(value=str(spec["unit_cells"])),
                "dinterface": tk.StringVar(value=str(spec.get("dinterface", 0.0))),
                "scale": tk.StringVar(value=str(spec.get("scale", 1.0))),
                "area_scale": tk.StringVar(value=str(spec.get("area_scale", 1.0))),
                **{
                    key: tk.StringVar(value=str(spec.get(key, 0.0)))
                    for key, _label, _lower, _upper in STACK_STRAIN_FIELDS
                },
            }
            for key in (
                "dinterface",
                "scale",
                "area_scale",
                *(field[0] for field in STACK_STRAIN_FIELDS),
            ):
                value = float(variables[key].get())
                parameter = fit_parameters.get(f"{prefix}.{key}")
                lower, upper = (
                    (float(parameter["min"]), float(parameter["max"]))
                    if parameter is not None
                    else default_bounds(key, value)
                )
                variables[f"{key}_min"] = tk.StringVar(value=f"{lower:g}")
                variables[f"{key}_max"] = tk.StringVar(value=f"{upper:g}")
                variables[f"{key}_fit"] = tk.StringVar(value="")
                variables[f"{key}_fit_enabled"] = tk.BooleanVar(value=parameter is not None)
            self.stack_row_vars.append(variables)
            ttk.Label(
                self.stack_layers_frame,
                text={"substrate": "Substrate", "repeat": "Repeat", "capping": "Capping"}[role],
            ).grid(row=row, column=0, sticky="w", padx=1, pady=1)
            ttk.Entry(
                self.stack_layers_frame,
                textvariable=variables["name"],
                width=8,
                state=tk.DISABLED if is_substrate else tk.NORMAL,
            ).grid(row=row, column=1, sticky="ew", padx=1, pady=1)
            structure_combo = ttk.Combobox(
                self.stack_layers_frame,
                textvariable=variables["filename"],
                values=POSCAR_FILES,
                state="readonly",
                width=18,
            )
            structure_combo.grid(row=row, column=2, sticky="ew", padx=1, pady=1)
            add_tooltip(structure_combo, "layer_structure")
            direction_combo = ttk.Combobox(
                self.stack_layers_frame,
                textvariable=variables["direction"],
                values=("1", "2", "3"),
                state="readonly",
                width=4,
            )
            direction_combo.grid(row=row, column=3, sticky="ew", padx=1, pady=1)
            add_tooltip(direction_combo, "layer_direction")
            unit_cells_entry = ttk.Entry(
                self.stack_layers_frame,
                textvariable=variables["unit_cells"],
                width=8,
            )
            unit_cells_entry.grid(row=row, column=4, sticky="ew", padx=1, pady=1)
            add_tooltip(unit_cells_entry, "unit_cells")

            for offset, (key, label, tooltip_key) in enumerate(
                (
                    ("dinterface", "Interface (A)", "interface_spacing"),
                    ("scale", "Lattice scale", "lattice_scale"),
                    ("area_scale", "Area scale", "area_scale"),
                ),
                start=1,
            ):
                parameter_row = row + offset
                add_tooltip(
                    ttk.Label(self.stack_layers_frame, text=label), tooltip_key
                ).grid(
                    row=parameter_row, column=1, sticky="w", padx=1, pady=1
                )
                add_tooltip(
                    ttk.Checkbutton(
                        self.stack_layers_frame,
                        variable=variables[f"{key}_fit_enabled"],
                    ),
                    "fit",
                ).grid(row=parameter_row, column=2, sticky="ew", padx=1, pady=1)
                for column, variable_key in (
                    (3, key),
                    (4, f"{key}_min"),
                    (5, f"{key}_max"),
                ):
                    ttk.Entry(
                        self.stack_layers_frame,
                        textvariable=variables[variable_key],
                        width=8,
                    ).grid(row=parameter_row, column=column, sticky="ew", padx=1, pady=1)
                self._add_range_indicator(
                    self.stack_layers_frame,
                    parameter_row,
                    f"{prefix} {label}",
                    variables[key],
                    variables[f"{key}_min"],
                    variables[f"{key}_max"],
                    variables[f"{key}_fit"],
                    variables[f"{key}_fit_enabled"],
                    column=6,
                )
                self.stack_structure_indicators.append(self.range_indicators[-1])
        self.stack_layers_frame.columnconfigure(2, weight=1)
        self.stack_layers_frame.columnconfigure(6, weight=1)

        for column, label in enumerate(
            ("Layer", "Parameter", "Fit", "Value", "Min", "Max", "Range")
        ):
            ttk.Label(self.stack_strain_frame, text=label).grid(
                row=0, column=column, sticky="ew", padx=1
            )
        strain_row = 1
        for role, variables in zip(self.stack_row_roles, self.stack_row_vars):
            if role == "substrate":
                continue
            prefix = role if role == "capping" else variables["name"].get().strip()
            for offset, (key, label, _lower, _upper) in enumerate(STACK_STRAIN_FIELDS):
                row = strain_row + offset
                if offset == 0:
                    ttk.Label(
                        self.stack_strain_frame,
                        textvariable=variables["name"],
                    ).grid(
                        row=row, column=0, rowspan=4, sticky="nw", padx=1, pady=2
                    )
                tooltip_key = {
                    "bottom_strain_amplitude": "bottom_strain_amplitude",
                    "bottom_strain_end": "bottom_strain_extent",
                    "top_strain_amplitude": "top_strain_amplitude",
                    "top_strain_end": "top_strain_extent",
                }[key]
                add_tooltip(
                    ttk.Label(self.stack_strain_frame, text=label), tooltip_key
                ).grid(
                    row=row, column=1, sticky="w", padx=1, pady=2
                )
                add_tooltip(
                    ttk.Checkbutton(
                        self.stack_strain_frame,
                        variable=variables[f"{key}_fit_enabled"],
                    ),
                    "fit",
                ).grid(row=row, column=2, sticky="ew", padx=1, pady=2)
                for column, variable_key in (
                    (3, key),
                    (4, f"{key}_min"),
                    (5, f"{key}_max"),
                ):
                    ttk.Entry(
                        self.stack_strain_frame,
                        textvariable=variables[variable_key],
                        width=8,
                    ).grid(row=row, column=column, sticky="ew", padx=1, pady=2)
                self._add_range_indicator(
                    self.stack_strain_frame,
                    row,
                    f"{prefix} {label}",
                    variables[key],
                    variables[f"{key}_min"],
                    variables[f"{key}_max"],
                    variables[f"{key}_fit"],
                    variables[f"{key}_fit_enabled"],
                    column=6,
                )
                self.stack_structure_indicators.append(self.range_indicators[-1])
            strain_row += len(STACK_STRAIN_FIELDS)
        self.stack_strain_frame.columnconfigure(1, weight=1)
        self.stack_strain_frame.columnconfigure(6, weight=1)

    def _stack_calculation_triplets(self) -> dict[str, tuple[float, float, float]]:
        model_variables = (
            (
                ("resolution", self.dynamic_resolution_start_var, self.dynamic_resolution_min_var, self.dynamic_resolution_max_var),
                ("intensity_scale", self.dynamic_intensity_start_var, self.dynamic_intensity_min_var, self.dynamic_intensity_max_var),
                ("background_a", self.dynamic_bkg_a_start_var, self.dynamic_bkg_a_min_var, self.dynamic_bkg_a_max_var),
                ("background_b", self.dynamic_bkg_b_start_var, self.dynamic_bkg_b_min_var, self.dynamic_bkg_b_max_var),
                ("background_c", self.dynamic_bkg_c_start_var, self.dynamic_bkg_c_min_var, self.dynamic_bkg_c_max_var),
            )
            if self.model_var.get() == "Dynamic"
            else (
                ("resolution", self.kin_resolution_start_var, self.kin_resolution_min_var, self.kin_resolution_max_var),
                ("intensity_scale", self.kin_scale_start_var, self.kin_scale_min_var, self.kin_scale_max_var),
                ("background_a", self.kin_bkg_a_start_var, self.kin_bkg_a_min_var, self.kin_bkg_a_max_var),
                ("background_b", self.kin_bkg_b_start_var, self.kin_bkg_b_min_var, self.kin_bkg_b_max_var),
                ("background_c", self.kin_bkg_c_start_var, self.kin_bkg_c_min_var, self.kin_bkg_c_max_var),
            )
        )
        return {
            name: (float(value.get()), float(minimum.get()), float(maximum.get()))
            for name, value, minimum, maximum in model_variables
        }

    def _stack_document_from_controls(
        self, allow_outside_start: bool = False
    ) -> dict[str, object]:
        if "repeat" not in self.stack_row_roles:
            raise ValueError("A superlattice requires a substrate and at least one repeated layer")
        document = copy.deepcopy(self.stack_document)
        document["name"] = self.stack_name_var.get().strip()
        document["model"] = self.model_var.get().lower()
        document["wavelength"] = self._parse_wavelength()
        fit_parameters: list[dict[str, object]] = []
        sequence = document["sequence"]
        sequence["repetitions"] = int(self.stack_repetitions_var.get())

        substrate_spec: dict[str, object] | None = None
        repeated_specs: list[dict[str, object]] = []
        capping_spec: dict[str, object] | None = None
        for role, base, variables in zip(
            self.stack_row_roles, self.stack_row_specs, self.stack_row_vars
        ):
            spec = copy.deepcopy(base)
            prefix = role if role in {"substrate", "capping"} else variables["name"].get().strip()
            if role != "substrate":
                spec["name"] = prefix
                if role == "capping":
                    spec["name"] = variables["name"].get().strip()
                spec.pop("monolayers", None)
            spec["filename"] = variables["filename"].get().strip()
            spec["direction"] = int(variables["direction"].get())
            spec["unit_cells"] = float(variables["unit_cells"].get())
            spec["dinterface"] = float(variables["dinterface"].get())
            spec["scale"] = float(variables["scale"].get())
            spec["area_scale"] = float(variables["area_scale"].get())
            for key, _label, _lower, _upper in STACK_STRAIN_FIELDS:
                spec[key] = float(variables[key].get())
            if spec["dinterface"] < 0 or spec["scale"] <= 0 or spec["area_scale"] <= 0:
                raise ValueError(
                    f"{prefix} interface must be non-negative; scale and area scale must be positive"
                )
            for key in ("dinterface", "scale", "area_scale"):
                if not variables[f"{key}_fit_enabled"].get():
                    continue
                values = (
                    float(spec[key]),
                    float(variables[f"{key}_min"].get()),
                    float(variables[f"{key}_max"].get()),
                )
                validate_start_min_max(
                    f"{prefix} {key}",
                    *values,
                    allow_outside_start=allow_outside_start,
                )
                if values[1] >= values[2]:
                    raise ValueError(f"{prefix} {key}: min must be less than max")
                if key == "dinterface" and values[1] < 0:
                    raise ValueError(f"{prefix} interface minimum must be non-negative")
                if key != "dinterface" and values[1] <= 0:
                    raise ValueError(f"{prefix} {key} minimum must be positive")
                if allow_outside_start:
                    spec[key] = float(np.clip(values[0], values[1], values[2]))
                fit_parameters.append(
                    {"target": f"{prefix}.{key}", "min": values[1], "max": values[2]}
                )
            if role != "substrate":
                for key, label, _lower, _upper in STACK_STRAIN_FIELDS:
                    value = float(spec[key])
                    minimum = float(variables[f"{key}_min"].get())
                    maximum = float(variables[f"{key}_max"].get())
                    if key.endswith("_end") and value < 0:
                        raise ValueError(f"{prefix} {label.lower()} must be non-negative")
                    if not variables[f"{key}_fit_enabled"].get():
                        continue
                    validate_start_min_max(
                        f"{prefix} {label.lower()}",
                        value,
                        minimum,
                        maximum,
                        allow_outside_start=allow_outside_start,
                    )
                    if minimum >= maximum:
                        raise ValueError(f"{prefix} {label.lower()}: min must be less than max")
                    if key.endswith("_end") and minimum < 0:
                        raise ValueError(
                            f"{prefix} {label.lower()} minimum must be non-negative"
                        )
                    if allow_outside_start:
                        spec[key] = float(np.clip(value, minimum, maximum))
                    fit_parameters.append(
                        {
                            "target": f"{prefix}.{key}",
                            "min": minimum,
                            "max": maximum,
                        }
                    )
            if role == "substrate":
                substrate_spec = spec
            elif role == "capping":
                capping_spec = spec
            else:
                repeated_specs.append(spec)
        document["substrate"] = substrate_spec
        sequence["layers"] = repeated_specs
        if capping_spec is None:
            document.pop("capping_layer", None)
        else:
            document["capping_layer"] = capping_spec

        calculation = dict(document.get("calculation", {}))
        if self.model_var.get() == "Dynamic":
            calculation.update(
                resolution=float(self.dynamic_resolution_start_var.get()),
                intensity_scale=float(self.dynamic_intensity_start_var.get()),
                background_a=float(self.dynamic_bkg_a_start_var.get()),
                background_b=float(self.dynamic_bkg_b_start_var.get()),
                background_c=float(self.dynamic_bkg_c_start_var.get()),
            )
        else:
            calculation.update(
                resolution=float(self.kin_resolution_start_var.get()),
                intensity_scale=float(self.kin_scale_start_var.get()),
                background_a=float(self.kin_bkg_a_start_var.get()),
                background_b=float(self.kin_bkg_b_start_var.get()),
                background_c=float(self.kin_bkg_c_start_var.get()),
            )
        calculation.update(
            density_slices=int(self.density_slices_var.get()),
            density_max_q0=float(self.density_max_q0_var.get()),
        )
        points_per_unit = float(self.stack_points_per_unit_var.get())
        if points_per_unit <= 0:
            raise ValueError("Superlattice simulation points per axis unit must be positive")
        calculation["q_step" if self._axis_is_q() else "twotheta_step"] = 1.0 / points_per_unit
        document["calculation"] = calculation
        calculation_triplets = self._stack_calculation_triplets()
        for name in ("resolution", "intensity_scale"):
            if calculation_triplets[name][0] <= 0:
                raise ValueError(
                    f"Superlattice {name.replace('_', ' ')} must be positive"
                )
        for name, values in calculation_triplets.items():
            if not self.stack_calculation_fit_enabled_vars[name].get():
                continue
            validate_start_min_max(
                name.replace("_", " "),
                *values,
                allow_outside_start=allow_outside_start,
            )
            if values[1] >= values[2]:
                raise ValueError(
                    f"Superlattice {name.replace('_', ' ')}: min must be less than max"
                )
            if name in {"resolution", "intensity_scale"} and values[1] <= 0:
                raise ValueError(
                    f"Superlattice {name.replace('_', ' ')} minimum must be positive"
                )
            if allow_outside_start:
                calculation[name] = float(np.clip(values[0], values[1], values[2]))
            fit_parameters.append(
                {
                    "target": f"calculation.{name}",
                    "min": values[1],
                    "max": values[2],
                }
            )
        document["calculation_ranges"] = {
            name: [values[1], values[2]] for name, values in calculation_triplets.items()
        }
        document["fit_parameters"] = fit_parameters
        StackDefinition(Path(self.stack_path_var.get()), document)
        return document

    def _set_stack_document(self, document: dict[str, object], sync_calculation: bool) -> None:
        self._populate_stack_rows(document)
        fit_parameters = {
            str(parameter["target"]): parameter
            for parameter in document.get("fit_parameters", [])
        }
        for name, variable in self.stack_calculation_fit_enabled_vars.items():
            variable.set(f"calculation.{name}" in fit_parameters)
        if not sync_calculation:
            return
        model = str(document.get("model", "kinematic")).title()
        self.model_var.set(model)
        self.wavelength_var.set(str(document.get("wavelength", CU_K_ALPHA_WAVELENGTH)))
        calculation = document.get("calculation", {})
        targets = (
            (
                self.dynamic_resolution_start_var,
                self.dynamic_intensity_start_var,
                self.dynamic_bkg_a_start_var,
                self.dynamic_bkg_b_start_var,
                self.dynamic_bkg_c_start_var,
            )
            if model == "Dynamic"
            else (
                self.kin_resolution_start_var,
                self.kin_scale_start_var,
                self.kin_bkg_a_start_var,
                self.kin_bkg_b_start_var,
                self.kin_bkg_c_start_var,
            )
        )
        for variable, key in zip(
            targets,
            (
                "resolution",
                "intensity_scale",
                "background_a",
                "background_b",
                "background_c",
            ),
        ):
            if key in calculation:
                variable.set(str(calculation[key]))
        if "background_c" not in calculation:
            targets[4].set("0.0")
        for variable, key in (
            (self.density_slices_var, "density_slices"),
            (self.density_max_q0_var, "density_max_q0"),
        ):
            if key in calculation:
                variable.set(str(calculation[key]))
        step = calculation.get("q_step" if self._axis_is_q() else "twotheta_step")
        if step is not None and float(step) > 0:
            self.stack_points_per_unit_var.set(f"{1.0 / float(step):.8g}")
        ranges = document.get("calculation_ranges", {})
        model_range_variables = (
            {
                "resolution": (self.dynamic_resolution_min_var, self.dynamic_resolution_max_var),
                "intensity_scale": (self.dynamic_intensity_min_var, self.dynamic_intensity_max_var),
                "background_a": (self.dynamic_bkg_a_min_var, self.dynamic_bkg_a_max_var),
                "background_b": (self.dynamic_bkg_b_min_var, self.dynamic_bkg_b_max_var),
                "background_c": (self.dynamic_bkg_c_min_var, self.dynamic_bkg_c_max_var),
            }
            if model == "Dynamic"
            else {
                "resolution": (self.kin_resolution_min_var, self.kin_resolution_max_var),
                "intensity_scale": (self.kin_scale_min_var, self.kin_scale_max_var),
                "background_a": (self.kin_bkg_a_min_var, self.kin_bkg_a_max_var),
                "background_b": (self.kin_bkg_b_min_var, self.kin_bkg_b_max_var),
                "background_c": (self.kin_bkg_c_min_var, self.kin_bkg_c_max_var),
            }
        )
        range_variables = model_range_variables
        range_value_variables = {
            "resolution": targets[0],
            "intensity_scale": targets[1],
            "background_a": targets[2],
            "background_b": targets[3],
            "background_c": targets[4],
        }
        for name, variables in range_variables.items():
            parameter = fit_parameters.get(f"calculation.{name}")
            values = ranges.get(name)
            if parameter is not None:
                values = [parameter["min"], parameter["max"]]
            if isinstance(values, list) and len(values) == 2:
                variables[0].set(str(values[0]))
                variables[1].set(str(values[1]))
            elif name in range_value_variables:
                value = float(range_value_variables[name].get())
                variables[0].set(str(min(float(variables[0].get()), value)))
                variables[1].set(str(max(float(variables[1].get()), value)))

    def _browse_stack_file(self) -> None:
        selected = filedialog.askopenfilename(
            title="Choose GenL superlattice file",
            initialdir=str(STACK_DIR),
            filetypes=[("GenL superlattice", "*.json"), ("All files", "*")],
        )
        if selected:
            self.stack_path_var.set(selected)
            self._load_stack_file()

    def _browse_superlattice_data_file(self) -> None:
        current_path = resolve_data_path(self.data_path_var.get())
        selected = filedialog.askopenfilename(
            title="Choose superlattice data file",
            initialdir=str(
                current_path.parent if current_path.parent.exists() else EXAMPLE_DATA_DIR
            ),
            filetypes=[("Text data", "*.txt *.dat *.csv"), ("All files", "*")],
        )
        if not selected:
            return
        try:
            twotheta, _observed = read_experimental_data(selected)
        except (OSError, ValueError) as exc:
            messagebox.showerror("Cannot load superlattice data", str(exc))
            return
        self.data_path_var.set(selected)
        self.superlattice_data_preview = True
        self._set_twotheta_window(twotheta)
        self.preview_data_path = Path(selected).resolve()
        if self.preview_after_id is not None:
            self.root.after_cancel(self.preview_after_id)
            self.preview_after_id = None
        self._draw_experimental_preview()
        self.status_var.set(
            f"Loaded superlattice data: {Path(selected).name} ({len(twotheta)} points)"
        )

    def _load_stack_file(self) -> None:
        try:
            definition = StackDefinition.load(self.stack_path_var.get())
            self.stack_path_var.set(str(definition.path))
            self._set_stack_document(definition.document, sync_calculation=True)
            self.status_var.set(f"Loaded superlattice: {definition.name}")
        except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
            messagebox.showerror("Cannot load superlattice", str(exc))

    def _save_stack_file(self) -> None:
        try:
            document = self._stack_document_from_controls()
            selected = filedialog.asksaveasfilename(
                title="Save GenL superlattice",
                initialdir=str(STACK_DIR),
                initialfile=Path(self.stack_path_var.get()).name,
                defaultextension=".json",
                filetypes=[("GenL superlattice", "*.json"), ("All files", "*")],
            )
            if not selected:
                return
            write_json_document(Path(selected), document)
        except Exception as exc:  # Save callbacks must not escape into Tk.
            self._report_save_error("Cannot save superlattice", exc)
            return
        self.stack_path_var.set(selected)
        self.stack_document = copy.deepcopy(document)
        self.status_var.set(f"Saved superlattice: {selected}")

    def _on_stack_capping_changed(self) -> None:
        try:
            document = self._stack_document_from_controls()
        except (TypeError, ValueError) as exc:
            self.stack_capping_enabled_var.set(not self.stack_capping_enabled_var.get())
            messagebox.showerror("Cannot update capping layer", str(exc))
            return
        if self.stack_capping_enabled_var.get():
            if document.get("capping_layer") is None:
                capping_layer = copy.deepcopy(document["sequence"]["layers"][-1])
                capping_layer.update(
                    name="Capping",
                    unit_cells=1.0,
                    dinterface=0.0,
                    scale=1.0,
                    area_scale=1.0,
                )
                document["capping_layer"] = capping_layer
        else:
            document.pop("capping_layer", None)
            document["fit_parameters"] = [
                parameter
                for parameter in document.get("fit_parameters", [])
                if not str(parameter["target"]).startswith("capping.")
            ]
        self._populate_stack_rows(document)

    def _add_stack_layer(self) -> None:
        try:
            document = self._stack_document_from_controls()
        except (TypeError, ValueError) as exc:
            messagebox.showerror("Cannot add layer", str(exc))
            return
        layers = document["sequence"]["layers"]
        existing = {str(layer["name"]) for layer in layers}
        number = 1
        while f"Layer {number}" in existing:
            number += 1
        layer = copy.deepcopy(layers[-1])
        layer["name"] = f"Layer {number}"
        layers.append(layer)
        self._populate_stack_rows(document)

    def _remove_stack_layer(self) -> None:
        try:
            document = self._stack_document_from_controls()
        except (TypeError, ValueError) as exc:
            messagebox.showerror("Cannot remove layer", str(exc))
            return
        layers = document["sequence"]["layers"]
        if len(layers) <= 1:
            messagebox.showinfo(
                "Superlattice", "A superlattice requires at least one repeated layer."
            )
            return
        removed_name = str(layers.pop()["name"])
        document["fit_parameters"] = [
            parameter
            for parameter in document.get("fit_parameters", [])
            if not str(parameter["target"]).startswith(f"{removed_name}.")
        ]
        self._populate_stack_rows(document)

    def _on_stack_enabled_changed(self, *_args: object) -> None:
        self._clear_fit_result_values()
        self._sync_workspace_panels()
        self._sync_model_controls()
        self._sync_optional_controls()
        state = tk.DISABLED if self.stack_enabled_var.get() else tk.NORMAL
        self.strain_checkbutton.configure(state=state)
        self.roughness_checkbutton.configure(state=state)
        self.run_button.configure(
            state=tk.DISABLED if self.running else tk.NORMAL
        )
        self.status_var.set(
            "Superlattice workspace selected"
            if self.stack_enabled_var.get()
            else "Film workspace selected"
        )
        if self.stack_enabled_var.get():
            if self.superlattice_data_preview:
                self._draw_experimental_preview()
            else:
                self._draw_empty_plot()
        else:
            self._schedule_data_preview()

    def _sync_workspace_panels(self) -> None:
        if self.stack_enabled_var.get():
            self.film_container.grid_remove()
            self.stack_container.grid()
        else:
            self.stack_container.grid_remove()
            self.film_container.grid()

    def _sync_model_controls(self) -> None:
        is_dynamic = self.model_var.get() == "Dynamic"
        self._set_widgets_enabled(self.dynamic_widgets, is_dynamic)
        self._set_widgets_enabled(self.kinematic_widgets, not is_dynamic)
        self._sync_stack_calculation_controls(is_dynamic)
        self._sync_kinematic_substrate_controls()
        self._sync_optional_controls()

    def _sync_stack_calculation_controls(self, is_dynamic: bool) -> None:
        variable_rows = (
            (
                (self.dynamic_resolution_start_var, self.dynamic_resolution_min_var, self.dynamic_resolution_max_var),
                (self.dynamic_intensity_start_var, self.dynamic_intensity_min_var, self.dynamic_intensity_max_var),
                (self.dynamic_bkg_a_start_var, self.dynamic_bkg_a_min_var, self.dynamic_bkg_a_max_var),
                (self.dynamic_bkg_b_start_var, self.dynamic_bkg_b_min_var, self.dynamic_bkg_b_max_var),
                (self.dynamic_bkg_c_start_var, self.dynamic_bkg_c_min_var, self.dynamic_bkg_c_max_var),
            )
            if is_dynamic
            else (
                (self.kin_resolution_start_var, self.kin_resolution_min_var, self.kin_resolution_max_var),
                (self.kin_scale_start_var, self.kin_scale_min_var, self.kin_scale_max_var),
                (self.kin_bkg_a_start_var, self.kin_bkg_a_min_var, self.kin_bkg_a_max_var),
                (self.kin_bkg_b_start_var, self.kin_bkg_b_min_var, self.kin_bkg_b_max_var),
                (self.kin_bkg_c_start_var, self.kin_bkg_c_min_var, self.kin_bkg_c_max_var),
            )
        )
        for value_entry, min_entry, max_entry, indicator, variables in zip(
            self.stack_calculation_entries,
            self.stack_calculation_min_entries,
            self.stack_calculation_max_entries,
            self.stack_model_calculation_indicators,
            variable_rows,
        ):
            value_var, min_var, max_var = variables
            value_entry.configure(textvariable=value_var)
            min_entry.configure(textvariable=min_var)
            max_entry.configure(textvariable=max_var)
            indicator.start_var = value_var
            indicator.min_var = min_var
            indicator.max_var = max_var
            self._draw_range_indicator(indicator)
        self._set_widgets_enabled(self.stack_dynamic_calculation_widgets, is_dynamic)

    def _redraw_stack_model_calculation_indicators(self, *_args: object) -> None:
        for indicator in self.stack_model_calculation_indicators:
            self._draw_range_indicator(indicator)

    def _sync_kinematic_substrate_controls(self) -> None:
        enabled = (
            not self.stack_enabled_var.get()
            and self.model_var.get() == "Kinematic"
            and bool(self.kin_substrate_var.get())
        )
        self._set_widgets_enabled(self.kinematic_substrate_widgets, enabled)

    def _sync_optional_controls(self) -> None:
        self._set_widgets_enabled(
            self.strain_widgets,
            not self.stack_enabled_var.get() and bool(self.strain_var.get()),
        )

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

    def _on_kinematic_substrate_changed(self, *_args: object) -> None:
        self._sync_kinematic_substrate_controls()
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

    def _project_document(self) -> dict[str, object]:
        setup = self._build_run_config()
        results = None
        if self.last_update is not None:
            results = fit_update_to_dict(self.last_update)
            results["progress_cost"] = list(self.history_y)
            results["progress_phase"] = list(self.history_phase)
            results["summary"] = self.summary_text.get("1.0", tk.END).strip()
        return {
            "format": PROJECT_FORMAT,
            "version": PROJECT_VERSION,
            "setup": setup,
            "display": {
                "axis": self.axis_var.get(),
                "minimum": self.min_var.get(),
                "maximum": self.max_var.get(),
            },
            "results": results,
        }

    def _report_save_error(self, title: str, exc: Exception) -> None:
        self.status_var.set("Save failed; current setup and results remain available")
        messagebox.showerror(
            title,
            f"{exc}\n\nThe current setup, results, and fit progress remain available in GenL.",
        )

    def _save_project(self) -> None:
        try:
            document = self._project_document()
            data_path = resolve_data_path(self.data_path_var.get())
            selected = filedialog.asksaveasfilename(
                title="Save GenL setup and results",
                initialdir=str(data_path.parent if data_path.parent.exists() else ROOT),
                initialfile=f"{data_path.stem}.genl.json",
                defaultextension=".genl.json",
                filetypes=[("GenL project", "*.genl.json"), ("JSON", "*.json")],
            )
            if not selected:
                return
            write_json_document(Path(selected), document)
        except Exception as exc:  # Save callbacks must not escape into Tk.
            self._report_save_error("Cannot save setup/results", exc)
            return
        self.status_var.set(f"Saved setup/results: {selected}")

    def _save_plot_image(self) -> None:
        if self.last_update is None:
            messagebox.showinfo("No results", "Run a simulation or fit before saving plots.")
            return
        has_density = (
            self.last_update.density_z is not None
            and self.last_update.density_rho_e is not None
        )
        if self.model_var.get() == "Dynamic" and not has_density:
            messagebox.showerror(
                "No electron-density profile",
                "Run a dynamic simulation or fit before saving the two plots.",
            )
            return
        try:
            data_path = resolve_data_path(self.data_path_var.get())
            selected = filedialog.asksaveasfilename(
                title="Save GenL plots",
                initialdir=str(data_path.parent if data_path.parent.exists() else ROOT),
                initialfile=f"{data_path.stem}_genl.png",
                defaultextension=".png",
                filetypes=[
                    ("PNG image", "*.png"),
                    ("JPEG image", "*.jpg *.jpeg"),
                    ("TIFF image", "*.tif *.tiff"),
                    ("SVG image", "*.svg"),
                    ("PDF document", "*.pdf"),
                ],
            )
            if not selected:
                return
            written = save_result_plots(
                self.last_update,
                self._parse_wavelength(),
                Path(selected),
            )
        except Exception as exc:  # Save callbacks must not escape into Tk.
            self._report_save_error("Cannot save plots", exc)
            return
        self.status_var.set("Saved plots: " + ", ".join(str(path) for path in written))
        messagebox.showinfo(
            "Plots saved",
            "\n".join(str(path) for path in written)
            + ("" if has_density else "\nElectron density is unavailable for kinematic results."),
        )

    def _export_graph_data(self) -> None:
        if self.last_update is None:
            messagebox.showinfo("No results", "Run a simulation or fit before exporting data.")
            return
        try:
            data_path = resolve_data_path(self.data_path_var.get())
            selected = filedialog.asksaveasfilename(
                title="Export GenL graph data",
                initialdir=str(data_path.parent if data_path.parent.exists() else ROOT),
                initialfile=f"{data_path.stem}_genl_results.csv",
                defaultextension=".csv",
                filetypes=[("CSV data", "*.csv"), ("Tab-delimited text", "*.txt")],
            )
            if not selected:
                return
            written = export_result_data(self.last_update, Path(selected))
        except Exception as exc:  # Save callbacks must not escape into Tk.
            self._report_save_error("Cannot export graph data", exc)
            return
        self.status_var.set("Exported graph data: " + ", ".join(str(path) for path in written))

    def _apply_project_setup(
        self,
        setup: dict[str, object],
        display: dict[str, object],
    ) -> None:
        def set_triplet(variables: tuple[tk.Variable, tk.Variable, tk.Variable], values: object) -> None:
            values_list = list(values)
            if len(values_list) != 3:
                raise ValueError("Saved parameter ranges must contain value, min, and max")
            for variable, value in zip(variables, values_list):
                variable.set(str(value))

        def set_flags(variables: tuple[tk.BooleanVar, ...], values: object) -> None:
            values_list = list(values)
            if len(values_list) == len(variables) - 1:
                values_list.append(False)
            if len(values_list) != len(variables):
                raise ValueError("Saved fit-selection list has the wrong length")
            for variable, value in zip(variables, values_list):
                variable.set(bool(value))

        sample_name = str(setup["sample_profile"])
        model_name = str(setup["model"])
        axis = str(display.get("axis", "2θ"))
        if sample_name not in SAMPLES:
            raise ValueError(f"Unknown saved sample profile: {sample_name}")
        if model_name not in {"Kinematic", "Dynamic"}:
            raise ValueError(f"Unknown saved scattering model: {model_name}")
        if axis not in {"2θ", "q"}:
            raise ValueError(f"Unknown saved horizontal axis: {axis}")

        self.updating_twotheta_window = True
        try:
            self.sample_var.set(sample_name)
            self.model_var.set(model_name)
            self.data_path_var.set(str(setup["data_path"]))
            self.wavelength_var.set(str(setup["wavelength"]))
            self.density_slices_var.set(str(setup.get("density_slices", 50)))
            self.density_max_q0_var.set(str(setup.get("density_max_q0", 30.0)))
            self.seed_var.set(str(setup["seed"]))
            self.maxiter_var.set(str(setup["maxiter"]))
            self.popsize_var.set(str(setup["popsize"]))
            self.local_var.set(str(setup["local_nfev"]))
            self.polish_var.set(str(setup["polish_iter"]))
            self.interval_var.set(str(setup["interval"]))
            self.strain_var.set(bool(setup["include_strain"]))
            self.roughness_var.set(bool(setup["include_roughness"]))
            self.kin_substrate_var.set(bool(setup.get("include_kinematic_substrate", False)))

            self.axis_mode = "q" if axis == "q" else "twotheta"
            self.axis_var.set(axis)
            min_label, max_label = self._axis_limit_labels()
            self.min_label_var.set(min_label)
            self.max_label_var.set(max_label)
            self.min_var.set(str(display.get("minimum", setup["twotheta_min"])))
            self.max_var.set(str(display.get("maximum", setup["twotheta_max"])))

            kinematic = setup["kinematic_settings"]
            for key, variables in (
                ("d_spacing", (self.kin_d_start_var, self.kin_d_min_var, self.kin_d_max_var)),
                ("planes", (self.kin_planes_start_var, self.kin_planes_min_var, self.kin_planes_max_var)),
                ("resolution", (self.kin_resolution_start_var, self.kin_resolution_min_var, self.kin_resolution_max_var)),
                ("scale", (self.kin_scale_start_var, self.kin_scale_min_var, self.kin_scale_max_var)),
                ("bkg_a", (self.kin_bkg_a_start_var, self.kin_bkg_a_min_var, self.kin_bkg_a_max_var)),
                ("bkg_b", (self.kin_bkg_b_start_var, self.kin_bkg_b_min_var, self.kin_bkg_b_max_var)),
            ):
                set_triplet(variables, kinematic[key])
            set_triplet(
                (
                    self.kin_bkg_c_start_var,
                    self.kin_bkg_c_min_var,
                    self.kin_bkg_c_max_var,
                ),
                kinematic.get("bkg_c", (0.0, -3.0, 3.0)),
            )
            self.kin_debye_var.set(str(kinematic["debye_waller_coeff"]))
            set_flags(
                (
                    self.kin_d_fit_enabled_var,
                    self.kin_planes_fit_enabled_var,
                    self.kin_resolution_fit_enabled_var,
                    self.kin_scale_fit_enabled_var,
                    self.kin_bkg_a_fit_enabled_var,
                    self.kin_bkg_b_fit_enabled_var,
                    self.kin_bkg_c_fit_enabled_var,
                ),
                setup["kinematic_fit_flags"],
            )
            kinematic_substrate = setup.get("kinematic_substrate_settings")
            if kinematic_substrate is not None:
                for key, variables in (
                    (
                        "intensity",
                        (
                            self.kin_substrate_intensity_start_var,
                            self.kin_substrate_intensity_min_var,
                            self.kin_substrate_intensity_max_var,
                        ),
                    ),
                    (
                        "width",
                        (
                            self.kin_substrate_width_start_var,
                            self.kin_substrate_width_min_var,
                            self.kin_substrate_width_max_var,
                        ),
                    ),
                    (
                        "d_spacing",
                        (
                            self.kin_substrate_d_start_var,
                            self.kin_substrate_d_min_var,
                            self.kin_substrate_d_max_var,
                        ),
                    ),
                ):
                    set_triplet(variables, kinematic_substrate[key])
            set_flags(
                (
                    self.kin_substrate_intensity_fit_enabled_var,
                    self.kin_substrate_width_fit_enabled_var,
                    self.kin_substrate_d_fit_enabled_var,
                ),
                setup.get("kinematic_substrate_fit_flags", (True, True, True)),
            )

            film = setup["film_settings"]
            self.film_filename_var.set(str(film["filename"]))
            self.film_direction_var.set(str(film["direction"]))
            for key, variables in (
                ("n", (self.film_n_start_var, self.film_n_min_var, self.film_n_max_var)),
                ("scale", (self.film_scale_start_var, self.film_scale_min_var, self.film_scale_max_var)),
                ("area_scale", (self.film_area_start_var, self.film_area_min_var, self.film_area_max_var)),
                ("dinterface", (self.film_interface_start_var, self.film_interface_min_var, self.film_interface_max_var)),
            ):
                set_triplet(variables, film[key])
            set_flags(
                (
                    self.film_n_fit_enabled_var,
                    self.film_scale_fit_enabled_var,
                    self.film_area_fit_enabled_var,
                    self.film_interface_fit_enabled_var,
                ),
                setup["film_fit_flags"],
            )

            dynamic = setup["dynamic_fit_settings"]
            for key, variables in (
                ("resolution", (self.dynamic_resolution_start_var, self.dynamic_resolution_min_var, self.dynamic_resolution_max_var)),
                ("intensity_scale", (self.dynamic_intensity_start_var, self.dynamic_intensity_min_var, self.dynamic_intensity_max_var)),
                ("bkg_a", (self.dynamic_bkg_a_start_var, self.dynamic_bkg_a_min_var, self.dynamic_bkg_a_max_var)),
                ("bkg_b", (self.dynamic_bkg_b_start_var, self.dynamic_bkg_b_min_var, self.dynamic_bkg_b_max_var)),
            ):
                set_triplet(variables, dynamic[key])
            set_triplet(
                (
                    self.dynamic_bkg_c_start_var,
                    self.dynamic_bkg_c_min_var,
                    self.dynamic_bkg_c_max_var,
                ),
                dynamic.get("bkg_c", (0.0, -0.1, 0.1)),
            )
            set_flags(
                (
                    self.dynamic_resolution_fit_enabled_var,
                    self.dynamic_intensity_fit_enabled_var,
                    self.dynamic_bkg_a_fit_enabled_var,
                    self.dynamic_bkg_b_fit_enabled_var,
                    self.dynamic_bkg_c_fit_enabled_var,
                ),
                setup["dynamic_fit_flags"],
            )

            substrate = setup["substrate_settings"]
            self.substrate_filename_var.set(str(substrate["filename"]))
            self.substrate_direction_var.set(str(substrate["direction"]))
            self.substrate_n_var.set(str(substrate["n"]))
            self.substrate_interface_var.set(str(substrate["dinterface"]))
            set_triplet(
                (self.substrate_scale_var, self.substrate_scale_min_var, self.substrate_scale_max_var),
                substrate["scale"],
            )
            self.substrate_area_var.set(str(substrate["area_scale"]))
            self.substrate_scale_fit_enabled_var.set(bool(setup["substrate_scale_fit_flag"]))

            strain = setup["strain_settings"]
            for index, key in enumerate(STRAIN_KEYS):
                set_triplet(
                    (self.strain_start_vars[index], self.strain_min_vars[index], self.strain_max_vars[index]),
                    strain[key],
                )
            set_flags(tuple(self.strain_fit_enabled_vars), setup["strain_fit_flags"])
            self._store_strain_values()

            roughness = setup["roughness_settings"]
            set_triplet(
                (self.film_rough_start_var, self.film_rough_min_var, self.film_rough_max_var),
                roughness["film"],
            )
            set_triplet(
                (self.substrate_rough_start_var, self.substrate_rough_min_var, self.substrate_rough_max_var),
                roughness["substrate"],
            )
            set_flags(
                (self.film_rough_fit_enabled_var, self.substrate_rough_fit_enabled_var),
                setup["roughness_fit_flags"],
            )
            stack_document = setup.get("stack_document")
            if isinstance(stack_document, dict):
                self.stack_path_var.set(str(setup.get("stack_path", self.stack_path_var.get())))
                self._set_stack_document(stack_document, sync_calculation=False)
            self.stack_points_per_unit_var.set(str(setup.get("stack_points_per_unit", 50)))
            stack_enabled = bool(setup.get("stack_enabled", False))
            self.superlattice_data_preview = stack_enabled
            self.stack_enabled_var.set(stack_enabled)
        finally:
            self.updating_twotheta_window = False

        self.preview_data_path = resolve_data_path(self.data_path_var.get())
        self._sync_model_controls()
        self._sync_optional_controls()
        self._clear_fit_result_values()
        self._redraw_range_indicators()

    def _load_project(self) -> None:
        if self.running:
            messagebox.showinfo("Fit running", "Stop the fit before loading another setup.")
            return
        selected = filedialog.askopenfilename(
            title="Load GenL setup and results",
            initialdir=str(ROOT),
            filetypes=[("GenL project", "*.genl.json"), ("JSON", "*.json"), ("All files", "*")],
        )
        if not selected:
            return
        try:
            with Path(selected).open("r", encoding="utf-8") as handle:
                document = json.load(handle)
            if not isinstance(document, dict) or document.get("format") != PROJECT_FORMAT:
                raise ValueError("The selected file is not a GenL GUI project")
            if document.get("version") != PROJECT_VERSION:
                raise ValueError(f"Unsupported GenL project version: {document.get('version')}")
            setup = document["setup"]
            display = document.get("display", {})
            if not isinstance(setup, dict) or not isinstance(display, dict):
                raise ValueError("The saved GenL setup is invalid")
            results = document.get("results")
            update = None
            saved_history: list[float] = []
            saved_phases: list[str] = []
            if results is not None:
                if not isinstance(results, dict):
                    raise ValueError("The saved GenL results are invalid")
                update = fit_update_from_dict(results)
                saved_history = [float(value) for value in results.get("progress_cost", [])]
                saved_phases = [str(value) for value in results.get("progress_phase", [])]

            if self.preview_after_id is not None:
                self.root.after_cancel(self.preview_after_id)
                self.preview_after_id = None
            self._apply_project_setup(setup, display)
            self.history_x.clear()
            self.history_y.clear()
            self.history_phase.clear()
            self.last_update = None
            self.fit_checkpoint = None
            self.pending_resume_config = None
            self.pause_button.configure(state=tk.DISABLED, text="Pause fit")
            self.summary_text.delete("1.0", tk.END)

            if update is not None:
                self.history_y[:] = saved_history[:-1] if saved_history else []
                self.history_x[:] = list(range(1, len(self.history_y) + 1))
                self.history_phase[:] = saved_phases[: len(self.history_y)]
                if len(self.history_phase) < len(self.history_y):
                    missing = len(self.history_y) - len(self.history_phase)
                    self.history_phase.extend([update.phase] * missing)
                self._apply_update(update)
                self.summary_text.insert(tk.END, str(results.get("summary", "")))
            else:
                self._draw_experimental_preview()
            self.status_var.set(f"Loaded setup/results: {selected}")
        except (OSError, TypeError, ValueError, KeyError, json.JSONDecodeError) as exc:
            messagebox.showerror("Cannot load setup/results", str(exc))

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
            self.stack_sampling_label_var.set(
                "Points per q unit" if self._axis_is_q() else "Points per degree"
            )
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
        self.kin_bkg_a_min_var.set("-3.0")
        self.kin_bkg_a_max_var.set("3.0")
        self.kin_bkg_b_start_var.set("0.1")
        self.kin_bkg_b_min_var.set("0.0")
        self.kin_bkg_b_max_var.set("3.0")
        self.kin_bkg_c_start_var.set("0.0")
        self.kin_bkg_c_min_var.set("-3.0")
        self.kin_bkg_c_max_var.set("3.0")
        self.kin_debye_var.set(f"{sample.debye_waller_coeff:.6g}")
        self.kin_substrate_d_start_var.set(f"{d0:.6g}")
        self.kin_substrate_d_min_var.set(f"{d0 * 0.995:.6g}")
        self.kin_substrate_d_max_var.set(f"{d0 * 1.005:.6g}")

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
        self.superlattice_data_preview = False
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
        if self.stack_enabled_var.get() and not self.superlattice_data_preview:
            self._draw_empty_plot()
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
            self.kin_bkg_c_fit_var,
            self.kin_debye_fit_var,
            self.kin_substrate_intensity_fit_var,
            self.kin_substrate_width_fit_var,
            self.kin_substrate_d_fit_var,
            self.film_n_fit_var,
            self.film_scale_fit_var,
            self.film_area_fit_var,
            self.film_interface_fit_var,
            self.dynamic_resolution_fit_var,
            self.dynamic_intensity_fit_var,
            self.dynamic_bkg_a_fit_var,
            self.dynamic_bkg_b_fit_var,
            self.dynamic_bkg_c_fit_var,
            self.substrate_scale_fit_var,
            *self.strain_fit_vars,
            self.film_rough_fit_var,
            self.substrate_rough_fit_var,
        ):
            variable.set("")
        for variables in self.stack_row_vars:
            for key in (
                "dinterface",
                "scale",
                "area_scale",
                *(field[0] for field in STACK_STRAIN_FIELDS),
            ):
                variables[f"{key}_fit"].set("")
        for variable in self.stack_calculation_fit_vars.values():
            variable.set("")

    def _set_fit_result_values(self, config: dict[str, object], params: np.ndarray) -> None:
        self._clear_fit_result_values()
        if bool(config.get("stack_enabled", False)):
            definition = StackDefinition(
                Path(str(config["stack_path"])),
                copy.deepcopy(config["stack_document"]),
            )
            values = definition.overrides(params)
            for role, variables in zip(self.stack_row_roles, self.stack_row_vars):
                prefix = (
                    role
                    if role in {"substrate", "capping"}
                    else variables["name"].get().strip()
                )
                keys = ["dinterface", "scale", "area_scale"]
                if role != "substrate":
                    keys.extend(field[0] for field in STACK_STRAIN_FIELDS)
                for key in keys:
                    target = f"{prefix}.{key}"
                    if target in values:
                        result = self._format_fit_result(values[target])
                        variables[key].set(result)
                        variables[f"{key}_fit"].set(result)
                    else:
                        variables[f"{key}_fit"].set("off")
            calculation_values = {
                key.removeprefix("calculation."): value
                for key, value in values.items()
                if key.startswith("calculation.")
            }
            calculation_start_vars = (
                {
                    "resolution": self.dynamic_resolution_start_var,
                    "intensity_scale": self.dynamic_intensity_start_var,
                    "background_a": self.dynamic_bkg_a_start_var,
                    "background_b": self.dynamic_bkg_b_start_var,
                    "background_c": self.dynamic_bkg_c_start_var,
                }
                if config["model"] == "Dynamic"
                else {
                    "resolution": self.kin_resolution_start_var,
                    "intensity_scale": self.kin_scale_start_var,
                    "background_a": self.kin_bkg_a_start_var,
                    "background_b": self.kin_bkg_b_start_var,
                    "background_c": self.kin_bkg_c_start_var,
                }
            )
            for key, variable in calculation_start_vars.items():
                if key in calculation_values:
                    result = self._format_fit_result(calculation_values[key])
                    variable.set(result)
                    self.stack_calculation_fit_vars[key].set(result)
                else:
                    self.stack_calculation_fit_vars[key].set("off")
            self._redraw_range_indicators()
            return
        if config["model"] == "Kinematic":
            self.kin_d_start_var.set(self._format_fit_result(float(params[0])))
            self.kin_planes_start_var.set(self._format_fit_result(float(params[1])))
            self.kin_resolution_start_var.set(self._format_fit_result(float(params[2])))
            self.kin_scale_start_var.set(self._format_fit_result(float(params[3])))
            self.kin_bkg_a_start_var.set(self._format_fit_result(float(params[4])))
            self.kin_bkg_b_start_var.set(self._format_fit_result(float(params[5])))
            self.kin_bkg_c_start_var.set(self._format_fit_result(float(params[6])))
            self.kin_d_fit_var.set(self.kin_d_start_var.get())
            self.kin_planes_fit_var.set(self.kin_planes_start_var.get())
            self.kin_resolution_fit_var.set(self.kin_resolution_start_var.get())
            self.kin_scale_fit_var.set(self.kin_scale_start_var.get())
            self.kin_bkg_a_fit_var.set(self.kin_bkg_a_start_var.get())
            self.kin_bkg_b_fit_var.set(self.kin_bkg_b_start_var.get())
            self.kin_bkg_c_fit_var.set(self.kin_bkg_c_start_var.get())
            self.kin_debye_fit_var.set(self.kin_debye_var.get())
            offset = 7
            if bool(config["include_kinematic_substrate"]):
                self.kin_substrate_intensity_start_var.set(
                    self._format_fit_result(float(params[offset]))
                )
                self.kin_substrate_width_start_var.set(
                    self._format_fit_result(float(params[offset + 1]))
                )
                self.kin_substrate_d_start_var.set(
                    self._format_fit_result(float(params[offset + 2]))
                )
                self.kin_substrate_intensity_fit_var.set(
                    self.kin_substrate_intensity_start_var.get()
                )
                self.kin_substrate_width_fit_var.set(
                    self.kin_substrate_width_start_var.get()
                )
                self.kin_substrate_d_fit_var.set(self.kin_substrate_d_start_var.get())
                offset += 3
            else:
                self.kin_substrate_intensity_fit_var.set("off")
                self.kin_substrate_width_fit_var.set("off")
                self.kin_substrate_d_fit_var.set("off")
        else:
            self.film_n_start_var.set(self._format_fit_result(float(params[0])))
            self.film_scale_start_var.set(self._format_fit_result(float(params[1])))
            self.film_area_start_var.set(self._format_fit_result(float(params[2])))
            self.film_interface_start_var.set(self._format_fit_result(float(params[3])))
            self.dynamic_resolution_start_var.set(self._format_fit_result(float(params[4])))
            self.dynamic_intensity_start_var.set(self._format_fit_result(float(params[5])))
            self.dynamic_bkg_a_start_var.set(self._format_fit_result(float(params[6])))
            self.dynamic_bkg_b_start_var.set(self._format_fit_result(float(params[7])))
            self.dynamic_bkg_c_start_var.set(self._format_fit_result(float(params[8])))
            self.film_n_fit_var.set(self.film_n_start_var.get())
            self.film_scale_fit_var.set(self.film_scale_start_var.get())
            self.film_area_fit_var.set(self.film_area_start_var.get())
            self.film_interface_fit_var.set(self.film_interface_start_var.get())
            self.dynamic_resolution_fit_var.set(self.dynamic_resolution_start_var.get())
            self.dynamic_intensity_fit_var.set(self.dynamic_intensity_start_var.get())
            self.dynamic_bkg_a_fit_var.set(self.dynamic_bkg_a_start_var.get())
            self.dynamic_bkg_b_fit_var.set(self.dynamic_bkg_b_start_var.get())
            self.dynamic_bkg_c_fit_var.set(self.dynamic_bkg_c_start_var.get())
            self.substrate_scale_var.set(self._format_fit_result(float(params[9])))
            self.substrate_scale_fit_var.set(self.substrate_scale_var.get())
            offset = 10

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
        if bool(config.get("stack_enabled", False)):
            definition = StackDefinition(
                Path(str(config["stack_path"])),
                copy.deepcopy(config["stack_document"]),
            )
            warnings = []
            for name, value, (lower, upper) in zip(
                definition.parameter_names, params, definition.bounds
            ):
                ratio = (float(value) - lower) / (upper - lower)
                if ratio <= 0.02:
                    warnings.append(f"{name} reached the lower bound region")
                elif ratio >= 0.98:
                    warnings.append(f"{name} reached the upper bound region")
                elif ratio <= 0.05:
                    warnings.append(f"{name} is close to the lower bound")
                elif ratio >= 0.95:
                    warnings.append(f"{name} is close to the upper bound")
            return warnings
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
                ("bkg c", params[6], settings["bkg_c"]),
            )
            rows.extend(
                (name, float(value), float(bounds[1]), float(bounds[2]))
                for (name, value, bounds), enabled in zip(base, config["kinematic_fit_flags"])
                if enabled
            )
            offset = 7
            if bool(config["include_kinematic_substrate"]):
                settings = config["kinematic_substrate_settings"]
                substrate_rows = (
                    ("substrate integrated intensity", params[offset], settings["intensity"]),
                    ("substrate FWHM", params[offset + 1], settings["width"]),
                    ("substrate d spacing", params[offset + 2], settings["d_spacing"]),
                )
                rows.extend(
                    (name, float(value), float(bounds[1]), float(bounds[2]))
                    for (name, value, bounds), enabled in zip(
                        substrate_rows, config["kinematic_substrate_fit_flags"]
                    )
                    if enabled
                )
                offset += 3
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
                ("dynamic bkg c", params[8], settings["bkg_c"]),
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
                        float(params[9]),
                        float(bounds[1]),
                        float(bounds[2]),
                    )
                )
            offset = 10

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

    def _build_run_config(self, allow_outside_start: bool = False) -> dict[str, object]:
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
            bool(self.kin_bkg_c_fit_enabled_var.get()),
        )
        kinematic_substrate_fit_flags = (
            bool(self.kin_substrate_intensity_fit_enabled_var.get()),
            bool(self.kin_substrate_width_fit_enabled_var.get()),
            bool(self.kin_substrate_d_fit_enabled_var.get()),
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
            bool(self.dynamic_bkg_c_fit_enabled_var.get()),
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
        kinematic_bkg_c = (
            float(self.kin_bkg_c_start_var.get()),
            float(self.kin_bkg_c_min_var.get()),
            float(self.kin_bkg_c_max_var.get()),
        )
        kinematic_debye = float(self.kin_debye_var.get())
        kinematic_substrate_intensity = (
            float(self.kin_substrate_intensity_start_var.get()),
            float(self.kin_substrate_intensity_min_var.get()),
            float(self.kin_substrate_intensity_max_var.get()),
        )
        kinematic_substrate_width = (
            float(self.kin_substrate_width_start_var.get()),
            float(self.kin_substrate_width_min_var.get()),
            float(self.kin_substrate_width_max_var.get()),
        )
        kinematic_substrate_d = (
            float(self.kin_substrate_d_start_var.get()),
            float(self.kin_substrate_d_min_var.get()),
            float(self.kin_substrate_d_max_var.get()),
        )
        for enabled, name, values in zip(
            kinematic_fit_flags,
            (
                "kinematic d spacing",
                "kinematic planes",
                "kinematic resolution",
                "kinematic scale",
                "kinematic bkg a",
                "kinematic bkg b",
                "kinematic bkg c",
            ),
            (
                kinematic_d,
                kinematic_planes,
                kinematic_resolution,
                kinematic_scale,
                kinematic_bkg_a,
                kinematic_bkg_b,
                kinematic_bkg_c,
            ),
        ):
            if enabled:
                validate_start_min_max(
                    name, *values, allow_outside_start=allow_outside_start
                )
        if kinematic_d[0] <= 0 or (kinematic_fit_flags[0] and kinematic_d[1] <= 0):
            raise ValueError("kinematic d spacing must be positive")
        if kinematic_planes[0] <= 0 or (kinematic_fit_flags[1] and kinematic_planes[1] <= 0):
            raise ValueError("kinematic planes must be positive")
        if kinematic_resolution[0] <= 0 or (kinematic_fit_flags[2] and kinematic_resolution[1] <= 0):
            raise ValueError("kinematic resolution must be positive")
        if kinematic_scale[0] <= 0 or (kinematic_fit_flags[3] and kinematic_scale[1] <= 0):
            raise ValueError("kinematic scale must be positive")
        include_kinematic_substrate = (
            self.model_var.get() == "Kinematic" and bool(self.kin_substrate_var.get())
        )
        if include_kinematic_substrate:
            for enabled, name, values in zip(
                kinematic_substrate_fit_flags,
                (
                    "kinematic substrate integrated intensity",
                    "kinematic substrate FWHM",
                    "kinematic substrate d spacing",
                ),
                (
                    kinematic_substrate_intensity,
                    kinematic_substrate_width,
                    kinematic_substrate_d,
                ),
            ):
                if enabled:
                    validate_start_min_max(
                        name, *values, allow_outside_start=allow_outside_start
                    )
            if kinematic_substrate_intensity[0] < 0 or (
                kinematic_substrate_fit_flags[0] and kinematic_substrate_intensity[1] < 0
            ):
                raise ValueError("kinematic substrate intensity must be non-negative")
            if kinematic_substrate_width[0] <= 0 or (
                kinematic_substrate_fit_flags[1] and kinematic_substrate_width[1] <= 0
            ):
                raise ValueError("kinematic substrate FWHM must be positive")
            minimum_d = wavelength / 2.0
            if kinematic_substrate_d[0] < minimum_d or (
                kinematic_substrate_fit_flags[2] and kinematic_substrate_d[1] < minimum_d
            ):
                raise ValueError(
                    "kinematic substrate d spacing must be at least half the X-ray wavelength"
                )

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
                validate_start_min_max(
                    name, *values, allow_outside_start=allow_outside_start
                )
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
        dynamic_bkg_c = (
            float(self.dynamic_bkg_c_start_var.get()),
            float(self.dynamic_bkg_c_min_var.get()),
            float(self.dynamic_bkg_c_max_var.get()),
        )
        for enabled, name, values in zip(
            dynamic_fit_flags,
            (
                "dynamic resolution",
                "dynamic intensity scale",
                "dynamic bkg a",
                "dynamic bkg b",
                "dynamic bkg c",
            ),
            (
                dynamic_resolution,
                dynamic_intensity_scale,
                dynamic_bkg_a,
                dynamic_bkg_b,
                dynamic_bkg_c,
            ),
        ):
            if enabled:
                validate_start_min_max(
                    name, *values, allow_outside_start=allow_outside_start
                )
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
            validate_start_min_max(
                "substrate lattice scale",
                *substrate_scale,
                allow_outside_start=allow_outside_start,
            )
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
                    validate_start_min_max(
                        key.replace("_", " "),
                        *strain_settings[key],
                        allow_outside_start=allow_outside_start,
                    )
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
            validate_start_min_max(
                "film roughness",
                *film_roughness,
                allow_outside_start=allow_outside_start,
            )
        if roughness_fit_flags[1]:
            validate_start_min_max(
                "substrate/interface roughness",
                *substrate_roughness,
                allow_outside_start=allow_outside_start,
            )
        twotheta_min, twotheta_max = self._window_twotheta_limits()
        density_slices = int(self.density_slices_var.get())
        density_max_q0 = float(self.density_max_q0_var.get())
        sampling_scale = (
            max(substrate_scale[0], substrate_scale[2])
            if substrate_scale_fit_flag
            else substrate_scale[0]
        )
        substrate_period = (
            propagation_period_for(substrate_filename, substrate_direction)
            * sampling_scale
        )
        validate_density_sampling(substrate_period, density_slices, density_max_q0)
        stack_enabled = bool(self.stack_enabled_var.get())
        try:
            stack_points_per_unit = float(self.stack_points_per_unit_var.get())
        except ValueError:
            if stack_enabled:
                raise ValueError(
                    "Superlattice simulation points per axis unit must be a number"
                ) from None
            stack_points_per_unit = 50.0
        stack_document = (
            self._stack_document_from_controls(allow_outside_start)
            if stack_enabled
            else copy.deepcopy(self.stack_document)
        )
        try:
            seed = int(self.seed_var.get())
        except ValueError:
            raise ValueError("optimizer seed must be an integer") from None
        if not 0 <= seed < 2**32:
            raise ValueError("optimizer seed must be between 0 and 4294967295")
        return {
            "sample_profile": self.sample_var.get(),
            "data_path": str(data_path),
            "model": self.model_var.get(),
            "wavelength": wavelength,
            "density_slices": density_slices,
            "density_max_q0": density_max_q0,
            "twotheta_min": twotheta_min,
            "twotheta_max": twotheta_max,
            "seed": seed,
            "maxiter": int(self.maxiter_var.get()),
            "popsize": int(self.popsize_var.get()),
            "local_nfev": int(self.local_var.get()),
            "polish_iter": int(self.polish_var.get()),
            "interval": max(1, int(self.interval_var.get())),
            "include_strain": include_strain,
            "include_roughness": bool(self.roughness_var.get()),
            "include_kinematic_substrate": include_kinematic_substrate,
            "kinematic_fit_flags": kinematic_fit_flags,
            "kinematic_substrate_fit_flags": kinematic_substrate_fit_flags,
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
                "bkg_c": kinematic_bkg_c,
                "debye_waller_coeff": kinematic_debye,
            },
            "kinematic_substrate_settings": {
                "intensity": kinematic_substrate_intensity,
                "width": kinematic_substrate_width,
                "d_spacing": kinematic_substrate_d,
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
                "bkg_c": dynamic_bkg_c,
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
            "stack_enabled": stack_enabled,
            "stack_path": self.stack_path_var.get(),
            "stack_document": stack_document,
            "stack_points_per_unit": stack_points_per_unit,
            "stack_axis_min": float(self.min_var.get()),
            "stack_axis_max": float(self.max_var.get()),
            "stack_q_axis": self._axis_is_q(),
        }

    def _make_model_and_start(
        self, config: dict[str, object]
    ) -> tuple[SampleConfig, KinematicModel | DynamicModel | StackModel, np.ndarray, np.ndarray, dict[str, object], dict[str, object]]:
        sample = SAMPLES[str(config["sample_profile"])]
        kinematic_settings = config["kinematic_settings"]
        film_settings = config["film_settings"]
        dynamic_fit_settings = config["dynamic_fit_settings"]
        substrate_settings = config["substrate_settings"]
        twotheta, observed = load_sample_data(
            float(config["twotheta_min"]), float(config["twotheta_max"]),
            str(config["data_path"]),
        )
        if bool(config.get("stack_enabled", False)):
            definition = StackDefinition(
                Path(str(config["stack_path"])),
                copy.deepcopy(config["stack_document"]),
            )
            model = StackModel(
                definition,
                twotheta,
                observed,
                model=str(config["model"]).lower(),
            )
            return sample, model, definition.bounds, definition.start, {}, {}
        if config["model"] == "Kinematic":
            model = KinematicModel(
                twotheta,
                observed,
                sample,
                debye_waller_coeff=float(kinematic_settings["debye_waller_coeff"]),
                include_substrate=bool(config["include_kinematic_substrate"]),
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
                include_substrate=bool(config["include_kinematic_substrate"]),
                substrate_peak_settings=config["kinematic_substrate_settings"],
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
                propagation_backend="reflection",
                density_slices=int(config["density_slices"]),
                density_max_q0=float(config["density_max_q0"]),
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
        if bool(config.get("stack_enabled", False)):
            if n_params == 0:
                raise ValueError("Select at least one superlattice parameter to fit")
            return np.ones(n_params, dtype=bool)
        if config["model"] == "Kinematic":
            flags = list(config["kinematic_fit_flags"])
            if bool(config["include_kinematic_substrate"]):
                flags.extend(list(config["kinematic_substrate_fit_flags"]))
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
        if bool(config.get("stack_enabled", False)):
            return bool(config["stack_document"].get("fit_parameters"))
        if config["model"] == "Kinematic":
            flags = list(config["kinematic_fit_flags"])
            if bool(config["include_kinematic_substrate"]):
                flags.extend(list(config["kinematic_substrate_fit_flags"]))
        else:
            flags = list(config["film_fit_flags"]) + list(config["dynamic_fit_flags"])
            flags.append(bool(config["substrate_scale_fit_flag"]))
        if bool(config["include_strain"]):
            flags.extend(list(config["strain_fit_flags"]))
        if bool(config["include_roughness"]):
            flags.extend(list(config["roughness_fit_flags"]))
        return any(flags)

    def _simulate_stack(self, config: dict[str, object]) -> None:
        if self.preview_after_id is not None:
            self.root.after_cancel(self.preview_after_id)
            self.preview_after_id = None
        document = copy.deepcopy(config["stack_document"])
        definition = StackDefinition(Path(str(config["stack_path"])), document)
        show_observed = self.superlattice_data_preview
        if show_observed:
            twotheta_min, twotheta_max = self._window_twotheta_limits()
            twotheta, observed = load_sample_data(
                twotheta_min, twotheta_max, str(config["data_path"])
            )
            q = np.asarray(q_from_twotheta(twotheta, float(config["wavelength"])))
        else:
            twotheta, q = stack_simulation_grid(
                float(config["stack_axis_min"]),
                float(config["stack_axis_max"]),
                float(config["stack_points_per_unit"]),
                bool(config["stack_q_axis"]),
                float(config["wavelength"]),
            )
            observed = np.zeros_like(twotheta)
        model = StackModel(definition, twotheta, model=str(config["model"]).lower())
        predicted = model.predict(definition.start)
        cost = 0.0
        if show_observed:
            positive = observed[observed > 0]
            if not len(positive):
                raise ValueError("Positive experimental intensities are required")
            floor = max(float(np.min(positive)) * 0.1, 1e-12)
            cost = float(
                np.mean(
                    np.abs(
                        np.log10(np.maximum(predicted, floor))
                        - np.log10(np.maximum(observed, floor))
                    )
                )
            )
        dynamic_result = model.last_dynamic_result
        self.stack_document = document
        self.history_x.clear()
        self.history_y.clear()
        self.history_phase.clear()
        self._clear_fit_result_values()
        self.summary_text.delete("1.0", tk.END)
        self._apply_update(
            FitUpdate(
                phase="simulation",
                cost=cost,
                twotheta=model.twotheta,
                q=q,
                observed=observed,
                predicted=predicted,
                params=definition.start,
                density_z=None if dynamic_result is None else dynamic_result.z,
                density_rho_e=None if dynamic_result is None else dynamic_result.rho_e,
                show_observed=show_observed,
            )
        )
        sequence = document["sequence"]
        layer_text = " / ".join(
            f"{layer['name']} ({float(layer['unit_cells']):g} cells)"
            for layer in sequence["layers"]
        )
        capping_layer = document.get("capping_layer")
        capping_text = (
            ""
            if capping_layer is None
            else f"Capping layer: {capping_layer['name']} "
            f"({float(capping_layer['unit_cells']):g} cells)\n"
        )
        summary = (
            "Superlattice simulation\n"
            f"Superlattice: {definition.name}\n"
            f"Model: {config['model']}\n"
            f"Sequence: [{layer_text}] x {sequence['repetitions']}\n"
            f"{capping_text}"
            f"Expanded layers including substrate: {len(definition.layers())}\n"
            f"X-ray wavelength: {float(config['wavelength']):.6g} A\n"
            f"Calculated points: {len(twotheta)} "
            + (
                "(experimental data grid)"
                if show_observed
                else f"({float(config['stack_points_per_unit']):g} per "
                f"{'q unit' if config['stack_q_axis'] else 'degree'})"
            )
        )
        self.summary_text.insert(tk.END, summary)
        self.status_var.set(
            f"Superlattice simulation complete: cost={cost:.5g}"
            if show_observed
            else f"Superlattice simulation complete: {len(twotheta)} points"
        )

    def simulate_pattern(self) -> None:
        if self.running:
            messagebox.showinfo("Fit running", "Stop the fit before simulating a new pattern.")
            return

        try:
            config = self._build_run_config()
            if config["stack_enabled"]:
                self._simulate_stack(config)
                return
            sample, model, _bounds_array, start, film_settings, substrate_settings = self._make_model_and_start(config)
            predicted = model.predict(start)
            cost = model.objective(start)
            residual = predicted - model.observed
            rmse = float(np.sqrt(np.mean(residual**2)))
            density_z = getattr(model, "last_z", None)
            density_rho_e = getattr(model, "last_rho_e", None)
            self.history_x.clear()
            self.history_y.clear()
            self.history_phase.clear()
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
                bool(config["include_kinematic_substrate"]),
                film_settings if config["model"] == "Dynamic" else None,
                substrate_settings if config["model"] == "Dynamic" else None,
            )
            backend_text = (
                f"Density sampling: {config['density_slices']} slices/cell, "
                f"Q max {float(config['density_max_q0']):.6g} 1/A\n"
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
        if self.stack_enabled_var.get() and not self.superlattice_data_preview:
            messagebox.showerror(
                "Superlattice data required",
                "Load experimental data in the superlattice Structure tab before fitting.",
            )
            return

        try:
            config = self._build_run_config()
        except ValueError as exc:
            messagebox.showerror("Invalid input", str(exc))
            return
        if not self._has_selected_fit_parameter(config):
            messagebox.showerror("Invalid input", "Select at least one parameter to fit.")
            return

        self.fit_checkpoint = None
        self.pending_resume_config = None
        self._start_fit(config, None, clear_history=True)

    def _start_fit(
        self,
        config: dict[str, object],
        checkpoint: FitCheckpoint | None,
        clear_history: bool,
    ) -> None:
        self.running = True
        self.paused = False
        self.active_fit_config = config
        self.stop_event.clear()
        self.pause_event.clear()
        self._clear_fit_result_values()
        self.simulate_button.configure(state=tk.DISABLED)
        self.run_button.configure(state=tk.DISABLED)
        self.pause_button.configure(state=tk.NORMAL, text="Pause fit")
        self.stop_button.configure(state=tk.NORMAL)
        if clear_history:
            self.history_x.clear()
            self.history_y.clear()
            self.history_phase.clear()
            self.summary_text.delete("1.0", tk.END)
            self._draw_empty_plot()
        self.status_var.set("Resuming fit from checkpoint..." if checkpoint else "Running fit...")

        thread = threading.Thread(
            target=self._fit_worker,
            args=(config, checkpoint),
            daemon=True,
        )
        thread.start()

    def _checkpoint_context(
        self,
        config: dict[str, object],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[object, ...]]:
        _sample, _model, bounds_array, start, _film, _substrate = self._make_model_and_start(config)
        fit_mask = self._fit_mask_for_config(config, len(start))
        return bounds_array, start, fit_mask, fit_resume_signature(config, start, fit_mask)

    def _resume_checkpoint(self, config: dict[str, object] | None = None) -> None:
        checkpoint = self.fit_checkpoint
        if checkpoint is None:
            messagebox.showinfo("No checkpoint", "No interrupted fit checkpoint is available.")
            return
        try:
            config = self._build_run_config(allow_outside_start=True) if config is None else config
            if not self._has_selected_fit_parameter(config):
                raise ValueError("Select at least one parameter to fit")
            _bounds, _start, fit_mask, signature = self._checkpoint_context(config)
            if signature != checkpoint.signature or not np.array_equal(fit_mask, checkpoint.fit_mask):
                raise ValueError(
                    "The scattering model, data window, fixed inputs, or fitted parameter selection changed. "
                    "Start a new fit for those changes."
                )
        except ValueError as exc:
            messagebox.showerror("Cannot resume checkpoint", str(exc))
            return
        self.pending_resume_config = None
        self._start_fit(config, checkpoint, clear_history=False)

    def pause_or_resume_fit(self) -> None:
        if not self.running:
            self._resume_checkpoint()
            return
        if not self.paused:
            self.paused = True
            self.pause_event.set()
            self.pause_button.configure(text="Resume fit")
            self.status_var.set("Pause requested; waiting for the current evaluation...")
            return

        try:
            current_config = self._build_run_config(allow_outside_start=True)
            current_bounds, _start, current_mask, signature = self._checkpoint_context(current_config)
            active_bounds, _active_start, active_mask, _signature = self._checkpoint_context(
                self.active_fit_config
            )
            checkpoint = self.fit_checkpoint
            if checkpoint is None or signature != checkpoint.signature:
                raise ValueError(
                    "Only parameter bounds and optimizer budgets can change during checkpoint resume"
                )
            if not np.array_equal(current_mask, active_mask):
                raise ValueError("The fitted parameter selection changed; start a new fit instead")
            bounds_changed = not np.array_equal(
                current_bounds[current_mask], active_bounds[active_mask]
            )
            budget_keys = ("maxiter", "local_nfev", "polish_iter")
            budget_changed = any(
                current_config[key] != self.active_fit_config[key] for key in budget_keys
            )
        except ValueError as exc:
            messagebox.showerror("Cannot resume fit", str(exc))
            return

        if bounds_changed or budget_changed:
            self.pending_resume_config = current_config
            self.paused = False
            self.pause_button.configure(state=tk.DISABLED, text="Resume fit")
            self.status_var.set("Restarting from checkpoint with updated bounds...")
            self.stop_event.set()
            self.pause_event.clear()
        else:
            self.paused = False
            self.pause_event.clear()
            self.pause_button.configure(text="Pause fit")
            self.status_var.set("Fit resumed")

    def stop_fit(self) -> None:
        if self.running:
            self.pending_resume_config = None
            self.paused = False
            self.stop_event.set()
            self.pause_event.clear()
            self.pause_button.configure(state=tk.DISABLED, text="Resume fit")
            self.status_var.set("Stopping after current optimizer evaluation...")

    def _fit_worker(
        self,
        config: dict[str, object],
        resume_checkpoint: FitCheckpoint | None = None,
    ) -> None:
        try:
            sample, model, bounds_array, start, film_settings, substrate_settings = self._make_model_and_start(config)
            fit_mask = self._fit_mask_for_config(config, len(start))
            fitted_bounds_array = bounds_array[fit_mask]
            bounds = [tuple(row) for row in fitted_bounds_array]
            counter = {"value": 0}
            signature = fit_resume_signature(config, start, fit_mask)
            rng = np.random.default_rng(int(config["seed"]))
            checkpoint_population = None
            if resume_checkpoint is not None:
                if signature != resume_checkpoint.signature or not np.array_equal(
                    fit_mask, resume_checkpoint.fit_mask
                ):
                    raise ValueError("Fit checkpoint is incompatible with the current setup")
                start[fit_mask] = np.clip(
                    resume_checkpoint.params[fit_mask],
                    fitted_bounds_array[:, 0],
                    fitted_bounds_array[:, 1],
                )
                checkpoint_population = clipped_checkpoint_population(
                    resume_checkpoint.population,
                    fitted_bounds_array,
                )
                rng.bit_generator.state = copy.deepcopy(resume_checkpoint.rng_state)
            fitted_start = start[fit_mask]
            best_full_params = start.copy()
            best_cost = float("inf")
            resume_phase = (
                resume_checkpoint.phase
                if resume_checkpoint is not None
                else "differential evolution"
            )
            phase = {"value": resume_phase}

            def full_params(fitted_params: np.ndarray) -> np.ndarray:
                params = start.copy()
                params[fit_mask] = fitted_params
                return params

            def check_stop() -> None:
                if self.stop_event.is_set():
                    raise FitCancelled("Fit stopped by user")
                while self.pause_event.is_set():
                    if self.stop_event.wait(0.1):
                        raise FitCancelled("Fit stopped by user")

            def record_checkpoint(
                params: np.ndarray,
                cost: float,
                population: np.ndarray | None = None,
            ) -> None:
                nonlocal best_full_params, best_cost, checkpoint_population
                changed = population is not None
                if population is not None:
                    checkpoint_population = np.asarray(population, dtype=float).copy()
                if np.isfinite(cost) and cost < best_cost:
                    best_cost = float(cost)
                    best_full_params = np.asarray(params, dtype=float).copy()
                    changed = True
                if changed:
                    self.fit_checkpoint = FitCheckpoint(
                        signature=signature,
                        fit_mask=fit_mask.copy(),
                        params=best_full_params.copy(),
                        population=None
                        if checkpoint_population is None
                        else checkpoint_population.copy(),
                        rng_state=copy.deepcopy(rng.bit_generator.state),
                        phase=phase["value"],
                        cost=best_cost,
                    )

            def objective(fitted_params: np.ndarray) -> float:
                check_stop()
                params = full_params(fitted_params)
                cost = model.objective(params)
                record_checkpoint(params, cost)
                return cost

            def residual_vector(fitted_params: np.ndarray) -> np.ndarray:
                check_stop()
                params = full_params(fitted_params)
                residual = model.residual_vector(params)
                record_checkpoint(params, float(np.mean(np.abs(residual))))
                return _least_squares_residual(residual)

            def send_update(phase: str, params: np.ndarray, force: bool = False) -> None:
                check_stop()
                counter["value"] += 1
                if not force and counter["value"] % int(config["interval"]) != 0:
                    return
                predicted = model.predict(params)
                cost = model.objective(params)
                record_checkpoint(params, cost)
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

            send_update("resumed checkpoint" if resume_checkpoint else "initial", start, force=True)

            def de_callback(intermediate_result) -> bool:
                params = full_params(np.asarray(intermediate_result.x, dtype=float))
                record_checkpoint(
                    params,
                    float(intermediate_result.fun),
                    np.asarray(intermediate_result.population, dtype=float),
                )
                send_update("differential evolution", params)
                return False

            if resume_phase in {"least squares", "Powell polish"}:
                global_x = fitted_start.copy()
                global_fun = objective(global_x)
            else:
                phase["value"] = "differential evolution"
                global_result = differential_evolution(
                    objective,
                    bounds,
                    rng=rng,
                    maxiter=int(config["maxiter"]),
                    popsize=int(config["popsize"]),
                    tol=2e-4,
                    polish=False,
                    updating="immediate",
                    workers=1,
                    callback=de_callback,
                    init="latinhypercube"
                    if checkpoint_population is None
                    else checkpoint_population,
                    x0=fitted_start,
                )
                check_stop()
                global_x = np.asarray(global_result.x, dtype=float)
                global_fun = float(global_result.fun)
                record_checkpoint(
                    full_params(global_x),
                    global_fun,
                    np.asarray(global_result.population, dtype=float),
                )

            if resume_phase == "Powell polish":
                local_x = fitted_start.copy()
            else:
                phase["value"] = "least squares"
                local_start = global_x if global_fun < objective(fitted_start) else fitted_start
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
                local_x = np.asarray(local_result.x, dtype=float)
                send_update("least squares", full_params(local_x), force=True)

            def polish_callback(xk: np.ndarray) -> None:
                send_update("Powell polish", full_params(np.asarray(xk)))

            phase["value"] = "Powell polish"
            polish_result = minimize(
                objective,
                local_x,
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
                full_params(global_x),
                full_params(local_x),
                full_params(polish_result.x),
            ]
            best_params = min(candidates, key=model.objective)
            predicted = model.predict(best_params)
            residual = predicted - model.observed
            cost = model.objective(best_params)
            rmse = float(np.sqrt(np.mean(residual**2)))
            phase["value"] = "complete"
            record_checkpoint(best_params, cost, checkpoint_population)
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
            if bool(config.get("stack_enabled", False)):
                definition = StackDefinition(
                    Path(str(config["stack_path"])),
                    copy.deepcopy(config["stack_document"]),
                )
                resolved_path = ROOT / "validation" / (
                    f"{sample_slug}_gui_{model_slug}_fit_superlattice.json"
                )
                write_json_document(
                    resolved_path, definition.resolved_document(best_params)
                )
                fitted_lines = "\n".join(
                    f"{name}: {value:.8g}"
                    for name, value in definition.overrides(best_params).items()
                )
                summary = (
                    "Superlattice fit\n"
                    f"Superlattice: {definition.name}\n"
                    f"Model: {config['model']}\n"
                    f"Data file: {resolve_data_path(str(config['data_path']))}\n"
                    f"Mean abs log10 error: {cost:.6e}\n"
                    f"RMSE: {rmse:.6e} cps\n"
                    f"{fitted_lines}\n"
                    f"Resolved superlattice: {resolved_path}"
                )
            else:
                summary = summarize_fit(
                    sample,
                    str(config["model"]),
                    best_params,
                    cost,
                    rmse,
                    bool(config["include_strain"]),
                    bool(config["include_roughness"]),
                    bool(config["include_kinematic_substrate"]),
                    film_settings if config["model"] == "Dynamic" else None,
                    substrate_settings if config["model"] == "Dynamic" else None,
                )
                backend_text = (
                    f"Density sampling: {config['density_slices']} slices/cell, "
                    f"Q max {float(config['density_max_q0']):.6g} 1/A\n"
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
                self.paused = False
                self.active_fit_config = None
                self.pending_resume_config = None
                self.pause_event.clear()
                self.simulate_button.configure(state=tk.NORMAL)
                self.run_button.configure(state=tk.NORMAL)
                self.pause_button.configure(
                    state=tk.NORMAL if self.fit_checkpoint is not None else tk.DISABLED,
                    text="Continue fit",
                )
                self.stop_button.configure(state=tk.DISABLED)
            elif message_type == "cancelled":
                pending_config = self.pending_resume_config
                self.running = False
                self.paused = False
                self.active_fit_config = None
                self.pending_resume_config = None
                self.pause_event.clear()
                self.simulate_button.configure(state=tk.NORMAL)
                self.run_button.configure(state=tk.NORMAL)
                self.pause_button.configure(
                    state=tk.NORMAL if self.fit_checkpoint is not None else tk.DISABLED,
                    text="Resume fit",
                )
                self.stop_button.configure(state=tk.DISABLED)
                if pending_config is not None:
                    self.status_var.set("Restarting from checkpoint with updated bounds...")
                    self.root.after(0, lambda config=pending_config: self._resume_checkpoint(config))
                elif self.fit_checkpoint is not None:
                    self.status_var.set("Fit stopped; checkpoint is ready to resume")
                else:
                    self.status_var.set(str(payload))
            elif message_type == "error":
                self.status_var.set("Fit failed")
                self.running = False
                self.paused = False
                self.active_fit_config = None
                self.pending_resume_config = None
                self.pause_event.clear()
                self.simulate_button.configure(state=tk.NORMAL)
                self.run_button.configure(state=tk.NORMAL)
                self.pause_button.configure(
                    state=tk.NORMAL if self.fit_checkpoint is not None else tk.DISABLED,
                    text="Resume fit",
                )
                self.stop_button.configure(state=tk.DISABLED)
                messagebox.showerror("Fit failed", str(payload))

        self.root.after(150, self._process_queue)

    def _apply_update(self, update: FitUpdate) -> None:
        self.last_update = update
        if self.active_fit_config is not None:
            self._set_fit_result_values(self.active_fit_config, update.params)
        if update.show_observed:
            self.history_x.append(len(self.history_x) + 1)
            self.history_y.append(update.cost)
            self.history_phase.append(update.phase)
            if self.paused:
                self.status_var.set(f"Fit paused at {update.phase}: cost={update.cost:.5g}")
            else:
                self.status_var.set(f"{update.phase}: cost={update.cost:.5g}")
        else:
            self.status_var.set(update.phase.title())

        self.loss_axis.clear()
        if self.history_x:
            self.loss_axis.plot(self.history_x, self.history_y, color="tab:blue", linewidth=1.5)
        for index in range(1, len(self.history_phase)):
            if self.history_phase[index] != self.history_phase[index - 1]:
                x = self.history_x[index]
                self.loss_axis.axvline(x, color="0.45", linewidth=0.8, linestyle="--")
                self.loss_axis.text(
                    x,
                    0.98,
                    self.history_phase[index],
                    transform=self.loss_axis.get_xaxis_transform(),
                    rotation=90,
                    va="top",
                    ha="right",
                    fontsize=8,
                    color="0.35",
                )
        self.loss_axis.set_xlabel("progress callback")
        self.loss_axis.set_ylabel("mean abs log10 error")
        self.loss_axis.set_title(
            f"{update.phase}: cost={update.cost:.5g}" if update.show_observed else update.phase.title()
        )
        self.loss_axis.grid(True, alpha=0.25)

        self.fit_axis.clear()
        x_values = self._x_values(update.twotheta, update.q, CU_K_ALPHA_WAVELENGTH)
        if update.show_observed:
            self.fit_axis.plot(
                x_values, update.observed, ".", color="black", markersize=3, label="data"
            )
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
