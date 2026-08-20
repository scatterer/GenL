from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import numpy as np

from .convolution import gauss_conv
from .dynamic import DynamicResult, DynamicWorkspace, calc_dynamic_density
from .kinematic import Control, Instrument, Layer, calc_kinematic
from .paths import FORM_FACTOR_DIR, STRUCTURE_DIR

STACK_FORMAT = "GenL sample stack"
STACK_VERSION = 1

_LAYER_DEFAULTS = {
    "dinterface": 0.0,
    "scale": 1.0,
    "area_scale": 1.0,
    "roughness": False,
    "sigma": 0.0,
    "bottom_strain_amplitude": 0.0,
    "bottom_strain_end": 0.0,
    "top_strain_amplitude": 0.0,
    "top_strain_end": 0.0,
}
_LAYER_TARGETS = {
    "unit_cells",
    "dinterface",
    "scale",
    "area_scale",
    "sigma",
    "bottom_strain_amplitude",
    "bottom_strain_end",
    "top_strain_amplitude",
    "top_strain_end",
}
_CALCULATION_DEFAULTS = {
    "resolution": 0.005,
    "intensity_scale": 1.0,
    "background_a": 0.0,
    "background_b": 0.0,
    "theta_m": 2.0,
    "polarization": 2,
    "vacuum_thickness": 5.0,
    "density_slices": 100,
    "density_max_q0": 30.0,
    "density_step_q0": 0.1,
    "dynamic_backend": "auto",
    "density_method": "sampled",
}


class StackDefinition:
    def __init__(self, path: Path, document: dict[str, object]) -> None:
        self.path = path
        self.document = document
        self._validate()

    @classmethod
    def load(cls, path: str | Path) -> StackDefinition:
        resolved = Path(path).expanduser().resolve()
        with resolved.open(encoding="utf-8") as handle:
            document = json.load(handle)
        if not isinstance(document, dict):
            raise ValueError("Superlattice definition must contain a JSON object")
        return cls(resolved, document)

    @property
    def name(self) -> str:
        return str(self.document["name"])

    @property
    def wavelength(self) -> float:
        return float(self.document.get("wavelength", 1.54056))

    @property
    def model(self) -> str:
        return str(self.document.get("model", "kinematic")).lower()

    @property
    def parameter_names(self) -> tuple[str, ...]:
        return tuple(str(item["target"]) for item in self._fit_parameters())

    @property
    def start(self) -> np.ndarray:
        values = []
        for item in self._fit_parameters():
            value = self._target_value(str(item["target"]))
            values.append(np.log10(value) if item.get("transform") == "log10" else value)
        return np.asarray(values, dtype=float)

    @property
    def bounds(self) -> np.ndarray:
        return np.asarray(
            [[float(item["min"]), float(item["max"])] for item in self._fit_parameters()],
            dtype=float,
        )

    def overrides(self, parameters: np.ndarray | None = None) -> dict[str, float]:
        values = self.start if parameters is None else np.asarray(parameters, dtype=float)
        specs = self._fit_parameters()
        if len(values) != len(specs):
            raise ValueError("Parameter vector does not match the superlattice fit definition")
        return {
            str(item["target"]): float(10.0**value if item.get("transform") == "log10" else value)
            for item, value in zip(specs, values)
        }

    def layers(self, parameters: np.ndarray | None = None) -> list[Layer]:
        overrides = self.overrides(parameters)
        substrate = self.document["substrate"]
        sequence = self.document["sequence"]
        layers = [self._make_layer("substrate", substrate, overrides)]
        for _ in range(int(sequence["repetitions"])):
            for layer in sequence["layers"]:
                layers.append(self._make_layer(str(layer["name"]), layer, overrides))
        capping_layer = self.document.get("capping_layer")
        if capping_layer is not None:
            layers.append(self._make_layer("capping", capping_layer, overrides))
        return layers

    def calculation(self, parameters: np.ndarray | None = None) -> dict[str, object]:
        values = dict(_CALCULATION_DEFAULTS)
        values.update(self.document.get("calculation", {}))
        for target, value in self.overrides(parameters).items():
            prefix, field = target.split(".", 1)
            if prefix == "calculation":
                values[field] = value
        return values

    def resolved_document(self, parameters: np.ndarray | None = None) -> dict[str, object]:
        layers = self.layers(parameters)
        sequence = self.document["sequence"]
        names = ["substrate"] + [
            str(layer["name"])
            for _ in range(int(sequence["repetitions"]))
            for layer in sequence["layers"]
        ]
        capping_layer = self.document.get("capping_layer")
        if capping_layer is not None:
            names.append(str(capping_layer["name"]))
        return {
            "format": "GenL resolved sample stack",
            "version": STACK_VERSION,
            "source": str(self.path),
            "name": self.name,
            "model": self.model,
            "wavelength": self.wavelength,
            "fit_values": self.overrides(parameters),
            "calculation": self.calculation(parameters),
            "layers": [
                {
                    "index": index,
                    "name": name,
                    "filename": str(layer.filename),
                    "direction": layer.direction,
                    "unit_cells": layer.n,
                    "dinterface": layer.dinterface,
                    "scale": layer.scale,
                    "area_scale": layer.area_scale,
                    "roughness": layer.roughness,
                    "sigma": layer.sigma,
                    "bottom_strain_amplitude": layer.bottom_strain_amplitude,
                    "bottom_strain_end": layer.bottom_strain_end,
                    "top_strain_amplitude": layer.top_strain_amplitude,
                    "top_strain_end": layer.top_strain_end,
                }
                for index, (name, layer) in enumerate(zip(names, layers))
            ],
        }

    def _fit_parameters(self) -> list[dict[str, object]]:
        return list(self.document.get("fit_parameters", []))

    def _layer_specs(self) -> dict[str, Mapping[str, object]]:
        sequence = self.document["sequence"]
        specs = {"substrate": self.document["substrate"]}
        specs.update({str(layer["name"]): layer for layer in sequence["layers"]})
        if self.document.get("capping_layer") is not None:
            specs["capping"] = self.document["capping_layer"]
        return specs

    def _target_value(self, target: str) -> float:
        prefix, field = target.split(".", 1)
        if prefix == "calculation":
            calculation = dict(_CALCULATION_DEFAULTS)
            calculation.update(self.document.get("calculation", {}))
            return float(calculation[field])
        spec = self._layer_specs()[prefix]
        return float(spec[field] if field in spec else _LAYER_DEFAULTS[field])

    def _make_layer(
        self,
        name: str,
        spec: Mapping[str, object],
        overrides: Mapping[str, float],
    ) -> Layer:
        def value(field: str) -> object:
            return overrides.get(
                f"{name}.{field}",
                spec.get(field, _LAYER_DEFAULTS.get(field)),
            )

        return Layer(
            direction=int(spec["direction"]),
            n=float(value("unit_cells")),
            filename=str(spec["filename"]),
            dinterface=float(value("dinterface")),
            scale=float(value("scale")),
            area_scale=float(value("area_scale")),
            roughness=bool(value("roughness")),
            sigma=float(value("sigma")),
            bottom_strain_amplitude=float(value("bottom_strain_amplitude")),
            bottom_strain_end=float(value("bottom_strain_end")),
            top_strain_amplitude=float(value("top_strain_amplitude")),
            top_strain_end=float(value("top_strain_end")),
        )

    def _validate(self) -> None:
        if self.document.get("format") != STACK_FORMAT:
            raise ValueError(f"Superlattice file format must be '{STACK_FORMAT}'")
        if self.document.get("version") != STACK_VERSION:
            raise ValueError(
                f"Unsupported superlattice version: {self.document.get('version')}"
            )
        if not self.document.get("name"):
            raise ValueError("Superlattice definition requires a name")
        if self.model not in {"kinematic", "dynamic"}:
            raise ValueError("Superlattice model must be kinematic or dynamic")
        sequence = self.document.get("sequence")
        if not isinstance(sequence, dict) or int(sequence.get("repetitions", 0)) < 1:
            raise ValueError("Superlattice sequence requires at least one repetition")
        sequence_layers = sequence.get("layers")
        if not isinstance(sequence_layers, list) or not sequence_layers:
            raise ValueError("Superlattice sequence requires at least one layer")
        names = [str(layer.get("name", "")) for layer in sequence_layers]
        if (
            not all(names)
            or len(names) != len(set(names))
            or {"substrate", "capping", "calculation"}.intersection(names)
        ):
            raise ValueError("Sequence layer names must be non-empty and unique")
        capping_layer = self.document.get("capping_layer")
        if capping_layer is not None:
            if not isinstance(capping_layer, dict) or not capping_layer.get("name"):
                raise ValueError("Capping layer must define a name")
            if str(capping_layer["name"]) in names:
                raise ValueError("Capping layer name must differ from repeated layer names")
        for name, spec in self._layer_specs().items():
            if "filename" not in spec or "direction" not in spec or "unit_cells" not in spec:
                raise ValueError(f"Layer '{name}' requires filename, direction, and unit_cells")
            if int(spec["direction"]) not in (1, 2, 3) or float(spec["unit_cells"]) <= 0:
                raise ValueError(f"Layer '{name}' has invalid direction or unit-cell count")
            poscar = Path(str(spec["filename"]))
            if not poscar.is_absolute():
                poscar = STRUCTURE_DIR / poscar
            if not poscar.is_file():
                raise ValueError(f"Layer '{name}' structure file not found: {poscar}")
        valid_targets = {
            *(f"{name}.{field}" for name in self._layer_specs() for field in _LAYER_TARGETS),
            *(f"calculation.{field}" for field in _CALCULATION_DEFAULTS),
        }
        for parameter in self._fit_parameters():
            target = str(parameter.get("target", ""))
            if target not in valid_targets:
                raise ValueError(f"Unknown fit target: {target}")
            if parameter.get("transform") not in (None, "log10"):
                raise ValueError(f"Unsupported transform for {target}")
            minimum, maximum = float(parameter["min"]), float(parameter["max"])
            start = self._target_value(target)
            if parameter.get("transform") == "log10":
                if start <= 0:
                    raise ValueError(f"Log-transformed fit target must be positive: {target}")
                start = np.log10(start)
            if minimum > maximum or not minimum <= start <= maximum:
                raise ValueError(f"Invalid bounds for {target}")


class StackModel:
    def __init__(
        self,
        definition: StackDefinition,
        twotheta: np.ndarray,
        observed: np.ndarray | None = None,
        model: str | None = None,
    ) -> None:
        self.definition = definition
        self.twotheta = np.asarray(twotheta, dtype=float)
        self.observed = None if observed is None else np.asarray(observed, dtype=float)
        self.model = (model or definition.model).lower()
        if self.model not in {"kinematic", "dynamic"}:
            raise ValueError("Superlattice model must be kinematic or dynamic")
        if self.observed is not None and self.observed.shape != self.twotheta.shape:
            raise ValueError("Observed data must match the 2theta array")
        self.q = 4.0 * np.pi / definition.wavelength * np.sin(np.deg2rad(self.twotheta / 2.0))
        self.workspace = DynamicWorkspace()
        self.last_dynamic_result: DynamicResult | None = None

    @property
    def last_z(self) -> np.ndarray | None:
        return None if self.last_dynamic_result is None else self.last_dynamic_result.z

    @property
    def last_rho_e(self) -> np.ndarray | None:
        return None if self.last_dynamic_result is None else self.last_dynamic_result.rho_e

    def predict(self, parameters: np.ndarray | None = None) -> np.ndarray:
        layers = self.definition.layers(parameters)
        calculation = self.definition.calculation(parameters)
        control = Control(pol=int(calculation["polarization"]), model=self.model)
        instrument = Instrument(theta_m=float(calculation["theta_m"]), theta=self.twotheta / 2.0)
        if self.model == "dynamic":
            self.last_dynamic_result = calc_dynamic_density(
                self.q,
                self.definition.wavelength,
                layers,
                control,
                instrument,
                poscar_dir=STRUCTURE_DIR,
                form_factor_dir=FORM_FACTOR_DIR,
                vacuum_thick=float(calculation["vacuum_thickness"]),
                slices=int(calculation["density_slices"]),
                max_q0=float(calculation["density_max_q0"]),
                step_q0=float(calculation["density_step_q0"]),
                propagation_backend=str(calculation["dynamic_backend"]),
                workspace=self.workspace,
                density_method=str(calculation["density_method"]),
            )
            scattering = self.last_dynamic_result.refl
        else:
            self.last_dynamic_result = None
            scattering = calc_kinematic(
                self.q,
                self.definition.wavelength,
                layers,
                control,
                instrument,
                poscar_dir=STRUCTURE_DIR,
                form_factor_dir=FORM_FACTOR_DIR,
            ).refl
        broadened = gauss_conv(self.q, scattering, float(calculation["resolution"]))
        return (
            float(calculation["intensity_scale"]) * broadened
            + float(calculation["background_a"]) * self.q
            + float(calculation["background_b"])
        )

    def residual_vector(self, parameters: np.ndarray) -> np.ndarray:
        if self.observed is None or not np.any(self.observed > 0):
            raise ValueError("Positive experimental intensities are required for fitting")
        predicted = self.predict(parameters)
        if not np.all(np.isfinite(predicted)):
            return np.full_like(predicted, np.inf)
        floor = max(float(np.min(self.observed[self.observed > 0])) * 0.1, 1e-12)
        return np.log10(np.maximum(predicted, floor)) - np.log10(
            np.maximum(self.observed, floor)
        )

    def objective(self, parameters: np.ndarray) -> float:
        return float(np.mean(np.abs(self.residual_vector(parameters))))
