from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class PoscarStructure:
    """Crystal structure normalized to fractional coordinates and Angstrom lattice vectors."""

    types: tuple[str, ...]
    type_counts: np.ndarray
    positions: np.ndarray
    a1: np.ndarray
    a2: np.ndarray
    a3: np.ndarray

    @property
    def lattice(self) -> np.ndarray:
        return np.vstack((self.a1, self.a2, self.a3))

    @property
    def volume(self) -> float:
        return float(abs(np.linalg.det(self.lattice)))


def _parse_scale(scale_line: str, raw_lattice: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.fromstring(scale_line, sep=" ", dtype=float)
    if values.size == 1:
        scale = float(values[0])
        if scale == 0.0:
            raise ValueError("POSCAR scale factor cannot be zero")
        if scale < 0.0:
            raw_volume = float(abs(np.linalg.det(raw_lattice)))
            if raw_volume == 0.0:
                raise ValueError("POSCAR lattice vectors are singular")
            scale = (abs(scale) / raw_volume) ** (1.0 / 3.0)
        component_scale = np.full(3, scale)
    elif values.size == 3:
        if np.any(values <= 0.0):
            raise ValueError("POSCAR three-component scale factors must be positive")
        component_scale = values
    else:
        raise ValueError("POSCAR scale line must contain one or three numbers")

    lattice = raw_lattice * component_scale[np.newaxis, :]
    return lattice, component_scale


def read_poscar(filename: str | Path, poscar_dir: str | Path | None = None) -> PoscarStructure:
    """Read a VASP POSCAR into GenL's fractional-coordinate structure model.

    Supports Direct/Cartesian coordinates, optional Selective dynamics, one
    universal scale factor (including negative target-volume form), and VASP's
    three Cartesian-component scale factors. Species names are required because
    GenL has no POTCAR from which to infer omitted names.
    """

    path = Path(filename)
    if poscar_dir is not None and not path.is_absolute():
        path = Path(poscar_dir) / path

    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 8:
        raise ValueError(f"POSCAR file is too short: {path}")

    raw_lattice = np.array(
        [np.fromstring(lines[i], sep=" ", dtype=float)[:3] for i in range(2, 5)],
        dtype=float,
    )
    if raw_lattice.shape != (3, 3):
        raise ValueError(f"Could not read lattice vectors from {path}")
    lattice, cartesian_scale = _parse_scale(lines[1], raw_lattice)

    species_tokens = lines[5].split()
    if species_tokens and all(token.lstrip("+-").isdigit() for token in species_tokens):
        raise ValueError(
            "POSCAR species names are omitted; GenL requires explicit element symbols"
        )
    types = tuple(species_tokens)
    type_counts = np.fromstring(lines[6], sep=" ", dtype=int)
    if len(types) != len(type_counts):
        raise ValueError(
            f"POSCAR element/count mismatch in {path}: {types} vs {type_counts}"
        )

    cursor = 7
    if lines[cursor].strip().lower().startswith("s"):
        cursor += 1
    if cursor >= len(lines):
        raise ValueError(f"POSCAR coordinate mode is missing: {path}")

    mode = lines[cursor].strip().lower()
    cartesian = mode.startswith(("c", "k"))
    cursor += 1

    n_positions = int(type_counts.sum())
    position_lines = lines[cursor : cursor + n_positions]
    coordinates = np.array(
        [np.asarray(line.split()[:3], dtype=float) for line in position_lines],
        dtype=float,
    )
    if coordinates.shape != (n_positions, 3):
        raise ValueError(f"Could not read {n_positions} POSCAR positions from {path}")

    if cartesian:
        cartesian_positions = coordinates * cartesian_scale[np.newaxis, :]
        positions = cartesian_positions @ np.linalg.inv(lattice)
    else:
        positions = coordinates

    return PoscarStructure(
        types,
        type_counts,
        positions,
        lattice[0],
        lattice[1],
        lattice[2],
    )
