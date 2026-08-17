from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class PoscarStructure:
    """Minimal POSCAR data used by the MATLAB GenL routines."""

    types: tuple[str, ...]
    type_counts: np.ndarray
    positions: np.ndarray
    a1: np.ndarray
    a2: np.ndarray
    a3: np.ndarray


def read_poscar(filename: str | Path, poscar_dir: str | Path | None = None) -> PoscarStructure:
    """Read the subset of POSCAR used by GenL.

    The original MATLAB reader treats the coordinate rows as fractional
    coordinates regardless of the POSCAR coordinate-system label. This port
    intentionally preserves that behavior because the bundled files rely on it.
    """

    path = Path(filename)
    if poscar_dir is not None and not path.is_absolute():
        path = Path(poscar_dir) / path

    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 8:
        raise ValueError(f"POSCAR file is too short: {path}")

    scale = float(lines[1].split()[0])
    a1 = np.fromstring(lines[2], sep=" ", dtype=float) * scale
    a2 = np.fromstring(lines[3], sep=" ", dtype=float) * scale
    a3 = np.fromstring(lines[4], sep=" ", dtype=float) * scale
    types = tuple(lines[5].split())
    type_counts = np.fromstring(lines[6], sep=" ", dtype=int)

    if len(types) != len(type_counts):
        raise ValueError(
            f"POSCAR element/count mismatch in {path}: {types} vs {type_counts}"
        )

    n_positions = int(type_counts.sum())
    position_lines = lines[8 : 8 + n_positions]
    positions = np.array(
        [np.fromstring(line, sep=" ", dtype=float)[:3] for line in position_lines],
        dtype=float,
    )
    if positions.shape != (n_positions, 3):
        raise ValueError(f"Could not read {n_positions} POSCAR positions from {path}")

    return PoscarStructure(types, type_counts, positions, a1, a2, a3)
