from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import reduce
from math import gcd
from pathlib import Path

import numpy as np

from .poscar import PoscarStructure, read_poscar


@dataclass(frozen=True)
class ProjectedStructure:
    """Exact 1D specular projection of a periodic crystal."""

    source: PoscarStructure
    hkl: tuple[int, int, int]
    period: float
    area: float
    z: np.ndarray


def _reduced_hkl(hkl: tuple[int, int, int]) -> tuple[int, int, int]:
    values = tuple(int(v) for v in hkl)
    if values == (0, 0, 0):
        raise ValueError("Miller indices cannot all be zero")
    divisor = reduce(gcd, (abs(v) for v in values if v != 0))
    return tuple(v // divisor for v in values)


def project_structure(
    structure: PoscarStructure,
    hkl: tuple[int, int, int],
) -> ProjectedStructure:
    """Project an arbitrary crystal onto the normal of the (hkl) planes.

    The returned area is V/d_hkl. Keeping all atoms from the input cell with
    this area makes the projected density invariant to primitive,
    conventional, and supercell representations of the same crystal.
    """

    reduced = _reduced_hkl(hkl)
    lattice = structure.lattice
    reciprocal = 2.0 * np.pi * np.linalg.inv(lattice).T
    normal_g = np.asarray(reduced, dtype=float) @ reciprocal
    period = float(2.0 * np.pi / np.linalg.norm(normal_g))
    area = structure.volume / period

    phase = structure.positions @ np.asarray(reduced, dtype=float)
    fractional_z = np.mod(phase, 1.0)
    fractional_z[np.isclose(fractional_z, 1.0, atol=1e-12)] = 0.0
    z = fractional_z * period

    return ProjectedStructure(structure, reduced, period, area, z)


def write_genl_projection(
    projected: ProjectedStructure,
    filename: str | Path,
) -> Path:
    """Write a GenL-compatible POSCAR representing the exact 1D projection.

    This is deliberately a projection cell, not a physical surface cell. The
    in-plane vectors are chosen only to reproduce the required projected area.
    """

    path = Path(filename)
    side = float(np.sqrt(projected.area))
    lines = [
        (
            f"GenL projection "
            f"({projected.hkl[0]} {projected.hkl[1]} {projected.hkl[2]})"
        ),
        "1.0",
        f"{side:.12f} 0.000000000000 0.000000000000",
        f"0.000000000000 {side:.12f} 0.000000000000",
        f"0.000000000000 0.000000000000 {projected.period:.12f}",
        " ".join(projected.source.types),
        " ".join(str(int(v)) for v in projected.source.type_counts),
        "Direct",
    ]
    for z in projected.z:
        lines.append(
            f"0.000000000000 0.000000000000 {z / projected.period:.12f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Project an arbitrary POSCAR onto a crystallographic (hkl) "
            "normal for GenL"
        )
    )
    parser.add_argument("poscar", type=Path)
    parser.add_argument("h", type=int)
    parser.add_argument("k", type=int)
    parser.add_argument("l", type=int)
    parser.add_argument("-o", "--output", type=Path)
    args = parser.parse_args()

    projected = project_structure(
        read_poscar(args.poscar),
        (args.h, args.k, args.l),
    )
    print(f"hkl: {projected.hkl}")
    print(f"normal repeat: {projected.period:.8f} A")
    print(f"effective in-plane area: {projected.area:.8f} A^2")
    for element, count in zip(
        projected.source.types,
        projected.source.type_counts,
    ):
        print(f"{element}: {int(count)} atoms")

    if args.output is not None:
        write_genl_projection(projected, args.output)
        print(f"wrote: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
