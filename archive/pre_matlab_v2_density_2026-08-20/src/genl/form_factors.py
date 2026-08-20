from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from .paths import FORM_FACTOR_DIR

ELEMENT_SYMBOLS = (
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
)


@dataclass(frozen=True)
class FormFactorCoefficients:
    a: np.ndarray
    b: np.ndarray
    c: np.ndarray
    f_1: np.ndarray
    f_2: np.ndarray


def default_data_dir() -> Path:
    if (FORM_FACTOR_DIR / "ASF.DAT").exists() and (
        FORM_FACTOR_DIR / "f1f2_BrennanCowanLong.dat"
    ).exists():
        return FORM_FACTOR_DIR
    raise FileNotFoundError(f"Could not locate GenL form-factor data in {FORM_FACTOR_DIR}")


def read_form_factor_coefficients(
    z: int | list[int] | tuple[int, ...] | np.ndarray,
    wavelength: float,
    data_dir: str | Path | None = None,
) -> FormFactorCoefficients:
    """Read Waasmaier/Kirfel and Brennan/Cowan coefficients for elements Z."""

    z_values = np.atleast_1d(np.asarray(z, dtype=int))
    source_dir = Path(data_dir) if data_dir is not None else default_data_dir()
    asf_rows = _read_asf(source_dir / "ASF.DAT")

    hc = 1.23984193e4
    energy_source = hc / wavelength

    a = np.zeros((len(z_values), 5), dtype=float)
    b = np.zeros((len(z_values), 5), dtype=float)
    c = np.zeros(len(z_values), dtype=float)
    f_1 = np.zeros(len(z_values), dtype=float)
    f_2 = np.zeros(len(z_values), dtype=float)

    for idx, atomic_number in enumerate(z_values):
        if atomic_number < 1 or atomic_number > len(ELEMENT_SYMBOLS):
            raise ValueError(f"Atomic number out of range: {atomic_number}")
        row = asf_rows[atomic_number]
        a[idx, :] = row[0:5]
        c[idx] = row[5]
        b[idx, :] = row[6:11]
        f_1[idx], f_2[idx] = _find_f1f2(
            atomic_number, energy_source, source_dir / "f1f2_BrennanCowanLong.dat"
        )

    return FormFactorCoefficients(a=a, b=b, c=c, f_1=f_1, f_2=f_2)


def form_factors(
    q: float | np.ndarray,
    coefficients: FormFactorCoefficients,
    composition: float | list[float] | tuple[float, ...] | np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate Q-dependent complex form factors.

    `composition` should be expressed as fractions, matching MATLAB calls such
    as `composition./100`.
    """

    q_values = np.atleast_1d(np.asarray(q, dtype=float))
    composition_values = np.atleast_1d(np.asarray(composition, dtype=float))
    if len(composition_values) != len(coefficients.c):
        raise ValueError("Composition length must match coefficient length")

    k_squared = (q_values / (4.0 * np.pi)) ** 2
    form_factor = (
        coefficients.c[:, np.newaxis]
        + coefficients.f_1[:, np.newaxis]
        - 1j * coefficients.f_2[:, np.newaxis]
    )
    form_factor = form_factor + np.sum(
        coefficients.a[:, :, np.newaxis]
        * np.exp(-coefficients.b[:, :, np.newaxis] * k_squared),
        axis=1,
    )

    weighted = composition_values[:, np.newaxis] * form_factor
    f = np.sum(weighted, axis=0)
    f_sqrd = np.sum(composition_values[:, np.newaxis] * form_factor * form_factor.conj(), axis=0)
    f_sqrd_real = np.real(f_sqrd)
    f_av_sqrd_real = np.real(f * f.conj())
    return f, f_sqrd_real, f_av_sqrd_real


@lru_cache(maxsize=8)
def _read_asf(path: Path) -> dict[int, np.ndarray]:
    rows: dict[int, np.ndarray] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("@"):
            continue
        parts = line.split()
        if len(parts) != 12:
            continue
        symbol = parts[0]
        if symbol not in ELEMENT_SYMBOLS:
            continue
        rows[ELEMENT_SYMBOLS.index(symbol) + 1] = np.array(parts[1:], dtype=float)

    if len(rows) < 90:
        raise ValueError(f"ASF table appears incomplete: {path}")
    return rows


@lru_cache(maxsize=256)
def _find_f1f2(z: int, energy_source: float, path: Path) -> tuple[float, float]:
    energies: list[float] = []
    f1_values: list[float] = []
    f2_values: list[float] = []
    in_section = False

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if stripped.startswith("#S"):
            if in_section:
                break
            parts = stripped.split()
            in_section = len(parts) >= 2 and int(parts[1]) == z
            continue
        if not in_section or not stripped or stripped.startswith("#"):
            continue

        parts = stripped.split()
        if len(parts) >= 3:
            energies.append(float(parts[0]))
            f1_values.append(float(parts[1]))
            f2_values.append(float(parts[2]))

    if not energies:
        return 0.0, 0.0

    idx = int(np.argmin(np.abs(np.asarray(energies) - energy_source)))
    return f1_values[idx], f2_values[idx]
