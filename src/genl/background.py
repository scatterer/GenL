from __future__ import annotations

import numpy as np


def centered_polynomial_background(
    q: np.ndarray,
    slope: float,
    offset: float,
    curvature: float = 0.0,
) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    half_range = float(np.ptp(q)) / 2.0
    x = (
        np.zeros_like(q)
        if half_range == 0.0
        else (q - (q.min() + q.max()) / 2.0) / half_range
    )
    return offset + slope * x + curvature * x**2
