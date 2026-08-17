from __future__ import annotations

import numpy as np


def pad_edge(values: np.ndarray, size: int = 154) -> np.ndarray:
    """Pad a vector with its edge values, matching MATLAB `padme`."""

    values = np.asarray(values)
    half = size // 2
    return np.pad(values, (half, half), mode="edge")


def gauss_conv(x: np.ndarray, y: np.ndarray, fwhm: float) -> np.ndarray:
    """Gaussian convolution used for instrumental resolution broadening."""

    x = np.asarray(x, dtype=float)
    y = np.asarray(y)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError("x and y must be one-dimensional arrays of equal length")
    if len(x) < 2:
        raise ValueError("at least two x points are required")
    if fwhm <= 0:
        return y.copy()

    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    dx = abs(x[0] - x[1])
    points = int(np.floor(4.0 * sigma / dx))

    padded = pad_edge(y, 154)
    offsets = np.arange(-points, points + 2) * dx
    kernel = 1.0 / (sigma * np.sqrt(2.0 * np.pi)) * np.exp(
        -((offsets / (np.sqrt(2.0) * sigma)) ** 2)
    )
    kernel = kernel / np.sum(kernel)

    convolved = np.convolve(padded, kernel, mode="same")
    half = 154 // 2
    return convolved[half:-half]
