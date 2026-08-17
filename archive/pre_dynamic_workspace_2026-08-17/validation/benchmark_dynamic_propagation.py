from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import genl.dynamic as dynamic  # noqa: E402
from genl import Control, Instrument, Layer  # noqa: E402


def benchmark_dynamic_propagation() -> dict[str, float]:
    data = np.loadtxt(ROOT / "kinematic_and_dynamic" / "examples" / "Example_data_10nmFe.txt")
    mask = (data[:, 0] >= 58.92) & (data[:, 0] <= 70.92)
    twotheta = data[mask, 0]
    wavelength = 1.54056
    q = 4.0 * np.pi / wavelength * np.sin(np.deg2rad(twotheta / 2.0))
    stack = [
        Layer(direction=1, n=1e6, filename="MgO_001_fractional.vasp"),
        Layer(
            direction=1,
            n=28.5,
            filename="Fe_fractional.vasp",
            dinterface=1.4,
            scale=1.04,
            area_scale=1.1927,
        ),
    ]
    kwargs = {
        "poscar_dir": ROOT / "kinematic_and_dynamic" / "POSCAR",
        "form_factor_dir": ROOT / "kinematic_and_dynamic" / "Form_Factor_and_Elemental_data",
        "slices": 50,
        "max_q0": 30.0,
        "step_q0": 0.1,
    }

    # Warm Numba compilation before timing.
    dynamic.calc_dynamic_density(
        q, wavelength, stack, Control(pol=2), Instrument(theta_m=2), propagation_backend="fused", **kwargs
    )

    start = time.perf_counter()
    fused_result = dynamic.calc_dynamic_density(
        q, wavelength, stack, Control(pol=2), Instrument(theta_m=2), propagation_backend="fused", **kwargs
    )
    fused_seconds = time.perf_counter() - start

    start = time.perf_counter()
    legacy_result = dynamic.calc_dynamic_density(
        q, wavelength, stack, Control(pol=2), Instrument(theta_m=2), propagation_backend="legacy", **kwargs
    )
    legacy_seconds = time.perf_counter() - start

    return {
        "fused_seconds": fused_seconds,
        "legacy_seconds": legacy_seconds,
        "speedup": legacy_seconds / fused_seconds,
        "max_abs_diff": float(np.max(np.abs(fused_result.refl - legacy_result.refl))),
    }


def main() -> int:
    metrics = benchmark_dynamic_propagation()
    print("Dynamic backend benchmark")
    print(f"fused seconds: {metrics['fused_seconds']:.6f}")
    print(f"legacy seconds: {metrics['legacy_seconds']:.6f}")
    print(f"speedup: {metrics['speedup']:.3f}x")
    print(f"max abs reflectivity diff: {metrics['max_abs_diff']:.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
