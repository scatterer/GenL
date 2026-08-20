from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import genl.dynamic as dynamic  # noqa: E402
from genl import Control, DynamicWorkspace, Instrument, Layer  # noqa: E402
from genl.paths import EXAMPLE_DATA_DIR, FORM_FACTOR_DIR, STRUCTURE_DIR  # noqa: E402


def benchmark_dynamic_propagation() -> dict[str, float]:
    data = np.loadtxt(EXAMPLE_DATA_DIR / "Example_data_10nmFe.txt")
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
        "poscar_dir": STRUCTURE_DIR,
        "form_factor_dir": FORM_FACTOR_DIR,
        "slices": 50,
        "max_q0": 30.0,
        "step_q0": 0.1,
    }

    reflection_workspace = DynamicWorkspace()
    dynamic.calc_dynamic_density(
        q,
        wavelength,
        stack,
        Control(pol=2),
        Instrument(theta_m=2),
        propagation_backend="reflection",
        workspace=reflection_workspace,
        **kwargs,
    )
    start = time.perf_counter()
    reflection_result = dynamic.calc_dynamic_density(
        q,
        wavelength,
        stack,
        Control(pol=2),
        Instrument(theta_m=2),
        propagation_backend="reflection",
        workspace=reflection_workspace,
        **kwargs,
    )
    reflection_seconds = time.perf_counter() - start

    workspace = DynamicWorkspace()
    # Warm Numba compilation and populate reusable material/substrate state.
    dynamic.calc_dynamic_density(
        q,
        wavelength,
        stack,
        Control(pol=2),
        Instrument(theta_m=2),
        propagation_backend="fused",
        workspace=workspace,
        **kwargs,
    )

    start = time.perf_counter()
    cached_result = dynamic.calc_dynamic_density(
        q,
        wavelength,
        stack,
        Control(pol=2),
        Instrument(theta_m=2),
        propagation_backend="fused",
        workspace=workspace,
        **kwargs,
    )
    cached_seconds = time.perf_counter() - start

    start = time.perf_counter()
    uncached_result = dynamic.calc_dynamic_density(
        q, wavelength, stack, Control(pol=2), Instrument(theta_m=2), propagation_backend="fused", **kwargs
    )
    uncached_seconds = time.perf_counter() - start

    start = time.perf_counter()
    legacy_result = dynamic.calc_dynamic_density(
        q, wavelength, stack, Control(pol=2), Instrument(theta_m=2), propagation_backend="legacy", **kwargs
    )
    legacy_seconds = time.perf_counter() - start

    analytic_workspace = DynamicWorkspace()
    start = time.perf_counter()
    analytic_result = dynamic.calc_dynamic_density(
        q,
        wavelength,
        stack,
        Control(pol=2),
        Instrument(theta_m=2),
        propagation_backend="reflection",
        density_method="analytic",
        workspace=analytic_workspace,
        **kwargs,
    )
    analytic_cold_seconds = time.perf_counter() - start
    start = time.perf_counter()
    dynamic.calc_dynamic_density(
        q,
        wavelength,
        stack,
        Control(pol=2),
        Instrument(theta_m=2),
        propagation_backend="reflection",
        density_method="analytic",
        workspace=analytic_workspace,
        **kwargs,
    )
    analytic_warm_seconds = time.perf_counter() - start

    return {
        "reflection_seconds": reflection_seconds,
        "cached_seconds": cached_seconds,
        "uncached_seconds": uncached_seconds,
        "legacy_seconds": legacy_seconds,
        "cache_speedup": uncached_seconds / cached_seconds,
        "legacy_speedup": legacy_seconds / cached_seconds,
        "reflection_speedup": cached_seconds / reflection_seconds,
        "reflection_diff": float(np.max(np.abs(reflection_result.refl - cached_result.refl))),
        "max_abs_diff": float(np.max(np.abs(cached_result.refl - legacy_result.refl))),
        "cached_uncached_diff": float(np.max(np.abs(cached_result.refl - uncached_result.refl))),
        "analytic_cold_seconds": analytic_cold_seconds,
        "analytic_warm_seconds": analytic_warm_seconds,
        "analytic_kernel_count": float(len(analytic_workspace.atomic_kernels)),
        "analytic_vacuum": float(analytic_result.diagnostics["vacuum_thickness_used"]),
    }


def main() -> int:
    metrics = benchmark_dynamic_propagation()
    print("Dynamic backend benchmark")
    print(f"cached reflection seconds: {metrics['reflection_seconds']:.6f}")
    print(f"cached fused seconds: {metrics['cached_seconds']:.6f}")
    print(f"uncached fused seconds: {metrics['uncached_seconds']:.6f}")
    print(f"legacy seconds: {metrics['legacy_seconds']:.6f}")
    print(f"cache speedup: {metrics['cache_speedup']:.3f}x")
    print(f"speedup over legacy: {metrics['legacy_speedup']:.3f}x")
    print(f"reflection speedup over fused: {metrics['reflection_speedup']:.3f}x")
    print(f"reflection/fused diff: {metrics['reflection_diff']:.6e}")
    print(f"max abs reflectivity diff: {metrics['max_abs_diff']:.6e}")
    print(f"cached/uncached diff: {metrics['cached_uncached_diff']:.6e}")
    print(f"analytic density cold seconds: {metrics['analytic_cold_seconds']:.6f}")
    print(f"analytic density warm seconds: {metrics['analytic_warm_seconds']:.6f}")
    print(f"analytic cached kernels: {int(metrics['analytic_kernel_count'])}")
    print(f"analytic vacuum used: {metrics['analytic_vacuum']:.6f} A")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
