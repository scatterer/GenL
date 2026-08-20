# Restore the pre-MATLAB-v2 density version

From the repository root:

```bash
cp archive/pre_matlab_v2_density_2026-08-20/src/genl/*.py src/genl/
cp archive/pre_matlab_v2_density_2026-08-20/apps/*.py apps/
cp archive/pre_matlab_v2_density_2026-08-20/examples/*.py examples/
cp archive/pre_matlab_v2_density_2026-08-20/tests/*.py tests/
cp archive/pre_matlab_v2_density_2026-08-20/data/stacks/*.json data/stacks/
cp archive/pre_matlab_v2_density_2026-08-20/validation/benchmark_dynamic_propagation.py validation/
cp archive/pre_matlab_v2_density_2026-08-20/README.md README.md
cp archive/pre_matlab_v2_density_2026-08-20/pyproject.toml pyproject.toml
```

This restores the active Python files but does not delete files added later.
