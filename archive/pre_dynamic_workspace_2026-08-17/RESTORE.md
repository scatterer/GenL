# Restore pre-dynamic-workspace version

This snapshot contains the calculation files from before material caching,
substrate reuse, and paired-polarization propagation were added on 2026-08-17.

From the repository root, restore them with:

```bash
cp archive/pre_dynamic_workspace_2026-08-17/src/genl/dynamic.py src/genl/dynamic.py
cp archive/pre_dynamic_workspace_2026-08-17/src/genl/fit_models.py src/genl/fit_models.py
cp archive/pre_dynamic_workspace_2026-08-17/src/genl/__init__.py src/genl/__init__.py
cp archive/pre_dynamic_workspace_2026-08-17/tests/test_core.py tests/test_core.py
cp archive/pre_dynamic_workspace_2026-08-17/validation/benchmark_dynamic_propagation.py validation/benchmark_dynamic_propagation.py
cp archive/pre_dynamic_workspace_2026-08-17/README.md README.md
```
