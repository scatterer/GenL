# Restore pre-reflection-recursion version

This snapshot contains the Python dynamic calculation, GUI, tests, benchmark,
and documentation immediately before the reflection-recursion backend was
added on 2026-08-17.

From the repository root, restore it with:

```bash
cp archive/pre_reflection_recursion_2026-08-17/src/genl/dynamic.py src/genl/dynamic.py
cp archive/pre_reflection_recursion_2026-08-17/src/genl/fit_models.py src/genl/fit_models.py
cp archive/pre_reflection_recursion_2026-08-17/src/genl/gui.py src/genl/gui.py
cp archive/pre_reflection_recursion_2026-08-17/src/genl/__init__.py src/genl/__init__.py
cp archive/pre_reflection_recursion_2026-08-17/tests/test_core.py tests/test_core.py
cp archive/pre_reflection_recursion_2026-08-17/validation/benchmark_dynamic_propagation.py validation/benchmark_dynamic_propagation.py
cp archive/pre_reflection_recursion_2026-08-17/README.md README.md
```
