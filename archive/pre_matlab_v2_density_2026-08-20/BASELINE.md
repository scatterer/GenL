# Pre-MATLAB-v2 density baseline

Snapshot date: 2026-08-20

Git HEAD: `a20e62ff82ef27e12d0cf85a36df364554f8ba46`

This snapshot captures the active Python source, launcher, examples, tests,
superlattice definitions, dynamic benchmark, README, and project metadata before
porting the analytic density calculation from `subroutines_v2`.

The working tree already contained uncommitted GUI and superlattice changes.
They are included in this snapshot rather than reconstructed from Git HEAD.

```text
 M .gitignore
 M README.md
 M src/genl/__init__.py
 M src/genl/dynamic.py
 M src/genl/gui.py
 M src/genl/kinematic.py
 M src/genl/paths.py
 M tests/test_core.py
?? data/examples/Fe-V.dat
?? data/stacks/
?? examples/superlattice.py
?? matlab/kinematic_and_dynamic/subroutines_v2/
?? src/genl/stack.py
?? validation/fe-v_fit_gui_dynamic_fit_superlattice.json
?? validation/fe-v_gui_dynamic_fit_superlattice.json
```

## Test baseline

Command:

```bash
.venv/bin/python -m unittest discover -s tests -v
```

Result: 28 tests run; 25 passed, 2 failed, and 1 errored. The three pre-existing
failures were:

- `test_fe_v_stack_expands_and_simulates`
- `test_superlattice_capping_layer_is_appended_and_fittable`
- `test_superlattice_fit_residual_uses_selected_parameters`

All dynamic density and propagation tests passed.

## Dynamic benchmark

```text
cached reflection seconds: 0.005847
cached fused seconds: 0.007627
uncached fused seconds: 0.010683
legacy seconds: 0.097354
reflection/fused diff: 1.145174e-16
max abs reflectivity diff: 2.114789e-16
```

Generated validation output and Python cache files are intentionally excluded.
