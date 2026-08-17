# Restore pre-fit-checkpoint version

This snapshot contains the GUI and tests from before pause/resume and optimizer
population checkpoints were added on 2026-08-17.

From the repository root, restore them with:

```bash
cp archive/pre_fit_checkpoint_2026-08-17/src/genl/gui.py src/genl/gui.py
cp archive/pre_fit_checkpoint_2026-08-17/tests/test_core.py tests/test_core.py
cp archive/pre_fit_checkpoint_2026-08-17/README.md README.md
```
