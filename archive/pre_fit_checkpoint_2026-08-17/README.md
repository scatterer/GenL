# GenL
X-ray reflectivity and diffraction fitting tool for single crystal films

The original MATLAB source code is kept under `matlab/`:

- `matlab/kinematic_only/`: original kinematic-only MATLAB GUI and command-line source.
- `matlab/kinematic_and_dynamic/`: original MATLAB command-line source for kinematic and dynamic diffraction, including fitting.

The shared example data, POSCAR files, form-factor tables, and bundled binaries
remain in the original `kinematic_only/` and `kinematic_and_dynamic/` data
folders so the MATLAB and Python ports can read the same inputs.

Folder `matlab/kinematic_only` includes a graphical user interface for fitting single layers and a command line version where more complicated structures can be fitted.

Binaries:

Binary is available for MacOS Ventura, Windows 10 and Windows 11 in `kinematic_only/binaries`.

When running the installer, the free MATLAB runtime environment will be installed and the application appears in Applications.

When running the application for the first time, it may take a couple of minutes to start.

Once the program starts, choose the working folder (on mac OS: default is /Applications/GenL/application/).


Folder `matlab/kinematic_and_dynamic` includes a command line version that can do either kinematic or dynamic diffraction including fitting.

We recommend starting with the kinematic version to get a rough fit and then switching to the dynamic for higher accuracy.

- Run and or modify `matlab/kinematic_and_dynamic/run_gaas_substrate.m` to simulate a pattern.
- Run and or modify `matlab/kinematic_and_dynamic/fit_laue_oscillations_fe.m` to fit diffraction using either kinematic or dynamic model.

Superlattice can be modeled using the engine without the GUI.

## Python translation

The Python port lives in `src/genl`. It currently includes:

- POSCAR reading
- tabulated form-factor coefficient loading
- Q-dependent form-factor calculation
- Debye-Waller prefactors
- Gaussian instrumental broadening
- kinematic scattering for POSCAR-backed layers
- density-based dynamic scattering with Numba-accelerated density generation
  and transfer-matrix propagation when `numba` is installed
- substrate-only and substrate + film dynamic stacks

Dynamic scattering supports three backends:

- `auto`: use the Numba density and fused propagation kernels when available,
  otherwise fall back to the legacy NumPy/vectorized path.
- `fused`: require the Numba density and fused propagation kernels.
- `legacy`: use the previous NumPy density and Python vectorized propagation
  path with the existing repeated-matrix combiner.

The GUI dynamic fit model automatically reuses a `DynamicWorkspace`. It caches
parsed POSCAR/form-factor material data, substrate density and transfer
matrices, computes sigma/pi polarization together in the fused Numba kernel,
and shares substrate work across the thickness distribution used for film
roughness. The one-entry substrate cache invalidates automatically when its
scale or other substrate inputs change. Direct API callers can get the same
reuse by passing one `DynamicWorkspace()` to repeated `calc_dynamic_density`
calls; omitting it retains the stateless behavior.

Set the backend in code with `propagation_backend=...`, in the GUI with the
`Dynamic backend` selector in `Dynamic fit parameters`, or globally with
`GENL_DYNAMIC_BACKEND=legacy` before launch.

If you use GenL, please cite the original GenL article:
Vassilios Kapaklis and Gunnar K. Palsson, J. Appl. Cryst. 59, 968-977 (2026),
https://doi.org/10.1107/S1600576726002566.

The Python project is organized as:

- `src/genl/`: reusable package code, including the core scattering routines,
  fit model helpers, and the GUI implementation.
- `apps/`: user-facing launchers, including the generic GUI launcher.
- `examples/`: small runnable examples.
- `validation/`: numerical validation scripts and compatibility wrappers.
- `tests/`: unit tests.
- `matlab/`: original MATLAB source, separated from the Python package.
- `archive/`: older or retired material kept for reference, not active code.

Run the bundled tests from the repository root:

```bash
.venv/bin/python -m unittest discover -s tests
```

Run a minimal Fe kinematic example:

```bash
.venv/bin/python examples/kinematic_fe.py
```

Run the Fe kinematic numerical validation:

```bash
.venv/bin/python validation/validate_fe_kinematic.py
```

This compares the translated closed-form single-layer kinematic calculation
with an independent brute-force finite sum over all atoms in the Fe layer.

## Command-line fits

Run a kinematic-only fit attempt for the bundled 100 Å Fe data:

```bash
.venv/bin/python validation/fit_fe_100a.py
```

Show fit progress while it runs and save a progress PNG:

```bash
.venv/bin/python validation/fit_fe_100a.py --plot-progress
```

Run a dynamic-density fit attempt for the bundled 100 Å Fe data:

```bash
.venv/bin/python validation/fit_fe_100a_dynamic.py --plot-progress
```

By default the dynamic Fe fit uses `58.92 <= 2theta <= 68.0` degrees to avoid
the sharp substrate contribution above 68 degrees. Adjust that with
`--twotheta-min` and `--twotheta-max`.

For a quick plotting smoke test, lower the optimizer budget:

```bash
.venv/bin/python validation/fit_fe_100a_dynamic.py --plot-progress --maxiter 2 --popsize 3 --local-max-nfev 20 --polish-maxiter 20
```

Benchmark the optional Numba density and transfer-matrix acceleration:

```bash
.venv/bin/python validation/benchmark_dynamic_propagation.py
```

The benchmark reports cached fused, uncached fused, and legacy timings plus
numerical differences. For the bundled 401-point Fe range on the development
machine, cached fused evaluation is about 0.008 seconds, roughly 1.4 times
faster than uncached fused and 13 times faster than legacy, with reflectivity
differences below `2e-14`.

Run the GaAs substrate-only dynamic example translated from
`matlab/kinematic_and_dynamic/run_gaas_substrate.m`:

```bash
.venv/bin/python validation/run_gaas_substrate.py
```

This writes `validation/gaas_substrate_dynamic.csv` and
`validation/gaas_substrate_dynamic.png`. For a faster coarse run:

```bash
.venv/bin/python validation/run_gaas_substrate.py --theta-step 0.5 --slices 100 --max-q0 30 --step-q0 0.05
```

## Fitting GUI

Run the GUI:

```bash
.venv/bin/python apps/genl_fit_gui.py
```

After installing the package in editable mode, the GUI can also be launched as:

```bash
.venv/bin/python -m pip install -e .
genl-fit-gui
```

The GUI includes bundled Fe, W, and V data/profile defaults:

- Fe: `Example_data_10nmFe.txt`, `Fe_fractional.vasp`, MgO substrate by default
- W: `Example_data_10nmW.txt`, `W_110_fractional.vasp`, Al2O3 substrate by default
- V: `Example_data_10nmV.txt`, `V_fractional.vasp`, MgO substrate by default

The `Experimental data file` field is the input loaded for plotting and
fitting; use `Browse...` to choose a two-column text, DAT, or CSV file. When a
known bundled Fe, W, or V data file is selected, the matching defaults are
applied automatically. The GUI plots the full 2theta range present in the data
file when it is loaded. The `X-ray wavelength (Å)` field defaults to Cu K alpha
(`1.5406 Å`) and is used by both kinematic and dynamic calculations. Use
`Horizontal axis` to switch the plotted and range-control axis between `2θ`
and `q`; input data files still use 2theta in the first column, and q limits
are converted internally. Edit the min/max fields to restrict the displayed
and fitted window. Experimental data is plotted as soon as the GUI opens and
refreshes when the data file, wavelength, horizontal axis, or range changes. Use
`Simulate` to calculate one pattern from the current parameter values and fixed
inputs without running the optimizer; use `Run fit` when the current setup is
ready to optimize. `Save setup/results...` writes a versioned `*.genl.json`
project containing the complete setup, fit selections, bounds, optimizer
settings, latest curves, progress history, electron-density profile, and fit
summary. `Load setup/results...` restores the controls and saved plots so the
simulation can be rerun or fitted again. The existing automatic CSV and PNG
outputs are unchanged. The left side keeps the run setup, Run/Stop buttons, status,
fit summary, and tabbed parameter controls visible while the monitoring plots
stay fixed on the right. Panels use restrained color accents to identify setup,
kinematic, dynamic film, dynamic fit, dynamic substrate, and roughness controls.
`Simulate`, `Run fit`, and `Stop fit` use blue, green, and red button styles,
and the status line changes color for running, complete, warning, and error
states. During fitting, the latest parameter values are written back into the
editable `value` boxes and range gauges; after completion they become the next
simulation or fit starting point. Each visible fitted parameter has a `Fit`
checkbox: checked rows are optimized, unchecked rows are held fixed at their
current value. Min/max fields are only used for checked rows; fixed rows show a
neutral `fixed` range marker. Bounded parameters also show a compact range
gauge whose marker turns blue before fitting, then green, yellow, or red
depending on whether the latest fitted value is comfortably inside, close to,
or at the selected min/max bounds. Boundary warnings are also added to the fit
summary.

The GUI can switch between kinematic and dynamic fitting. Model-specific
controls are enabled according to the selected model:

- `Kinematic parameters` tab: plane spacing, coherent planes, resolution,
  intensity scale, linear background, and Debye-Waller coefficient.
- `Dynamic film parameters` tab: dynamic-model film structure file, layer
  direction, number of layers, lattice scale, area scale, and interface
  spacing. The structure-file dropdown is populated from all `*.vasp`
  files currently present in `kinematic_and_dynamic/POSCAR`.
- `Dynamic fit parameters` tab: dynamic-model resolution, intensity scale, and
  linear background controls, plus the `Dynamic backend` selector. Use `auto`
  for normal work, `fused` to require the Numba density and fused propagation
  kernels, and `legacy` to force the previous NumPy density and vectorized
  propagation implementation.
- `Dynamic substrate setup` tab: dynamic-model substrate structure file, layer
  direction, number of layers, interface spacing, lattice scale, and area scale.
  Lattice scale has `Fit`, `Value`, `Min`, and `Max` controls; it is unchecked by
  default with limits at +/-0.5% of the loaded sample value.
- `Strain parameters` tab: bottom/top strain amplitudes and affected extents,
  each with `Fit`, `Value`, `Min`, `Max`, and range controls. Extents are in
  planes for the kinematic model and atomic positions for the dynamic model.
- `Roughness parameters` tab: film roughness and substrate/interface roughness
  value/min/max controls.

Use `Include strain` to add bottom/top strain amplitude and affected-depth
parameters to either model. The strain tab retains separate values and limits
for the kinematic and dynamic models, and its `Fit` checkboxes select which
strain parameters are optimized. Use `Include roughness` to add film
roughness and substrate/interface roughness parameters with editable
start/min/max controls.
Film roughness is modeled as a Gaussian distribution of coherent film
thicknesses, with sigma measured in planes for the kinematic model and repeated
layers for the dynamic model; substrate/interface roughness is measured in
angstrom and damps the coherent signal. Dynamic roughness follows the MATLAB
path by averaging the complex amplitudes before calculating intensity. Its
electron-density monitor shows the probability-weighted density distribution.

The GUI shows fit progress live, plots the current fit curve, and plots the
evolving electron-density profile for dynamic fits. The `Stop fit` button requests
cancellation and interrupts the fit between optimizer evaluations. GUI outputs
are written under `validation/` using sample-specific names, for example
`w_10_nm_gui_dynamic_fit.csv` and `w_10_nm_gui_dynamic_fit_progress.png`.
The plot area includes a citation watermark with the authors and DOI, and
clicking it opens the article DOI.

## Current limitations

- Substrate layer count, interface spacing, and area scale remain fixed model
  inputs; substrate lattice scale can be fitted.
- Kinematic roughness currently averages the intensities of the discrete
  thickness components; dynamic roughness uses the MATLAB coherent-amplitude
  average.
- Older validation launchers are kept as compatibility wrappers where useful.

The pre-project-export GUI, tests, and README are preserved in
`archive/pre_project_export_2026-08-16/` for rollback.
The pre-dynamic-workspace calculation files are preserved in
`archive/pre_dynamic_workspace_2026-08-17/`.
