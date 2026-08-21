# GenL

GenL simulates and fits X-ray reflectivity and Laue-oscillation patterns from
single-crystal films and superlattices using kinematic or dynamic scattering
models.

## Quick start

GenL requires Python 3.10 or newer. The commands below install the optional
Numba acceleration used by the dynamic model.

### macOS and Linux

```bash
git clone https://github.com/scatterer/GenL.git
cd GenL
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[fast]"
.venv/bin/python apps/genl_fit_gui.py
```

### Windows PowerShell

```powershell
git clone https://github.com/scatterer/GenL.git
cd GenL
py -m venv .venv
.venv\Scripts\python.exe -m pip install -e ".[fast]"
.venv\Scripts\python.exe apps\genl_fit_gui.py
```

If the repository is already downloaded, start with the `cd GenL` command.
Keep the editable installation (`-e`) because GenL reads the shared `data/`
folder from the repository. On Linux, install your distribution's
`python3-tk` package if Python reports that `tkinter` is unavailable.

To launch GenL again later, open a terminal in the repository and run:

```bash
.venv/bin/python apps/genl_fit_gui.py
```

On Windows, use:

```powershell
.venv\Scripts\python.exe apps\genl_fit_gui.py
```

## First simulation or fit

1. For films, use **Browse** beside **Film data** to load a two-column text,
   DAT, or CSV file. The columns must contain 2theta in degrees and intensity.
2. Choose **Kinematic** or **Dynamic** in the scattering-model selector.
3. Enter the film, substrate, wavelength, and calculation parameters.
4. Press **Simulate** to check the current values without fitting.
5. Select the **Fit** checkbox only for parameters that should be optimized,
   then set their minimum and maximum values.
6. Press **Run fit**. Use **Pause fit**, **Resume fit**, or **Stop fit** while
   monitoring the diffraction, cost, and electron-density plots.

The complete experimental range is displayed immediately after loading a
file. Use the range controls to restrict the simulation or fit, and switch the
horizontal axis between 2theta and q when needed. The default wavelength is
Cu K alpha (`1.5406 Å`).

## Python translation

The Python implementation is in `src/genl/`. It provides:

- POSCAR/VASP crystal-structure and atomic form-factor loading with
  element-specific Debye-Waller damping
- kinematic scattering for crystalline films
- density-based dynamic scattering for substrates, multilayers, and superlattices
- versioned superlattice files for repeated crystalline layers
- strain, film roughness, and substrate/interface roughness controls
- Gaussian instrumental broadening, scale, and centered quadratic background fitting
- Numba-accelerated reflection recursion and transfer-matrix propagation
- a Tk GUI for simulation, fitting, progress monitoring, and result export

Dynamic calculations use the MATLAB-v2 smooth, tail-aware analytic density
method and Numba reflection-recursion propagation. For converged Fe
calculations, start with 200 density slices per cell. The previous sampled
density implementation is retained only in the dated archive.

The dynamic fitter caches material, density, and substrate calculations during
optimization. Code using the Python API can obtain the same reuse by passing a
single `DynamicWorkspace()` to repeated `calc_dynamic_density()` calls.
The simple kinematic GUI uses its displayed Debye-Waller coefficient; the
POSCAR-based kinematic and dynamic engines use the element-specific values from
the original GenL table automatically.

## GUI controls

The film-parameter notebook is divided into four groups:

- **Kinematic**: film plane spacing, coherent planes, resolution, intensity
  scale, background, Debye-Waller coefficient, and an optional substrate peak.
- **Dynamic**: film and substrate structure files, crystal directions, layer
  counts, lattice and area scales, interface spacing, density sampling,
  resolution, intensity scale, and background.
- **Strain / roughness**: model-specific strain and film or substrate roughness
  parameters.
- **Optimization settings**: random seed, progress update interval,
  differential-evolution limits, local evaluations, and polish iterations.

GenL initializes the optimizer seed from the operating system's random source.
Use **New seed** to generate an independent differential-evolution population;
the displayed seed is retained in saved setups for reproducible reruns.
Hover over specialized scientific and optimizer controls to see a short
description. Tooltips are intentionally omitted from the self-explanatory
simulation buttons and **Value**, **Min**, and **Max** fields.

Use **Sample configuration** in the setup panel to switch between the mutually
exclusive film and superlattice workspaces. The superlattice workspace has
**Structure**, **Strain**, and **Calculation** tabs. It provides a reusable superlattice file,
substrate and repeated-layer rows, repetitions, unit-cell counts, interface
spacing, lattice scales, and axis-aware sampling. An optional capping layer is
appended after the repeated sequence and has the same editable structure and fit
controls. Bottom and top strain amplitudes and extents can be simulated or fitted
for each repeated layer and the capping layer; repeated-layer values are shared by
every occurrence in the sequence. The layer and strain parameter lists have their
own vertical scrollbars. Its calculation values share
the same resolution, intensity-scale, background, and density controls
as the main model tabs. Interface spacing, lattice scale, area scale, resolution,
intensity scale, and all three background coefficients have fit selectors, editable
bounds, and boundary gauges. Sampling density, density slices, and density Q
maximum remain value-only controls.

Film and superlattice backgrounds use a q coordinate centered and normalized to
the selected scan range. Background curvature is fixed at zero by default and is
included in fitting only when its **Fit** checkbox is selected.

Use **Load data...** in the superlattice **Structure** tab to select a two-column
`2theta/intensity` file or a three-column `2theta/q/intensity` file such as
`data/examples/Fe-V_fit.txt`. Loading sets the full available axis range; the
data are previewed in the center plot. Pressing **Simulate** overlays the
calculated superlattice pattern on the selected experimental data range. **Run
fit** uses the selected superlattice parameters and the shared optimization,
pause/resume, progress, and export workflow.

Structure-file menus are populated from the `*.vasp` files in
`data/structures/`. Each optimizable parameter has a **Fit** checkbox, current
value, bounds, and a range gauge. The gauge changes color when a fitted value
approaches or reaches a boundary. Latest fit values are written back into the
editable value fields so they can be simulated or used as the next starting
point.

Use **Include strain** or **Include roughness** to enable those models. Strain
depths are measured in planes for kinematic calculations and atomic positions
for dynamic calculations. Film roughness is a Gaussian thickness distribution;
dynamic roughness averages complex amplitudes before calculating intensity.

The GUI updates the fitted curve, objective cost, parameter values, and dynamic
electron-density profile while fitting. Pausing preserves the optimizer
checkpoint. Bounds and iteration budgets can be changed before resuming;
changing the model, data window, fixed inputs, seed, population size, or fitted
parameter selection starts a new fit.

## Saving results

- **Save setup/results** writes a versioned `*.genl.json` project containing
  the setup, selected fit parameters, bounds, optimizer settings, curves,
  progress, density profile, and fit summary. It can be loaded again later.
- **Save plots** exports separate diffraction and electron-density figures as
  PNG, JPEG, TIFF, SVG, or PDF.
- **Export graph data** writes diffraction data, fit, and residuals as CSV or
  tab-delimited text, plus a separate density table when available.

Project and superlattice JSON files are replaced atomically. If any save or
plot backend fails, the GUI reports the error and retains the current fit and
progress in memory so the user can retry. Save dialogs use single-extension
filters to avoid a native Tk 8.6 crash on macOS.

Automatic fit outputs are written under `validation/`. Exported diffraction
figures use logarithmic intensity with 2theta and q axes. The GUI also contains
a link to the GenL article in its citation watermark.

## Examples and validation

Simulate the nominal `[Fe(4 ML)/V(28 ML)] x 11` superlattice from the GenL
paper using the reusable structure definition in
`data/stacks/fe_v_4_28_x11.json`:

```bash
.venv/bin/python examples/superlattice.py
```

Use the dynamic model for the same superlattice:

```bash
.venv/bin/python examples/superlattice.py --model dynamic
```

Fit the parameters selected in the superlattice file after supplying a two-column
experimental data file:

```bash
.venv/bin/python examples/superlattice.py --data path/to/fe_v_data.txt --fit --maxiter 20 --popsize 8
```

Each run saves the calculated pattern and a fully expanded
`*_resolved_superlattice.json` record of every layer and effective fitted value.
Repeated layers with the same name share fitted parameters, so all 11 Fe/V
bilayers remain tied during optimization.

Run a minimal Fe kinematic simulation:

```bash
.venv/bin/python examples/kinematic_fe.py
```

Run the independent Fe kinematic numerical validation:

```bash
.venv/bin/python validation/validate_fe_kinematic.py
```

Run a kinematic fit of the bundled 100 Å Fe data:

```bash
.venv/bin/python validation/fit_fe_100a.py --plot-progress
```

Run the corresponding dynamic fit:

```bash
.venv/bin/python validation/fit_fe_100a_dynamic.py --plot-progress
```

The dynamic Fe fit defaults to `58.92 <= 2theta <= 68.0` degrees to exclude
the sharp substrate feature above 68 degrees. For a quick plotting test:

```bash
.venv/bin/python validation/fit_fe_100a_dynamic.py --plot-progress --maxiter 2 --popsize 3 --local-max-nfev 20 --polish-maxiter 20
```

Run the GaAs substrate-only dynamic example:

```bash
.venv/bin/python validation/run_gaas_substrate.py
```

This writes `validation/gaas_substrate_dynamic.csv` and
`validation/gaas_substrate_dynamic.png`. The Python density slicing preserves
the shared unit-cell boundary samples and prevents nonphysical low-angle
reflectivity above one. A faster coarse run is available with:

```bash
.venv/bin/python validation/run_gaas_substrate.py --theta-step 0.5 --slices 100 --max-q0 30
```

Benchmark the dynamic propagation backends:

```bash
.venv/bin/python validation/benchmark_dynamic_propagation.py
```

Check analytic-density slice convergence:

```bash
.venv/bin/python validation/validate_dynamic_density_convergence.py
```

MATLAB-v2 parity references can be exported with
`matlab/kinematic_and_dynamic/export_subroutines_v2_reference.m`, then checked
with `.venv/bin/python validation/validate_matlab_v2_density.py`.

Run all tests:

```bash
.venv/bin/python -m unittest discover -s tests
```

Windows users can replace `.venv/bin/python` in these commands with
`.venv\Scripts\python.exe`.

## Project layout

- `src/genl/`: Python scattering, fitting, and GUI package
- `apps/`: user-facing launchers
- `examples/`: small runnable Python examples
- `validation/`: numerical validation, fitting, and benchmark scripts
- `tests/`: automated tests
- `data/examples/`: bundled experimental and reference patterns
- `data/form_factors/`: atomic and scattering-factor tables
- `data/structures/`: POSCAR/VASP crystal structures
- `data/stacks/`: reusable superlattice, substrate, and fit definitions
- `archive/`: older material retained for reference, not active code

## Current limitations

- Substrate layer count, interface spacing, and area scale are fixed inputs;
  substrate lattice scale can be fitted.
- Kinematic roughness averages intensities of discrete thickness components;
  dynamic roughness averages complex amplitudes.
- Checkpoints are retained only for the current GUI session.
- Superlattice strain and roughness are not yet exposed in the GUI fitter.

Rollback snapshots are stored under `archive/`, including versions before
project export, dynamic workspaces, fit checkpoints, reflection recursion, and
the MATLAB-v2 density port.

## Citation

If you use GenL, please cite:

Anna L. Ravensburg, Johan Bylin, Vassilios Kapaklis and Gunnar K. Pálsson,
*Journal of Applied Crystallography* **59**, 968-977 (2026),
[https://doi.org/10.1107/S1600576726002566](https://doi.org/10.1107/S1600576726002566).

## License

GenL is distributed under the GNU General Public License version 3. See
`LICENSE`.

## MATLAB version

The original MATLAB implementation is retained under `matlab/`:

- `matlab/kinematic_only/`: kinematic GUI and command-line source
- `matlab/kinematic_and_dynamic/`: kinematic and dynamic command-line source,
  including fitting
- `matlab/kinematic_and_dynamic/subroutines_updated/`: updated dynamic
  reflection and transmission propagation routines
- `matlab/kinematic_and_dynamic/subroutines_v2/`: analytic density and shared
  sigma/pi propagation updates
- `matlab/binaries/`: installers for macOS Ventura, Windows 10, and Windows 11

The installers add the free MATLAB Runtime and place GenL with the installed
applications. First launch can take several minutes. When prompted, select the
GenL working folder; the default macOS application folder is
`/Applications/GenL/application/`.

For source-based use:

- Modify `matlab/kinematic_and_dynamic/run_gaas_substrate.m` to simulate the
  GaAs substrate pattern.
- Modify `matlab/kinematic_and_dynamic/fit_laue_oscillations_fe.m` to fit with
  the kinematic or dynamic model.

Starting with the kinematic model is recommended for a quick initial estimate,
followed by the dynamic model for higher accuracy. Superlattices can be edited,
simulated, and fitted in the GUI or with the command-line runner.
