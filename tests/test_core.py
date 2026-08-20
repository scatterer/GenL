import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genl import (
    Control,
    DynamicWorkspace,
    Instrument,
    Layer,
    StackDefinition,
    StackModel,
    calc_kinematic,
    calc_dynamic_density,
    form_factors,
    gauss_conv,
    read_form_factor_coefficients,
    read_poscar,
    validate_density_sampling,
)
from genl.gui import (
    FitUpdate,
    KinematicModel,
    SAMPLES,
    _least_squares_residual,
    clipped_checkpoint_population,
    dynamic_bounds_and_start,
    export_result_data,
    fit_resume_signature,
    fit_update_from_dict,
    fit_update_to_dict,
    kinematic_bounds_and_start,
    kinematic_substrate_peak,
    load_sample_data,
    read_experimental_data,
    save_result_plots,
    stack_simulation_grid,
)
from genl.fit_models import DynamicModel, roughness_distribution
import genl.dynamic as dynamic
import genl.kinematic as kinematic
from genl.paths import EXAMPLE_DATA_DIR, FORM_FACTOR_DIR, STACK_DIR, STRUCTURE_DIR
from validation.validate_fe_kinematic import validate_fe_kinematic

DATA = FORM_FACTOR_DIR
POSCAR = STRUCTURE_DIR
FE_DATA = EXAMPLE_DATA_DIR / "Example_data_10nmFe.txt"


class CoreTests(unittest.TestCase):
    def test_stack_simulation_grid_uses_selected_axis_density(self):
        twotheta, q = stack_simulation_grid(60.0, 62.0, 10.0, False, 1.5406)
        self.assertEqual(len(twotheta), 21)
        self.assertAlmostEqual(twotheta[0], 60.0)
        q_twotheta, q_values = stack_simulation_grid(q[0], q[-1], 100.0, True, 1.5406)
        self.assertEqual(len(q_values), int(np.ceil((q[-1] - q[0]) * 100.0)) + 1)
        self.assertAlmostEqual(q_twotheta[0], 60.0)

    def test_density_sampling_rejects_aliasing(self):
        self.assertAlmostEqual(validate_density_sampling(4.0, 100, 30.0), 25.0 * np.pi)
        with self.assertRaisesRegex(ValueError, "Nyquist"):
            validate_density_sampling(4.0, 20, 20.0)

    def test_dynamic_form_factors_include_matlab_debye_waller_factor(self):
        q0 = np.array([0.0, 4.0])
        prepared = dynamic._prepare_layer(
            Layer(direction=1, n=1, filename="Fe_fractional.vasp"),
            1.54056,
            q0,
            POSCAR,
            DATA,
        )
        raw, _, _ = form_factors(
            q0,
            read_form_factor_coefficients(26, 1.54056, DATA),
            1.0,
        )
        expected = raw * np.exp(-0.35 * (q0 / (4.0 * np.pi)) ** 2)

        for atom_form_factor in np.asarray(prepared["ff"]).T:
            np.testing.assert_allclose(atom_form_factor, expected)

    def test_fe_v_stack_expands_and_simulates(self):
        definition = StackDefinition.load(STACK_DIR / "fe_v_4_28_x11.json")
        layers = definition.layers()

        sequence = definition.document["sequence"]
        cap_count = int(definition.document.get("capping_layer") is not None)
        self.assertEqual(
            len(layers),
            1 + int(sequence["repetitions"]) * len(sequence["layers"]) + cap_count,
        )
        self.assertEqual([layers[1].n, layers[2].n], [13, 2])
        self.assertEqual(
            [Path(layer.filename).stem for layer in layers[1:5]],
            ["V_fractional", "Fe_fractional", "V_fractional", "Fe_fractional"],
        )
        model = StackModel(definition, np.linspace(64.0, 67.0, 31))
        predicted = model.predict()
        self.assertTrue(np.all(np.isfinite(predicted)))
        self.assertTrue(np.all(predicted >= 0))

    def test_superlattice_fit_residual_uses_selected_parameters(self):
        base = StackDefinition.load(STACK_DIR / "fe_v_4_28_x11.json")
        document = copy.deepcopy(base.document)
        document["model"] = "kinematic"
        fe_layer = next(
            layer for layer in document["sequence"]["layers"] if layer["name"] == "Fe"
        )
        interface = float(fe_layer["dinterface"])
        intensity_scale = float(document["calculation"]["intensity_scale"])
        document["fit_parameters"] = [
            {
                "target": "Fe.dinterface",
                "min": max(0.0, interface - 0.5),
                "max": interface + 0.5,
            },
            {
                "target": "calculation.intensity_scale",
                "min": intensity_scale * 0.5,
                "max": intensity_scale * 1.5,
            },
        ]
        definition = StackDefinition(base.path, document)
        twotheta = np.linspace(64.0, 64.2, 5)
        observed = StackModel(definition, twotheta, model="kinematic").predict(
            definition.start
        )
        model = StackModel(
            definition, twotheta, observed=observed, model="kinematic"
        )

        self.assertEqual(
            definition.parameter_names,
            ("Fe.dinterface", "calculation.intensity_scale"),
        )
        np.testing.assert_allclose(model.residual_vector(definition.start), 0.0)
        self.assertEqual(model.objective(definition.start), 0.0)

    def test_superlattice_capping_layer_is_appended_and_fittable(self):
        base = StackDefinition.load(STACK_DIR / "fe_v_4_28_x11.json")
        document = copy.deepcopy(base.document)
        document.pop("capping_layer", None)
        document["fit_parameters"] = []
        uncapped = StackDefinition(base.path, copy.deepcopy(document))
        document["capping_layer"] = {
            "name": "Fe cap",
            "filename": "Fe_fractional.vasp",
            "direction": 1,
            "unit_cells": 1.0,
            "dinterface": 0.5,
            "scale": 1.0,
            "area_scale": 1.0,
        }
        document["fit_parameters"] = [
            {"target": "capping.scale", "min": 0.9, "max": 1.1}
        ]
        definition = StackDefinition(base.path, document)

        layers = definition.layers(np.array([1.02]))
        self.assertEqual(len(layers), len(uncapped.layers()) + 1)
        self.assertEqual(Path(layers[-1].filename).stem, "Fe_fractional")
        self.assertAlmostEqual(layers[-1].scale, 1.02)
        self.assertEqual(
            definition.resolved_document(np.array([1.02]))["layers"][-1]["name"],
            "Fe cap",
        )

    def test_kinematic_multilayer_does_not_repeat_interface_phase(self):
        q = np.array([1.0])
        stack = [
            Layer(direction=1, n=1, filename="Fe_fractional.vasp"),
            Layer(direction=1, n=1, filename="Fe_fractional.vasp", dinterface=0.2),
            Layer(direction=1, n=1, filename="Fe_fractional.vasp", dinterface=0.7),
        ]
        film_amplitude = np.ones_like(q, dtype=complex)
        with patch.object(
            kinematic,
            "_substrate_like_amplitude",
            return_value=np.zeros_like(q, dtype=complex),
        ), patch.object(
            kinematic,
            "_layer_amplitude",
            side_effect=[
                (film_amplitude, np.array([0.0, 1.0])),
                (film_amplitude, np.array([2.0, 3.0])),
            ],
        ):
            result = calc_kinematic(
                q,
                1.54056,
                stack,
                Control(pol=0),
                poscar_dir=POSCAR,
                form_factor_dir=DATA,
            )

        np.testing.assert_allclose(result.refl, [4.0])

    def test_save_result_plots_writes_separate_nonblank_images(self):
        update = FitUpdate(
            phase="best fit",
            cost=0.1,
            twotheta=np.linspace(58.0, 62.0, 20),
            q=np.linspace(3.9, 4.2, 20),
            observed=np.geomspace(1.0, 100.0, 20),
            predicted=np.geomspace(1.2, 90.0, 20),
            params=np.array([1.0]),
            density_z=np.linspace(0.0, 100.0, 30),
            density_rho_e=np.linspace(0.0, 2.0, 30).astype(complex),
        )
        with tempfile.TemporaryDirectory() as directory:
            paths = save_result_plots(update, 1.5406, Path(directory) / "result.png")

            self.assertEqual(
                [path.name for path in paths],
                ["result_diffraction.png", "result_density.png"],
            )
            self.assertTrue(all(path.stat().st_size > 10_000 for path in paths))

    def test_export_result_data_writes_diffraction_and_density_tables(self):
        update = FitUpdate(
            phase="simulation",
            cost=0.1,
            twotheta=np.array([60.0, 61.0]),
            q=np.array([4.0, 4.1]),
            observed=np.array([10.0, 20.0]),
            predicted=np.array([11.0, 19.0]),
            params=np.array([1.0]),
            density_z=np.array([0.0, 1.0, 2.0]),
            density_rho_e=np.array([2.0 + 0.1j, 3.0 + 0.2j, 4.0 + 0.3j]),
        )
        with tempfile.TemporaryDirectory() as directory:
            paths = export_result_data(update, Path(directory) / "result.csv")

            self.assertEqual([path.name for path in paths], ["result.csv", "result_density.csv"])
            self.assertEqual(np.loadtxt(paths[0], delimiter=",", skiprows=1).shape, (2, 5))
            self.assertEqual(np.loadtxt(paths[1], delimiter=",", skiprows=1).shape, (3, 3))

    def test_least_squares_residual_matches_displayed_cost(self):
        residual = np.array([-2.0, -0.25, 0.0, 1.0])
        transformed = _least_squares_residual(residual)

        self.assertAlmostEqual(
            float(np.sum(transformed**2)),
            float(np.mean(np.abs(residual))),
        )

    def test_fit_checkpoint_clips_bounds_and_rejects_fixed_input_changes(self):
        bounds = np.array([[0.0, 1.0], [2.0, 4.0]])
        population = np.array([[-1.0, 3.0], [0.5, 5.0]])
        np.testing.assert_allclose(
            clipped_checkpoint_population(population, bounds),
            [[0.0, 3.0], [0.5, 4.0]],
        )

        config = {
            "sample_profile": "Fe 10 nm",
            "data_path": FE_DATA,
            "model": "Dynamic",
            "wavelength": 1.5406,
            "dynamic_backend": "auto",
            "density_slices": 100,
            "density_max_q0": 30.0,
            "twotheta_min": 58.92,
            "twotheta_max": 68.0,
            "include_strain": False,
            "include_roughness": False,
            "seed": 1,
            "popsize": 6,
            "film_settings": {"filename": "Fe_fractional.vasp", "direction": 1},
            "substrate_settings": {
                "filename": "MgO_001_fractional.vasp",
                "direction": 1,
                "n": 1e6,
                "dinterface": 0.0,
                "area_scale": 1.0,
            },
            "kinematic_settings": {"debye_waller_coeff": -0.34},
        }
        mask = np.array([True, False])
        signature = fit_resume_signature(config, np.array([1.0, 2.0]), mask)
        self.assertEqual(
            signature,
            fit_resume_signature(config, np.array([1.5, 2.0]), mask),
        )
        self.assertNotEqual(
            signature,
            fit_resume_signature(config, np.array([1.0, 2.1]), mask),
        )
        changed_sampling = dict(config, density_slices=120)
        self.assertNotEqual(
            signature,
            fit_resume_signature(changed_sampling, np.array([1.0, 2.0]), mask),
        )

    def test_gui_result_json_round_trip(self):
        update = FitUpdate(
            phase="best fit",
            cost=0.125,
            twotheta=np.array([60.0, 61.0]),
            q=np.array([4.0, 4.1]),
            observed=np.array([10.0, 20.0]),
            predicted=np.array([11.0, 19.0]),
            params=np.array([1.0, 2.0]),
            density_z=np.array([0.0, 1.0]),
            density_rho_e=np.array([2.0 + 0.5j, 3.0 - 0.25j]),
        )

        restored = fit_update_from_dict(json.loads(json.dumps(fit_update_to_dict(update))))

        self.assertEqual(restored.phase, update.phase)
        self.assertEqual(restored.cost, update.cost)
        np.testing.assert_allclose(restored.predicted, update.predicted)
        np.testing.assert_allclose(restored.params, update.params)
        np.testing.assert_allclose(restored.density_rho_e, update.density_rho_e)

    def test_read_poscar_fe(self):
        structure = read_poscar("Fe_fractional.vasp", POSCAR)

        self.assertEqual(structure.types, ("Fe",))
        self.assertEqual(structure.type_counts.tolist(), [2])
        np.testing.assert_allclose(structure.a1, [2.866, 0.0, 0.0])
        self.assertEqual(structure.positions.shape, (2, 3))

    def test_form_factors_fe_are_complex_and_vectorized(self):
        q = np.array([1.0, 2.0, 3.0])
        coeffs = read_form_factor_coefficients(26, 1.54056, DATA)
        f, f_sqrd_real, f_av_sqrd_real = form_factors(q, coeffs, 1.0)

        self.assertEqual(f.shape, q.shape)
        self.assertTrue(np.iscomplexobj(f))
        self.assertTrue(np.all(f_sqrd_real > 0))
        self.assertTrue(np.all(f_av_sqrd_real > 0))

    def test_gauss_conv_preserves_constant_signal(self):
        x = np.linspace(0.0, 1.0, 101)
        y = np.ones_like(x) * 7.0

        convolved = gauss_conv(x, y, 0.05)

        np.testing.assert_allclose(convolved, y)

    def test_calc_kinematic_fe_returns_finite_intensity(self):
        theta = np.arange(35.0, 45.0, 0.1)
        wavelength = 1.54056
        q = 4.0 * np.pi / wavelength * np.sin(np.deg2rad(theta))
        layer = Layer(direction=3, n=10, filename="Fe_fractional.vasp")

        result = calc_kinematic(
            q,
            wavelength,
            [layer],
            Control(pol=0),
            Instrument(theta=theta),
            poscar_dir=POSCAR,
            form_factor_dir=DATA,
        )

        self.assertEqual(result.refl.shape, q.shape)
        self.assertTrue(np.all(np.isfinite(result.refl)))
        self.assertTrue(np.all(result.refl >= 0))

    def test_calc_dynamic_density_returns_finite_reflectivity(self):
        wavelength = 1.54056
        q = np.linspace(4.0, 4.8, 20)
        stack = [
            Layer(direction=1, n=1e4, filename="MgO_001_fractional.vasp"),
            Layer(
                direction=1,
                n=12,
                filename="Fe_fractional.vasp",
                dinterface=1.4,
                scale=1.04,
                area_scale=1.19,
            ),
        ]

        result = calc_dynamic_density(
            q,
            wavelength,
            stack,
            Control(pol=2),
            Instrument(theta_m=2),
            poscar_dir=POSCAR,
            form_factor_dir=DATA,
            slices=20,
            max_q0=14,
            step_q0=0.2,
        )

        self.assertEqual(result.refl.shape, q.shape)
        self.assertTrue(np.all(np.isfinite(result.refl)))
        self.assertTrue(np.all(result.refl >= 0))

    def test_dynamic_backends_match(self):
        if (
            dynamic._propagate_vectorized_fused_numba is None
            or dynamic._prepare_substrate_reflection_pair_numba is None
        ):
            self.skipTest("numba dynamic backends are unavailable")

        wavelength = 1.54056
        q = np.linspace(4.0, 4.8, 8)
        stack = [
            Layer(direction=1, n=1e4, filename="MgO_001_fractional.vasp"),
            Layer(
                direction=1,
                n=8,
                filename="Fe_fractional.vasp",
                dinterface=1.4,
                scale=1.04,
                area_scale=1.19,
            ),
        ]
        kwargs = {
            "control": Control(pol=2),
            "instrument": Instrument(theta_m=2),
            "poscar_dir": POSCAR,
            "form_factor_dir": DATA,
            "slices": 12,
            "max_q0": 8,
            "step_q0": 0.5,
        }

        fused = calc_dynamic_density(q, wavelength, stack, propagation_backend="fused", **kwargs)
        reflection = calc_dynamic_density(
            q, wavelength, stack, propagation_backend="reflection", **kwargs
        )
        legacy = calc_dynamic_density(q, wavelength, stack, propagation_backend="legacy", **kwargs)

        np.testing.assert_allclose(fused.refl, legacy.refl, rtol=1e-8, atol=1e-13)
        np.testing.assert_allclose(reflection.refl, fused.refl, rtol=1e-10, atol=1e-13)
        np.testing.assert_allclose(reflection.amplitude_s, fused.amplitude_s, rtol=1e-10, atol=1e-13)
        np.testing.assert_allclose(reflection.amplitude_p, fused.amplitude_p, rtol=1e-10, atol=1e-13)

    def test_analytic_density_matches_backends_and_reuses_kernels(self):
        if dynamic._prepare_substrate_reflection_pair_numba is None:
            self.skipTest("numba reflection backend is unavailable")
        q = np.linspace(2.9, 3.1, 31)
        stack = [
            Layer(direction=1, n=1e4, filename="MgO_001_fractional.vasp"),
            Layer(
                direction=1,
                n=8,
                filename="Fe_fractional.vasp",
                dinterface=1.4,
                bottom_strain_amplitude=0.02,
                bottom_strain_end=3,
            ),
        ]
        kwargs = {
            "poscar_dir": POSCAR,
            "form_factor_dir": DATA,
            "slices": 50,
            "max_q0": 30.0,
            "density_method": "analytic",
        }
        workspace = DynamicWorkspace()
        reflection = calc_dynamic_density(
            q,
            1.54056,
            stack,
            propagation_backend="reflection",
            workspace=workspace,
            **kwargs,
        )
        material_builds = workspace.material_builds
        kernel_count = len(workspace.atomic_kernels)
        cached = calc_dynamic_density(
            q,
            1.54056,
            stack,
            propagation_backend="reflection",
            workspace=workspace,
            **kwargs,
        )
        fused = calc_dynamic_density(
            q, 1.54056, stack, propagation_backend="fused", **kwargs
        )
        legacy = calc_dynamic_density(
            q, 1.54056, stack, propagation_backend="legacy", **kwargs
        )

        self.assertEqual(workspace.material_builds, material_builds)
        self.assertEqual(len(workspace.atomic_kernels), kernel_count)
        self.assertEqual(reflection.diagnostics["density_method"], "analytic")
        self.assertTrue(reflection.diagnostics["vacuum_was_extended"])
        self.assertTrue(np.all(np.isfinite(reflection.rho_e)))
        np.testing.assert_allclose(cached.refl, reflection.refl, rtol=0, atol=0)
        np.testing.assert_allclose(fused.refl, reflection.refl, rtol=1e-9, atol=1e-11)
        np.testing.assert_allclose(legacy.refl, reflection.refl, rtol=1e-9, atol=1e-11)

    def test_gaas_low_angle_reflectivity_is_passive(self):
        angle = 0.69
        q = 4.0 * np.pi / 1.54056 * np.sin(np.deg2rad(angle))
        result = calc_dynamic_density(
            np.array([q]),
            1.54056,
            [
                Layer(
                    direction=3,
                    n=1e8,
                    filename="GaAs_alt_fractional.vasp",
                    scale=1.001,
                    area_scale=1.001,
                )
            ],
            control=Control(pol=0, model="density"),
            instrument=Instrument(theta_m=2),
            poscar_dir=POSCAR,
            form_factor_dir=DATA,
            vacuum_thick=20,
            slices=400,
            max_q0=75,
            step_q0=0.01,
            propagation_backend="reflection",
        )

        self.assertLessEqual(result.refl[0], 1.0)

    def test_reflection_matches_fused_for_bundled_materials_and_polarizations(self):
        if dynamic._prepare_substrate_reflection_pair_numba is None:
            self.skipTest("numba reflection backend is unavailable")
        cases = [
            [
                Layer(direction=1, n=1e4, filename="MgO_001_fractional.vasp"),
                Layer(direction=1, n=8, filename="V_fractional.vasp", dinterface=1.4),
            ],
            [
                Layer(direction=1, n=1e4, filename="Al2O3_11-20_fractional.vasp"),
                Layer(direction=3, n=8, filename="W_110_fractional.vasp", dinterface=1.4),
            ],
            [Layer(direction=3, n=1e4, filename="GaAs_alt_fractional.vasp", scale=1.001)],
        ]
        kwargs = {
            "instrument": Instrument(theta_m=2),
            "poscar_dir": POSCAR,
            "form_factor_dir": DATA,
            "slices": 16,
            "max_q0": 8,
            "step_q0": 1.0,
        }
        q = np.linspace(3.5, 6.0, 5)

        for stack in cases:
            for pol in (0, 1, 2):
                reflection = calc_dynamic_density(
                    q,
                    1.54056,
                    stack,
                    control=Control(pol=pol),
                    propagation_backend="reflection",
                    **kwargs,
                )
                fused = calc_dynamic_density(
                    q,
                    1.54056,
                    stack,
                    control=Control(pol=pol),
                    propagation_backend="fused",
                    **kwargs,
                )
                np.testing.assert_allclose(reflection.refl, fused.refl, rtol=1e-10, atol=1e-13)

    def test_reflection_backend_reuses_substrate_reflection(self):
        if dynamic._prepare_substrate_reflection_pair_numba is None:
            self.skipTest("numba reflection backend is unavailable")
        q = np.linspace(4.0, 4.8, 8)
        workspace = DynamicWorkspace()
        kwargs = {
            "control": Control(pol=2),
            "instrument": Instrument(theta_m=2),
            "poscar_dir": POSCAR,
            "form_factor_dir": DATA,
            "slices": 12,
            "max_q0": 8,
            "step_q0": 0.5,
            "propagation_backend": "reflection",
            "workspace": workspace,
        }
        stack = [
            Layer(direction=1, n=1e4, filename="MgO_001_fractional.vasp"),
            Layer(direction=1, n=8, filename="Fe_fractional.vasp", dinterface=1.4),
        ]

        first = calc_dynamic_density(q, 1.54056, stack, **kwargs)
        second = calc_dynamic_density(q, 1.54056, stack, **kwargs)

        np.testing.assert_allclose(second.refl, first.refl, rtol=0.0, atol=0.0)
        self.assertEqual(workspace.substrate_reflection_builds, 1)

    def test_dynamic_workspace_reuses_and_invalidates_substrate(self):
        if dynamic._prepare_substrate_transfer_pair_numba is None:
            self.skipTest("numba fused backend is unavailable")
        wavelength = 1.54056
        q = np.linspace(4.0, 4.8, 8)
        workspace = DynamicWorkspace()
        kwargs = {
            "control": Control(pol=2),
            "instrument": Instrument(theta_m=2),
            "poscar_dir": POSCAR,
            "form_factor_dir": DATA,
            "slices": 12,
            "max_q0": 8,
            "step_q0": 0.5,
            "propagation_backend": "fused",
            "workspace": workspace,
        }

        def stack(substrate_scale=1.0, film_n=8):
            return [
                Layer(
                    direction=1,
                    n=1e4,
                    filename="MgO_001_fractional.vasp",
                    scale=substrate_scale,
                ),
                Layer(
                    direction=1,
                    n=film_n,
                    filename="Fe_fractional.vasp",
                    dinterface=1.4,
                    scale=1.04,
                    area_scale=1.19,
                ),
            ]

        first = calc_dynamic_density(q, wavelength, stack(), **kwargs)
        cached = calc_dynamic_density(q, wavelength, stack(), **kwargs)
        np.testing.assert_allclose(cached.refl, first.refl, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(cached.rho_e, first.rho_e, rtol=0.0, atol=0.0)
        calc_dynamic_density(q, wavelength, stack(film_n=10), **kwargs)
        self.assertEqual(workspace.material_builds, 2)
        self.assertEqual(workspace.substrate_builds, 1)
        self.assertEqual(workspace.substrate_transfer_builds, 1)

        calc_dynamic_density(q, wavelength, stack(substrate_scale=1.001), **kwargs)
        self.assertEqual(workspace.material_builds, 2)
        self.assertEqual(workspace.substrate_builds, 2)
        self.assertEqual(workspace.substrate_transfer_builds, 2)

    def test_dynamic_roughness_averages_amplitudes_and_density(self):
        model = DynamicModel(
            np.linspace(64.0, 65.0, 4),
            np.ones(4),
            include_roughness=True,
        )
        params = np.array(
            [28.5, 1.04, 1.19, 1.4, 0.005, 5000.0, 0.0, 0.0, 1.0, 0.5, 0.0]
        )
        n_values, weights = roughness_distribution(params[0], params[9])
        self.assertEqual(n_values.tolist(), [28, 29, 30, 31])

        results = [
            model._single_result(params, n, 0.0, 0.0, 0.0, 0.0)
            for n in n_values
        ]
        expected_s = sum(
            weight * result.amplitude_s
            for weight, result in zip(weights, results)
        )
        expected_p = sum(
            weight * result.amplitude_p
            for weight, result in zip(weights, results)
        )
        mono = np.cos(np.deg2rad(4.0)) ** 2
        expected_refl = (
            np.abs(expected_s) ** 2 + mono * np.abs(expected_p) ** 2
        ) / (1.0 + mono)

        actual_refl = model.reflectivity(params)
        np.testing.assert_allclose(actual_refl, expected_refl, rtol=1e-12)
        self.assertEqual(model.workspace.material_builds, 2)
        self.assertEqual(model.workspace.substrate_builds, 1)
        self.assertEqual(model.workspace.substrate_reflection_builds, 1)

        longest = max(results, key=lambda result: len(result.z))
        expected_density = np.zeros_like(longest.rho_e)
        for weight, result in zip(weights, results):
            expected_density[: len(result.rho_e)] += weight * result.rho_e
        np.testing.assert_allclose(model.last_z, longest.z)
        np.testing.assert_allclose(model.last_rho_e, expected_density, rtol=1e-12, atol=1e-12)

    def test_cached_paired_propagation_matches_legacy_with_strain_and_roughness(self):
        if (
            dynamic._prepare_substrate_transfer_pair_numba is None
            or dynamic._prepare_substrate_reflection_pair_numba is None
        ):
            self.skipTest("numba dynamic backends are unavailable")
        twotheta = np.linspace(64.0, 65.0, 4)
        fused = DynamicModel(
            twotheta,
            np.ones(4),
            include_strain=True,
            include_roughness=True,
            propagation_backend="fused",
        )
        legacy = DynamicModel(
            twotheta,
            np.ones(4),
            include_strain=True,
            include_roughness=True,
            propagation_backend="legacy",
        )
        reflection = DynamicModel(
            twotheta,
            np.ones(4),
            include_strain=True,
            include_roughness=True,
            propagation_backend="reflection",
        )
        params = np.array(
            [
                28.5,
                1.04,
                1.19,
                1.4,
                0.005,
                5000.0,
                0.0,
                0.0,
                1.0,
                0.02,
                3.0,
                -0.01,
                2.0,
                0.5,
                0.0,
            ]
        )

        np.testing.assert_allclose(
            fused.reflectivity(params),
            legacy.reflectivity(params),
            rtol=1e-8,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            reflection.reflectivity(params),
            fused.reflectivity(params),
            rtol=1e-9,
            atol=1e-12,
        )
        np.testing.assert_allclose(fused.last_rho_e, legacy.last_rho_e, rtol=1e-10, atol=1e-10)
        self.assertEqual(reflection.workspace.substrate_reflection_builds, 1)

    def test_gui_strain_settings_feed_both_model_parameter_vectors(self):
        strain = {
            "bottom_amplitude": (0.03, -0.1, 0.1),
            "bottom_extent": (4.0, 1.0, 10.0),
            "top_amplitude": (-0.02, -0.08, 0.08),
            "top_extent": (3.0, 0.0, 8.0),
        }
        expected_start = [0.03, 4.0, -0.02, 3.0]
        expected_bounds = [[-0.1, 0.1], [1.0, 10.0], [-0.08, 0.08], [0.0, 8.0]]
        substrate = {"scale": (1.002, 0.997, 1.007)}

        kin_bounds, kin_start = kinematic_bounds_and_start(
            SAMPLES["Fe 10 nm"], True, False, strain_settings=strain
        )
        dyn_bounds, dyn_start = dynamic_bounds_and_start(
            SAMPLES["Fe 10 nm"], True, False,
            strain_settings=strain,
            substrate_settings=substrate,
        )

        np.testing.assert_allclose(kin_start[6:10], expected_start)
        np.testing.assert_allclose(kin_bounds[6:10], expected_bounds)
        self.assertEqual(dyn_start[8], 1.002)
        np.testing.assert_allclose(dyn_bounds[8], [0.997, 1.007])
        np.testing.assert_allclose(dyn_start[9:13], expected_start)
        np.testing.assert_allclose(dyn_bounds[9:13], expected_bounds)

    def test_kinematic_substrate_peak_is_in_parameter_vector_and_prediction(self):
        sample = SAMPLES["Fe 10 nm"]
        center = 65.0
        d_spacing = 1.5406 / (2.0 * np.sin(np.deg2rad(center / 2.0)))
        substrate = {
            "intensity": (25.0, 0.0, 100.0),
            "width": (0.08, 0.01, 0.2),
            "d_spacing": (d_spacing, d_spacing * 0.999, d_spacing * 1.001),
        }
        bounds, start = kinematic_bounds_and_start(
            sample,
            False,
            False,
            include_substrate=True,
            substrate_peak_settings=substrate,
        )
        twotheta = np.linspace(64.5, 65.5, 101)
        observed = np.ones_like(twotheta)
        with_substrate = KinematicModel(
            twotheta, observed, sample, include_substrate=True
        ).predict(start)
        film_only = KinematicModel(twotheta, observed, sample).predict(start[:6])

        np.testing.assert_allclose(start[6:9], [25.0, 0.08, d_spacing])
        np.testing.assert_allclose(
            bounds[6:9],
            [[0.0, 100.0], [0.01, 0.2], [d_spacing * 0.999, d_spacing * 1.001]],
        )
        np.testing.assert_allclose(
            with_substrate - film_only,
            kinematic_substrate_peak(twotheta, 25.0, 0.08, d_spacing, 1.5406),
        )
        self.assertEqual(int(np.argmax(with_substrate - film_only)), 50)

    def test_gui_data_loader_reads_full_file_before_windowing(self):
        twotheta, observed = read_experimental_data(FE_DATA)
        window_twotheta, window_observed = load_sample_data(
            float(np.min(twotheta)),
            float(np.max(twotheta)),
            FE_DATA,
        )

        self.assertEqual(len(twotheta), len(window_twotheta))
        np.testing.assert_allclose(observed, window_observed)

        restricted_twotheta, _ = load_sample_data(58.92, 68.0, FE_DATA)
        self.assertLess(len(restricted_twotheta), len(twotheta))
        self.assertGreaterEqual(float(np.min(restricted_twotheta)), 58.92)
        self.assertLessEqual(float(np.max(restricted_twotheta)), 68.0)

        with tempfile.TemporaryDirectory() as directory:
            three_column_path = Path(directory) / "superlattice.txt"
            np.savetxt(
                three_column_path,
                [[50.0, 3.44, 0.036], [50.1, 3.45, 0.042]],
            )
            three_column_twotheta, three_column_intensity = read_experimental_data(
                three_column_path
            )
        np.testing.assert_allclose(three_column_twotheta, [50.0, 50.1])
        np.testing.assert_allclose(three_column_intensity, [0.036, 0.042])

    def test_fe_kinematic_matches_brute_force_sum(self):
        metrics = validate_fe_kinematic()

        self.assertLess(metrics["max_abs_error"], 1e-14)
        self.assertLess(metrics["max_rel_error"], 1e-9)
        self.assertLess(metrics["rms_rel_error"], 1e-10)


if __name__ == "__main__":
    unittest.main()
