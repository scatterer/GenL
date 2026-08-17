import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genl import (
    Control,
    Instrument,
    Layer,
    calc_kinematic,
    calc_dynamic_density,
    form_factors,
    gauss_conv,
    read_form_factor_coefficients,
    read_poscar,
)
from genl.gui import (
    FitUpdate,
    SAMPLES,
    dynamic_bounds_and_start,
    fit_update_from_dict,
    fit_update_to_dict,
    kinematic_bounds_and_start,
    load_sample_data,
    read_experimental_data,
)
from genl.fit_models import DynamicModel, roughness_distribution
import genl.dynamic as dynamic
from validation.validate_fe_kinematic import validate_fe_kinematic

DATA = ROOT / "kinematic_and_dynamic" / "Form_Factor_and_Elemental_data"
POSCAR = ROOT / "kinematic_and_dynamic" / "POSCAR"
FE_DATA = ROOT / "kinematic_and_dynamic" / "examples" / "Example_data_10nmFe.txt"


class CoreTests(unittest.TestCase):
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
            max_q0=15,
            step_q0=0.2,
        )

        self.assertEqual(result.refl.shape, q.shape)
        self.assertTrue(np.all(np.isfinite(result.refl)))
        self.assertTrue(np.all(result.refl >= 0))

    def test_fused_dynamic_backend_matches_legacy_backend(self):
        if dynamic._propagate_vectorized_fused_numba is None:
            self.skipTest("numba fused backend is unavailable")

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
            "control": Control(pol=0),
            "instrument": Instrument(theta_m=2),
            "poscar_dir": POSCAR,
            "form_factor_dir": DATA,
            "slices": 12,
            "max_q0": 10,
            "step_q0": 0.5,
        }

        fused = calc_dynamic_density(q, wavelength, stack, propagation_backend="fused", **kwargs)
        legacy = calc_dynamic_density(q, wavelength, stack, propagation_backend="legacy", **kwargs)

        np.testing.assert_allclose(fused.refl, legacy.refl, rtol=1e-8, atol=1e-13)

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
        amplitude_weights = np.sqrt(weights)
        amplitude_weights /= np.sum(amplitude_weights)
        expected_s = sum(
            weight * result.amplitude_s
            for weight, result in zip(amplitude_weights, results)
        )
        expected_p = sum(
            weight * result.amplitude_p
            for weight, result in zip(amplitude_weights, results)
        )
        mono = np.cos(np.deg2rad(4.0)) ** 2
        expected_refl = (
            np.abs(expected_s) ** 2 + mono * np.abs(expected_p) ** 2
        ) / (1.0 + mono)

        actual_refl = model.reflectivity(params)
        np.testing.assert_allclose(actual_refl, expected_refl, rtol=1e-12)

        longest = max(results, key=lambda result: len(result.z))
        expected_density = np.zeros_like(longest.rho_e)
        for weight, result in zip(weights, results):
            expected_density[: len(result.rho_e)] += weight * result.rho_e
        np.testing.assert_allclose(model.last_z, longest.z)
        np.testing.assert_allclose(model.last_rho_e, expected_density, rtol=1e-12, atol=1e-12)

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

    def test_fe_kinematic_matches_brute_force_sum(self):
        metrics = validate_fe_kinematic()

        self.assertLess(metrics["max_abs_error"], 1e-14)
        self.assertLess(metrics["max_rel_error"], 1e-9)
        self.assertLess(metrics["rms_rel_error"], 1e-10)


if __name__ == "__main__":
    unittest.main()
