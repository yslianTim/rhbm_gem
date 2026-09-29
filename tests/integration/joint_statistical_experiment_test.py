"""Contract tests for the fixed-seed Joint statistical research experiment."""
import hashlib
import unittest

import numpy as np
import joint_statistical_experiment as experiment


class StatisticalExperimentTest(unittest.TestCase):
    def test_fixed_matrix_and_paired_seed_order(self):
        cases = list(experiment.conditions())
        self.assertEqual(len(cases), 326)
        self.assertEqual(sum(len(case[-1]) for case in cases), 450)
        self.assertEqual(experiment.SEEDS, list(range(20260921, 20260941)))
        for case in cases:
            if case[4] is not None:
                self.assertIn(case[4], experiment.SEEDS)
            if 'production' in case[-1]:
                self.assertTrue(case[3] in (0., .05))

    def test_seeded_noise_matches_the_retained_generator_stream(self):
        expected = {
            'iid': 'b4975f62fded3e787d44433e81f6b29234e3df69e5150d30dbf50952099f6375',
            'correlated': '86780b5eca52108530d0778ecdcd266f65c03860356b7c50a75f2e4bc6debf76',
        }
        for kind, digest in expected.items():
            values = experiment.noise_field(experiment.SEEDS[0], kind)
            self.assertEqual(hashlib.sha256(values.astype('<f8').tobytes()).hexdigest(), digest)
            np.testing.assert_array_equal(values, experiment.noise_field(experiment.SEEDS[0], kind))

    def test_historical_generated_observation_hash_is_unchanged(self):
        clean = experiment.clean_map(1.2, 0.)
        scale = float(np.sqrt(np.mean(clean[experiment.TARGET] ** 2)))
        values = clean + .01 * scale * experiment.noise_field(experiment.SEEDS[0], 'iid')
        self.assertEqual(hashlib.sha256(values.astype('<f8').tobytes()).hexdigest(),
                         'c5147ebda1fe72d1012b49d91838743ce4e8ca00abf4b3f8d3dcfa233fba02e1')

    def test_summary_separates_execution_and_numerical_denominators(self):
        available_pass = dict(state_available=True, runtime_convergence='passed', errors=[[1.]],
                              failure_reasons=[], prediction_rmse=None, residual_rmse=None,
                              residual_neighbor_correlation=None)
        available_fail = dict(state_available=True, runtime_convergence='failed', errors=[[3.]],
                              failure_reasons=['rank:failed'], prediction_rmse=None, residual_rmse=None,
                              residual_neighbor_correlation=None)
        unavailable = dict(state_available=False, runtime_convergence='unavailable', errors=None,
                           failure_reasons=[], prediction_rmse=None, residual_rmse=None,
                           residual_neighbor_correlation=None)
        records = [
            dict(distance=1.2, shift=0., noise='iid', sigma_fraction=.05, modes=['fixed'],
                 status='completed', results={'fixed': available_pass}),
            dict(distance=1.2, shift=0., noise='iid', sigma_fraction=.05, modes=['fixed'],
                 status='completed', results={'fixed': available_fail}),
            dict(distance=1.2, shift=0., noise='iid', sigma_fraction=.05, modes=['fixed'],
                 status='completed', results={'fixed': unavailable}),
            dict(distance=1.2, shift=0., noise='iid', sigma_fraction=.05, modes=['fixed'],
                 status='process-failure'),
        ]
        row = experiment.statistical_summary(records)[0]
        self.assertEqual((row['attempted_cases'], row['completed_cases'], row['qualified_cases'],
                          row['failed_cases'], row['unavailable_cases']), (4, 3, 1, 1, 1))
        self.assertEqual(row['execution_status_counts'], {'completed': 3, 'process-failure': 1})
        self.assertEqual(row['process_incomplete'], 1)
        self.assertEqual(row['all_available']['n'], 2)
        self.assertEqual(row['converged_only']['n'], 1)
        self.assertEqual(row['nonconverged_only']['n'], 1)

    def test_shift_changes_data_without_changing_target_selection(self):
        unshifted = experiment.clean_map(1.2, 0.)
        shifted = experiment.clean_map(1.2, .15)
        self.assertGreater(float(np.linalg.norm(unshifted - shifted)), 0)
        self.assertEqual(int(experiment.TARGET.sum()), 515)
        self.assertTrue(np.isfinite(shifted).all())

    def test_missing_state_stays_missing(self):
        outcome = dict(assembled_state=None, runtime_convergence='unavailable',
                       initialization={'valid': False}, components=[])
        metric = experiment.metrics(outcome, None, np.zeros(1), np.zeros(1), np.array([0]))
        self.assertIsNone(metric['errors'])
        self.assertIsNone(metric['parameters'])
        self.assertFalse(metric['state_available'])


if __name__ == '__main__':
    unittest.main()
