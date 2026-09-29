"""Unit tests for neutral experiment support helpers."""
import math
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

import experiment_io as io
import experiment_process as process
import experiment_provenance as provenance
import joint_validation_checks as checks


class ExperimentIoTest(unittest.TestCase):
    def test_json_round_trip_and_hashes_are_stable(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'nested/result.json'
            expected = {'z': [1, 2], 'a': {'value': True}}
            io.write(path, expected)
            self.assertEqual(io.read(path), expected)
            self.assertEqual(io.sha(path), io.sha256_file(path))
            self.assertEqual(io.digest({'a': 1, 'b': 2}), io.digest({'b': 2, 'a': 1}))

    def test_missing_file_hash_keeps_file_error(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                io.sha256_file(Path(directory) / 'missing')


class ExperimentProcessTest(unittest.TestCase):
    class FakeProcess:
        pid = 1234

        def __init__(self, exit_code):
            self.returncode = exit_code

        def poll(self):
            return self.returncode

        def wait(self, timeout=None):
            if self.returncode is None:
                self.returncode = -15
            return self.returncode

    def launch(self, exit_code):
        def start(command, **kwargs):
            kwargs['stdout'].write('out\n')
            kwargs['stderr'].write('err\n123 maximum resident set size\n')
            return self.FakeProcess(exit_code)
        return start

    def run_process(self, directory, exit_code=0, **kwargs):
        with patch.object(process.subprocess, 'Popen', side_effect=self.launch(exit_code)):
            return process.monitored([sys.executable, '-c', 'pass'], directory,
                                     time.monotonic() + 5, **kwargs)

    def test_completed_process_captures_output(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_process(directory)
            self.assertEqual(result['status'], 'completed')
            self.assertEqual(result['exit_code'], 0)
            self.assertIn('out', (Path(directory) / 'stdout.txt').read_text())
            self.assertIn('err', (Path(directory) / 'stderr.txt').read_text())

    def test_nonzero_exit_remains_process_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_process(directory, exit_code=7)
            self.assertEqual(result['status'], 'process-failure')
            self.assertEqual(result['exit_code'], 7)

    def test_timeout_remains_time_limit(self):
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(process.subprocess, 'Popen', side_effect=self.launch(None)), \
                patch.object(process, 'process_tree_rss', return_value=0), \
                patch.object(process.os, 'killpg'):
            result = process.monitored([sys.executable, '-c', 'pass'], directory,
                                       time.monotonic() + 5, seconds=.15)
            self.assertEqual(result['status'], 'time-limit')

    def test_rss_limit_remains_distinct(self):
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(process.subprocess, 'Popen', side_effect=self.launch(None)), \
                patch.object(process, 'process_tree_rss', return_value=1024 * 1024), \
                patch.object(process.os, 'killpg'):
            result = process.monitored([sys.executable, '-c', 'pass'], directory,
                                       time.monotonic() + 5, rss_limit=1)
            self.assertEqual(result['status'], 'rss-limit')

    def test_exhausted_budget_is_not_run(self):
        with tempfile.TemporaryDirectory() as directory:
            result = process.monitored([sys.executable, '-c', 'pass'], directory,
                                       time.monotonic() - 1)
            self.assertEqual(result['status'], 'not-run-budget')


class ProvenanceTest(unittest.TestCase):
    def test_source_hash_is_repeatable_and_tracks_file_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'src').mkdir()
            (root / 'CMakeLists.txt').write_text('project(test)')
            (root / 'src/model.cpp').write_text('first')
            original = provenance.source_hash(root)
            self.assertEqual(provenance.source_hash(root), original)
            (root / 'src/model.cpp').write_text('changed')
            self.assertNotEqual(provenance.source_hash(root), original)


class ValidationPrimitiveTest(unittest.TestCase):
    def test_scaled_comparison_preserves_shape_and_finite_requirements(self):
        self.assertEqual(checks.scaled([1, 2], [1, 2]), 0.0)
        self.assertIsNotNone(checks.scaled([1, 2], [1, 2 + 1e-12]))
        self.assertIsNone(checks.scaled([1], [1, 2]))
        self.assertIsNone(checks.scaled([math.nan], [0]))
        self.assertIsNone(checks.scaled([math.inf], [0]))

    def test_vector_comparison_keeps_explicit_tolerances(self):
        self.assertTrue(checks.vector_close([1.0], [1.0]))
        self.assertTrue(checks.vector_close([1.0], [1.0 + 1e-9], relative=1e-8))
        self.assertFalse(checks.vector_close([1.0], [1.0 + 1e-5], relative=1e-8))
        self.assertFalse(checks.vector_close([math.nan], [math.nan]))
        self.assertFalse(checks.vector_close([1.0], [1.0, 2.0]))


class JointComparisonTest(unittest.TestCase):
    @staticmethod
    def endpoint():
        return dict(valid=True, reason='qualified-inner', feasible=True, free_rank=2,
                    kkt_passed=True, active_atoms=[], beta=[2., .2], objective=.1,
                    relative_residual=.2, b_gradient=[1e-13])

    @staticmethod
    def fixed_action():
        return dict(valid=True, reason='full-profile-operator', mode='normal', normal_q_actions=2,
                    state_control=dict(eta=[.2], beta=[1., .1], free_columns=[0, 1], scale=2.,
                                       rank_rows=10, residual=[.1, .2], gradient=[.03], objective=.025),
                    rank=[dict(valid=True, rank=2, rows=2, columns=2, relative_threshold=1e-12,
                               absolute_override=-1, threshold=1e-12, singular_values=[1., .5])],
                    apply=[.1, .2], adjoint=[.3], normal=[.4], operator_gradient=[.03],
                    work=dict(reference_solves=0, derivative_preparations=0),
                    steps=[dict(kind='schwarz', valid=True, step=[.2], predicted=.001,
                                true_residual=1e-11)])

    @staticmethod
    def search_result():
        return dict(stage='complete', search_completed=True, search=dict(residual_scale=2),
                    returned_assessment=dict(runtime_convergence='passed', runtime_failure='none',
                        runtime_checks=dict(inner=True), design_spectrum=dict(rank=2),
                        width_spectrum=dict(rank=1), profile_jacobian_spectrum=dict(rank=1)),
                    returned_state=dict(objective=.04, beta=[2., .2], b=[.5], active_atoms=[]))

    def test_sparse_and_compact_parity_keep_rank_spectrum_and_unavailability(self):
        endpoint = self.endpoint()
        row = dict(primary=endpoint, reference=endpoint)
        self.assertTrue(all(check['passed'] for check in checks.sparse_parity(row, row).values()))
        changed = dict(primary={**endpoint, 'free_rank': 1}, reference=endpoint)
        self.assertFalse(checks.sparse_parity(row, changed)['primary']['passed'])
        spectrum = dict(valid=True, rank=2, threshold=1e-9,
                        singular_values=[1., 2e-9], solution=[1., 2.])
        self.assertTrue(checks.svd_parity(spectrum, spectrum)['passed'])
        self.assertFalse(checks.svd_parity(spectrum, {**spectrum, 'rank': 1})['passed'])
        self.assertFalse(checks.svd_parity(spectrum, {**spectrum, 'threshold': 1e-8})['passed'])

    def test_audit_comparison_requires_trust_and_complete_derivative(self):
        endpoint = self.endpoint()
        spectrum = dict(valid=True, rank=2, threshold=1e-9, singular_values=[1., .5], solution=[])
        derivative = dict(valid=True, **{name: [1., 2.] for name in
                           ('coefficients', 'correction', 'projected', 'jacobian', 'response')})
        report = dict(primary=endpoint, reference=endpoint, trust=dict(passed=True), initial_b=[.5],
                      derivative_audit=derivative, derivative_reason='full-profile-derivative',
                      svd_records=[spectrum])
        self.assertTrue(checks.audit_parity(report, report)['passed'])
        self.assertFalse(checks.audit_parity(report, {**report, 'trust': dict(passed=False)})['passed'])
        broken = {**report, 'derivative_audit': {**derivative, 'correction': []}}
        self.assertFalse(checks.audit_parity(report, broken)['passed'])

    def test_search_and_fixed_action_comparison_reject_changed_state_or_extra_work(self):
        search = self.search_result()
        self.assertTrue(checks.scientific_parity(search, search)['passed'])
        changed = {**search, 'returned_state': {**search['returned_state'], 'beta': [2., .21]}}
        self.assertFalse(checks.scientific_parity(search, changed)['passed'])
        fixed = self.fixed_action()
        self.assertTrue(checks.fixed_action_parity(fixed, fixed, ('schwarz',))['passed'])
        changed_fixed = {**fixed, 'work': {'reference_solves': 1, 'derivative_preparations': 0}}
        self.assertFalse(checks.fixed_action_parity(fixed, changed_fixed, ('schwarz',))['passed'])


if __name__ == '__main__':
    unittest.main()
