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


if __name__ == '__main__':
    unittest.main()
