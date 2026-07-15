#!/usr/bin/env python3

from __future__ import annotations

import json
import sqlite3
import subprocess
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest import mock

import fold_168_regression as regression


class Fold168RegressionTest(unittest.TestCase):
    def test_reads_atoms_and_rejects_duplicate_or_nonfinite_results(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            database = Path(temp_dir) / "results.sqlite"
            with sqlite3.connect(database) as connection:
                connection.execute(
                    "CREATE TABLE model_atom_local_potential ("
                    "key_tag TEXT, serial_id INTEGER, amplitude_estimate_mdpde REAL, "
                    "width_estimate_mdpde REAL, intercept_estimate_mdpde REAL, alpha_r REAL)")
                connection.executemany(
                    "INSERT INTO model_atom_local_potential VALUES (?, ?, ?, ?, ?, ?)",
                    [
                        (regression.SAVED_KEY, 1, 6.0, 0.5, 0.1, 0.2),
                        (regression.SAVED_KEY, 2, 7.0, 0.6, -0.1, 0.3),
                    ])
            self.assertEqual(len(regression.read_atom_results(database)), 2)

            with sqlite3.connect(database) as connection:
                connection.execute(
                    "INSERT INTO model_atom_local_potential VALUES (?, ?, ?, ?, ?, ?)",
                    (regression.SAVED_KEY, 2, 7.0, 0.6, -0.1, 0.3))
            with self.assertRaises(regression.RegressionError):
                regression.read_atom_results(database)

    def test_applies_atom_tolerance_and_rejects_nonfinite_values(self) -> None:
        expected = {"atoms": [{"serial_id": 1, "width_mdpde": 1.0}]}
        within_tolerance = {"atoms": [{"serial_id": 1, "width_mdpde": 1.0 + 5.0e-7}]}
        outside_tolerance = {"atoms": [{"serial_id": 1, "width_mdpde": 1.0 + 2.0e-6}]}
        nonfinite = {"atoms": [{"serial_id": 1, "width_mdpde": float("nan")}]}
        self.assertEqual(regression.compare_values(expected, within_tolerance), [])
        self.assertTrue(regression.compare_values(expected, outside_tolerance))
        self.assertTrue(regression.compare_values(expected, nonfinite))

    def test_reports_missing_atom(self) -> None:
        expected = {"atoms": [{"serial_id": 1}]}
        actual = {"atoms": []}
        differences = "\n".join(regression.compare_values(expected, actual))
        self.assertIn("root.atoms: expected 1 entries, got 0", differences)

    def test_hash_failure_does_not_execute_and_preserves_reports(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            executable = root / "RHBM-GEM"
            executable.write_text("not executed", encoding="utf-8")
            model = root / "model.cif"
            map_path = root / "map.map"
            for path in (model, map_path):
                path.write_bytes(b"wrong fixture")
            baseline = root / "baseline.json"
            baseline.write_text(json.dumps({
                "schema_version": regression.SCHEMA_VERSION,
                "input_hashes": regression.EXPECTED_INPUT_HASHES,
                "command_arguments": regression.COMMAND_ARGUMENT_TEMPLATE,
            }), encoding="utf-8")
            output_dir = root / "output"
            arguments = [
                "--executable", str(executable),
                "--model", str(model),
                "--map", str(map_path),
                "--baseline", str(baseline),
                "--output-dir", str(output_dir),
            ]

            with mock.patch.object(regression.subprocess, "run") as run_mock:
                with redirect_stdout(StringIO()):
                    self.assertEqual(regression.run(arguments), 1)
                run_mock.assert_not_called()
            self.assertTrue((output_dir / "run.log").exists())
            self.assertTrue((output_dir / "actual.json").exists())
            report = json.loads((output_dir / "report.json").read_text(encoding="utf-8"))
            self.assertFalse(report["passed"])
            self.assertIn("SHA-256 mismatch", "\n".join(report["errors"]))

    def test_uses_nonexistent_temporary_output_database(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            executable = root / "RHBM-GEM"
            executable.write_text("mock executable", encoding="utf-8")
            model = root / "model.cif"
            map_path = root / "map.map"
            model.write_bytes(b"model")
            map_path.write_bytes(b"map")
            baseline = root / "baseline.json"
            baseline.write_text(json.dumps({
                "schema_version": regression.SCHEMA_VERSION,
                "input_hashes": regression.EXPECTED_INPUT_HASHES,
                "command_arguments": regression.COMMAND_ARGUMENT_TEMPLATE,
            }), encoding="utf-8")
            output_dir = root / "output"
            arguments = [
                "--executable", str(executable),
                "--model", str(model),
                "--map", str(map_path),
                "--baseline", str(baseline),
                "--output-dir", str(output_dir),
            ]

            def reject_after_inspection(
                command: list[str],
                **_: object,
            ) -> subprocess.CompletedProcess[bytes]:
                database_path = Path(command[command.index("-d") + 1])
                self.assertEqual(database_path.name, "database.sqlite")
                self.assertFalse(database_path.exists())
                return subprocess.CompletedProcess(command, 1, stdout=b"mock failure")

            with mock.patch.object(
                regression,
                "validate_input_hashes",
                return_value=regression.EXPECTED_INPUT_HASHES,
            ), mock.patch.object(
                regression.subprocess,
                "run",
                side_effect=reject_after_inspection,
            ), redirect_stdout(StringIO()):
                self.assertEqual(regression.run(arguments), 1)

    def test_baseline_schema_version_three_is_required(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            baseline_path = Path(temp_dir) / "baseline.json"
            baseline = {
                "schema_version": regression.SCHEMA_VERSION,
                "input_hashes": regression.EXPECTED_INPUT_HASHES,
                "command_arguments": regression.COMMAND_ARGUMENT_TEMPLATE,
            }
            baseline_path.write_text(json.dumps(baseline), encoding="utf-8")
            self.assertEqual(regression.load_baseline(baseline_path), baseline)

            baseline["schema_version"] = 2
            baseline_path.write_text(json.dumps(baseline), encoding="utf-8")
            with self.assertRaises(regression.RegressionError):
                regression.load_baseline(baseline_path)


if __name__ == "__main__":
    unittest.main()
