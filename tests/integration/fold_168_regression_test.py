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


def make_baseline() -> dict[str, object]:
    return {
        "schema_version": regression.SCHEMA_VERSION,
        "input_hashes": regression.EXPECTED_INPUT_HASHES,
        "command_arguments": regression.COMMAND_ARGUMENT_TEMPLATE,
        "expected_residue_count": 20,
        "maximum_residues_per_cluster": 10,
        "minimum_cluster_count": 2,
        "serial_ids": list(range(1, regression.EXPECTED_ATOM_COUNT + 1)),
        "reference_quality_metrics": {
            "amplitude_rmse": 1.0,
            "width_rmse": 1.0,
            "offset_rmse": 1.0,
            "maximum_absolute_offset": 1.0,
        },
    }


class Fold168RegressionTest(unittest.TestCase):
    def test_reads_atoms_and_rejects_duplicate_or_nonfinite_results(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            database = Path(temp_dir) / "results.sqlite"
            with sqlite3.connect(database) as connection:
                connection.execute(
                    "CREATE TABLE model_atom_local_potential ("
                    "key_tag TEXT, serial_id INTEGER, amplitude_estimate_mdpde_3rd REAL, "
                    "width_estimate_mdpde_3rd REAL, intercept_estimate_mdpde_3rd REAL, "
                    "alpha_r_3rd REAL)")
                connection.execute(
                    "CREATE TABLE model_atom ("
                    "key_tag TEXT, serial_id INTEGER, element INTEGER)")
                connection.executemany(
                    "INSERT INTO model_atom VALUES (?, ?, ?)",
                    [
                        (regression.SAVED_KEY, 1, 6),
                        (regression.SAVED_KEY, 2, 7),
                    ])
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

    def test_calculates_truth_quality_metrics(self) -> None:
        atoms = [
            {
                "serial_id": 1,
                "element": 6,
                "amplitude_mdpde": 6.0,
                "width_mdpde": 0.5,
                "intercept_mdpde": 1.0,
                "alpha_r": 0.2,
            },
            {
                "serial_id": 2,
                "element": 8,
                "amplitude_mdpde": 10.0,
                "width_mdpde": 0.7,
                "intercept_mdpde": -1.0,
                "alpha_r": 0.3,
            },
        ]
        metrics = regression.calculate_quality_metrics(atoms)
        self.assertAlmostEqual(metrics["amplitude_rmse"], 2.0 ** 0.5)
        self.assertAlmostEqual(metrics["width_rmse"], (0.04 / 2.0) ** 0.5)
        self.assertAlmostEqual(metrics["offset_rmse"], (4.0 / 2.0) ** 0.5)
        self.assertEqual(metrics["maximum_absolute_offset"], 1.0)

    def test_quality_gate_checks_validity_metrics_and_iteration_limit(self) -> None:
        baseline = make_baseline()
        atoms = [
            {
                "serial_id": serial_id,
                "element": 6,
                "amplitude_mdpde": 6.0,
                "width_mdpde": 0.5,
                "intercept_mdpde": 1.0,
                "alpha_r": 0.2,
            }
            for serial_id in range(1, regression.EXPECTED_ATOM_COUNT + 1)
        ]
        actual = {
            "atoms": atoms,
            "quality_metrics": {
                "amplitude_rmse": 1.05,
                "width_rmse": 1.05,
                "offset_rmse": 1.05,
                "maximum_absolute_offset": 1.05,
            },
            "second_stage_summary": {"accepted_iterations": 25},
            "residue_cutoff_summary": {
                "residue_count": 20,
                "limit": 10,
                "cluster_count": 2,
                "maximum_residue_count": 10,
                "cut_edge_count": 12,
            },
        }
        self.assertEqual(regression.validate_quality_gate(baseline, actual), [])
        actual["quality_metrics"]["width_rmse"] = 1.051
        actual["second_stage_summary"]["accepted_iterations"] = 26
        differences = "\n".join(regression.validate_quality_gate(baseline, actual))
        self.assertIn("quality_metrics.width_rmse", differences)
        self.assertIn("accepted_iterations", differences)

        actual["second_stage_summary"]["accepted_iterations"] = 25
        actual["residue_cutoff_summary"]["cluster_count"] = 1
        actual["residue_cutoff_summary"]["maximum_residue_count"] = 11
        cutoff_differences = "\n".join(
            regression.validate_quality_gate(baseline, actual))
        self.assertIn("residue_cutoff_summary.cluster_count", cutoff_differences)
        self.assertIn(
            "residue_cutoff_summary.maximum_residue_count",
            cutoff_differences)

    def test_parses_stable_second_stage_summary(self) -> None:
        summary = regression.parse_second_stage_summary(
            "Second-stage local fitting summary: accepted_iterations=4, "
            "best_iteration=2, stop_reason=audit-patience, "
            "best_audit_objective=1.25000000e-03, final_uses_polish=yes, "
            "final_state_source=best-audit.\n")
        self.assertEqual(summary["accepted_iterations"], 4)
        self.assertEqual(summary["best_iteration"], "2")
        self.assertEqual(summary["stop_reason"], "audit-patience")
        self.assertIs(summary["final_uses_polish"], True)
        self.assertEqual(summary["final_state_source"], "best-audit")

        no_polish = regression.parse_second_stage_summary(
            "Second-stage local fitting summary: accepted_iterations=0, "
            "best_iteration=initial, stop_reason=all-rejected-minimum-radius, "
            "best_audit_objective=2.50000000e-03, final_uses_polish=no, "
            "final_state_source=latest-validated.\n")
        self.assertIs(no_polish["final_uses_polish"], False)
        self.assertEqual(no_polish["final_state_source"], "latest-validated")

        no_retry_progress = regression.parse_second_stage_summary(
            "Second-stage local fitting summary: accepted_iterations=2, "
            "best_iteration=1, stop_reason=all-rejected-no-retry-progress, "
            "best_audit_objective=2.00000000e-03, final_uses_polish=no, "
            "final_state_source=best-audit.\n")
        self.assertEqual(
            no_retry_progress["stop_reason"],
            "all-rejected-no-retry-progress")

        unavailable = regression.parse_second_stage_summary(
            "Second-stage local fitting summary: accepted_iterations=0, "
            "best_iteration=unavailable, stop_reason=no-valid-seed, "
            "best_audit_objective=unavailable, final_uses_polish=unavailable, "
            "final_state_source=unavailable.\n")
        self.assertIsNone(unavailable["final_uses_polish"])
        self.assertEqual(unavailable["final_state_source"], "unavailable")

        with self.assertRaises(regression.RegressionError):
            regression.parse_second_stage_summary(
                "Second-stage local fitting summary: accepted_iterations=4, "
                "best_iteration=2, stop_reason=audit-patience, "
                "best_audit_objective=1.25000000e-03, final_uses_polish=yes.\n")

    def test_parses_stable_residue_cutoff_summary(self) -> None:
        summary = regression.parse_residue_cutoff_summary(
            "Local-fitting residue cutoff: residues=20, limit=10, clusters=2, "
            "max-residues=10, cutoff-edges=42.\n")
        self.assertEqual(summary, {
            "residue_count": 20,
            "limit": 10,
            "cluster_count": 2,
            "maximum_residue_count": 10,
            "cut_edge_count": 42,
        })
        with self.assertRaises(regression.RegressionError):
            regression.parse_residue_cutoff_summary("missing")

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
            baseline.write_text(json.dumps(make_baseline()), encoding="utf-8")
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
            baseline.write_text(json.dumps(make_baseline()), encoding="utf-8")
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

    def test_baseline_schema_version_five_is_required(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            baseline_path = Path(temp_dir) / "baseline.json"
            baseline = make_baseline()
            baseline_path.write_text(json.dumps(baseline), encoding="utf-8")
            self.assertEqual(regression.load_baseline(baseline_path), baseline)

            baseline["schema_version"] = 3
            baseline_path.write_text(json.dumps(baseline), encoding="utf-8")
            with self.assertRaises(regression.RegressionError):
                regression.load_baseline(baseline_path)


if __name__ == "__main__":
    unittest.main()
