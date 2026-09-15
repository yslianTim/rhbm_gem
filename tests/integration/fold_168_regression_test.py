#!/usr/bin/env python3

from __future__ import annotations

import copy
import json
import math
import sqlite3
import subprocess
import tempfile
import unittest
from contextlib import closing, redirect_stdout
from io import StringIO
from pathlib import Path
from unittest import mock

import fold_168_regression as regression


BASELINE_PATH = Path(__file__).parents[1] / "benchmarks" / "fold_168_simulation_baseline.json"
SUMMARY_LOG = (
    "Second-stage local fitting summary: accepted_iterations=25, "
    "best_iteration=20, stop_reason=audit-patience, best_audit_objective=0.1, "
    "final_uses_polish=no, final_state_source=best-audit.\n"
    "Local-fitting atom cutoff: atoms=168, limit=100, clusters=2, max-atoms=100, cutoff-edges=12.\n"
)


def make_manifest(count: int = 4) -> dict:
    contract = json.loads(BASELINE_PATH.read_text())["generation_contract"]
    atoms = [{
        "preparation_index": i, "serial_id": i + 1, "sequence_id": i // 2,
        "chain_id": "A", "component_id": "ALA", "atom_id": f"C{i}",
        "alternate_indicator": ".", "element": 6, "residue": 1, "spot": 6200,
        "structure": 1, "position": [i / 4.0, 0.0, 1.0], "table": "buried",
        "lookup_status": "found", "charge_used": (-0.5, 0.25, 0.0, 0.75)[i % 4],
        "lookup_charge": (-0.5, 0.25, 0.0, 0.75)[i % 4],
    } for i in range(count)]
    atoms[-1].update(lookup_status="unsupported_spot", charge_used=0.0, lookup_charge=None)
    return {
        "schema_version": 1,
        "generator": {"version": "2.0.0", "source_sha256": "0" * 64,
                      "configuration_sha256": "1" * 64, "build_sha256": "2" * 64},
        "source": {"model_path": "/original/model.cif", "pdb_id": "", "model_sha256": "a" * 64},
        "output": {"map_file": "original__map.map", "map_sha256": "b" * 64,
                   "format": "ccp4", "calculation_precision": "float64", "storage_precision": "float32"},
        **contract, "atom_count": count, "fallback_charge_count": 1,
        "execution": {"openmp_enabled": True, "requested_job_count": 4, "actual_job_count": 4,
                      "accumulation": "z_planes_in_preparation_order"},
        "atoms": atoms,
    }


def make_atoms(manifest: dict) -> list[dict]:
    return [{**copy.deepcopy(atom), "amplitude_mdpde": float(atom["element"]),
             "width_mdpde": manifest["settings"]["blurring_width"],
             "intercept_mdpde": atom["charge_used"], "alpha_r": 0.2}
            for atom in manifest["atoms"]]


def write_database(path: Path, atoms: list[dict]) -> None:
    with closing(sqlite3.connect(path)) as db, db:
        db.execute("CREATE TABLE model_atom (key_tag TEXT, serial_id INTEGER, chain_id TEXT, "
                   "sequence_id INTEGER, component_id TEXT, atom_id TEXT, indicator TEXT, "
                   "element INTEGER, structure INTEGER, position_x REAL, position_y REAL, position_z REAL)")
        db.execute("CREATE TABLE model_atom_local_potential (key_tag TEXT, serial_id INTEGER, "
                   "amplitude_estimate_mdpde_2nd REAL, width_estimate_mdpde_2nd REAL, "
                   "intercept_estimate_mdpde_2nd REAL, alpha_r_2nd REAL)")
        for atom in atoms:
            db.execute("INSERT INTO model_atom VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", (
                regression.SAVED_KEY, *(atom[name] for name in regression.IDENTITY_FIELDS),
                atom["element"], atom["structure"], *atom["position"]))
            db.execute("INSERT INTO model_atom_local_potential VALUES (?,?,?,?,?,?)", (
                regression.SAVED_KEY, atom["serial_id"], atom["amplitude_mdpde"], atom["width_mdpde"],
                atom["intercept_mdpde"], atom["alpha_r"]))


def make_fixture(root: Path) -> tuple[list[str], dict]:
    manifest = make_manifest(168)
    for name in ("model.cif", "map.map", "RHBM-GEM"):
        (root / name).write_text(name)
    manifest["source"]["model_sha256"] = regression.sha256_file(root / "model.cif")
    manifest["output"]["map_sha256"] = regression.sha256_file(root / "map.map")
    save_fixture(root, manifest)
    args = ["--executable", str(root / "RHBM-GEM"), "--model", str(root / "model.cif"),
            "--map", str(root / "map.map"), "--baseline", str(root / "baseline.json"),
            "--output-dir", str(root / "output")]
    return args, manifest


def save_fixture(root: Path, manifest: dict) -> None:
    regression.write_json(root / "map.map.simulation.json", manifest)
    baseline = json.loads(BASELINE_PATH.read_text())
    baseline["input_hashes"] = {name: regression.sha256_file(root / path) for name, path in (
        ("model", "model.cif"), ("map", "map.map"), ("manifest", "map.map.simulation.json"))}
    baseline["generation_contract"] = {key: manifest[key] for key in baseline["generation_contract"]}
    regression.write_json(root / "baseline.json", baseline)


class Fold168RegressionTest(unittest.TestCase):
    def load_manifest(self, root: Path, manifest: dict) -> dict:
        path = root / "truth.json"
        regression.write_json(path, manifest)
        return regression.load_simulation_manifest(path, {
            "model": manifest["source"]["model_sha256"], "map": manifest["output"]["map_sha256"]})

    def test_truth_scores_zero_and_preserves_zero_provenance(self) -> None:
        manifest = make_manifest()
        pairs = regression.pair_atom_truth(make_atoms(manifest), manifest)
        metrics = regression.calculate_quality_metrics(pairs)
        self.assertTrue(all(value == 0.0 for value in metrics.values()))
        self.assertNotIn("maximum_absolute_offset", metrics)
        self.assertEqual(pairs[2]["lookup_status"], "found")
        self.assertEqual(pairs[3]["lookup_status"], "unsupported_spot")
        self.assertEqual(pairs[2]["truth"]["offset"], pairs[3]["truth"]["offset"])

    def test_fixed_offset_shift_and_independent_parameter_errors(self) -> None:
        manifest = make_manifest()
        atoms = make_atoms(manifest)
        for atom in atoms:
            atom["intercept_mdpde"] -= 0.125
        atoms[0]["amplitude_mdpde"] += 2.0
        atoms[0]["width_mdpde"] += 0.2
        metrics = regression.calculate_quality_metrics(regression.pair_atom_truth(atoms, manifest))
        self.assertEqual(metrics["offset_rmse"], 0.125)
        self.assertEqual(metrics["offset_bias"], -0.125)
        self.assertEqual(metrics["offset_max_absolute_error"], 0.125)
        self.assertEqual(metrics["amplitude_rmse"], 1.0)
        self.assertAlmostEqual(metrics["width_rmse"], 0.1)

    def test_pairing_is_order_independent_and_does_not_modify_inputs(self) -> None:
        manifest = make_manifest()
        atoms = make_atoms(manifest)
        original = copy.deepcopy((atoms, manifest))
        expected = regression.pair_atom_truth(atoms, manifest)
        self.assertEqual((atoms, manifest), original)
        manifest["atoms"].reverse()
        self.assertEqual(regression.pair_atom_truth(atoms[::-1], manifest), expected)
        self.assertEqual([atom["preparation_index"] for atom in expected], list(range(4)))

    def test_rejects_missing_extra_duplicate_and_mismatched_atoms(self) -> None:
        manifest = make_manifest()
        original = make_atoms(manifest)
        cases = [original[:-1], original + [original[0]]]
        for field, value in (("chain_id", "B"), ("sequence_id", 999), ("component_id", "VAL"),
                             ("atom_id", "wrong"), ("alternate_indicator", "A"), ("serial_id", 99),
                             ("element", 8), ("structure", 2), ("position", [0.0, 0.0, 2.0])):
            changed = copy.deepcopy(original)
            changed[0][field] = value
            cases.append(changed)
        for atoms in cases:
            with self.subTest(atoms=atoms):
                with self.assertRaises(regression.RegressionError):
                    regression.pair_atom_truth(atoms, manifest)

    def test_same_element_uses_individual_charge_and_full_width(self) -> None:
        manifest = make_manifest()
        manifest["settings"]["blurring_width"] = 0.5123456789012345
        manifest["settings"]["blurring_width_list"] = [0.5123456789012345]
        manifest["kernel"]["effective_charge_width"] = 0.5123456789012345
        manifest["atoms"][0]["chain_id"] = 'A"\\line\n'
        with tempfile.TemporaryDirectory() as temp_dir:
            loaded = self.load_manifest(Path(temp_dir), manifest)
        atoms = make_atoms(loaded)
        self.assertEqual(regression.pair_atom_truth(atoms, loaded)[0]["truth"]["width"], 0.5123456789012345)
        atoms[0]["intercept_mdpde"], atoms[1]["intercept_mdpde"] = atoms[1]["intercept_mdpde"], atoms[0]["intercept_mdpde"]
        metrics = regression.calculate_quality_metrics(regression.pair_atom_truth(atoms, loaded))
        self.assertAlmostEqual(metrics["offset_rmse"], math.sqrt(2 * 0.75 ** 2 / 4))

    def test_neutral_amber_and_all_fallback_reasons(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            for status in sorted(regression.CHARGE_FAILURE_STATUSES):
                manifest = make_manifest()
                manifest["atoms"][-1]["lookup_status"] = status
                manifest["atoms"][-1]["table"] = None
                self.load_manifest(root, manifest)
            for mode, code in (("neutral", 0), ("amber", 2)):
                manifest = make_manifest()
                manifest["settings"].update(charge_mode=mode, charge_mode_code=code)
                if mode == "neutral":
                    manifest["fallback_charge_count"] = 0
                    for atom in manifest["atoms"]:
                        atom.update(lookup_status="neutral_mode", table=None, lookup_charge=None, charge_used=0.0)
                else:
                    for atom in manifest["atoms"]:
                        atom["table"] = "amber95"
                loaded = self.load_manifest(root, manifest)
                self.assertEqual(regression.calculate_quality_metrics(
                    regression.pair_atom_truth(make_atoms(loaded), loaded))["offset_rmse"], 0.0)

    def test_rejects_bad_manifest_schema_evidence_and_nonfinite_values(self) -> None:
        mutations = [
            lambda d: d.update(schema_version=2),
            lambda d: d.update(atom_count=5),
            lambda d: d.update(fallback_charge_count=0),
            lambda d: d["atoms"][0].update(serial_id=d["atoms"][1]["serial_id"]),
            lambda d: d["atoms"][0].update(preparation_index=9),
            lambda d: d["atoms"][0].update(preparation_index=True),
            lambda d: d["atoms"][0].update(charge_used=math.nan),
            lambda d: d["atoms"][0].update(charge_used=0.1),
            lambda d: d["atoms"][-1].update(charge_used=0.1),
            lambda d: d["atoms"][0].update(table=None),
            lambda d: d["atoms"][0].pop("lookup_charge"),
            lambda d: d["atoms"][0].update(position=[0.0, math.inf, 0.0]),
            lambda d: d["settings"].update(potential_model="five_gaus_charge", potential_model_code=1),
            lambda d: d["settings"].update(normalization_applied=True),
            lambda d: d["settings"].update(occupancy_applied=True),
            lambda d: d["settings"].update(temperature_factor_applied=True),
            lambda d: d["settings"].update(blurring_width=0.0),
            lambda d: d["settings"].update(blurring_width_list=[0.6]),
            lambda d: d["settings"].update(grid_size=[1, 2, True]),
            lambda d: d["kernel"].pop("minimum_charge_width"),
            lambda d: d["kernel"].update(effective_charge_width=0.6),
        ]
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "truth.json"
            for i, mutate in enumerate(mutations):
                manifest = make_manifest()
                mutate(manifest)
                path.write_text(json.dumps(manifest))
                with self.subTest(case=i), self.assertRaises(regression.RegressionError):
                    regression.load_simulation_manifest(path, {
                        "model": manifest["source"]["model_sha256"], "map": manifest["output"]["map_sha256"]})
            for text in ('{"x":1,"x":2}', '{"x":1e999}', '{"x":NaN}'):
                path.write_text(text)
                with self.assertRaises(regression.RegressionError):
                    regression.read_json(path)

    def test_reads_identity_and_physical_offset_without_writing_database(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "database.sqlite"
            manifest = make_manifest()
            write_database(path, make_atoms(manifest))
            before = regression.sha256_file(path)
            atoms = regression.read_atom_results(path)
            pairs = regression.pair_atom_truth(atoms, manifest)
            self.assertEqual(regression.calculate_quality_metrics(pairs)["offset_rmse"], 0.0)
            self.assertEqual(before, regression.sha256_file(path))
            self.assertEqual(atoms[0]["intercept_mdpde"], -0.5)

    def test_rejects_invalid_sqlite_results_and_orphan_rows(self) -> None:
        changes = [
            "UPDATE model_atom_local_potential SET amplitude_estimate_mdpde_2nd=0 WHERE serial_id=1",
            "UPDATE model_atom_local_potential SET width_estimate_mdpde_2nd=-1 WHERE serial_id=1",
            "UPDATE model_atom_local_potential SET intercept_estimate_mdpde_2nd=NULL WHERE serial_id=1",
            "UPDATE model_atom_local_potential SET intercept_estimate_mdpde_2nd=1e999 WHERE serial_id=1",
            "UPDATE model_atom_local_potential SET alpha_r_2nd=2 WHERE serial_id=1",
            "DELETE FROM model_atom WHERE serial_id=1",
            "INSERT INTO model_atom_local_potential SELECT * FROM model_atom_local_potential WHERE serial_id=1",
        ]
        for sql in changes:
            with self.subTest(sql=sql), tempfile.TemporaryDirectory() as temp_dir:
                path = Path(temp_dir) / "database.sqlite"
                write_database(path, make_atoms(make_manifest()))
                with closing(sqlite3.connect(path)) as db, db:
                    db.execute(sql)
                with self.assertRaises(regression.RegressionError):
                    regression.read_atom_results(path)

    def test_large_finite_errors_do_not_overflow_squares(self) -> None:
        metrics = regression.calculate_quality_metrics([
            {"error": {"amplitude": 1e200, "width": -1e200, "offset": 1e200}},
            {"error": {"amplitude": 1e200, "width": -1e200, "offset": -1e200}},
        ])
        self.assertEqual(metrics["offset_rmse"], 1e200)
        self.assertEqual(metrics["offset_bias"], 0.0)

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
            "final_state_source=best-audit.\n")
        self.assertIs(no_polish["final_uses_polish"], False)
        self.assertEqual(no_polish["final_state_source"], "best-audit")

        converged = regression.parse_second_stage_summary(
            "Second-stage local fitting summary: accepted_iterations=3, "
            "best_iteration=2, stop_reason=converged, "
            "best_audit_objective=1.75000000e-03, final_uses_polish=no, "
            "final_state_source=latest-validated.\n")
        self.assertIs(converged["final_uses_polish"], False)
        self.assertEqual(converged["final_state_source"], "latest-validated")

        multiline = regression.parse_second_stage_summary(
            " Second-Stage Local Fitting Summary : \n"
            " - accepted_iterations = 7\n"
            " - best_iteration = 4\n"
            " - stop_reason = audit-patience\n"
            " - best_audit_objective = 5.19e-01\n"
            " - final_uses_polish = yes\n"
            " - final_state_source = best-audit\n")
        self.assertEqual(multiline["accepted_iterations"], 7)
        self.assertEqual(multiline["best_iteration"], "4")
        self.assertEqual(multiline["stop_reason"], "audit-patience")
        self.assertIs(multiline["final_uses_polish"], True)

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

    def test_parses_stable_atom_cutoff_summary(self) -> None:
        log = (
            "Local-fitting atom cutoff: atoms=168, limit=100, clusters=2, "
            "max-atoms=100, cutoff-edges=42.\n")
        summary = regression.parse_atom_cutoff_summary(log)
        self.assertEqual(summary, {
            "atom_count": 168,
            "limit": 100,
            "cluster_count": 2,
            "maximum_atom_count": 100,
            "cut_edge_count": 42,
        })
        for invalid_log in (
            "missing", log + log,
            "Local-fitting residue cutoff: residues=20, limit=10, clusters=2, "
            "max-residues=10, cutoff-edges=42.\n",
        ):
            with self.subTest(log=invalid_log):
                with self.assertRaises(regression.RegressionError):
                    regression.parse_atom_cutoff_summary(invalid_log)

    def test_final_certificate_requires_persisted_state_and_finite_evidence(self) -> None:
        record = {"final_polish_applied": False, "attempts": 12,
                  "recovery_operator_evaluations": 1, "certificate_operator_evaluations": 1,
                  "background_reference": "last_frozen_background",
                  "certificate": {"reference": "persisted_state", "status": "evaluated",
                                  "qualified": False, "complete": True,
                                  "operator_nominal_p99": [0.0, 0.0, 0.0]}}
        prefix = "Second-stage final state: schema=1, payload="
        line = prefix + json.dumps(record)
        self.assertEqual(regression.parse_final_state_certificate(line), record)
        self.assertIsNone(regression.parse_final_state_certificate(SUMMARY_LOG))
        for field, value in (("reference", "iteration_previous"), ("qualified", 1),
                             ("operator_nominal_p99", [float("nan"), 0.0, 0.0]),
                             ("operator_nominal_p99", [-1.0, 0.0, 0.0])):
            changed = copy.deepcopy(record)
            changed["certificate"][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(regression.RegressionError):
                regression.parse_final_state_certificate(prefix + json.dumps(changed))
        with self.assertRaises(regression.RegressionError):
            regression.parse_final_state_certificate(line + "\n" + line)

    def test_small_residual_and_zero_truth_error_cannot_bypass_solver_or_budget(self) -> None:
        for qualified, complete, attempts, reason, expected in (
                (False, True, 12, "converged", False),
                (True, False, 12, "converged", False),
                (True, True, 12, "recovery-failed", False),
                (True, True, 26, "converged", False),
                (True, True, 12, "converged", True)):
            with self.subTest(qualified=qualified, complete=complete, attempts=attempts, reason=reason), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                args, manifest = make_fixture(root)
                record = {"final_polish_applied": False, "attempts": attempts,
                          "recovery_operator_evaluations": 1, "certificate_operator_evaluations": 1,
                          "background_reference": "last_frozen_background",
                          "certificate": {"reference": "persisted_state", "status": "evaluated",
                                          "qualified": qualified, "complete": complete,
                                          "operator_nominal_p99": [0.0, 0.0, 0.0]}}
                log = SUMMARY_LOG.replace("audit-patience", reason)
                log += "Second-stage final state: schema=1, payload=" + json.dumps(record)
                def execute(command, **kwargs):
                    write_database(Path(command[command.index("-d") + 1]), make_atoms(manifest))
                    return subprocess.CompletedProcess(command, 0, log.encode())
                with mock.patch.object(regression.subprocess, "run", side_effect=execute), redirect_stdout(StringIO()):
                    self.assertEqual(regression.run(args), 1)
                report = json.loads((root / "output/report.json").read_text())
                self.assertEqual(report["errors"], [])
                self.assertEqual(report["convergence_acceptance"]["passed"], expected)
                self.assertEqual(report["quality_gate"]["status"], "uncalibrated")

    def test_separate_uncalibrated_quality_and_existing_gates(self) -> None:
        baseline = regression.load_baseline(BASELINE_PATH)
        actual = regression.make_empty_actual({})
        actual.update(atoms=make_atoms(make_manifest(168)),
                      second_stage_summary=regression.parse_second_stage_summary(SUMMARY_LOG),
                      atom_cutoff_summary=regression.parse_atom_cutoff_summary(SUMMARY_LOG))
        gates = regression.evaluate_gates(baseline, actual)
        self.assertEqual(gates["quality_gate"]["status"], "uncalibrated")
        self.assertFalse(gates["quality_gate"]["passed"])
        self.assertTrue(gates["iteration_gate"]["passed"])
        self.assertTrue(gates["atom_cutoff_gate"]["passed"])
        actual["second_stage_summary"]["accepted_iterations"] = 26
        actual["atom_cutoff_summary"].update(cluster_count=1, maximum_atom_count=101)
        gates = regression.evaluate_gates(baseline, actual)
        self.assertFalse(gates["iteration_gate"]["passed"])
        self.assertFalse(gates["atom_cutoff_gate"]["passed"])
        self.assertIn("maximum_atom_count", " ".join(gates["atom_cutoff_gate"]["differences"]))

    def test_input_failures_never_execute_or_publish_scores(self) -> None:
        for case in ("model", "map", "manifest", "missing_manifest", "unsupported_model", "selection", "manifest_hash"):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as temp_dir:
                root = Path(temp_dir)
                args, manifest = make_fixture(root)
                if case == "missing_manifest":
                    (root / "map.map.simulation.json").unlink()
                elif case == "unsupported_model":
                    manifest["settings"].update(potential_model="five_gaus_charge", potential_model_code=1)
                    save_fixture(root, manifest)
                elif case == "selection":
                    manifest["settings"]["only_backbone"] = True
                    regression.write_json(root / "map.map.simulation.json", manifest)
                    baseline = json.loads((root / "baseline.json").read_text())
                    baseline["input_hashes"]["manifest"] = regression.sha256_file(root / "map.map.simulation.json")
                    regression.write_json(root / "baseline.json", baseline)
                elif case == "manifest_hash":
                    manifest["output"]["map_sha256"] = "0" * 64
                    save_fixture(root, manifest)
                else:
                    path = root / {"model": "model.cif", "map": "map.map", "manifest": "map.map.simulation.json"}[case]
                    path.write_bytes(b"changed")
                with mock.patch.object(regression.subprocess, "run") as run_mock, redirect_stdout(StringIO()):
                    self.assertEqual(regression.run(args), 1)
                    run_mock.assert_not_called()
                report = json.loads((root / "output/report.json").read_text())
                actual = json.loads((root / "output/actual.json").read_text())
                self.assertEqual(report["truth_scoring"]["status"], "failed")
                self.assertFalse(report["passed"])
                self.assertTrue(report["errors"])
                self.assertIsNone(actual["quality_metrics"])
                self.assertEqual(actual["atoms"], [])
                self.assertTrue((root / "output/run.log").exists())

    def test_successful_scoring_remains_uncalibrated_and_temporary_database_is_removed(self) -> None:
        for explicit_manifest in (False, True):
            with self.subTest(explicit=explicit_manifest), tempfile.TemporaryDirectory() as temp_dir:
                root = Path(temp_dir)
                args, manifest = make_fixture(root)
                if explicit_manifest:
                    (root / "map.map.simulation.json").rename(root / "relocated.json")
                    args.extend(["--simulation-manifest", str(root / "relocated.json")])
                captured = []

                def execute(command, **kwargs):
                    database = Path(command[command.index("-d") + 1])
                    self.assertFalse(database.exists())
                    self.assertEqual(Path(kwargs["cwd"]), database.parent)
                    self.assertNotIn("--simulation-manifest", command)
                    self.assertEqual(command[-6:], ["--simulation", "true", "-r", "0.50", "--exclude-hydrogen", "true"])
                    captured.append(database)
                    write_database(database, make_atoms(manifest))
                    return subprocess.CompletedProcess(command, 0, SUMMARY_LOG.encode())

                with mock.patch.object(regression.subprocess, "run", side_effect=execute), redirect_stdout(StringIO()):
                    self.assertEqual(regression.run(args), 1)
                self.assertEqual(len(captured), 1)
                self.assertFalse(captured[0].parent.exists())
                report = json.loads((root / "output/report.json").read_text())
                actual = json.loads((root / "output/actual.json").read_text())
                self.assertEqual(report["truth_scoring"]["status"], "complete")
                self.assertEqual(report["errors"], [])
                self.assertEqual(report["quality_gate"]["status"], "uncalibrated")
                self.assertFalse(report["passed"])
                self.assertEqual(report["stop_reason"], "audit-patience")
                self.assertFalse(report["convergence_acceptance"]["passed"])
                self.assertTrue(report["iteration_gate"]["passed"])
                self.assertTrue(report["atom_cutoff_gate"]["passed"])
                self.assertEqual(actual["quality_metrics"]["offset_rmse"], 0.0)
                self.assertEqual(actual["diagnostics"]["maximum_absolute_offset"], 0.75)
                self.assertEqual(actual["atoms"][-1]["lookup_status"], "unsupported_spot")
                self.assertEqual(actual["generation_record"]["output"]["map_file"], "original__map.map")

    def test_inputs_changed_during_execution_are_not_scored(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            args, manifest = make_fixture(root)

            def execute(command, **kwargs):
                write_database(Path(command[command.index("-d") + 1]), make_atoms(manifest))
                (root / "map.map").write_text("changed during fit")
                return subprocess.CompletedProcess(command, 0, SUMMARY_LOG.encode())

            with mock.patch.object(regression.subprocess, "run", side_effect=execute), redirect_stdout(StringIO()):
                self.assertEqual(regression.run(args), 1)
            report = json.loads((root / "output/report.json").read_text())
            self.assertIn("SHA-256 mismatch", " ".join(report["errors"]))
            self.assertIsNone(json.loads((root / "output/actual.json").read_text())["quality_metrics"])

    def test_baseline_requires_schema_seven_and_no_old_thresholds(self) -> None:
        baseline = regression.load_baseline(BASELINE_PATH)
        self.assertEqual(baseline["schema_version"], 7)
        self.assertNotIn("--fit-max", regression.COMMAND_ARGUMENT_TEMPLATE)
        self.assertNotIn("--fit-min", regression.COMMAND_ARGUMENT_TEMPLATE)
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "baseline.json"
            for version in (5, 6):
                regression.write_json(path, {**baseline, "schema_version": version})
                with self.assertRaises(regression.RegressionError):
                    regression.load_baseline(path)
            regression.write_json(path, {**baseline, "reference_quality_metrics": {"offset_rmse": 1.0}})
            with self.assertRaises(regression.RegressionError):
                regression.load_baseline(path)


if __name__ == "__main__":
    unittest.main()
