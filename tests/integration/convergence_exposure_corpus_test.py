#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RUNNER = load_module(
    "exposure_corpus_runner",
    PROJECT_ROOT / "resources" / "tools" / "developer" /
    "run_convergence_exposure_corpus.py")
ANALYZER = load_module(
    "exposure_corpus_analyzer",
    PROJECT_ROOT / "resources" / "tools" / "developer" /
    "analyze_convergence_exposure_corpus.py")
COUNTERFACTUAL = load_module(
    "counterfactual_convergence_analyzer_for_exposure_test",
    PROJECT_ROOT / "resources" / "tools" / "developer" /
    "analyze_counterfactual_convergence.py")


def policy(
    objective: float,
    truth_error: float,
    extra_attempts: int = 0,
) -> dict[str, object]:
    return {
        "reached": True,
        "extra_attempts": extra_attempts,
        "objective": objective,
        "material_objective_improvement": objective < 0.999,
        "raw_median": [1.0e-5, 1.0e-5, 1.0e-5],
        "truth_metrics": {"transformed_aggregate_rmse": truth_error},
    }


def case_summary(
    case_id: str,
    legacy_population: bool,
    maximum_gate: bool,
    solver_qualification: bool,
    candidate_objective: float = 0.998,
    candidate_truth: float = 0.08,
) -> dict[str, object]:
    policies = {
        "production": policy(1.0, 0.1),
        "legacy-population": policy(candidate_objective, candidate_truth, 2),
        "legacy-maximum": policy(candidate_objective, candidate_truth, 2),
        "solver-qualified": policy(candidate_objective, candidate_truth, 2),
    }
    return {
        "status": "complete",
        "case": {
            "case_id": case_id, "family": "fixture",
            "topology": "fixture", "level": 0,
        },
        "safety_regression": False,
        "audit": {
            "status": "counterfactual_records",
            "experiments": [{
                "trigger_try": 4,
                "termination_reason": "all-policies-reached",
                "exposures": {
                    "legacy_population": legacy_population,
                    "maximum_gate": maximum_gate,
                    "solver_qualification": solver_qualification,
                },
                "policies": policies,
            }],
            "trajectory_audit": {
                "p99_implies_max": {
                    "production": {"accepted:N<=91:samples": 2},
                },
                "unique_blocker_count": {
                    "legacy_maximum": {"accepted_max": 1, "raw_max": 0},
                },
            },
        },
    }


class ConvergenceExposureCorpusTest(unittest.TestCase):
    def test_manifest_expands_to_reproducible_six_hundred_cases(self) -> None:
        manifest = json.loads((
            PROJECT_ROOT / "tests" / "benchmarks" /
            "convergence_exposure_manifest.json").read_text(encoding="utf-8"))
        cases = RUNNER.expand_manifest(manifest)

        self.assertEqual(len(cases), 600)
        self.assertEqual(len({case["case_id"] for case in cases}), 600)
        self.assertEqual(cases[0]["seed"], 410000)
        self.assertEqual(cases[200]["seed"], 420000)
        self.assertEqual(cases[400]["seed"], 430000)
        self.assertEqual(RUNNER.expand_manifest(manifest), cases)

    def test_manifest_rejects_incomplete_axis(self) -> None:
        manifest = json.loads((
            PROJECT_ROOT / "tests" / "benchmarks" /
            "convergence_exposure_manifest.json").read_text(encoding="utf-8"))
        manifest["families"][0]["levels"].pop()
        with self.assertRaisesRegex(ValueError, "five topologies and levels"):
            RUNNER.expand_manifest(manifest)

    def test_truth_parser_and_transformed_truth_metrics(self) -> None:
        text = "\n".join((
            "Convergence exposure truth: schema=1, case=x, serial=1, "
            "amplitude=6, width=0.5, offset=0.1.",
            "[Debug] Counterfactual convergence atom: schema=3, "
            "experiment=1-1, policy=production, serial=1, "
            "amplitude=6, width=0.5, offset=0.1",
        ))
        truth = RUNNER.parse_truth(text)
        atom_rows = COUNTERFACTUAL.parse_log(text)["atoms"][("1-1", "production")]
        metrics = COUNTERFACTUAL._truth_metrics(atom_rows, truth)

        self.assertEqual(len(truth), 1)
        self.assertEqual(metrics["transformed_aggregate_rmse"], 0.0)

    def test_analyzer_classifies_exposures_outcomes_and_replay_order(self) -> None:
        summaries = []
        for index in range(5):
            summaries.append(case_summary(f"p-{index}", True, False, False))
            summaries.append(case_summary(f"m-{index}", False, True, False))
            summaries.append(case_summary(f"s-{index}", False, False, True))
        report = ANALYZER.analyze(summaries)

        self.assertEqual(report["genuine_exposure_count"], 15)
        self.assertEqual(
            report["comparator_exposure_counts"]["maximum_gate"], 5)
        self.assertEqual(len(report["replay_case_ids"]), 15)
        self.assertFalse(report["corpus_target_met"])
        self.assertEqual(
            report["policy_decisions"]["solver-qualified"]["benefit_ratio"], 1.0)
        self.assertEqual(
            report["maximum_evidence"]["production"]
            ["accepted:N<=91:samples"], 30)
        self.assertEqual(
            report["maximum_evidence"]["production"]
            ["accepted_max_unique_catches"], 15)

    def test_analyzer_reports_harm_and_safety_regression(self) -> None:
        harmful = case_summary("harm", False, False, True, 1.01, 0.12)
        harmful["safety_regression"] = True
        report = ANALYZER.analyze([harmful])
        outcome = report["cases"][0]["outcomes"]["solver-qualified"]

        self.assertEqual(outcome["category"], "material-harm")
        self.assertTrue(outcome["safety_regression"])
        self.assertFalse(
            report["policy_decisions"]["solver-qualified"][
                "redesign_candidate"])

    def test_analyzer_keeps_agreement_and_no_trigger_as_controls(self) -> None:
        agreement = case_summary("agreement", False, False, False)
        no_trigger = {
            "status": "complete",
            "case": {
                "case_id": "no-trigger", "family": "fixture",
                "topology": "fixture", "level": 0,
            },
            "audit": {
                "status": "no_convergence_trigger", "experiments": [],
                "trajectory_audit": {},
            },
        }
        report = ANALYZER.analyze([agreement, no_trigger])

        self.assertEqual(report["exposure_counts"]["policy-agreement"], 1)
        self.assertEqual(report["exposure_counts"]["no-trigger"], 1)
        self.assertEqual(report["genuine_exposure_count"], 0)

    def test_analyzer_separates_budget_and_existing_safeguard_termination(self) -> None:
        budget = case_summary("budget", False, False, True)
        budget_experiment = budget["audit"]["experiments"][0]
        budget_experiment["termination_reason"] = "budget-exhausted"
        budget_experiment["policies"]["solver-qualified"] = {"reached": False}
        safeguarded = case_summary("safeguarded", False, False, True)
        safeguard_experiment = safeguarded["audit"]["experiments"][0]
        safeguard_experiment["termination_reason"] = "audit-patience"
        safeguard_experiment["policies"]["solver-qualified"] = {"reached": False}

        self.assertEqual(
            ANALYZER.classify_policy_outcome(
                budget, "solver-qualified")["category"],
            "unresolved-budget-exhausted")
        self.assertEqual(
            ANALYZER.classify_policy_outcome(
                safeguarded, "solver-qualified")["category"],
            "terminated-existing-safeguard")

    def test_runner_resumes_completed_case(self) -> None:
        case = {
            "case_id": "natural-v00-r0", "family": "natural",
            "topology": "unk-c", "level": 0, "replica": 0,
            "seed": 410000, "noise_sigma": 0.0,
        }
        with tempfile.TemporaryDirectory() as temp_directory:
            case_directory = (
                Path(temp_directory) / "cases" / case["case_id"])
            case_directory.mkdir(parents=True)
            legacy = {
                "schema_version": 5,
                "status": "complete", "case": case, "thread_count": 1,
            }
            RUNNER.write_json(case_directory / "case-summary.json", legacy)
            with mock.patch.object(
                    RUNNER.subprocess, "run",
                    return_value=subprocess.CompletedProcess(
                        [], 0, stdout=b"")) as initial_run_mock:
                expected = RUNNER.run_case(
                    Path("unused"), case, Path(temp_directory), 1)
            with mock.patch.object(RUNNER.subprocess, "run") as resume_run_mock:
                actual = RUNNER.run_case(
                    Path("unused"), case, Path(temp_directory), 1)

        self.assertEqual(actual, expected)
        self.assertEqual(
            actual["schema_version"], RUNNER.CASE_SUMMARY_SCHEMA_VERSION)
        initial_run_mock.assert_called_once()
        resume_run_mock.assert_not_called()

    def test_runner_writes_complete_case_artifacts(self) -> None:
        case = {
            "case_id": "natural-v00-r0", "family": "natural",
            "topology": "unk-c", "level": 0, "replica": 0,
            "seed": 410000, "noise_sigma": 0.0,
        }
        log = (
            "Convergence exposure truth: schema=1, case=natural-v00-r0, "
            "serial=1, amplitude=6, width=0.5, offset=0.1.\n")
        with tempfile.TemporaryDirectory() as temp_directory:
            output_directory = Path(temp_directory)
            with mock.patch.object(
                    RUNNER.subprocess, "run",
                    return_value=subprocess.CompletedProcess(
                        [], 0, stdout=log.encode("utf-8"))):
                summary = RUNNER.run_case(
                    Path("unused"), case, output_directory, 1)
            case_directory = output_directory / "cases" / case["case_id"]
            artifact_names = {
                path.name for path in case_directory.iterdir()
            }

        self.assertEqual(summary["status"], "complete")
        self.assertTrue({
            "run.log", "scenario-truth.json", "trajectory-schema-5.json",
            "counterfactual-schema-3.json", "case-summary.json",
        }.issubset(artifact_names))

    def test_runner_isolates_failed_case(self) -> None:
        case = {
            "case_id": "failed", "family": "natural",
            "topology": "unk-c", "level": 0, "replica": 0,
            "seed": 410000, "noise_sigma": 0.0,
        }
        with tempfile.TemporaryDirectory() as temp_directory:
            output_directory = Path(temp_directory)
            with mock.patch.object(
                    RUNNER.subprocess, "run",
                    return_value=subprocess.CompletedProcess(
                        [], 7, stdout=b"isolated failure\n")):
                summary = RUNNER.run_case(
                    Path("unused"), case, output_directory, 1)
            log = (output_directory / "cases" / "failed" / "run.log").read_text()

        self.assertEqual(summary["status"], "failed")
        self.assertEqual(summary["return_code"], 7)
        self.assertEqual(log, "isolated failure\n")


if __name__ == "__main__":
    unittest.main()
