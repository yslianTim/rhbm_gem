#!/usr/bin/env python3
"""Contract tests for production-only convergence corpus tooling."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest


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


def case_summary(
    case_id: str,
    *,
    objective: float = 1.0,
    truth_rmse: float = 0.1,
    accepted_iteration: int = 4,
    elapsed_seconds: float = 1.0,
    trajectory_digest: str = "trajectory",
    terminal_digest: str = "terminal",
    stop_reason: str = "converged",
) -> dict[str, object]:
    return {
        "schema_version": 13,
        "status": "complete",
        "case": {
            "case_id": case_id,
            "family": "natural",
            "topology": "c-c",
            "seed": 410000,
        },
        "elapsed_seconds": elapsed_seconds,
        "production_converged": stop_reason == "converged",
        "semantic_trajectory_sha256": trajectory_digest,
        "terminal_state_sha256": terminal_digest,
        "frozen_truth_sha256": "truth",
        "safety_regression": False,
        "terminal": {
            "reason": stop_reason,
            "try": accepted_iteration,
            "accepted_iteration": accepted_iteration,
            "objective": objective,
            "truth_metrics": {
                "transformed_aggregate_rmse": truth_rmse,
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

    def test_compact_baseline_is_production_only(self) -> None:
        summary = case_summary("case-a")
        aggregate = ANALYZER.analyze([summary])
        baseline = RUNNER.build_compact_baseline(
            b'{"schema_version":1}', [summary], aggregate)
        self.assertEqual(baseline["schema_version"], 3)
        self.assertEqual(baseline["schema_contract"]["trajectory"], 10)
        self.assertEqual(baseline["schema_contract"]["terminal"], 2)
        self.assertEqual(baseline["schema_contract"]["case_summary"], 13)
        self.assertEqual(baseline["schema_contract"]["aggregate"], 8)
        self.assertEqual(baseline["schema_contract"]["comparison"], 4)
        self.assertNotIn("comparator_set", baseline["schema_contract"])
        self.assertNotIn("comparator_definitions", baseline)
        self.assertEqual(baseline["cases"][0][4], "converged")

    def test_truth_parser_and_transformed_truth_metrics(self) -> None:
        truth = RUNNER.parse_truth(
            "Convergence exposure truth: schema=1, case=x, serial=1, "
            "amplitude=6, width=0.5, offset=0.1.")
        metrics = RUNNER._truth_metrics([{
            "serial": "1", "group": "1", "amplitude": "6",
            "width": "0.5", "offset": "0.1",
        }], truth)
        self.assertIsNotNone(metrics)
        assert metrics is not None
        self.assertEqual(metrics["transformed_aggregate_rmse"], 0.0)

    def test_analyzer_reports_only_production_metrics(self) -> None:
        report = ANALYZER.analyze([
            case_summary("a"),
            case_summary("b", stop_reason="audit-patience"),
        ])
        self.assertEqual(report["schema_version"], 8)
        self.assertEqual(report["case_count"], 2)
        self.assertEqual(report["production_convergence_count"], 1)
        self.assertEqual(
            report["termination_counts"],
            {"audit-patience": 1, "converged": 1})
        serialized = json.dumps(report)
        for retired in ("comparator", "exposure", "accepted_only", "rho"):
            self.assertNotIn(retired, serialized)

    def test_paired_gate_checks_digests_deltas_safety_and_cost(self) -> None:
        before = ANALYZER.analyze([
            case_summary("a", elapsed_seconds=2.0),
            case_summary("b", elapsed_seconds=4.0),
        ])
        after = ANALYZER.analyze([
            case_summary("a", elapsed_seconds=1.0),
            case_summary("b", elapsed_seconds=3.0),
        ])
        comparison = ANALYZER.compare(before, after)
        self.assertEqual(comparison["schema_version"], 4)
        self.assertEqual(comparison["production_semantic_match_count"], 2)
        self.assertEqual(comparison["terminal_state_match_count"], 2)
        self.assertEqual(comparison["objective_delta"]["median"], 0.0)
        self.assertEqual(comparison["truth_rmse_delta"]["p90"], 0.0)
        self.assertEqual(
            comparison["accepted_iteration_delta"]["median"], 0.0)
        self.assertTrue(comparison["elapsed_seconds"]["strictly_lower"])
        self.assertTrue(comparison["blocking_gate"]["passed"])

    def test_timing_is_excluded_from_semantic_digest(self) -> None:
        value = [{"try": "1", "certificate": "1/1/1/1/1/1"}]
        self.assertEqual(
            RUNNER.semantic_digest(value), RUNNER.semantic_digest(value))
        before = case_summary("a", elapsed_seconds=10.0)
        after = case_summary("a", elapsed_seconds=1.0)
        self.assertEqual(
            before["semantic_trajectory_sha256"],
            after["semantic_trajectory_sha256"])
        self.assertEqual(
            before["terminal_state_sha256"], after["terminal_state_sha256"])


if __name__ == "__main__":
    unittest.main()
