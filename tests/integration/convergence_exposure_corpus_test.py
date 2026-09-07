#!/usr/bin/env python3
"""Contract tests for production-only convergence corpus tooling."""

from __future__ import annotations

import importlib.util
import json
import copy
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
    trajectory_digest: str = "a" * 64,
    terminal_digest: str = "b" * 64,
    stop_reason: str = "converged",
) -> dict[str, object]:
    return {
        "schema_version": 14,
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
        "frozen_truth_sha256": "c" * 64,
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
        self.assertEqual(baseline["schema_version"], 4)
        self.assertEqual(baseline["schema_contract"]["trajectory"], 10)
        self.assertEqual(baseline["schema_contract"]["terminal"], 2)
        self.assertEqual(baseline["schema_contract"]["case_summary"], 14)
        self.assertEqual(baseline["schema_contract"]["aggregate"], 9)
        self.assertEqual(baseline["schema_contract"]["comparison"], 5)
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
        target = {"serial": "1", "amplitude": "6", "width": "0.5", "offset": "0.1"}
        neighbor = {"serial": "2", "amplitude": "1", "width": "0.5", "offset": "0"}
        self.assertEqual(RUNNER._truth_metrics([target, neighbor], truth)["transformed_aggregate_rmse"], 0.0)
        self.assertIsNone(RUNNER._truth_metrics([target, target], truth))
        self.assertIsNone(RUNNER._truth_metrics([neighbor], truth))
        parser = load_module(
            "independent_offset_audit_parser",
            PROJECT_ROOT / "resources" / "tools" / "developer" /
            "analyze_convergence_audit.py")
        for schema, group_field in (("1", ", group=7"), ("2", "")):
            record = (
                f"Second-stage audit terminal atom: schema={schema}, serial=1"
                f"{group_field}, amplitude=6, width=0.5, offset=0.1.")
            parsed = parser.parse_log(record)["terminal_atoms"]
            self.assertEqual(len(parsed), 1)
            self.assertEqual(parsed[0]["offset"], "0.1")
            self.assertEqual("group" in parsed[0], schema == "1")
            terminal = RUNNER.normalized_terminal({
                "terminal": {"reason": "converged", "try": "1", "acc": "1",
                             "objective": "0/0/0/0"},
                "terminal_atoms": parsed,
            })
            expected_atom = {name: parsed[0][name] for name in
                             ("serial", "amplitude", "width", "offset")}
            if schema == "1":
                expected_atom["group"] = "7"
            self.assertEqual(terminal["atoms"], [expected_atom])
        with self.assertRaises(ValueError):
            parser.parse_log(
                "Second-stage audit terminal atom: schema=2, serial=1, width=0.5.")
        with self.assertRaises(ValueError):
            parser.parse_log(
                "Second-stage audit terminal atom: schema=1, serial=1, "
                "amplitude=6, width=0.5, offset=0.1.")

    def test_analyzer_reports_only_production_metrics(self) -> None:
        report = ANALYZER.analyze([
            case_summary("a"),
            case_summary("b", stop_reason="audit-patience"),
        ])
        self.assertEqual(report["schema_version"], 9)
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
        self.assertEqual(comparison["schema_version"], 5)
        self.assertEqual(comparison["production_semantic_match_count"], 2)
        self.assertEqual(comparison["terminal_state_match_count"], 2)
        self.assertEqual(comparison["objective_delta"]["median"], 0.0)
        self.assertEqual(comparison["truth_rmse_delta"]["p90"], 0.0)
        self.assertEqual(
            comparison["accepted_iteration_delta"]["median"], 0.0)
        self.assertTrue(comparison["elapsed_seconds"]["strictly_lower"])
        self.assertFalse(comparison["safety_gate"]["passed"])  # Only two cases.
        self.assertFalse(comparison["quality_gate"]["passed"])

    def test_full_pair_separates_safety_quality_and_efficiency(self) -> None:
        before = ANALYZER.analyze([case_summary(str(i)) for i in range(600)])
        after = copy.deepcopy(before)
        report = ANALYZER.compare(before, after)
        self.assertTrue(report["safety_gate"]["passed"])
        self.assertTrue(report["quality_gate"]["passed"])
        self.assertFalse(report["efficiency_gate"]["passed"])
        after["cases"][0]["objective"] += 0.01
        self.assertFalse(ANALYZER.compare(before, after)["quality_gate"]["passed"])

    def test_missing_and_nonfinite_evidence_fails_closed(self) -> None:
        before = ANALYZER.analyze([case_summary(str(i)) for i in range(600)])
        for field, value in (("objective", None), ("truth_rmse", float("nan")),
                             ("accepted_iteration", float("inf")),
                             ("semantic_trajectory_sha256", None), ("terminal_state_sha256", None)):
            after = copy.deepcopy(before)
            after["cases"][0][field] = value
            self.assertFalse(ANALYZER.compare(before, after)["quality_gate"]["passed"])
        after = copy.deepcopy(before)
        after["cases"][0]["safety_regression"] = None
        self.assertFalse(ANALYZER.compare(before, after)["safety_gate"]["passed"])
        after = copy.deepcopy(before)
        after["cases"][0]["elapsed_seconds"] = None
        self.assertFalse(ANALYZER.compare(before, after)["efficiency_gate"]["passed"])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            ANALYZER.analyze([case_summary("same"), case_summary("same")])
        with self.assertRaisesRegex(ValueError, "identities differ"):
            ANALYZER.compare(before, ANALYZER.analyze([case_summary("missing")]))
        self.assertFalse(ANALYZER.compare(ANALYZER.analyze([]), ANALYZER.analyze([]))["safety_gate"]["passed"])

    def test_safety_requires_certificate_and_polish_evidence(self) -> None:
        record = {"certificate": "1/1/1/1/1/1", "blockers": "0/0/0/0",
                  "accepted-active-p99": "0/0/0", "operator-nominal-residual-p99": "0/0/0"}
        parsed = {"terminal": {"objective": "1/0/0/1", "reason": "converged", "acc": "1"},
                  "terminal_atoms": [{"amplitude": "1", "width": "1", "offset": "0"}],
                  "trajectory_records": [record],
                  "polish": {"accepted": "no", "applied": "no", "residual-safety": "not-evaluated"}}
        self.assertFalse(RUNNER._safety_regression(parsed))
        rejected = copy.deepcopy(parsed)
        rejected["terminal"].update({"reason": "all-rejected-backtracking-exhausted", "acc": "0"})
        rejected["trajectory_records"] = []
        self.assertFalse(RUNNER._safety_regression(rejected))
        rejected["terminal"]["acc"] = "1"
        self.assertTrue(RUNNER._safety_regression(rejected))
        for key in ("polish", "trajectory_records", "terminal_atoms"):
            bad = copy.deepcopy(parsed)
            bad.pop(key)
            self.assertTrue(RUNNER._safety_regression(bad))
        bad = copy.deepcopy(parsed)
        bad["trajectory_records"][0]["operator-nominal-residual-p99"] = "0/inf/0"
        self.assertTrue(RUNNER._safety_regression(bad))
        bad = copy.deepcopy(parsed)
        bad["polish"].update({"accepted": "yes", "applied": "yes", "residual-safety": "failed"})
        self.assertTrue(RUNNER._safety_regression(bad))

    def test_diagnostics_are_validated_and_excluded_from_production_digest(self) -> None:
        trajectory = ("Convergence safeguard audit: schema=10, try=1, acc=1, atoms=1, "
                      "quarantine=0, accepted-active-population=1/1/1, operator-nominal-population=1/1/1, "
                      "certificate=1/1/1/1/1/1, accepted-active-p99=0/0/0, accepted-active-max=0/0/0, "
                      "operator-nominal-residual-p99=0/0/0, operator-nominal-residual-max=0/0/0, blockers=0/0/0/0.")
        diagnostic = ("Second-stage conditioning: schema=1, solve=joint-offset, phase=outer-operator, "
                      "columns=1, pivot-ratio=1, guard=0, ridge-min=1, ridge-max=1.")
        baseline = RUNNER.CONVERGENCE_ANALYZER.parse_log(trajectory)
        candidate = RUNNER.CONVERGENCE_ANALYZER.parse_log(trajectory + "\n" + diagnostic)
        self.assertEqual(RUNNER.semantic_digest(RUNNER.semantic_trajectory(baseline["trajectory_records"])),
                         RUNNER.semantic_digest(RUNNER.semantic_trajectory(candidate["trajectory_records"])))
        self.assertEqual(len(candidate["diagnostics"]["conditioning"]), 1)
        for invalid in (diagnostic.replace("pivot-ratio=1", "pivot-ratio=nan"),
                        diagnostic.replace(", columns=1", ""), diagnostic.replace("schema=1", "schema=9")):
            with self.assertRaises(ValueError):
                RUNNER.CONVERGENCE_ANALYZER.parse_log(invalid)
        summary = case_summary("a")
        summary["diagnostics"] = candidate["diagnostics"]
        stats = ANALYZER.analyze([summary])["diagnostics"]
        self.assertEqual(stats["conditioning_case_count"], 1)
        self.assertEqual(stats["conditioning"]["all/outer-operator/joint-offset"]["pivot_ratio"]["minimum"], 1.0)

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
