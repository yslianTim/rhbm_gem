#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = (
    PROJECT_ROOT / "resources" / "tools" / "developer"
    / "analyze_counterfactual_convergence.py"
)
SPEC = importlib.util.spec_from_file_location("counterfactual_audit", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def checkpoint(policy: str, attempt: int, objective: float) -> str:
    return (
        "[Debug] Counterfactual convergence checkpoint: schema=4, comparator-set=1, "
        f"experiment=4-4, policy={policy}, try={attempt}, acc={attempt}, "
        f"extra-try={attempt - 4}, extra-acc={attempt - 4}, extra-ms=0, "
        "final-polish=0, "
        "solver-qualified=1, restricted=0, all-fixed=0, population=2/2/2, "
        f"objective=1/0/0/{objective}, accepted-median=5e-6/5e-6/5e-6, "
        "accepted-p99=1e-5/1e-5/1e-5, "
        "accepted-max=2e-5/2e-5/2e-5, raw-p99=1e-5/1e-5/1e-5, "
        "raw-median=5e-6/5e-6/5e-6, raw-max=2e-5/2e-5/2e-5")


def atom(policy: str, serial: int, amplitude: float) -> str:
    return (
        "[Debug] Counterfactual convergence atom: schema=4, comparator-set=1, "
        f"experiment=4-4, policy={policy}, serial={serial}, "
        "group=1, shape-active=1, offset-active=1, quarantined=0, "
        f"amplitude={amplitude}, width=0.5, offset=1.0")


class CounterfactualConvergenceAnalyzerTest(unittest.TestCase):
    def test_reports_exposures_objective_and_truth_separately(self) -> None:
        text = "\n".join((
            checkpoint("production", 4, 1.0),
            atom("production", 1, 5.0),
            atom("production", 2, 7.0),
            checkpoint("historical-all-selected", 4, 1.0),
            checkpoint("historical-cluster-active-proposal-maximum", 6, 0.998),
            checkpoint("historical-active-proposal", 6, 0.998),
            atom("historical-active-proposal", 1, 6.0),
            atom("historical-active-proposal", 2, 7.0),
            checkpoint("production-maximum", 6, 0.998),
            "[Debug] Counterfactual convergence termination: schema=4, comparator-set=1, "
            "experiment=4-4, reason=all-policies-reached, try=6, acc=6, "
            "extra-try=2, extra-acc=2, extra-ms=0, checkpoints=1/1/1/1/1",
        ))
        parsed = MODULE.parse_log(text)
        truth = {
            1: {"amplitude": 6.0, "width": 0.5, "offset": 1.0},
            2: {"amplitude": 7.0, "width": 0.5, "offset": 1.0},
        }

        report = MODULE.analyze(parsed, truth)
        experiment = report["experiments"][0]

        self.assertFalse(experiment["exposures"]["historical-all-selected"])
        self.assertTrue(experiment["exposures"][
            "historical-cluster-active-proposal-maximum"])
        self.assertTrue(experiment["exposures"]["historical-active-proposal"])
        self.assertTrue(experiment["actual_continuation"])
        solver_qualified = experiment["policies"]["historical-active-proposal"]
        self.assertTrue(solver_qualified["material_objective_improvement"])
        self.assertEqual(
            solver_qualified["truth_metrics"]["amplitude_rmse"], 0.0)
        self.assertEqual(report["exposure_counts"]["historical-active-proposal"], 1)
        self.assertEqual(report["exposure_counts"]["historical-all-selected"], 0)
        self.assertEqual(
            report["exposure_overlap"],
            {"historical-cluster-active-proposal-maximum+historical-active-proposal+production-maximum": 1})
        self.assertEqual(
            report["material_objective_improvement_counts"][
                "historical-active-proposal"], 1)
        self.assertEqual(report["actual_continuation_count"], 1)

    def test_reports_unresolved_policy_and_budget_termination(self) -> None:
        text = "\n".join((
            checkpoint("production", 4, 1.0),
            "[Debug] Counterfactual convergence termination: schema=4, comparator-set=1, "
            "experiment=4-4, reason=budget-exhausted, try=29, acc=10, "
            "extra-try=25, extra-acc=6, extra-ms=0, checkpoints=1/0/0/0/0",
        ))
        report = MODULE.analyze(MODULE.parse_log(text), {})
        experiment = report["experiments"][0]
        self.assertEqual(experiment["termination_reason"], "budget-exhausted")
        self.assertFalse(experiment["policies"][
            "historical-active-proposal"]["reached"])
        self.assertEqual(report["termination_counts"], {"budget-exhausted": 1})
        self.assertEqual(
            report["termination_category_counts"], {"budget_exhaustion": 1})
        self.assertEqual(
            report["unresolved_policy_counts"]["historical-active-proposal"], 1)

    def test_no_checkpoint_is_negative_control(self) -> None:
        report = MODULE.analyze(MODULE.parse_log("ordinary log"), {})
        self.assertEqual(report["status"], "no_convergence_trigger")
        self.assertEqual(report["experiment_count"], 0)

    def test_schema4_rejects_legacy_policy_and_missing_fields(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unknown counterfactual policy"):
            MODULE.parse_log(checkpoint("legacy-population", 4, 1.0))
        with self.assertRaisesRegex(ValueError, "record is missing"):
            MODULE.parse_log(
                "Counterfactual convergence checkpoint: schema=4, "
                "comparator-set=1, experiment=4-4, policy=production")


if __name__ == "__main__":
    unittest.main()
