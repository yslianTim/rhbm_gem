#!/usr/bin/env python3
"""Contract tests for the developer-only trust-model analyzer."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PATH = (PROJECT_ROOT / "resources" / "tools" / "developer" /
        "analyze_trust_model_experiment.py")
SPEC = importlib.util.spec_from_file_location("trust_model_analyzer", PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class TrustModelExperimentAnalyzerTest(unittest.TestCase):
    def test_trial_diagnostics_are_experiment_only(self) -> None:
        text = "\n".join((
            "[Debug] Trust-model funnel: schema=1, try=2, acc=1, atoms=2, "
            "key-first=1, key-last=2, disposition=accepted, generated=1, "
            "invalid=0, trust-skipped=0, guard-rejected=0, nonmaterial=0, "
            "objective-evaluated=1, polish-objective-evaluated=0",
            "[Debug] Trust-model shadow: schema=2, try=2, acc=1, atoms=2, "
            "key-first=1, key-last=2, disposition=accepted, "
            "boundary-touched=1, boundary-rescued=0, readiness-eligible=1, "
            "final-local-candidate=1, status=available, source=base, "
            "search-pass=1, trial=1, factor=1, trial-disposition=accepted, "
            "rejected-by-previous=0, rejected-by-best=0, "
            "rejected-by-strict-polish=0, step-norm=1, "
            "actual-reduction=1, polish-reduction=-, "
            "predicted-residual-reduction=1, predicted-penalty-reduction=0, "
            "predicted-reduction=1, rho=0.5, boundary-utilization=1, "
            "current-action=keep, shadow-action=keep, "
            "objective-backtracked=0, unselected-dependencies=0, "
            "elapsed-ms=0.25",
        ))
        report = MODULE.analyze(MODULE.parse_log(text))
        self.assertTrue(report["diagnostic_only"])
        self.assertEqual(report["funnel_record_count"], 1)
        self.assertEqual(report["trial_record_count"], 1)
        self.assertEqual(report["rho_bins"], {"mid": 1})
        self.assertEqual(report["action_confusion"], {"keep->keep": 1})
        self.assertEqual(report["elapsed_milliseconds"], 0.25)


if __name__ == "__main__":
    unittest.main()
