#!/usr/bin/env python3
"""Contract tests for the production-only convergence audit parser."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
ANALYZER_PATH = (
    PROJECT_ROOT / "resources" / "tools" / "developer" /
    "analyze_convergence_audit.py")
SPEC = importlib.util.spec_from_file_location("convergence_audit", ANALYZER_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def schema9_record(certificate: str = "1/1/1/1/1/1/1") -> str:
    return (
        "[Debug] Convergence safeguard audit: schema=9, "
        "certificate-definition=1, try=4, acc=4, atoms=2, quarantine=0, "
        "accepted-active-population=2/2/1, operator-nominal-population=2/2/1, "
        f"certificate[solver/accepted-p99/operator-complete/operator-p99/invariants/orthogonal/production]={certificate}, "
        "accepted-active-p99=5e-5/5e-5/5e-5, accepted-active-max=5e-4/5e-4/5e-4, "
        "operator-nominal-residual-p99=5e-5/5e-5/5e-5, "
        "operator-nominal-residual-max=5e-4/5e-4/5e-4, "
        "operator-nominal-unavailable[height/width/offset]=0/0/0, "
        "operator-nominal-unavailable-reasons[offset-solver/invalid-offset/shape-refit]=0/0/0, "
        "operator-nominal-tail[height/width/offset]=0/0/0, "
        "residual-state=fixed-point-converged, operator-shadow-refit=0, "
        "accepted-equals-operator=1, "
        "unified-search[trials/invalid/trust-skipped/guard-rejected/objective-rejected/accepted-limited/terminal]=1/0/0/0/0/0/0, "
        "path[limited/polish/boundary/rescue]=0/0/0/0, "
        "joint-status[converged/system-build/empty/initial-solve/irls-solve/objective-deteriorated/max-iter]=1/0/0/0/0/0/0, "
        "qualification[production/solver/restricted/all-fixed/active-shape/solver-shape/soft-shape/hard-shape/fixed-shape/quarantine-shape/active-offset/solver-offset/soft-offset/hard-offset/fixed-offset/quarantine-offset/mixed-offset]=1/1/0/0/2/2/0/0/0/0/1/1/0/0/0/0/0, "
        "local-status[success/max-iter/single/insufficient/numerical/unavailable]=2/0/0/0/0/0, "
        "offset-groups[total/active/fixed/quarantine/mixed]=1/1/0/0/0, "
        "ratios[shape-active/offset-active/quarantine]=1/1/0, "
        "certificate-blockers[objective-domain/quarantine-transition/suspicious-offset/rejected-cluster]=0/0/0/0, "
        "limiters[guard/fixed/quarantine/trust/objective/reject/polish/boundary/rescue]=0/0/0/0/0/0/0/0/0, "
        "fixed[shape/offset/hard]=0/0/0, "
        "blockers[suspicious/rejected/quarantine-transition/domain-change]=0/0/0/0."
    )


class ConvergenceAuditAnalyzerTest(unittest.TestCase):
    def test_schema9_accepts_production_certificate(self) -> None:
        record = MODULE.parse_record(schema9_record())
        self.assertIsNotNone(record)
        assert record is not None
        self.assertEqual(record["schema"], "9")
        self.assertEqual(record["certificate"], "1/1/1/1/1/1/1")
        self.assertEqual(MODULE.analyze_records([record])["record_count"], 1)

    def test_schema9_rejects_nonconjunctive_or_retired_fields(self) -> None:
        with self.assertRaises(ValueError):
            MODULE.parse_record(schema9_record("1/1/1/1/1/1/0"))
        with self.assertRaises(ValueError):
            MODULE.parse_record(
                schema9_record().replace(
                    ", try=4", ", comparator-set=1, try=4"))

    def test_schema8_has_one_frozen_compatibility_fixture(self) -> None:
        fixture = (
            PROJECT_ROOT / "tests" / "fixtures" / "convergence-audit" /
            "schema-8-record.log").read_text(encoding="utf-8")
        record = MODULE.parse_record(fixture)
        self.assertIsNotNone(record)
        assert record is not None
        normalized = MODULE.normalize_record(record)
        self.assertEqual(normalized["source-schema"], "8")
        self.assertEqual(normalized["certificate"], "1/1/1/1/1/1/1")
        self.assertNotIn("comparator-set", normalized)

if __name__ == "__main__":
    unittest.main()
