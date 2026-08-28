#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = (
    PROJECT_ROOT
    / "resources"
    / "tools"
    / "developer"
    / "analyze_convergence_audit.py"
)
SPEC = importlib.util.spec_from_file_location("convergence_audit", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def make_record(
    *,
    production: str,
    legacy: str,
    solver_qualified: str,
    stops: str,
    exposures: str,
    qualification: str,
    production_population: str = "10/10/2",
    legacy_population: str = "10/10/10",
    production_values: str = "5e-5/5e-5/5e-5",
    legacy_values: str = "5e-4/5e-4/5e-4",
    legacy_maximum: str = "1/1/1/1/1",
) -> dict[str, str]:
    line = (
        "[Debug] Convergence safeguard audit: schema=7, try=1, acc=1, "
        f"production-population={production_population}, "
        f"legacy-population={legacy_population}, "
        f"operator-population={legacy_population}, "
        f"production-predicates[q/a99/r99]={production}, "
        f"legacy-population-predicates[q/a99/r99]={legacy}, "
        f"legacy-maximum-predicates[q/a99/amax/r99/rmax]={legacy_maximum}, "
        f"solver-qualified-predicates[q/a99/r99]={solver_qualified}, "
        f"operator-predicates[q/a99/fp99/complete/fpmax]={production}/1/1, "
        f"production-accepted-p99={production_values}, "
        "production-accepted-max=5e-4/5e-4/5e-4, "
        f"legacy-accepted-p99={legacy_values}, "
        "legacy-accepted-max=5e-4/5e-4/5e-4, "
        f"fixed-point-residual-p99={production_values}, "
        "fixed-point-residual-max=5e-4/5e-4/5e-4, "
        "operator-unavailable[height/width/offset]=0/0/0, "
        "operator-unavailable-reasons[offset-solver/invalid-offset/shape-refit]=0/0/0, "
        "operator-tail[height/width/offset]=0/0/0, "
        "residual-state=fixed-point-converged, "
        "unified-search[trials/invalid/trust-skipped/guard-rejected/objective-rejected/accepted-limited/terminal]=1/0/0/0/0/0/0, "
        "path[limited/polish/boundary/rescue]=0/0/0/0, "
        "limiters[guard/fixed/quarantine/trust/objective/reject/polish/boundary/rescue]=0/0/0/0/0/0/0/0/0, "
        "joint-status[converged/system-build/empty/initial-solve/irls-solve/objective-deteriorated/max-iter]=1/0/0/0/0/1/0, "
        f"qualification[production/solver/restricted/all-fixed/active-shape/solver-shape/soft-shape/hard-shape/fixed-shape/quarantine-shape/active-offset/solver-offset/soft-offset/hard-offset/fixed-offset/quarantine-offset/mixed-offset]={qualification}, "
        "offset-groups[total/active/fixed/quarantine/mixed/min/p50/p99/max]=2/2/0/0/0/1/5/9/9, "
        "ratios[shape-active/offset-active/quarantine]=0.01/0.01/0, "
        "blockers[suspicious/rejected/quarantine-transition/domain-change]=0/0/0/0, "
        f"stop-candidates[orthogonal-clear/production/legacy-population/legacy-maximum/solver-qualified/operator/operator-maximum]={stops}/1/1, "
        f"exposures[legacy-population/maximum-gate/solver-qualification/operator/operator-maximum]={exposures}/0/0."
    )
    record = MODULE.parse_record(line)
    assert record is not None
    return record


def make_schema8_record(certificate: str = "1/1/1/1/1/1/1") -> dict[str, str]:
    line = (
        "[Debug] Convergence safeguard audit: schema=8, "
        "certificate-definition=1, comparator-set=1, try=1, acc=1, "
        "accepted-active-population=10/10/2, "
        "historical-all-selected-population=10/10/10, "
        "operator-nominal-population=10/10/2, "
        f"certificate[solver/accepted-p99/operator-complete/operator-p99/invariants/orthogonal/production]={certificate}, "
        "historical-all-selected-predicates[q/a99/op99]=1/1/1, "
        "historical-cluster-active-proposal-maximum-predicates[q/a99/amax/op99/opmax]=1/1/1/1/1, "
        "historical-active-proposal-predicates[q/a99/op99]=1/1/1, "
        "production-maximum-predicates[production/accepted-max/operator-max]=1/1/1, "
        "accepted-active-p99=5e-5/5e-5/5e-5, "
        "accepted-active-max=5e-4/5e-4/5e-4, "
        "historical-all-selected-accepted-p99=5e-5/5e-5/5e-5, "
        "historical-all-selected-accepted-max=5e-4/5e-4/5e-4, "
        "historical-active-proposal-p99=5e-5/5e-5/5e-5, "
        "historical-active-proposal-max=5e-4/5e-4/5e-4, "
        "operator-nominal-residual-p99=5e-5/5e-5/5e-5, "
        "operator-nominal-residual-max=5e-4/5e-4/5e-4, "
        "operator-nominal-unavailable[height/width/offset]=0/0/0, "
        "operator-nominal-unavailable-reasons[offset-solver/invalid-offset/shape-refit]=0/0/0, "
        "operator-nominal-tail[height/width/offset]=0/0/0, "
        "residual-state=fixed-point-converged, "
        "unified-search[trials/invalid/trust-skipped/guard-rejected/objective-rejected/accepted-limited/terminal]=1/0/0/0/0/0/0, "
        "path[limited/polish/boundary/rescue]=0/0/0/0, "
        "limiters[guard/fixed/quarantine/trust/objective/reject/polish/boundary/rescue]=0/0/0/0/0/0/0/0/0, "
        "joint-status[converged/system-build/empty/initial-solve/irls-solve/objective-deteriorated/max-iter]=1/0/0/0/0/0/0, "
        "qualification[production/solver/restricted/all-fixed/active-shape/solver-shape/soft-shape/hard-shape/fixed-shape/quarantine-shape/active-offset/solver-offset/soft-offset/hard-offset/fixed-offset/quarantine-offset/mixed-offset]=1/1/0/0/10/10/0/0/0/0/2/2/0/0/0/0/0, "
        "offset-groups[total/active/fixed/quarantine/mixed/min/p50/p99/max]=2/2/0/0/0/1/1/1/1, "
        "ratios[shape-active/offset-active/quarantine]=1/1/0, "
        "certificate-blockers[objective-domain/quarantine-transition/suspicious-offset/rejected-cluster]=0/0/0/0, "
        "stop-candidates[production/historical-all-selected/historical-cluster-active-proposal-maximum/historical-active-proposal/production-maximum]=1/1/1/1/1, "
        "exposures[historical-all-selected/historical-cluster-active-proposal-maximum/historical-active-proposal/production-maximum]=0/0/0/0."
    )
    record = MODULE.parse_record(line)
    assert record is not None
    return record


class ConvergenceAuditAnalyzerTest(unittest.TestCase):
    def test_aggregates_exposures_implications_and_strata(self) -> None:
        exposed = make_record(
            production="1/1/1",
            legacy="1/0/0",
            solver_qualified="0/1/1",
            stops="1/1/0/1/0",
            exposures="1/0/1",
            qualification="1/0/0/0/10/10/0/0/0/0/2/1/1/0/0/0/0",
        )
        stable = make_record(
            production="1/1/1",
            legacy="1/1/1",
            solver_qualified="1/1/1",
            stops="1/1/1/1/1",
            exposures="0/0/0",
            qualification="1/1/0/0/10/10/0/0/0/0/2/2/0/0/0/0/0",
            legacy_values="5e-5/5e-5/5e-5",
        )

        report = MODULE.analyze_records([exposed, stable])

        self.assertEqual(report["record_count"], 2)
        self.assertEqual(report["actual_stop_exposure_count"], 1)
        self.assertNotIn(
            "production_qualification=>solver_qualified",
            report["implication_counterexamples"],
        )
        self.assertEqual(
            report["implication_counterexamples"][
                "production_accepted_p99=>historical_all_selected_accepted_p99"
            ],
            1,
        )
        self.assertEqual(report["exposure_count"]["historical-all-selected"], 1)
        self.assertEqual(
            report["strata"]["solver_status"][
                "converged+objective_deteriorated:actual_exposures"
            ],
            1,
        )
        self.assertEqual(
            report["truth_tables"]["historical_all_selected"][
                "qualification|accepted_p99"
            ]["10"],
            1,
        )

    def test_counts_unique_blocker(self) -> None:
        record = make_record(
            production="1/1/1",
            legacy="1/0/1",
            solver_qualified="1/1/1",
            stops="0/0/0/0/0",
            exposures="0/0/0",
            qualification="1/1/0/0/10/10/0/0/0/0/2/2/0/0/0/0/0",
        )
        report = MODULE.analyze_records([record])
        self.assertEqual(
            report["unique_blocker_count"]["historical_all_selected"]["accepted_p99"],
            1,
        )

    def test_schema8_contract_and_serialized_decision(self) -> None:
        record = make_schema8_record()
        report = MODULE.analyze_records([record])
        self.assertEqual(record["certificate-definition"], "1")
        self.assertEqual(record["comparator-set"], "1")
        self.assertEqual(report["record_count"], 1)
        with self.assertRaisesRegex(ValueError, "contradicts certificate"):
            MODULE.analyze_records([make_schema8_record("1/1/1/1/1/1/0")])

    def test_schema8_rejects_mixed_or_reordered_definitions(self) -> None:
        line = (
            "Convergence safeguard audit: schema=8, certificate-definition=1, "
            "comparator-set=1, production-population=1/1/1")
        with self.assertRaisesRegex(ValueError, "record is missing|mixes"):
            MODULE.parse_record(line)

        record_line = (
            "[Debug] Convergence safeguard audit: schema=8, "
            "certificate-definition=1, comparator-set=1, try=1, acc=1, "
            "accepted-active-population=1/1/1, "
            "historical-all-selected-population=1/1/1, "
            "operator-nominal-population=1/1/1, "
            "certificate[solver/accepted-p99/operator-complete/operator-p99/invariants/orthogonal/production]=1/1/1/1/1/1/1, "
            "historical-all-selected-predicates[q/a99/op99]=1/1/1, "
            "historical-cluster-active-proposal-maximum-predicates[q/a99/amax/op99/opmax]=1/1/1/1/1, "
            "historical-active-proposal-predicates[q/a99/op99]=1/1/1, "
            "production-maximum-predicates[production/accepted-max/operator-max]=1/1/1, "
            "accepted-active-p99=0/0/0, accepted-active-max=0/0/0, "
            "historical-all-selected-accepted-p99=0/0/0, "
            "historical-all-selected-accepted-max=0/0/0, "
            "historical-active-proposal-p99=0/0/0, "
            "historical-active-proposal-max=0/0/0, "
            "operator-nominal-residual-p99=0/0/0, "
            "operator-nominal-residual-max=0/0/0, "
            "operator-nominal-unavailable=0/0/0, "
            "operator-nominal-unavailable-reasons=0/0/0, "
            "operator-nominal-tail=0/0/0, residual-state=fixed-point-converged, "
            "qualification=1/1/0/0/1/1/0/0/0/0/1/1/0/0/0/0/0, "
            "path=0/0/0/0, ratios=1/1/0, offset-groups=1/1/0/0/0/1/1/1/1, "
            "joint-status=1/0/0/0/0/0/0, limiters=0/0/0/0/0/0/0/0/0, "
            "certificate-blockers=0/0/0/0, unified-search=1/0/0/0/0/0/0, "
            "stop-candidates[production/historical-all-selected/historical-active-proposal/historical-cluster-active-proposal-maximum/production-maximum]=1/1/1/1/1, "
            "exposures[historical-all-selected/historical-cluster-active-proposal-maximum/historical-active-proposal/production-maximum]=0/0/0/0")
        with self.assertRaisesRegex(ValueError, "names or order"):
            MODULE.parse_record(record_line)

    def test_schema7_normalizes_without_mixing_field_definitions(self) -> None:
        record = make_record(
            production="1/1/1",
            legacy="1/1/1",
            solver_qualified="1/1/1",
            stops="1/1/1/1/1",
            exposures="0/0/0",
            qualification="1/1/0/0/10/10/0/0/0/0/2/2/0/0/0/0/0",
        )
        self.assertEqual(record["source-schema"], "7")
        self.assertIn("accepted-active-population", record)
        self.assertIn("operator-nominal-population", record)

    def test_ignores_older_schema(self) -> None:
        self.assertIsNone(
            MODULE.parse_record("Convergence safeguard audit: schema=5, try=1")
        )


if __name__ == "__main__":
    unittest.main()
