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
    strict_dof: str,
    stops: str,
    exposures: str,
    strict: str,
    production_population: str = "10/10/2",
    legacy_population: str = "10/10/10",
    member_population: str = "10/10/10",
    production_values: str = "5e-5/5e-5/5e-5",
    legacy_values: str = "5e-4/5e-4/5e-4",
    member_values: str = "5e-4/5e-4/5e-4",
    legacy_maximum: str = "1/1/1/1/1",
) -> dict[str, str]:
    line = (
        "[Debug] Convergence safeguard audit: schema=4, try=1, acc=1, "
        f"production-population={production_population}, "
        f"legacy-population={legacy_population}, "
        f"member-population={member_population}, "
        f"production-predicates[s/a99/r99]={production}, "
        f"legacy-population-predicates[s/a99/r99]={legacy}, "
        f"legacy-maximum-predicates[s/a99/amax/r99/rmax]={legacy_maximum}, "
        f"strict-dof-predicates[s/a99/r99]={strict_dof}, "
        "member-diagnostic-predicates[s/a99/r99]=1/0/0, "
        f"production-accepted-p99={production_values}, "
        "production-accepted-max=5e-4/5e-4/5e-4, "
        f"production-raw-p99={production_values}, "
        "production-raw-max=5e-4/5e-4/5e-4, "
        f"legacy-accepted-p99={legacy_values}, "
        "legacy-accepted-max=5e-4/5e-4/5e-4, "
        f"legacy-raw-p99={legacy_values}, legacy-raw-max=5e-4/5e-4/5e-4, "
        f"member-accepted-p99={member_values}, "
        "member-accepted-max=5e-4/5e-4/5e-4, "
        f"member-raw-p99={member_values}, member-raw-max=5e-4/5e-4/5e-4, "
        "path[trust/backtrack/polish/boundary/rescue]=0/0/0/0/0, "
        "stationarity[current/full/active-ineligible/refit-ineligible/soft-joint/hard-joint]=1/0/0/0/1/0, "
        "joint-status[converged/system-build/empty/initial-solve/irls-solve/objective-deteriorated/max-iter]=1/0/0/0/0/1/0, "
        f"strict-stationarity[current/strict/restricted/all-fixed/active-shape/qualified-shape/soft-shape/hard-shape/fixed-shape/quarantine-shape/active-offset/qualified-offset/soft-offset/hard-offset/fixed-offset/quarantine-offset/mixed-offset]={strict}, "
        "offset-groups[total/active/fixed/quarantine/mixed/member-count/min/p50/p99/max]=2/2/0/0/0/10/1/5/9/9, "
        "ratios[shape-active/offset-member-active/quarantine]=0.01/0.01/0, "
        f"stop-candidates[orthogonal-clear/production/legacy-population/legacy-maximum/strict-dof]={stops}, "
        f"exposures[legacy-population/maximum-gate/strict-stationarity]={exposures}."
    )
    record = MODULE.parse_record(line)
    assert record is not None
    return record


class ConvergenceAuditAnalyzerTest(unittest.TestCase):
    def test_aggregates_exposures_implications_and_strata(self) -> None:
        exposed = make_record(
            production="1/1/1",
            legacy="1/0/0",
            strict_dof="0/1/1",
            stops="1/1/0/1/0",
            exposures="1/0/1",
            strict="1/0/0/0/10/10/0/0/0/0/2/1/1/0/0/0/0",
        )
        stable = make_record(
            production="1/1/1",
            legacy="1/1/1",
            strict_dof="1/1/1",
            stops="1/1/1/1/1",
            exposures="0/0/0",
            strict="1/1/0/0/10/10/0/0/0/0/2/2/0/0/0/0/0",
            legacy_values="5e-5/5e-5/5e-5",
            member_values="5e-5/5e-5/5e-5",
        )

        report = MODULE.analyze_records([exposed, stable])

        self.assertEqual(report["record_count"], 2)
        self.assertEqual(report["actual_stop_exposure_count"], 1)
        self.assertEqual(
            report["implication_counterexamples"][
                "production_stationarity=>strict"],
            1,
        )
        self.assertEqual(
            report["implication_counterexamples"][
                "production_accepted_p99=>legacy_population_accepted_p99"
            ],
            1,
        )
        self.assertEqual(report["exposure_count"]["legacy_population"], 1)
        self.assertEqual(
            report["strata"]["solver_status"][
                "converged+objective_deteriorated:actual_exposures"
            ],
            1,
        )
        self.assertEqual(
            report["truth_tables"]["legacy_population"][
                "stationarity|accepted_p99"
            ]["10"],
            1,
        )

    def test_counts_unique_blocker(self) -> None:
        record = make_record(
            production="1/1/1",
            legacy="1/0/1",
            strict_dof="1/1/1",
            stops="0/0/0/0/0",
            exposures="0/0/0",
            strict="1/1/0/0/10/10/0/0/0/0/2/2/0/0/0/0/0",
        )
        report = MODULE.analyze_records([record])
        self.assertEqual(
            report["unique_blocker_count"]["legacy_population"]["accepted_p99"],
            1,
        )

    def test_ignores_older_schema(self) -> None:
        self.assertIsNone(
            MODULE.parse_record("Convergence safeguard audit: schema=3, try=1")
        )


if __name__ == "__main__":
    unittest.main()
