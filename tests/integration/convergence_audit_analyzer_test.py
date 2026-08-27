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
    current: str,
    member: str,
    dof: str,
    stops: str,
    strict: str,
    current_population: str = "1000/1000/1000",
    member_population: str = "10/10/10",
    dof_population: str = "10/10/2",
    current_values: str = "5e-5/5e-5/5e-5",
    member_values: str = "5e-4/5e-4/5e-4",
    dof_values: str = "5e-4/5e-4/5e-4",
) -> dict[str, str]:
    line = (
        "[Debug] Convergence safeguard audit: schema=2, try=1, acc=1, "
        f"population={current_population}, shadow-population={member_population}, "
        f"active-dof-population={dof_population}, "
        f"predicates[s/a99/amax/r99/rmax]={current}, "
        "shadow-predicates[s/a99/amax/r99/rmax]=1/0/1/0/1, "
        f"member-strict-predicates[s/a99/amax/r99/rmax]={member}, "
        f"dof-strict-predicates[s/a99/amax/r99/rmax]={dof}, "
        f"accepted-p99={current_values}, accepted-max=5e-4/5e-4/5e-4, "
        f"raw-p99={current_values}, raw-max=5e-4/5e-4/5e-4, "
        f"shadow-accepted-p99={member_values}, shadow-accepted-max=5e-4/5e-4/5e-4, "
        f"shadow-raw-p99={member_values}, shadow-raw-max=5e-4/5e-4/5e-4, "
        f"dof-accepted-p99={dof_values}, dof-accepted-max=5e-4/5e-4/5e-4, "
        f"dof-raw-p99={dof_values}, dof-raw-max=5e-4/5e-4/5e-4, "
        "path[trust/backtrack/polish/boundary/rescue]=0/0/0/0/0, "
        "stationarity[current/full/active-ineligible/refit-ineligible/soft-joint/hard-joint]=1/0/0/0/1/0, "
        "joint-status[converged/system-build/empty/initial-solve/irls-solve/objective-deteriorated/max-iter]=1/0/0/0/0/1/0, "
        f"strict-stationarity[current/strict/restricted/all-fixed/active-shape/qualified-shape/soft-shape/hard-shape/fixed-shape/quarantine-shape/active-offset/qualified-offset/soft-offset/hard-offset/fixed-offset/quarantine-offset/mixed-offset]={strict}, "
        "offset-groups[total/active/fixed/quarantine/mixed/member-count/min/p50/p99/max]=2/2/0/0/0/10/1/5/9/9, "
        "ratios[shape-active/offset-member-active/quarantine]=0.01/0.01/0, "
        f"stop-candidates[orthogonal-clear/production/member/dof/stationarity-exposure/member-population-exposure/dof-population-exposure]={stops}."
    )
    record = MODULE.parse_record(line)
    assert record is not None
    return record


class ConvergenceAuditAnalyzerTest(unittest.TestCase):
    def test_aggregates_exposures_implications_and_strata(self) -> None:
        exposed = make_record(
            current="1/1/1/1/1",
            member="0/0/1/0/1",
            dof="0/0/1/0/1",
            stops="1/1/0/0/1/1/1",
            strict="1/0/0/0/10/10/0/0/0/0/2/1/1/0/0/0/0",
        )
        stable = make_record(
            current="1/1/1/1/1",
            member="1/1/1/1/1",
            dof="1/1/1/1/1",
            stops="1/1/1/1/0/0/0",
            strict="1/1/0/0/10/10/0/0/0/0/2/2/0/0/0/0/0",
            member_values="5e-5/5e-5/5e-5",
            dof_values="5e-5/5e-5/5e-5",
        )

        report = MODULE.analyze_records([exposed, stable])

        self.assertEqual(report["record_count"], 2)
        self.assertEqual(report["actual_stop_exposure_count"], 1)
        self.assertEqual(
            report["implication_counterexamples"]["current_stationarity=>strict"],
            1,
        )
        self.assertEqual(
            report["implication_counterexamples"][
                "current_accepted_p99=>active_member_accepted_p99"
            ],
            1,
        )
        self.assertEqual(
            report["strata"]["solver_status"][
                "converged+objective_deteriorated:actual_exposures"
            ],
            1,
        )
        self.assertEqual(
            report["truth_tables"]["active_member"][
                "stationarity|accepted_p99"
            ]["00"],
            1,
        )

    def test_counts_unique_blocker(self) -> None:
        record = make_record(
            current="1/1/1/1/1",
            member="1/1/0/1/1",
            dof="1/1/1/1/1",
            stops="0/0/0/0/0/0/0",
            strict="1/1/0/0/10/10/0/0/0/0/2/2/0/0/0/0/0",
            member_values="5e-5/5e-5/5e-5",
            dof_values="5e-5/5e-5/5e-5",
        )
        report = MODULE.analyze_records([record])
        self.assertEqual(
            report["unique_blocker_count"]["active_member"]["accepted_max"],
            1,
        )

    def test_ignores_older_schema(self) -> None:
        self.assertIsNone(
            MODULE.parse_record("Convergence safeguard audit: schema=1, try=1")
        )


if __name__ == "__main__":
    unittest.main()
