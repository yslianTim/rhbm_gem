import unittest

from joint_fixed_neighbor_prepared_block import analyze, summarize


def row(topology, atoms, *, symbolic=4, reuses=12, numeric=16,
        preparations=2, converged=True):
    fixed = {
        "search_converged": converged,
        "sweeps": 4,
        "confirmed_stationarity_sweep": 4 if converged else None,
        "blocks_per_sweep": preparations,
        "block_solves": 16,
        "profile_evaluations": 32,
        "accepted_local_updates": 8,
        "objective": 1e-8,
        "final_global_ac_kkt": 1e-12,
        "final_raw_width_gradient_inf_norm": 1e-14,
        "search_seconds": 10.0,
        "total_elapsed_seconds": 11.0,
        "prepared_block_count": preparations,
        "block_preparations": preparations,
        "domain_preparations": preparations,
        "mapping_preparations": preparations,
        "symbolic_factorizations": symbolic,
        "symbolic_reuses": reuses,
        "numeric_factorizations": numeric,
        "fresh_workspace_symbolic_factorizations": numeric,
    }
    return {"topology": topology, "atoms": atoms, "status": "completed",
            "measurement_scope": "fixed-neighbor-search-only", "core_atoms": 12,
            "fixed_neighbor": fixed}


class PreparedBlockQualificationTest(unittest.TestCase):
    def test_all_matched_cases_require_preparation_and_reuse_attribution(self):
        cases = [row(topology, 24) for topology in ("chain", "cube")]
        report = analyze([summarize(case) for case in cases])
        self.assertEqual(report["qualification_gate"], "passed")
        self.assertEqual(len(report["cases"]), 2)
        self.assertEqual(report["wall_time_gate"], "not-run")

    def test_missing_or_unreused_case_fails_without_relaxing_numerical_gate(self):
        cases = [row("chain", 24)]
        broken = row("cube", 24, symbolic=16, reuses=0, numeric=16)
        report = analyze([summarize(case) for case in cases] + [summarize(broken)])
        self.assertEqual(report["qualification_gate"], "failed")
        reasons = report["cases"][-1]["reasons"]
        self.assertIn("symbolic-reuse-gate", reasons)


if __name__ == "__main__":
    unittest.main()
