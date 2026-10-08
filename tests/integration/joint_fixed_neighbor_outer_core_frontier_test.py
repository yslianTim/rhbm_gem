import unittest

from joint_fixed_neighbor_outer_core_frontier import analyze


def row(topology, atoms, core, repetition, search, total, rss, correct=True):
    return {
        "topology": topology,
        "atoms": atoms,
        "core_atoms": core,
        "warmup": False,
        "repetition": repetition,
        "status": "completed",
        "search_converged": correct,
        "search_reason": "block-stationary" if correct else "maximum-sweeps",
        "final_global_ac_kkt": 1e-12 if correct else None,
        "final_raw_width_gradient_inf_norm": 1e-14 if correct else None,
        "cache_replay_within_limit": correct,
        "objective_replay_within_limit": correct,
        "search_seconds": search,
        "total_seconds": total,
        "peak_rss_mb": rss,
        "sweeps": 18,
        "block_solves": 100,
        "local_work_seconds": search * 0.7,
        "maximum_local_rows": core * 2,
        "maximum_local_columns": core,
    }


def frontier_rows(searches, totals=None, rss=100.0, incorrect=None):
    totals = totals or {}
    incorrect = incorrect or set()
    rows = []
    for topology, atoms, cores in (
        ("chain", 512, (12, 16, 64)),
        ("cube", 512, (12, 16, 64)),
        ("chain", 1024, (12, 64)),
        ("cube", 1024, (12, 64)),
    ):
        for core in cores:
            rows.extend([
                row(topology, atoms, core, repetition,
                    searches[(topology, atoms, core)],
                    totals.get((topology, atoms, core), searches[(topology, atoms, core)] * 2),
                    rss, (topology, atoms, core, repetition) not in incorrect)
                for repetition in (1, 2, 3)
            ])
    return rows


class FrontierTest(unittest.TestCase):
    def test_selects_512_winner_and_passes_frontier_gate(self):
        searches = {
            ("chain", 512, 12): 8.0, ("chain", 512, 16): 9.0, ("chain", 512, 64): 12.0,
            ("cube", 512, 12): 8.5, ("cube", 512, 16): 9.5, ("cube", 512, 64): 13.0,
            ("chain", 1024, 12): 10.0, ("chain", 1024, 64): 15.0,
            ("cube", 1024, 12): 11.0, ("cube", 1024, 64): 16.0,
        }
        report = analyze(frontier_rows(searches))
        self.assertEqual(report["selected_core_size"], 12)
        self.assertEqual(report["qualification_gate"], "passed")
        self.assertEqual(len(report["comparisons"]), 4)

    def test_fails_when_large_frontier_improvement_is_below_gate(self):
        searches = {
            ("chain", 512, 12): 8.0, ("chain", 512, 16): 9.0, ("chain", 512, 64): 12.0,
            ("cube", 512, 12): 8.5, ("cube", 512, 16): 9.5, ("cube", 512, 64): 13.0,
            ("chain", 1024, 12): 14.0, ("chain", 1024, 64): 15.0,
            ("cube", 1024, 12): 11.0, ("cube", 1024, 64): 16.0,
        }
        report = analyze(frontier_rows(searches))
        self.assertEqual(report["selected_core_size"], 12)
        self.assertEqual(report["qualification_gate"], "failed")
        self.assertTrue(any("chain-1024" in reason for reason in report["gate_reasons"]))


if __name__ == "__main__":
    unittest.main()
