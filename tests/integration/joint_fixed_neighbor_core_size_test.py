import unittest

from joint_fixed_neighbor_core_size import _parse_core_sizes, analyze, summarize


def row(topology, atoms, core, search, converged=True, kkt=1e-12, width=1e-14):
    return {"topology": topology, "atoms": atoms, "core_atoms": core,
            "policy": "OneAccepted", "status": "completed", "search_converged": converged,
            "search_reason": "block-stationary", "maximum_cache_replay_error": 1e-15,
            "maximum_objective_replay_error": 1e-15, "cache_replay_within_limit": True,
            "objective_replay_within_limit": True, "objective": 1.0,
            "final_global_ac_kkt": kkt, "final_raw_width_gradient_inf_norm": width,
            "search_seconds": search, "sweeps": 7, "confirmed_stationarity_sweep": 7}


class CoreSizeTest(unittest.TestCase):
    def test_selects_fastest_core_when_all_points_are_correct(self):
        rows = [row("chain", 512, 64, 30.), row("chain", 512, 128, 20.), row("chain", 512, 256, 25.),
                row("cube", 512, 64, 40.), row("cube", 512, 128, 35.), row("cube", 512, 256, 45.)]
        rows += [row("cube", 1024, 64, 80.), row("cube", 1024, 128, 70.), row("cube", 1024, 256, 90.)]
        report = analyze(rows)
        self.assertEqual(report["correctness_gate"], "passed")
        self.assertEqual(report["selected_core_size"], 128)

    def test_failed_core_is_not_selected(self):
        rows = [row("chain", 512, 64, 30.), row("chain", 512, 128, 20.),
                row("chain", 512, 256, 25., converged=False)]
        report = analyze(rows)
        self.assertEqual(report["correctness_gate"], "failed")
        self.assertIsNone(report["selected_core_size"])

    def test_accepts_explicit_core_size_list(self):
        rows = [row("chain", 512, core, float(core)) for core in (8, 16, 32)]
        report = analyze(rows, core_sizes=(8, 16, 32))
        self.assertEqual(report["correctness_gate"], "passed")
        self.assertEqual(report["selected_core_size"], 8)
        self.assertEqual(_parse_core_sizes("8, 16,32"), (8, 16, 32))

    def test_summarize_reports_tail_aware_geometry_and_normalized_costs(self):
        result = {"topology": "chain", "atoms": 10, "core_atoms": 6,
                  "peak_rss_mb": 12.0, "fixed_neighbor": {
                      "outer_core_atoms": 6, "sweeps": 2, "search_seconds": 10.0,
                      "search_converged": True, "search_reason": "block-stationary",
                      "objective": 1.0, "block_solves": 4, "profile_evaluations": 8,
                      "accepted_local_updates": 3, "final_global_ac_kkt": 1e-12,
                      "final_raw_width_gradient_inf_norm": 1e-14,
                      "sweep_telemetry": [{"cache_replay_error": 1e-15,
                                           "objective_replay_error": 1e-15}],
                      "block_telemetry": [{"block": 1, "atoms": 6},
                                           {"block": 2, "atoms": 4}],
                      "fixed_neighbor_work": {
                          "local_search_seconds": 6.0, "candidate_replay_seconds": 2.0,
                          "sweep_replay_seconds": 1.0, "sweep_global_state_seconds": 1.0,
                          "local_profile_work": {"total": {
                              "derivative_prepare_seconds": 0.5,
                              "derivative_reduce_seconds": 0.25}}}}}
        summary = summarize({"result": result, "status": "completed",
                             "topology": "chain", "atoms": 10, "core_atoms": 6})
        self.assertEqual(summary["realized_outer_block_count"], 2)
        self.assertEqual(summary["realized_outer_core_atoms_minimum"], 4)
        self.assertEqual(summary["realized_outer_core_atoms_mean"], 5.0)
        self.assertEqual(summary["realized_outer_core_atoms_maximum"], 6)
        self.assertEqual(summary["local_search_fraction"], 0.6)
        self.assertEqual(summary["candidate_replay_seconds_per_block_solve"], 0.5)
        self.assertEqual(summary["derivative_prepare_seconds"], 0.5)


if __name__ == "__main__":
    unittest.main()
