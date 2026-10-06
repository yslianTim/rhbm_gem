import unittest

from joint_fixed_neighbor_inexact_qualification import analyze, summarize


def row(topology, policy, local_work, search, converged=True, kkt=1e-12, width=1e-14):
    return {"topology": topology, "atoms": 1024, "policy": policy, "status": "completed",
            "search_converged": converged, "final_global_ac_kkt": kkt,
            "final_raw_width_gradient_inf_norm": width, "local_work_seconds": local_work,
            "search_seconds": search, "sweeps": 7, "confirmed_stationarity_sweep": 7,
            "block_solves": 56, "profile_evaluations": 100, "accepted_local_updates": 20,
            "total_seconds": search, "peak_rss_mb": 100.0}


class QualificationTest(unittest.TestCase):
    def test_positive_reduction_passes(self):
        report = analyze([row("chain", "Full", 100.0, 110.0),
                          row("chain", "OneAccepted", 40.0, 50.0)])
        self.assertEqual(report["qualification_gate"], "passed")
        self.assertEqual(report["selected_policy"], "OneAccepted")
        self.assertGreater(report["cases"][0]["candidate"]["search_speedup"], 2.0)

    def test_nonconverged_candidate_fails(self):
        report = analyze([row("cube", "Full", 100.0, 110.0),
                          row("cube", "OneAccepted", 40.0, 50.0, converged=False)])
        self.assertEqual(report["qualification_gate"], "failed")
        self.assertIsNone(report["selected_policy"])

    def test_scaling_summary_preserves_search_fields(self):
        result = {"topology": "chain", "atoms": 1024,
                  "measurement_scope": "fixed-neighbor-search-only",
                  "fixed_neighbor": {"search_converged": True, "search_reason": "block-stationary",
                                     "sweeps": 7, "confirmed_stationarity_sweep": 7,
                                     "block_solves": 56, "profile_evaluations": 100,
                                     "accepted_local_updates": 20, "local_factor_seconds": 40.0,
                                     "profile_factor_seconds": 3.0, "search_seconds": 50.0,
                                     "total_elapsed_seconds": 51.0,
                                     "final_global_ac_kkt": 1e-12,
                                     "final_raw_width_gradient_inf_norm": 1e-14},
                  "peak_rss_mb": 123.0}
        summary = summarize({"topology": "chain", "atoms": 1024,
                             "policy": "OneAccepted", "status": "completed", "result": result})
        self.assertEqual(summary["local_work_seconds"], 40.0)
        self.assertEqual(summary["profile_factor_seconds"], 3.0)
        self.assertEqual(summary["search_seconds"], 50.0)


if __name__ == "__main__":
    unittest.main()
