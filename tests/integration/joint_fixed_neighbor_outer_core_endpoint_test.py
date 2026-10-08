import unittest

from joint_fixed_neighbor_outer_core_endpoint import analyze


def row(topology, atoms, core, eta_offset=0.0, passed=True):
    return {
        "topology": topology, "atoms": atoms, "core_atoms": core,
        "status": "completed", "search_converged": True,
        "search_reason": "block-stationary", "sweeps": 8,
        "confirmed_stationarity_sweep": 8, "global_ac_kkt": 1e-12,
        "width_gradient_inf_norm": 1e-14,
        "cache_replay_within_limit": True, "objective_replay_within_limit": True,
        "endpoint_certified": passed, "runtime_convergence": "Passed" if passed else "Failed",
        "assessment_inner": "Passed" if passed else "Failed",
        "assessment_gradient": "Passed" if passed else "Failed",
        "assessment_local": "Passed" if passed else "Failed",
        "assessment_identified": "Passed" if passed else "Failed",
        "endpoint_trust": passed, "objective": 1.0,
        "final_eta": [0.1 + eta_offset, 0.2 + eta_offset],
        "final_beta": [1.0 + eta_offset, 2.0 + eta_offset],
        "final_ac_scaling_weights": [1.0, 1.0],
        "rank_evidence": {"projected_width_rank": 2, "corrected_jacobian_rank": 2,
                          "normalized_width_rank": 2},
    }


class EndpointTest(unittest.TestCase):
    def test_requires_every_core_on_every_case_and_reports_parity(self):
        rows = []
        for topology in ("chain", "cube"):
            rows.extend(row(topology, 256, core, eta_offset=0.0 if core == 64 else 1e-12)
                        for core in (12, 16, 64))
        report = analyze(rows)
        self.assertEqual(report["qualification_gate"], "passed")
        self.assertEqual(report["qualified_core_sizes"], [12, 16])
        self.assertTrue(report["cases"][0]["parity_vs_core64"]["12"]["parameter_parity"])

    def test_endpoint_failure_blocks_qualification(self):
        rows = [row("chain", 256, 12), row("chain", 256, 16), row("chain", 256, 64, passed=False)]
        report = analyze(rows)
        self.assertEqual(report["qualification_gate"], "failed")
        self.assertEqual(report["qualified_core_sizes"], [])


if __name__ == "__main__":
    unittest.main()
