import unittest

from joint_fixed_neighbor_outer_core_qualification import analyze, build_parser


def screen_row(topology, atoms, core, search, converged=True):
    return {
        "topology": topology, "atoms": atoms, "core_atoms": core,
        "policy": "OneAccepted", "status": "completed",
        "search_converged": converged,
        "search_reason": "block-stationary" if converged else "maximum-sweeps",
        "maximum_cache_replay_error": 1e-15,
        "maximum_objective_replay_error": 1e-15,
        "cache_replay_within_limit": converged,
        "objective_replay_within_limit": converged,
        "final_global_ac_kkt": 1e-12 if converged else None,
        "final_raw_width_gradient_inf_norm": 1e-14 if converged else None,
        "search_seconds": search,
    }


def endpoint_row(topology, atoms, core, passed=True):
    return {
        "topology": topology, "atoms": atoms, "core_atoms": core,
        "status": "completed", "search_converged": True,
        "search_reason": "block-stationary", "sweeps": 8,
        "confirmed_stationarity_sweep": 8,
        "global_ac_kkt": 1e-12, "width_gradient_inf_norm": 1e-14,
        "cache_replay_within_limit": True, "objective_replay_within_limit": True,
        "endpoint_certified": passed, "runtime_convergence": "Passed" if passed else "Failed",
        "assessment_inner": "Passed" if passed else "Failed",
        "assessment_gradient": "Passed" if passed else "Failed",
        "assessment_local": "Passed" if passed else "Failed",
        "assessment_identified": "Passed" if passed else "Failed",
        "endpoint_trust": passed, "objective": 1.0,
        "final_eta": [0.1, 0.2], "final_beta": [1.0, 2.0],
        "final_ac_scaling_weights": [1.0, 1.0],
        "rank_evidence": {"projected_width_rank": 2,
                          "corrected_jacobian_rank": 2,
                          "normalized_width_rank": 2},
    }


def frontier_row(topology, atoms, core, repetition, search, correct=True):
    return {
        "topology": topology, "atoms": atoms, "core_atoms": core,
        "warmup": False, "repetition": repetition, "status": "completed",
        "search_converged": correct,
        "search_reason": "block-stationary" if correct else "maximum-sweeps",
        "final_global_ac_kkt": 1e-12 if correct else None,
        "final_raw_width_gradient_inf_norm": 1e-14 if correct else None,
        "cache_replay_within_limit": correct,
        "objective_replay_within_limit": correct,
        "search_seconds": search, "total_seconds": search * 2,
        "peak_rss_mb": 100.0, "sweeps": 18, "block_solves": 100,
        "local_work_seconds": search * 0.7,
        "maximum_local_rows": core * 2, "maximum_local_columns": core,
    }


class OuterCoreQualificationTest(unittest.TestCase):
    def test_one_parser_dispatches_all_phases_and_parses_generic_core_sizes(self):
        args = build_parser().parse_args([
            "--phase", "screen", "--build-dir", "build", "--output-dir", "out",
            "--core-sizes", "64, 128,256",
        ])
        self.assertEqual(args.core_sizes, (64, 128, 256))
        self.assertEqual(analyze("screen", [screen_row("chain", 512, 64, 3.0)],
                                 core_sizes=(64,)),
                         analyze("screen", [screen_row("chain", 512, 64, 3.0)],
                                 core_sizes=(64,)))

    def test_screen_gate_selects_fastest_correct_core(self):
        rows = [screen_row("chain", 512, core, search)
                for core, search in ((64, 30.0), (128, 20.0), (256, 25.0))]
        rows += [screen_row("cube", 512, core, search)
                 for core, search in ((64, 40.0), (128, 35.0), (256, 45.0))]
        report = analyze("screen", rows, core_sizes=(64, 128, 256))
        self.assertEqual(report["correctness_gate"], "passed")
        self.assertEqual(report["selected_core_size"], 128)

        broken = list(rows)
        broken[-1] = screen_row("cube", 512, 256, 45.0, converged=False)
        self.assertEqual(analyze("screen", broken, core_sizes=(64, 128, 256))[
            "correctness_gate"], "failed")

    def test_endpoint_gate_requires_all_finalists_and_control(self):
        rows = [endpoint_row(topology, 256, core)
                for topology in ("chain", "cube") for core in (12, 16, 64)]
        report = analyze("endpoint", rows, core_sizes=(12, 16, 64))
        self.assertEqual(report["qualification_gate"], "passed")
        self.assertEqual(report["qualified_core_sizes"], [12, 16])

        rows[-1] = endpoint_row("cube", 256, 64, passed=False)
        self.assertEqual(analyze("endpoint", rows, core_sizes=(12, 16, 64))[
            "qualification_gate"], "failed")

    def test_frontier_gate_selects_512_winner_and_checks_large_cases(self):
        searches = {
            ("chain", 512, 12): 8.0, ("chain", 512, 16): 9.0, ("chain", 512, 64): 12.0,
            ("cube", 512, 12): 8.5, ("cube", 512, 16): 9.5, ("cube", 512, 64): 13.0,
            ("chain", 1024, 12): 10.0, ("chain", 1024, 64): 15.0,
            ("cube", 1024, 12): 11.0, ("cube", 1024, 64): 16.0,
        }
        rows = [frontier_row(topology, atoms, core, repetition,
                             searches[(topology, atoms, core)])
                for topology, atoms, cores in (
                    ("chain", 512, (12, 16, 64)), ("cube", 512, (12, 16, 64)),
                    ("chain", 1024, (12, 64)), ("cube", 1024, (12, 64)))
                for core in cores for repetition in (1, 2, 3)]
        report = analyze("frontier", rows, finalist_core_sizes=(12, 16),
                         control_core_size=64)
        self.assertEqual(report["selected_core_size"], 12)
        self.assertEqual(report["qualification_gate"], "passed")
        self.assertEqual(len(report["comparisons"]), 4)


if __name__ == "__main__":
    unittest.main()
