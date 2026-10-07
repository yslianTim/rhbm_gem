from joint_fixed_neighbor_workspace_timing import analyze


def _report(variant, repetition, search_seconds, symbolic, reuses):
    fixed = {
        "search_converged": True, "search_reason": "block-stationary", "sweeps": 2,
        "first_order_stationarity_sweep": 1, "confirmed_stationarity_sweep": 2,
        "block_solves": 4, "profile_evaluations": 8, "accepted_local_updates": 4,
        "objective": 1.0, "final_global_ac_kkt": 1e-12,
        "final_raw_width_gradient_inf_norm": 1e-14, "search_seconds": search_seconds,
        "total_elapsed_seconds": search_seconds + 0.1, "symbolic_factorizations": symbolic,
        "symbolic_reuses": reuses, "numeric_factorizations": 8,
        "matrix_preparation_seconds": 0.2, "symbolic_seconds": 0.3, "numeric_seconds": 0.5,
        "final_eta": [0.1, 0.2], "final_beta": [0.3, 0.4, 0.5, 0.6],
        "block_telemetry": [{"sweep": 1, "block": 1, "status": "accepted",
                             "accepted": True, "accepted_updates": 1}],
        "fixed_neighbor_work": {"candidate_replay_seconds": 0.1,
                                 "candidate_copy_seconds": 0.01,
                                 "full_candidate_replays": 3,
                                 "candidate_state_full_copies": 3},
    }
    return {
        "topology": "chain", "atoms": 1024, "variant": variant,
        "repetition": repetition, "warmup": False, "status": "completed",
        "result": {"workspace_mode": "persistent" if variant == "persistent-workspace"
                    else "fresh-per-block-visit", "measurement_scope": "fixed-neighbor-search-only",
                    "fixed_neighbor": fixed, "peak_rss_mb": 100.0},
    }


def test_matched_workspace_analysis_checks_trajectory_and_reports_phases():
    reports = []
    for repetition in range(3):
        reports.append(_report("fresh-workspace", repetition, 10.0, 8, 0))
        reports.append(_report("persistent-workspace", repetition, 8.0, 2, 6))
    report = analyze(reports, [("chain", 1024)], warmup=1, measurements=3)
    case = report["cases"][0]
    assert report["qualification_gate"] == "passed"
    assert case["numerical_gate"] == "passed"
    assert case["performance_class"] == "material"
    assert case["treatment"]["candidate_replay_copy_fraction_median"] == 0.11 / 8.0
    assert all(comparison["passed"] for comparison in case["comparisons"])


def test_matched_workspace_analysis_rejects_endpoint_mismatch():
    control = _report("fresh-workspace", 0, 10.0, 8, 0)
    treatment = _report("persistent-workspace", 0, 8.0, 2, 6)
    treatment["result"]["fixed_neighbor"]["final_eta"][0] = 0.3
    report = analyze([control, treatment], [("chain", 1024)], warmup=1, measurements=1)
    assert report["qualification_gate"] == "failed"
    assert "eta-endpoint-mismatch" in report["cases"][0]["comparisons"][0]["reasons"]


if __name__ == "__main__":
    test_matched_workspace_analysis_checks_trajectory_and_reports_phases()
    test_matched_workspace_analysis_rejects_endpoint_mismatch()
