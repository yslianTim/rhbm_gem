from joint_fixed_neighbor_local_certification import compare


def baseline(passed, state=0.0):
    endpoint = {"assessment": {"inner": passed, "gradient": passed, "local": passed,
        "identified": passed, "endpoint_trust": {"passed": passed}}}
    return {"topology": "cube", "atoms": 256,
        "fixed_neighbor": {"first_order_stationarity_sweep": 11, "search_converged": True,
            "sweep_telemetry": [{"wall_seconds": 1.0, "objective_before": 1.0, "objective_after": 0.5}],
            "block_telemetry": [{"sweep": 11, "profile_evaluations": 2}]},
        "endpoint_decomposition": [{"sweep": 11, "raw": endpoint,
            "eta": [state], "beta": [state]}], "peak_rss_mb": 100.0}


def certified(passed, state=0.0, rejected=False):
    return {"topology": "cube", "atoms": 256, "peak_rss_mb": 110.0,
        "fixed_neighbor": {"sweeps": 12, "local_assessment_count": 1,
            "local_assessment_seconds": 0.25, "search_seconds": 1.25,
            "endpoint_certified": passed, "runtime_convergence": "Passed" if passed else "Failed",
            "final_eta": [state], "final_beta": [state],
            "maximum_local_assessment_rows": 500, "maximum_local_assessment_columns": 256,
            "block_telemetry": [{"sweep": 12, "block": 1, "profile_evaluations": 4, "objective_before": 1.0,
                "objective_after": 0.5, "local_assessment_attempted": True,
                "local_assessment_passed": not rejected, "accepted": not rejected,
                "local_assessment_failure": "b-not-stationary" if rejected else "none",
                "local_trust_reason": "trusted", "local_correction_inf_norm": 2e-10,
                "local_profile_gradient_inf_norm": 1e-13}]}}


def test_comparison_reports_improvement_rejections_and_resource_bounds():
    result = compare(baseline(False), certified(True, 1e-6, rejected=True))
    candidate = result["CertifiedLocal"]
    assert result["certified_local_improves_global_endpoint"]
    assert result["local_endpoint_differs"]
    assert candidate["rejected_local_candidates"] == 1
    assert candidate["objective_monotone"]
    assert candidate["maximum_local_assessment_columns"] == 256


def test_added_cost_without_global_endpoint_improvement_is_visible():
    result = compare(baseline(False), certified(False))
    assert not result["certified_local_improves_global_endpoint"]
    assert result["CertifiedLocal"]["local_assessment_count"] == 1
    assert result["CertifiedLocal"]["local_assessment_seconds"] > 0


if __name__ == "__main__":
    test_comparison_reports_improvement_rejections_and_resource_bounds()
    test_added_cost_without_global_endpoint_improvement_is_visible()
