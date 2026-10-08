import joint_fixed_neighbor_operator_diagnostic as diagnostic


def _case(method, correction, objective, offset=0.0):
    return {
        "topology": "chain", "atoms": 4,
        "core_atoms": 2,
        "fixed_neighbor": {
            "core_atoms": 2, "block_order": "forward", "local_work_policy": "OneAcceptedLocalUpdate",
            "local_search_method": method, "local_preconditioner": "Schwarz",
            "search_reason": "block-stationary", "endpoint_certified": method == "LegacyCompact",
            "runtime_convergence": "Passed" if method == "LegacyCompact" else "Failed",
            "objective": objective, "final_eta": [0.1 + offset] * 4,
            "final_beta": [1.0 + offset] * 8, "final_ac_scaling_weights": [2.0] * 8,
            "endpoint_assessment": {
                "coefficient_difference": correction / 2,
                "primary_gradient_inf_norm": 1e-13,
                "reference_gradient_inf_norm": 2e-13,
                "correction": {"inf_norm": correction, "vector": [correction] * 4,
                               "max_coordinate": {"atom": 0, "value": correction}},
            },
            "sweep_telemetry": [
                {"sweep": 1, "objective_after": 1.0, "global_ac_kkt": 1e-3,
                 "global_width_gradient_inf_norm": 1e-5, "eta_change_inf": 1e-3,
                 "beta_scaled_change": 1e-3},
                {"sweep": 2, "objective_after": objective, "global_ac_kkt": 1e-13,
                 "global_width_gradient_inf_norm": 1e-13, "eta_change_inf": 1e-12,
                 "beta_scaled_change": 1e-12},
            ],
        },
        "endpoint_assessment_by_sweep": [
            {"sweep": 1, "search_stop_reason": "block-stationary", "assessment": {
                "correction": {"inf_norm": correction * 2, "vector": [correction * 2] * 4},
                "projected_width": {"rank": 4}, "corrected_jacobian": {"rank": 4},
            }},
            {"sweep": 2, "search_stop_reason": "block-stationary", "assessment": {
                "coefficient_difference": correction / 2,
                "primary_gradient_inf_norm": 1e-13, "reference_gradient_inf_norm": 2e-13,
                "correction": {"inf_norm": correction, "vector": [correction] * 4,
                               "max_coordinate": {"atom": 0, "value": correction}},
                "projected_width": {"rank": 4}, "corrected_jacobian": {"rank": 4},
            }},
        ],
    }


def test_analyze_case_preserves_trajectory_and_comparison():
    operator = _case("OperatorPcg", 1.4e-10, 1.0e-4)
    legacy = _case("LegacyCompact", 3.0e-13, 1.0e-4, offset=1e-12)
    report = diagnostic.analyze_case(operator, legacy)

    assert report["correction_diagnosis"]["classification"] == "still decreasing"
    assert len(report["correction_diagnosis"]["trajectory"]) == 2
    assert abs(report["comparison"]["eta_inf_difference"] - 1e-12) < 1e-15
    assert abs(report["comparison"]["correction_vector_inf_difference"] - 1.397e-10) < 1e-15
    assert report["operator"]["reference_coefficient_difference"] == 7e-11


def test_rank_change_takes_precedence_over_trajectory_shape():
    operator = _case("OperatorPcg", 1e-10, 1.0e-4)
    operator["endpoint_assessment_by_sweep"][0]["assessment"]["projected_width"]["rank"] = 3
    legacy = _case("LegacyCompact", 3e-13, 1.0e-4)
    assert diagnostic.analyze_case(operator, legacy)["correction_diagnosis"]["classification"] == \
        "rank/active-face driven"


if __name__ == "__main__":
    test_analyze_case_preserves_trajectory_and_comparison()
    test_rank_change_takes_precedence_over_trajectory_shape()
