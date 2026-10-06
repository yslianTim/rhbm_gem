from joint_fixed_neighbor_attribution import analyze, analyze_case


def _case():
    return {
        "topology": "chain", "atoms": 256,
        "fixed_neighbor": {
            "local_trajectory_telemetry": True, "sweeps": 1,
            "search_seconds": 3.0, "profile_factor_seconds": 2.0,
            "block_telemetry": [{
                "sweep": 1, "block": 1, "profile_evaluations": 4,
                "profile_trials": [
                    {"accepted": True, "accepted_update": 0, "objective_reduction": 0.0,
                     "factor_seconds": 0.2},
                    {"accepted": False, "accepted_update": None, "objective_reduction": 0.4,
                     "factor_seconds": 0.3},
                    {"accepted": True, "accepted_update": 1, "objective_reduction": 0.5,
                     "factor_seconds": 0.8},
                    {"accepted": True, "accepted_update": 2, "objective_reduction": 0.25,
                     "factor_seconds": 0.7},
                ],
            }],
        },
    }


def test_attribution_metrics_preserve_accepted_update_semantics():
    result = analyze_case(_case())
    assert result["trajectory_available"]
    assert result["accepted_updates_per_visit"]["mean"] == 2
    assert result["profile_evaluations_per_visit"]["mean"] == 4
    assert result["first_update_reduction_fraction"]["mean"] == 2 / 3
    assert result["first_two_update_reduction_fraction"]["mean"] == 1
    assert result["factor_time_fraction_after_first_update"]["mean"] == 0.7 / 2
    assert result["factor_time_fraction_after_second_update"]["mean"] == 0


def test_decision_is_explicit_and_not_hidden_in_thresholds():
    report = analyze([_case()], "yes")
    assert report["over_solving_decision"] == "yes"
    assert "no percentage gate" in report["decision_basis"]


if __name__ == "__main__":
    test_attribution_metrics_preserve_accepted_update_semantics()
    test_decision_is_explicit_and_not_hidden_in_thresholds()
