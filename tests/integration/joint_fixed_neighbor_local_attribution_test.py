from joint_fixed_neighbor_local_attribution import analyze_case


def test_phase_attribution_reports_both_denominators_and_roles():
    result = analyze_case({
        "topology": "chain", "atoms": 32,
        "fixed_neighbor": {
            "search_seconds": 10.0,
            "local_profile_work": {
                "total_seconds": 8.0,
                "total": {"evaluations": 6, "profile_basis_seconds": 2.0,
                           "linear_numeric_seconds": 3.0, "replay_trust_seconds": 1.0},
                "initial_profile": {"evaluations": 2},
                "trial_profile": {"evaluations": 4},
                "accepted_endpoint": {"evaluations": 2},
                "reference_evaluation": {"evaluations": 1},
            },
        },
    })
    assert result["profile_evaluations"] == 6
    assert result["initial_profile_evaluations"] == 2
    assert result["trial_profile_evaluations"] == 4
    phases = {row["phase"]: row for row in result["phases"]}
    assert phases["Basis construction"]["percent_of_local_search"] == 25.0
    assert phases["Numeric QR"]["percent_of_total_search"] == 30.0


if __name__ == "__main__":
    test_phase_attribution_reports_both_denominators_and_roles()
