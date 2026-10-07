from joint_fixed_neighbor_numeric_reuse import analyze_case


def test_exact_rate_uses_numeric_factor_requests_and_preserves_roles():
    result = analyze_case({
        "topology": "chain",
        "atoms": 32,
        "fixed_neighbor": {
            "numeric_factor_requests": 10,
            "numeric_factor_exact_reuse_opportunities": 2,
            "numeric_factor_pattern_only_matches": 3,
            "numeric_factor_value_mismatches": 3,
            "numeric_factor_column_mismatches": 1,
            "numeric_factor_policy_mismatches": 1,
            "numeric_factor_pattern_mismatches": 1,
            "initial_profile_exact_reuse_opportunities": 1,
            "trial_profile_exact_reuse_opportunities": 1,
            "reference_exact_reuse_opportunities": 0,
            "accepted_endpoint_exact_reuse_opportunities": 0,
        },
    })
    assert result["numeric_factor_requests"] == 10
    assert result["numeric_factor_exact_reuse_opportunities"] == 2
    assert result["exact_reuse_rate"] == 0.2
    assert result["initial_profile_exact_reuse_opportunities"] == 1
    assert result["trial_profile_exact_reuse_opportunities"] == 1


if __name__ == "__main__":
    test_exact_rate_uses_numeric_factor_requests_and_preserves_roles()
