from joint_fixed_neighbor_endpoint_decomposition import classify


def state(inner=True, gradient=True, local=True, identified=True, trust=True,
          correction=1e-12, width_gradient=1e-13, singular=1.0, threshold=1e-8):
    assessment = {"inner": inner, "gradient": gradient, "local": local,
        "identified": identified, "local_correction_inf_norm": correction,
        "endpoint_trust": {"passed": trust, "reason": "trusted"},
        "projected_width": {"minimum_singular_value": singular, "rank_threshold": threshold},
        "corrected_jacobian": {"minimum_singular_value": singular, "rank_threshold": threshold},
        "normalized_width": {"minimum_singular_value": singular, "rank_threshold": threshold}}
    return {"assessment": assessment, "width_gradient_inf_norm": width_gradient}


def sample(raw, primary, reference, block_primary=1e-9, primary_reference=1e-12):
    return {"sweep": 11, "raw": raw, "same_eta_primary": primary, "same_eta_reference": reference,
        "coefficient_difference_block_primary": block_primary,
        "coefficient_difference_primary_reference": primary_reference}


def test_same_eta_profiles_passing_classifies_ac_profile_gap():
    result = classify([sample(state(inner=False, local=False, correction=1e-9), state(), state())])
    assert result["classification"] == "ac-profile-manifold-gap"


def test_profiled_width_failure_with_material_beta_gap_classifies_mixed():
    result = classify([sample(state(inner=False), state(gradient=False, local=False,
        width_gradient=2e-12), state(gradient=False, local=False, width_gradient=2e-12))])
    assert result["classification"] == "mixed"


def test_profiled_width_failure_without_beta_gap_classifies_width_gap():
    result = classify([sample(state(inner=False), state(gradient=False, local=False,
        width_gradient=2e-12), state(gradient=False, local=False, width_gradient=2e-12),
        block_primary=1e-12)])
    assert result["classification"] == "width-stationarity-gap"


def test_only_near_rank_decision_boundary_classifies_rank_related():
    result = classify([sample(state(identified=False, singular=1.5e-8, threshold=1e-8),
        state(identified=False, singular=1.5e-8, threshold=1e-8),
        state(identified=False, singular=1.5e-8, threshold=1e-8), block_primary=1e-12)])
    assert result["classification"] == "rank-boundary-related"


def test_missing_positive_evidence_is_unidentified():
    assert classify([])["classification"] == "unidentified"


if __name__ == "__main__":
    test_same_eta_profiles_passing_classifies_ac_profile_gap()
    test_profiled_width_failure_with_material_beta_gap_classifies_mixed()
    test_profiled_width_failure_without_beta_gap_classifies_width_gap()
    test_only_near_rank_decision_boundary_classifies_rank_related()
    test_missing_positive_evidence_is_unidentified()
