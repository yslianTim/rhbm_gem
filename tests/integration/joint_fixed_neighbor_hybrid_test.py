import joint_fixed_neighbor_hybrid as hybrid


def _case(topology, policy, converged=True, runtime="Passed", local="Passed", eta=1e-12):
    fixed = {
        "outer_core_atoms": 64, "local_schwarz_core_atoms": 32,
        "local_schwarz_overlap_hops": 0, "local_schwarz_max_block_atoms": 64,
        "search_converged": converged, "search_reason": "block-stationary",
        "search_seconds": 10.0, "total_elapsed_seconds": 12.0,
        "global_kkt": {"value": 1e-13}, "width_gradient_inf_norm": 1e-13,
        "assessment_inner": {"status": "Passed"}, "assessment_gradient": {"status": "Passed"},
        "assessment_local": {"status": local}, "assessment_identified": {"status": "Passed"},
        "endpoint_assessment": {"local_correction_inf_norm": 1e-12,
                                 "endpoint_trust": {"passed": runtime == "Passed"}},
        "endpoint_certified": runtime == "Passed", "runtime_convergence": runtime,
        "sweeps": [{"eta_change_inf": eta, "beta_scaled_change": 1e-12}],
    }
    if policy.startswith("+"):
        fixed.update({"operator_phase_seconds": 7.0, "legacy_polish_seconds": 5.0,
                      "total_search_seconds": 10.0, "polish_sweeps": 1})
    return hybrid.summarize({"topology": topology, "atoms": 256, "fixed_neighbor": fixed}, policy)


def test_two_polish_is_selected_when_one_lacks_confirmation():
    rows = [_case("chain", "Operator", converged=False, runtime="Failed", local="Failed"),
            _case("chain", "+1 Legacy polish", converged=False),
            _case("chain", "+2 Legacy polish"), _case("chain", "Legacy control"),
            _case("cube", "Operator", converged=False, runtime="Failed", local="Failed"),
            _case("cube", "+1 Legacy polish", converged=False),
            _case("cube", "+2 Legacy polish"), _case("cube", "Legacy control")]
    report = hybrid.analyze(rows)
    assert report["correctness_gate"] == "passed"
    assert report["selected_policy"] == "+2 Legacy polish"


def test_runtime_pass_without_search_confirmation_is_not_enough():
    rows = [_case("chain", policy, converged=False, runtime="Passed", local="Passed")
            for policy in ("+1 Legacy polish",)]
    assert hybrid.analyze(rows)["correctness_gate"] == "failed"


if __name__ == "__main__":
    test_two_polish_is_selected_when_one_lacks_confirmation()
    test_runtime_pass_without_search_confirmation_is_not_enough()
