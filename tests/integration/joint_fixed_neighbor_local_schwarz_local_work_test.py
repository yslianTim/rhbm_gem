import joint_fixed_neighbor_local_schwarz_local_work as qualification


def _case(topology, policy, local_work, runtime="Passed", local="Passed"):
    return qualification.summarize({
        "topology": topology, "atoms": 256,
        "fixed_neighbor": {
            "outer_core_atoms": 64, "local_work_policy": local_work,
            "local_search_method": "OperatorPcg", "local_preconditioner": "Schwarz",
            "local_schwarz_core_atoms": 32, "local_schwarz_overlap_hops": 0,
            "local_schwarz_max_block_atoms": 64, "search_converged": True,
            "sweeps": 7, "search_seconds": 10.0, "assessment_seconds": 2.0,
            "total_elapsed_seconds": 12.0, "global_kkt": {"value": 1e-13},
            "width_gradient_inf_norm": 1e-13,
            "assessment_inner": {"status": "Passed"},
            "assessment_gradient": {"status": "Passed"},
            "assessment_local": {"status": local},
            "assessment_identified": {"status": "Passed"},
            "endpoint_assessment": {"endpoint_trust": {"passed": runtime == "Passed"}},
            "endpoint_certified": runtime == "Passed", "runtime_convergence": runtime,
            "sweep_telemetry": [], "block_telemetry": [],
            "local_operator_work": {"pcg_mean_iterations": 2.0, "pcg_max_iterations": 3},
        },
    }, policy)


def test_two_accepted_is_selected_before_full_when_all_pass():
    rows = [_case("chain", "OneAccepted", "OneAcceptedLocalUpdate", runtime="Failed", local="Failed"),
            _case("chain", "TwoAccepted", "TwoAcceptedLocalUpdates"),
            _case("chain", "Full", "FullLocalSearch"),
            _case("cube", "OneAccepted", "OneAcceptedLocalUpdate", runtime="Failed", local="Failed"),
            _case("cube", "TwoAccepted", "TwoAcceptedLocalUpdates"),
            _case("cube", "Full", "FullLocalSearch")]
    report = qualification.analyze(rows)
    assert report["correctness_gate"] == "passed"
    assert report["selected_policy"] == "TwoAccepted"


def test_failed_endpoint_blocks_local_work_selection():
    rows = [_case("chain", policy, local_work, runtime="Failed", local="Failed")
            for policy, local_work in (("OneAccepted", "OneAcceptedLocalUpdate"),
                                       ("TwoAccepted", "TwoAcceptedLocalUpdates"),
                                       ("Full", "FullLocalSearch"))]
    assert qualification.analyze(rows)["correctness_gate"] == "failed"


if __name__ == "__main__":
    test_two_accepted_is_selected_before_full_when_all_pass()
    test_failed_endpoint_blocks_local_work_selection()
