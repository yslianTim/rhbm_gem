import joint_fixed_neighbor_operator_route as route


def _row(topology, atoms, converged=True):
    return {
        "topology": topology, "atoms": atoms, "policy": route.SELECTED_POLICY,
        "status": "completed", "search_converged": converged,
        "global_ac_kkt": 1e-13, "global_width_gradient_inf_norm": 1e-13,
        "assessment_inner": "Passed", "assessment_gradient": "Passed",
        "assessment_local": "Passed", "assessment_identified": "Passed",
        "endpoint_trust": True, "endpoint_certified": True,
        "runtime_convergence": "Passed", "eta_change_inf": 1e-12,
    }


def test_route_requires_every_256_and_512_case_to_pass():
    rows = [_row(topology, atoms) for topology, atoms in route.EXPECTED]
    assert route.analyze(rows)["qualification"] == "passed"


def test_confirmation_failure_blocks_route():
    rows = [_row(topology, atoms, converged=not (topology == "cube" and atoms == 512))
            for topology, atoms in route.EXPECTED]
    report = route.analyze(rows)
    assert report["qualification"] == "failed"
    assert report["promotion_blocked"] is True


if __name__ == "__main__":
    test_route_requires_every_256_and_512_case_to_pass()
    test_confirmation_failure_blocks_route()
