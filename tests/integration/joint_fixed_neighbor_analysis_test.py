from joint_fixed_neighbor_analysis import analyze


def result():
    return {"runtime_convergence": "Passed", "objective": 1e-10,
        "global_kkt": {"status": "Passed", "value": 1e-12},
        "width_gradient_inf_norm": 1e-13, "search_seconds": 1.0,
        "assessment_inner": {"status": "Passed"}, "assessment_gradient": {"status": "Passed"},
        "assessment_local": {"status": "Passed"}, "assessment_identified": {"status": "Passed"}}


def case():
    fixed = result() | {"search_converged": True, "endpoint_certified": True, "sweeps": 3}
    return {"topology": "chain", "atoms": 32, "global_legacy_compact": result(),
        "global_operator_pcg": result(), "fixed_neighbor": fixed, "peak_rss_mb": 100.0,
        "comparisons": [{"method": "FixedNeighbor", "objective_difference": 1e-15,
            "scaled_ac_inf_difference": 1e-9}]}


def test_passes_only_for_certified_parity_case():
    expected = {("chain", 32)}
    report = analyze([case()], expected)
    assert report["fixed_neighbor_f2_gate"] == "passed"
    broken = case()
    broken["fixed_neighbor"]["endpoint_certified"] = False
    assert analyze([broken], expected)["fixed_neighbor_f2_gate"] == "failed"


def test_missing_expected_case_fails_gate():
    report = analyze([], {("cube", 256)})
    assert report["fixed_neighbor_f2_gate"] == "failed"
    assert report["missing_cases"] == [{"topology": "cube", "atoms": 256}]


if __name__ == "__main__":
    test_passes_only_for_certified_parity_case()
    test_missing_expected_case_fails_gate()
