from joint_fixed_neighbor_stationarity import analyze


def make_case(topology="chain", atoms=32, confirmed_sweep=3, certified_sweep=3):
    sweeps = []
    assessments = []
    for sweep in range(1, 5):
        sweeps.append({"sweep": sweep, "global_ac_kkt": 1e-12 if sweep >= 1 else 1.0,
            "global_width_gradient_inf_norm": 1e-13, "eta_change_inf": 1e-11 if sweep >= confirmed_sweep else 1e-8,
            "beta_scaled_change": 1e-11 if sweep >= confirmed_sweep else 1e-8,
            "coordinate_confirmation_available": sweep > 1})
        assessments.append({"sweep": sweep, "inner": sweep >= certified_sweep,
            "gradient": sweep >= certified_sweep, "local": sweep >= certified_sweep,
            "identified": sweep >= certified_sweep, "endpoint_trust": sweep >= certified_sweep,
            "runtime_convergence": "Passed" if sweep >= certified_sweep else "Failed"})
    return {"topology": topology, "atoms": atoms,
        "fixed_neighbor": {"sweeps": 4, "runtime_convergence": "Passed", "sweep_telemetry": sweeps},
        "endpoint_assessment_by_sweep": assessments}


def test_analyzer_classification_is_deterministic():
    case = make_case()
    first = analyze([case], {("chain", 32)})
    second = analyze([case], {("chain", 32)})
    assert first == second
    assert first["qualification"] == "passed"
    assert first["cases"][0]["cheap_stationarity_sweep"] == 1
    assert first["cases"][0]["first_certified_sweep"] == 3
    assert first["candidate_rule_safety"]["eta_only"]["safe_across_cases"]


def test_confirmation_requires_a_previous_complete_sweep_and_cannot_precede_certification():
    early = make_case(confirmed_sweep=2, certified_sweep=3)
    late = make_case("cube", 32, confirmed_sweep=4, certified_sweep=3)
    report = analyze([early, late], {("chain", 32), ("cube", 32)})
    chain = next(row for row in report["cases"] if row["topology"] == "chain")
    cube = next(row for row in report["cases"] if row["topology"] == "cube")
    assert chain["eta_only_confirmation_sweep"] == 2
    assert chain["candidate_results"]["eta_only"]["confirmation_precedes_certification"]
    assert cube["candidate_results"]["eta_only"]["confirmation_delay"] == 1
    assert not report["candidate_rule_safety"]["eta_only"]["safe_across_cases"]


if __name__ == "__main__":
    test_analyzer_classification_is_deterministic()
    test_confirmation_requires_a_previous_complete_sweep_and_cannot_precede_certification()
