from joint_fixed_neighbor_order import analyze


def make_case(topology="chain", atoms=256, order="forward", eta=0.0, beta=0.0,
        endpoint=True, runtime="Passed", converged=True):
    fixed = {
        "block_order": order,
        "search_converged": converged,
        "search_reason": "block-stationary" if converged else "block-stationarity-unconfirmed",
        "sweeps": 7,
        "first_order_stationarity_sweep": 6,
        "confirmed_stationarity_sweep": 7 if converged else None,
        "objective": 1e-10,
        "global_kkt": {"value": 1e-13},
        "width_gradient_inf_norm": 1e-14,
        "endpoint_certified": endpoint if atoms <= 512 else None,
        "runtime_convergence": runtime if atoms <= 512 else None,
        "final_eta": [eta, eta],
        "final_beta": [beta, beta],
        "final_ac_scaling_weights": [1.0, 0.5],
    }
    return {"topology": topology, "atoms": atoms, "fixed_neighbor": fixed}


def test_compares_only_confirmed_endpoints_and_allows_search_only_scope():
    cases = [make_case("chain", 256, order) for order in ("forward", "reverse")]
    report = analyze(cases, {("chain", 256)})
    row = report["cases"][0]
    assert report["order_sensitivity_gate"] == "passed"
    assert row["comparison_status"] == "compared"
    assert row["measurement_scope"] == "endpoint-assessed"
    assert row["endpoint_certification_parity"] is True
    assert row["runtime_convergence_parity"] is True
    assert row["endpoint_and_runtime_passed_both"] is True
    search_only = [make_case("cube", 1024, order) for order in ("forward", "reverse")]
    for case in search_only:
        fixed = case["fixed_neighbor"]
        case["block_order"] = fixed.pop("block_order")
        fixed["final_global_ac_kkt"] = fixed.pop("global_kkt")["value"]
        fixed["final_raw_width_gradient_inf_norm"] = fixed.pop("width_gradient_inf_norm")
    row = analyze(search_only, {("cube", 1024)})["cases"][0]
    assert row["measurement_scope"] == "fixed-neighbor-search-only"
    assert row["endpoint_certification_parity"] is None
    assert row["endpoint_and_runtime_passed_both"] is None


def test_does_not_compare_unconfirmed_searches_or_accept_material_divergence():
    asymmetric = [make_case("chain", 256, "forward", converged=True),
        make_case("chain", 256, "reverse", converged=False)]
    report = analyze(asymmetric, {("chain", 256)})
    assert report["cases"][0]["comparison_status"] == "convergence-asymmetry"
    assert report["order_sensitivity_gate"] == "failed"
    unconfirmed = [make_case("chain", 256, order, converged=False) for order in ("forward", "reverse")]
    report = analyze(unconfirmed, {("chain", 256)})
    assert report["cases"][0]["comparison_status"] == "unconfirmed"
    assert report["order_sensitivity_gate"] == "failed"
    divergent = [make_case("cube", 1024, "forward"), make_case("cube", 1024, "reverse", eta=2e-10)]
    report = analyze(divergent, {("cube", 1024)})
    assert report["cases"][0]["eta_inf_difference"] == 2e-10
    assert report["order_sensitivity_gate"] == "failed"


def test_assessed_sizes_require_endpoint_and_runtime_success():
    cases = [make_case("cube", 512, "forward"), make_case("cube", 512, "reverse", endpoint=False)]
    report = analyze(cases, {("cube", 512)})
    assert report["cases"][0]["endpoint_and_runtime_passed_both"] is False
    assert report["order_sensitivity_gate"] == "failed"


if __name__ == "__main__":
    test_compares_only_confirmed_endpoints_and_allows_search_only_scope()
    test_does_not_compare_unconfirmed_searches_or_accept_material_divergence()
    test_assessed_sizes_require_endpoint_and_runtime_success()
