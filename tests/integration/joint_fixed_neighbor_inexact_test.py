from joint_fixed_neighbor_inexact import analyze, summarize


def _case(policy, endpoint=True, search=10.0):
    return summarize({
        "topology": "chain", "atoms": 256, "peak_rss_mb": 100,
        "fixed_neighbor": {
            "sweeps": 3, "confirmed_stationarity_sweep": 3,
            "block_solves": 6, "profile_evaluations": 12,
            "accepted_local_updates": 6, "profile_factor_seconds": 4.0,
            "search_seconds": search, "local_assessment_seconds": 2.0,
            "total_elapsed_seconds": search + 2.0, "objective": 1.0,
            "endpoint_certified": endpoint,
            "runtime_convergence": "Passed" if endpoint else "Failed",
            "search_converged": True,
        },
    }, policy)


def test_inexact_analyzer_requires_full_endpoint_correctness():
    report = analyze([_case("Full"), _case("OneAccepted", endpoint=False),
                      _case("TwoAccepted")])
    assert report["correctness_gate"] == "failed"
    assert report["cases"][0]["correct_policies"] == ["Full", "TwoAccepted"]


def test_inexact_analyzer_selects_correct_faster_candidate():
    rows = [_case("Full", search=20.0), _case("OneAccepted", search=8.0),
            _case("TwoAccepted", search=12.0)]
    report = analyze(rows)
    assert report["correctness_gate"] == "passed"
    assert report["selected_policy"] == "OneAccepted"


if __name__ == "__main__":
    test_inexact_analyzer_requires_full_endpoint_correctness()
    test_inexact_analyzer_selects_correct_faster_candidate()
