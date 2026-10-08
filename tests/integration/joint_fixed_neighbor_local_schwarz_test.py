from joint_fixed_neighbor_local_schwarz import analyze


def row(topology, candidate, preconditioner="Schwarz", search=10.0, correct=True):
    return {
        "topology": topology, "atoms": 256, "candidate": candidate,
        "preconditioner": preconditioner, "status": "completed",
        "search_converged": correct, "final_global_ac_kkt": 1e-12 if correct else 2e-10,
        "final_raw_width_gradient_inf_norm": 1e-14 if correct else 2e-12,
        "maximum_cache_replay_error": 1e-15, "maximum_objective_replay_error": 1e-15,
        "search_seconds": search,
    }


def test_retains_fastest_three_common_correct_schwarz_candidates():
    names = ["Schwarz-128/1", "Schwarz-64/0", "Schwarz-64/1", "Schwarz-32/0",
             "Schwarz-32/1", "Schwarz-16/0", "Schwarz-16/1", "Schwarz-8/0", "Schwarz-8/1"]
    rows = []
    for topology in ("chain", "cube"):
        rows.append(row(topology, "Identity", "Identity", 30.0))
        rows.append(row(topology, "Diagonal", "Diagonal", 20.0))
        for index, name in enumerate(names):
            rows.append(row(topology, name, search=float(index + 1)))
    report = analyze(rows)
    assert report["correctness_gate"] == "passed"
    assert report["retained_schwarz_candidates"] == ["Schwarz-128/1", "Schwarz-64/0", "Schwarz-64/1"]
    assert report["gate_s_status"] == "pending-full-endpoint-qualification"


def test_failed_geometry_is_not_retained():
    rows = [row("chain", "Schwarz-32/1", search=3.0), row("cube", "Schwarz-32/1", search=4.0),
            row("chain", "Schwarz-16/1", search=2.0, correct=False),
            row("cube", "Schwarz-16/1", search=2.0, correct=False)]
    report = analyze(rows)
    assert report["correctness_gate"] == "failed"
    assert report["common_correct_schwarz_candidates"] == ["Schwarz-32/1"]
    assert report["retained_schwarz_candidates"] == ["Schwarz-32/1"]


if __name__ == "__main__":
    test_retains_fastest_three_common_correct_schwarz_candidates()
    test_failed_geometry_is_not_retained()
