"""Compare raw block search with bounded, fully certified local candidates."""
import argparse
import json
from pathlib import Path


def _endpoint_passed(state):
    assessment = state["assessment"]
    return all(assessment[name] for name in ("inner", "gradient", "local", "identified")) and \
        assessment["endpoint_trust"]["passed"]


def _baseline(case):
    fixed = case["fixed_neighbor"]
    stationarity = fixed["first_order_stationarity_sweep"]
    snapshot = next(row for row in case["endpoint_decomposition"] if row["sweep"] == stationarity)
    sweeps = [row for row in fixed["sweep_telemetry"][:stationarity]]
    blocks = [row for row in fixed["block_telemetry"] if row["sweep"] <= stationarity]
    passed = _endpoint_passed(snapshot["raw"])
    return {"sweeps": stationarity,
        "local_assessment_count": 0,
        "local_assessment_seconds": 0.0,
        "profile_evaluations": sum(row["profile_evaluations"] for row in blocks),
        "search_seconds": sum(row["wall_seconds"] for row in sweeps),
        "endpoint_certification": passed,
        "runtime_convergence": "Passed" if passed and fixed["search_converged"] else "Failed",
        "final_eta": snapshot["eta"], "final_beta": snapshot["beta"],
        "rejected_local_candidates": 0, "objective_monotone": all(
            row["objective_after"] <= row["objective_before"] + 1e-12 +
            2e-12 * abs(row["objective_before"]) for row in sweeps),
        "maximum_local_assessment_rows": 0, "maximum_local_assessment_columns": 0,
        "process_peak_rss_mb": case["peak_rss_mb"],
        "process_peak_rss_scope": "Q1 decomposition run includes global profile and endpoint assessments"}


def _certified(case):
    fixed = case["fixed_neighbor"]
    blocks = fixed["block_telemetry"]
    attempts = [row for row in blocks if row["local_assessment_attempted"]]
    return {"sweeps": fixed["sweeps"],
        "local_assessment_count": fixed["local_assessment_count"],
        "local_assessment_seconds": fixed["local_assessment_seconds"],
        "profile_evaluations": sum(row["profile_evaluations"] for row in blocks),
        "search_seconds": fixed["search_seconds"],
        "endpoint_certification": fixed["endpoint_certified"],
        "runtime_convergence": fixed["runtime_convergence"],
        "final_eta": fixed["final_eta"], "final_beta": fixed["final_beta"],
        "rejected_local_candidates": sum(row["local_assessment_attempted"] and
            not row["local_assessment_passed"] and not row["accepted"] for row in blocks),
        "local_assessment_failures": [{"sweep": row["sweep"], "block": row["block"],
            "failure": row["local_assessment_failure"], "trust_reason": row["local_trust_reason"],
            "correction_inf_norm": row["local_correction_inf_norm"],
            "profile_gradient_inf_norm": row["local_profile_gradient_inf_norm"]}
            for row in attempts if not row["local_assessment_passed"]],
        "objective_monotone": all(row["objective_after"] <= row["objective_before"] + 1e-12 +
            2e-12 * abs(row["objective_before"]) for row in blocks),
        "maximum_local_assessment_rows": fixed["maximum_local_assessment_rows"],
        "maximum_local_assessment_columns": fixed["maximum_local_assessment_columns"],
        "process_peak_rss_mb": case["peak_rss_mb"],
        "process_peak_rss_scope": "CertifiedLocal run includes existing final global endpoint assessment"}


def compare(search_only, certified):
    baseline = _baseline(search_only)
    candidate = _certified(certified)
    eta_difference = max((abs(a - b) for a, b in zip(baseline["final_eta"], candidate["final_eta"])), default=0.0)
    beta_difference = max((abs(a - b) for a, b in zip(baseline["final_beta"], candidate["final_beta"])), default=0.0)
    for result in (baseline, candidate):
        result.pop("final_eta")
        result.pop("final_beta")
    return {"SearchOnlyLocal": baseline, "CertifiedLocal": candidate,
        "eta_inf_difference": eta_difference, "beta_inf_difference": beta_difference,
        "local_endpoint_differs": eta_difference != 0.0 or beta_difference != 0.0,
        "certified_local_improves_global_endpoint": candidate["endpoint_certification"] and
            candidate["runtime_convergence"] == "Passed" and not baseline["endpoint_certification"]}


def analyze(chain_baseline, chain_certified, cube_baseline, cube_certified):
    return {"strategy": "bounded-local-assessment", "cases": {
        "chain-256": compare(chain_baseline, chain_certified),
        "cube-256": compare(cube_baseline, cube_certified)}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--chain-baseline", required=True)
    parser.add_argument("--chain-certified", required=True)
    parser.add_argument("--cube-baseline", required=True)
    parser.add_argument("--cube-certified", required=True)
    args = parser.parse_args()
    paths = (args.chain_baseline, args.chain_certified, args.cube_baseline, args.cube_certified)
    cases = [json.loads(Path(path).read_text()) for path in paths]
    report = analyze(*cases)
    Path(args.output).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value["certified_local_improves_global_endpoint"]
        for key, value in report["cases"].items()}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
