"""Analyze bounded tuned OperatorPcg -> LegacyCompact polish experiments."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


POLICIES = ("Operator", "+1 Legacy polish", "+2 Legacy polish", "Legacy control")
KKT_TOLERANCE = 1e-10
WIDTH_TOLERANCE = 1e-12
LOCAL_TOLERANCE = 1e-10


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _value(check):
    return check.get("value") if isinstance(check, dict) else check


def _status(check):
    return check.get("status") if isinstance(check, dict) else None


def summarize(case, policy):
    fixed = case["fixed_neighbor"]
    assessment = fixed.get("endpoint_assessment", {})
    trust = assessment.get("endpoint_trust", {})
    if policy == "Operator":
        operator_seconds = fixed.get("search_seconds")
        polish_seconds = 0.0
        search_seconds = fixed.get("search_seconds")
    elif policy == "Legacy control":
        operator_seconds = 0.0
        polish_seconds = fixed.get("search_seconds")
        search_seconds = fixed.get("search_seconds")
    else:
        operator_seconds = fixed.get("operator_phase_seconds")
        polish_seconds = fixed.get("legacy_polish_seconds")
        search_seconds = fixed.get("total_search_seconds")
    sweeps = fixed.get("sweeps", [])
    if not isinstance(sweeps, list):
        sweeps = fixed.get("sweep_telemetry", [])
    last_sweep = sweeps[-1] if sweeps else {}
    return {
        "topology": case.get("topology"), "atoms": case.get("atoms"), "policy": policy,
        "status": "completed", "measurement_scope": case.get("measurement_scope"),
        "outer_core_atoms": fixed.get("outer_core_atoms", fixed.get("core_atoms")),
        "local_schwarz_core_atoms": fixed.get("local_schwarz_core_atoms"),
        "local_schwarz_overlap_hops": fixed.get("local_schwarz_overlap_hops"),
        "local_schwarz_max_block_atoms": fixed.get("local_schwarz_max_block_atoms"),
        "search_converged": fixed.get("search_converged"),
        "search_reason": fixed.get("search_reason"),
        "polish_sweeps": fixed.get("polish_sweeps", 0),
        "operator_seconds": operator_seconds,
        "polish_seconds": polish_seconds,
        "search_seconds": search_seconds,
        "total_elapsed_seconds": fixed.get("total_elapsed_seconds"),
        "eta_change_inf": last_sweep.get("eta_change_inf"),
        "beta_scaled_change": last_sweep.get("beta_scaled_change"),
        "objective": fixed.get("objective"),
        "global_ac_kkt": _value(fixed.get("global_kkt")),
        "global_width_gradient_inf_norm": fixed.get("width_gradient_inf_norm"),
        "assessment_inner": _status(fixed.get("assessment_inner")),
        "assessment_gradient": _status(fixed.get("assessment_gradient")),
        "assessment_local": _status(fixed.get("assessment_local")),
        "assessment_identified": _status(fixed.get("assessment_identified")),
        "assessment_local_correction_inf_norm": assessment.get("local_correction_inf_norm"),
        "endpoint_trust": trust.get("passed"),
        "endpoint_certified": fixed.get("endpoint_certified"),
        "runtime_convergence": fixed.get("runtime_convergence"),
        "peak_rss_mb": case.get("peak_rss_mb"),
    }


def _correct(row):
    return row["status"] == "completed" and row["search_converged"] is True and \
        _finite(row["global_ac_kkt"]) and row["global_ac_kkt"] <= KKT_TOLERANCE and \
        _finite(row["global_width_gradient_inf_norm"]) and \
        row["global_width_gradient_inf_norm"] <= WIDTH_TOLERANCE and \
        row["assessment_inner"] == "Passed" and row["assessment_gradient"] == "Passed" and \
        row["assessment_local"] == "Passed" and row["assessment_identified"] == "Passed" and \
        row["endpoint_trust"] is True and row["endpoint_certified"] is True and \
        row["runtime_convergence"] == "Passed" and \
        (_finite(row["eta_change_inf"]) is False or row["eta_change_inf"] <= 1e-10)


def analyze(rows):
    groups = {}
    for row in rows:
        groups.setdefault((row["topology"], row["atoms"]), []).append(row)
    cases = []
    for key, entries in sorted(groups.items()):
        entries = sorted(entries, key=lambda row: POLICIES.index(row["policy"]))
        cases.append({"topology": key[0], "atoms": key[1], "policies": entries,
                      "correct_policies": [row["policy"] for row in entries if _correct(row)]})
    hybrid_candidates = [policy for policy in ("+1 Legacy polish", "+2 Legacy polish")
                        if cases and all(policy in case["correct_policies"] for case in cases)]
    return {
        "phase": "P5 bounded tuned FixedNeighbor OperatorPcg + LegacyCompact polish",
        "geometry": {"outer_core_atoms": 64, "local_schwarz_core_atoms": 32,
                     "local_schwarz_overlap_hops": 0, "local_schwarz_max_block_atoms": 64},
        "cases": cases,
        "correctness_gate": "passed" if hybrid_candidates else "failed",
        "selected_policy": hybrid_candidates[0] if hybrid_candidates else None,
        "hybrid_candidates": hybrid_candidates,
        "selection_basis": "bounded polish only after pure Operator failure; strict search, confirmation, endpoint and RuntimeConvergence gates",
        "thresholds": {"global_ac_kkt": KKT_TOLERANCE, "width_gradient_inf_norm": WIDTH_TOLERANCE,
                       "assessment_local": LOCAL_TOLERANCE, "eta_confirmation": 1e-10},
    }


def write_outputs(report, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = ["topology", "atoms", "policy", "operator_seconds", "polish_seconds",
               "search_seconds", "total_elapsed_seconds", "polish_sweeps", "eta_change_inf",
               "beta_scaled_change", "assessment_local_correction_inf_norm", "global_ac_kkt",
               "global_width_gradient_inf_norm", "endpoint_trust", "endpoint_certified",
               "runtime_convergence", "peak_rss_mb"]
    with output.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            for row in case["policies"]:
                writer.writerow({column: row.get(column) for column in columns})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--case", action="append", required=True,
                        help="topology:policy:raw-json; may be repeated")
    args = parser.parse_args(argv)
    rows = []
    for specification in args.case:
        topology, policy, path = specification.split(":", 2)
        case = json.loads(Path(path).read_text())
        if case.get("topology") != topology or policy not in POLICIES:
            parser.error(f"invalid case specification: {specification}")
        rows.append(summarize(case, policy))
    report = analyze(rows)
    write_outputs(report, args.output)
    print(json.dumps({"output": str(args.output),
                      "correctness_gate": report["correctness_gate"],
                      "selected_policy": report["selected_policy"]}))
    return 0 if report["correctness_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
