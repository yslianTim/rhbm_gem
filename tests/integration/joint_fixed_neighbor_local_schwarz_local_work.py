"""Qualify local-work budgets with a fixed tuned Schwarz geometry."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


POLICIES = ("OneAccepted", "TwoAccepted", "Full")
KKT_TOLERANCE = 1e-10
WIDTH_TOLERANCE = 1e-12


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _value(check):
    return check.get("value") if isinstance(check, dict) else check


def _status(check):
    return check.get("status") if isinstance(check, dict) else None


def summarize(case, policy):
    fixed = case["fixed_neighbor"]
    work = fixed.get("local_operator_work", {})
    blocks = fixed.get("block_telemetry", [])
    sweeps = fixed.get("sweep_telemetry", [])
    assessment = fixed.get("endpoint_assessment", {})
    trust = assessment.get("endpoint_trust", {})
    return {
        "topology": case.get("topology"), "atoms": case.get("atoms"), "policy": policy,
        "status": "completed", "measurement_scope": "fixed-neighbor-full-endpoint",
        "outer_core_atoms": fixed.get("outer_core_atoms", fixed.get("core_atoms")),
        "local_work_policy": fixed.get("local_work_policy"),
        "local_search_method": fixed.get("local_search_method"),
        "local_preconditioner": fixed.get("local_preconditioner"),
        "local_schwarz_core_atoms": fixed.get("local_schwarz_core_atoms"),
        "local_schwarz_overlap_hops": fixed.get("local_schwarz_overlap_hops"),
        "local_schwarz_max_block_atoms": fixed.get("local_schwarz_max_block_atoms"),
        "search_converged": fixed.get("search_converged"),
        "search_reason": fixed.get("search_reason"),
        "sweeps": fixed.get("sweeps", len(sweeps)),
        "block_solves": sum(item.get("block_solves", 0) for item in sweeps),
        "profile_evaluations": sum(item.get("profile_evaluations", 0) for item in sweeps),
        "accepted_local_updates": sum(item.get("accepted_updates", 0) for item in blocks),
        "search_seconds": fixed.get("search_seconds"),
        "assessment_seconds": fixed.get("assessment_seconds"),
        "total_seconds": fixed.get("total_elapsed_seconds"),
        "pcg_solves": work.get("pcg_solves"),
        "pcg_iterations": work.get("pcg_iterations"),
        "pcg_mean_iterations": work.get("pcg_mean_iterations"),
        "pcg_max_iterations": work.get("pcg_max_iterations"),
        "preconditioner_setup_seconds": work.get("preconditioner_setup_seconds"),
        "assessment_inner": _status(fixed.get("assessment_inner")),
        "assessment_gradient": _status(fixed.get("assessment_gradient")),
        "assessment_local": _status(fixed.get("assessment_local")),
        "assessment_identified": _status(fixed.get("assessment_identified")),
        "endpoint_trust": trust.get("passed"),
        "endpoint_certified": fixed.get("endpoint_certified"),
        "runtime_convergence": fixed.get("runtime_convergence"),
        "global_ac_kkt": _value(fixed.get("global_kkt")),
        "global_width_gradient_inf_norm": fixed.get("width_gradient_inf_norm"),
        "assessment_local_correction_inf_norm":
            (assessment.get("local_correction_inf_norm") if assessment else None),
        "peak_rss_mb": case.get("peak_rss_mb"),
        "objective": fixed.get("objective"),
    }


def _correct(row):
    return row["status"] == "completed" and row["search_converged"] is True and \
        _finite(row["global_ac_kkt"]) and row["global_ac_kkt"] <= KKT_TOLERANCE and \
        _finite(row["global_width_gradient_inf_norm"]) and \
        row["global_width_gradient_inf_norm"] <= WIDTH_TOLERANCE and \
        row["assessment_inner"] == "Passed" and row["assessment_gradient"] == "Passed" and \
        row["assessment_local"] == "Passed" and row["assessment_identified"] == "Passed" and \
        row["endpoint_trust"] is True and row["endpoint_certified"] is True and \
        row["runtime_convergence"] == "Passed"


def analyze(rows):
    groups = {}
    for row in rows:
        groups.setdefault((row["topology"], row["atoms"]), []).append(row)
    cases = []
    for key, entries in sorted(groups.items()):
        entries = sorted(entries, key=lambda row: POLICIES.index(row["policy"]))
        cases.append({"topology": key[0], "atoms": key[1], "policies": entries,
                      "correct_policies": [row["policy"] for row in entries if _correct(row)]})
    common = [policy for policy in POLICIES if cases and all(
        policy in case["correct_policies"] for case in cases)]
    return {
        "phase": "P4 tuned FixedNeighbor local-work qualification",
        "geometry": {"outer_core_atoms": 64, "local_schwarz_core_atoms": 32,
                     "local_schwarz_overlap_hops": 0, "local_schwarz_max_block_atoms": 64},
        "cases": cases,
        "correctness_gate": "passed" if common else "failed",
        "common_correct_policies": common,
        "selected_policy": common[0] if common else None,
        "selection_basis": "first smallest accepted-update budget passing all endpoint gates; no threshold relaxation",
        "thresholds": {"global_ac_kkt": KKT_TOLERANCE, "width_gradient_inf_norm": WIDTH_TOLERANCE,
                       "assessment_local": 1e-10, "eta_confirmation": 1e-10},
    }


def write_outputs(report, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = ["topology", "atoms", "policy", "local_work_policy", "search_seconds",
               "assessment_seconds", "total_seconds", "sweeps", "block_solves",
               "profile_evaluations", "accepted_local_updates", "pcg_mean_iterations",
               "pcg_max_iterations", "preconditioner_setup_seconds", "assessment_local",
               "global_ac_kkt", "global_width_gradient_inf_norm", "endpoint_trust",
               "endpoint_certified", "runtime_convergence", "peak_rss_mb"]
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
