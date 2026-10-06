"""Summarize Full, OneAccepted, and TwoAccepted FixedNeighbor policies."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


POLICIES = ("Full", "OneAccepted", "TwoAccepted")


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _factor_seconds(fixed):
    if _finite(fixed.get("profile_factor_seconds")):
        return fixed["profile_factor_seconds"]
    return sum(block.get("profile_factor_seconds", 0.0)
               for block in fixed.get("block_telemetry", []))


def summarize(case, policy, process=None):
    fixed = case.get("fixed_neighbor", case)
    process = process or {}
    search_seconds = fixed.get("search_seconds")
    total_seconds = fixed.get("total_elapsed_seconds")
    assessment_seconds = fixed.get("assessment_seconds")
    if assessment_seconds is None and _finite(search_seconds) and _finite(total_seconds):
        assessment_seconds = max(0.0, total_seconds - search_seconds)
    return {
        "topology": case.get("topology"), "atoms": case.get("atoms"), "policy": policy,
        "measurement_scope": case.get("measurement_scope", "fixed-neighbor-full-endpoint"),
        "status": process.get("status", "completed"),
        "sweeps": fixed.get("sweeps"),
        "confirmed_stationarity_sweep": fixed.get("confirmed_stationarity_sweep"),
        "block_solves": fixed.get("block_solves", sum(s.get("block_solves", 0)
                                                        for s in fixed.get("sweep_telemetry", []))),
        "profile_evaluations": fixed.get("profile_evaluations", sum(
            s.get("profile_evaluations", 0) for s in fixed.get("sweep_telemetry", []))),
        "accepted_local_updates": fixed.get("accepted_local_updates", sum(
            block.get("accepted_updates", 0) for block in fixed.get("block_telemetry", []))),
        "factor_seconds": _factor_seconds(fixed),
        "search_seconds": search_seconds,
        "assessment_seconds": assessment_seconds,
        "total_seconds": total_seconds,
        "peak_rss_mb": case.get("peak_rss_mb", process.get("sampled_tree_peak_rss_bytes", 0) / 1024**2
                               if process.get("sampled_tree_peak_rss_bytes") is not None else None),
        "objective": fixed.get("objective"),
        "endpoint_certified": fixed.get("endpoint_certified"),
        "runtime_convergence": fixed.get("runtime_convergence"),
        "search_converged": fixed.get("search_converged"),
        "failure_reason": fixed.get("search_reason"),
    }


def _correct(row):
    return row["status"] == "completed" and row["endpoint_certified"] is True and \
        row["runtime_convergence"] == "Passed"


def analyze(records):
    groups = {}
    for row in records:
        groups.setdefault((row["topology"], row["atoms"]), []).append(row)
    cases = []
    for key, rows in sorted(groups.items()):
        rows = sorted(rows, key=lambda row: POLICIES.index(row["policy"]))
        correct = [row["policy"] for row in rows if _correct(row)]
        fastest = min((row for row in rows if _finite(row.get("search_seconds"))),
                      key=lambda row: row["search_seconds"], default=None)
        cases.append({"topology": key[0], "atoms": key[1],
                      "correct_policies": correct,
                      "fastest_correct_policy": (fastest["policy"] if fastest and _correct(fastest) else None),
                      "policies": rows})
    all_cases_correct = bool(cases) and all(
        set(case["correct_policies"]) == set(POLICIES) for case in cases)
    candidates = [case["fastest_correct_policy"] for case in cases if case["fastest_correct_policy"]]
    selected = min(candidates, key=lambda policy: sum(
        row["search_seconds"] for case in cases for row in case["policies"]
        if row["policy"] == policy and _finite(row.get("search_seconds")))) if candidates else None
    return {"phase": "P3 bounded FixedNeighbor local-work comparison", "cases": cases,
            "correctness_gate": "passed" if all_cases_correct else "failed",
            "selected_policy": selected,
            "selection_basis": "smallest accepted-update budget among correct measured candidates; "
                               "search wall time is reported without a preset percentage gate"}


def write_outputs(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = ["topology", "atoms", "policy", "sweeps", "block_solves", "profile_evaluations",
               "accepted_local_updates", "factor_seconds", "search_seconds", "assessment_seconds",
               "total_seconds", "peak_rss_mb", "endpoint_certified", "runtime_convergence"]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            for row in case["policies"]:
                writer.writerow({column: row.get(column) for column in columns})


def _read_runs(paths):
    rows = []
    for path in paths:
        wrapper = json.loads(Path(path).read_text())
        result = wrapper.get("result", wrapper)
        rows.append(summarize(result, wrapper["policy"], wrapper.get("process")))
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("runs", nargs="+")
    args = parser.parse_args(argv)
    report = analyze(_read_runs(args.runs))
    write_outputs(report, args.output)
    print(report["correctness_gate"])
    return 0 if report["correctness_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
