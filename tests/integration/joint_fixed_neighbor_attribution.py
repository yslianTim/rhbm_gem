"""Analyze accepted-update and factor-time telemetry from FixedNeighbor block visits."""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _distribution(values):
    values = [value for value in values if _finite(value)]
    if not values:
        return {"count": 0, "minimum": None, "maximum": None, "mean": None, "median": None}
    return {"count": len(values), "minimum": min(values), "maximum": max(values),
            "mean": statistics.fmean(values), "median": statistics.median(values)}


def _accepted_update(trial):
    return trial.get("accepted") is True and (trial.get("accepted_update") or 0) > 0


def _visit(block):
    trials = block.get("profile_trials") or []
    accepted = [(index, trial) for index, trial in enumerate(trials) if _accepted_update(trial)]
    reductions = [trial.get("objective_reduction") for _, trial in accepted
                  if _finite(trial.get("objective_reduction"))]
    total_reduction = sum(reductions)
    first_fraction = reductions[0] / total_reduction if reductions and total_reduction > 0 else None
    first_two_fraction = (sum(reductions[:2]) / total_reduction
                          if reductions and total_reduction > 0 else None)
    factor_values = [trial.get("factor_seconds", 0.0) for trial in trials]
    factor_seconds = sum(value for value in factor_values if _finite(value))
    first_index = accepted[0][0] if accepted else None
    second_index = accepted[1][0] if len(accepted) > 1 else None
    after_first = (sum(value for index, value in enumerate(factor_values)
                       if _finite(value) and first_index is not None and index > first_index)
                   if first_index is not None else None)
    after_second = (sum(value for index, value in enumerate(factor_values)
                        if _finite(value) and second_index is not None and index > second_index)
                    if second_index is not None else None)
    return {
        "sweep": block.get("sweep"), "block": block.get("block"),
        "profile_evaluations": block.get("profile_evaluations", len(trials)),
        "accepted_updates": len(accepted), "search_seconds": block.get("search_seconds"),
        "factor_seconds": factor_seconds,
        "total_objective_reduction": total_reduction,
        "first_update_reduction_fraction": first_fraction,
        "first_two_update_reduction_fraction": first_two_fraction,
        "factor_time_fraction_after_first_update": (
            after_first / factor_seconds if after_first is not None and factor_seconds > 0 else None),
        "factor_time_fraction_after_second_update": (
            after_second / factor_seconds if after_second is not None and factor_seconds > 0 else None),
    }


def analyze_case(case):
    fixed = case.get("fixed_neighbor", case)
    blocks = fixed.get("block_telemetry") or []
    visits = [_visit(block) for block in blocks]
    return {
        "topology": case.get("topology"), "atoms": case.get("atoms"),
        "measurement_scope": case.get("measurement_scope", "fixed-neighbor-search-only"),
        "trajectory_available": fixed.get("local_trajectory_telemetry") is True and
            bool(visits) and all("profile_trials" in block for block in blocks),
        "block_visits": len(visits),
        "sweeps": fixed.get("sweeps"),
        "search_seconds": fixed.get("search_seconds"),
        "profile_factor_seconds": fixed.get("profile_factor_seconds"),
        "accepted_updates_per_visit": _distribution([visit["accepted_updates"] for visit in visits]),
        "profile_evaluations_per_visit": _distribution([visit["profile_evaluations"] for visit in visits]),
        "search_seconds_per_visit": _distribution([visit["search_seconds"] for visit in visits]),
        "factor_seconds_per_visit": _distribution([visit["factor_seconds"] for visit in visits]),
        "first_update_reduction_fraction": _distribution(
            [visit["first_update_reduction_fraction"] for visit in visits]),
        "first_two_update_reduction_fraction": _distribution(
            [visit["first_two_update_reduction_fraction"] for visit in visits]),
        "factor_time_fraction_after_first_update": _distribution(
            [visit["factor_time_fraction_after_first_update"] for visit in visits]),
        "factor_time_fraction_after_second_update": _distribution(
            [visit["factor_time_fraction_after_second_update"] for visit in visits]),
        "visits": visits,
    }


def analyze(cases, decision=None):
    results = [analyze_case(case) for case in cases]
    return {
        "phase": "P2 FixedNeighbor local-work attribution",
        "cases": sorted(results, key=lambda row: (row.get("topology") or "", row.get("atoms") or 0)),
        "over_solving_decision": decision or "insufficient evidence",
        "decision_basis": (
            "The analyzer reports first-update progress and post-update factor fractions. "
            "The over-solving label is supplied after reviewing these measured distributions; "
            "no percentage gate is encoded here."
        ),
    }


def write_outputs(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = ["topology", "atoms", "block_visits", "accepted_updates_mean",
               "profile_evaluations_mean", "search_seconds_mean", "factor_seconds_mean",
               "first_update_fraction_median",
               "first_two_update_fraction_median", "tail_factor_fraction_after_first_median",
               "tail_factor_fraction_after_second_median"]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            def median(name):
                return case[name]["median"]
            writer.writerow({
                "topology": case["topology"], "atoms": case["atoms"],
                "block_visits": case["block_visits"],
                "accepted_updates_mean": case["accepted_updates_per_visit"]["mean"],
                "profile_evaluations_mean": case["profile_evaluations_per_visit"]["mean"],
                "search_seconds_mean": case["search_seconds_per_visit"]["mean"],
                "factor_seconds_mean": case["factor_seconds_per_visit"]["mean"],
                "first_update_fraction_median": median("first_update_reduction_fraction"),
                "first_two_update_fraction_median": median("first_two_update_reduction_fraction"),
                "tail_factor_fraction_after_first_median": median("factor_time_fraction_after_first_update"),
                "tail_factor_fraction_after_second_median": median("factor_time_fraction_after_second_update"),
            })


def _read_cases(paths):
    cases = []
    for path in paths:
        value = json.loads(Path(path).read_text())
        cases.append(value.get("result", value))
    return cases


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--decision", choices=("yes", "no", "insufficient evidence"))
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args(argv)
    report = analyze(_read_cases(args.cases), args.decision)
    write_outputs(report, args.output)
    print(report["over_solving_decision"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
