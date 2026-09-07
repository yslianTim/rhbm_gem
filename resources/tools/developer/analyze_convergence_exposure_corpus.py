#!/usr/bin/env python3
"""Aggregate and compare production-only convergence corpus summaries."""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import statistics
from typing import Any, Iterable, Sequence


SCHEMA_VERSION = 9
COMPARISON_SCHEMA_VERSION = 5
TRUTH_METRIC = "transformed_aggregate_rmse"


def _percentile(values: Iterable[float], percentile: float) -> float | None:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return None
    index = (len(ordered) - 1) * percentile
    lower = int(index)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = index - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _distribution(values: Iterable[float]) -> dict[str, float | None]:
    rows = [v for v in values if finite(v)]
    return {
        "median": statistics.median(rows) if rows else None,
        "p90": _percentile(rows, 0.90),
        "p99": _percentile(rows, 0.99),
    }


def _terminal(summary: dict[str, Any]) -> dict[str, Any]:
    terminal = summary.get("terminal")
    return terminal if isinstance(terminal, dict) else {}


def _case_row(summary: dict[str, Any]) -> dict[str, Any]:
    terminal = _terminal(summary)
    truth_metrics = terminal.get("truth_metrics", {})
    return {
        "case_id": summary["case"]["case_id"],
        "family": summary["case"]["family"],
        "topology": summary["case"]["topology"],
        "seed": summary["case"]["seed"],
        "production_converged": bool(summary.get("production_converged", False)),
        "stop_reason": terminal.get("reason"),
        "objective": terminal.get("objective"),
        "truth_rmse": truth_metrics.get(TRUTH_METRIC),
        "accepted_iteration": terminal.get("accepted_iteration"),
        "elapsed_seconds": summary.get("elapsed_seconds"),
        "semantic_trajectory_sha256": summary.get(
            "semantic_trajectory_sha256"),
        "terminal_state_sha256": summary.get("terminal_state_sha256"),
        "safety_regression": summary.get("safety_regression") is not False,
        "frozen_truth_sha256": summary.get("frozen_truth_sha256"),
        "diagnostics": summary.get("diagnostics", {}),
    }


def analyze(summaries: Iterable[dict[str, Any]]) -> dict[str, Any]:
    summaries = list(summaries)
    complete = [
        summary for summary in summaries if summary.get("status") == "complete"
    ]
    cases = [_case_row(summary) for summary in complete]
    if len({row["case_id"] for row in cases}) != len(cases):
        raise ValueError("Duplicate corpus case identity")
    return {
        "schema_version": SCHEMA_VERSION,
        "case_count": len(cases),
        "failed_case_count": len(summaries) - len(complete),
        "production_convergence_count": sum(
            row["production_converged"] for row in cases),
        "termination_counts": dict(sorted(Counter(
            row["stop_reason"] for row in cases).items())),
        "safety_regression_count": sum(
            row["safety_regression"] for row in cases),
        "elapsed_seconds": _distribution(
            row["elapsed_seconds"] for row in cases
            if row["elapsed_seconds"] is not None),
        "cases": cases,
        "diagnostics": diagnostic_distribution(cases),
    }


def _paired_rows(
    before: dict[str, Any],
    after: dict[str, Any],
) -> list[tuple[dict[str, Any], dict[str, Any]]]:
    for corpus in (before, after):
        if len({row["case_id"] for row in corpus["cases"]}) != len(corpus["cases"]):
            raise ValueError("Duplicate corpus case identity")
    before_by_id = {row["case_id"]: row for row in before["cases"]}
    after_by_id = {row["case_id"]: row for row in after["cases"]}
    if before_by_id.keys() != after_by_id.keys():
        raise ValueError("Before/after corpus case identities differ")
    return [
        (before_by_id[case_id], after_by_id[case_id])
        for case_id in sorted(before_by_id)
    ]


def _numeric_delta(
    pairs: list[tuple[dict[str, Any], dict[str, Any]]],
    field: str,
) -> dict[str, float | None]:
    deltas = [
        float(after[field]) - float(before[field])
        for before, after in pairs
        if finite(before.get(field)) and finite(after.get(field))
    ]
    return {"count": len(deltas), **_distribution(deltas)}


def finite(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def digest(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def diagnostic_distribution(cases: list[dict[str, Any]]) -> dict[str, Any]:
    groups: dict[str, list[dict[str, str]]] = {}
    for case in cases:
        for event in case.get("diagnostics", {}).get("conditioning", []):
            for key in ("all/" + event["phase"] + "/" + event["solve"],
                        case["family"] + "/" + case["topology"] + "/" + event["phase"] + "/" + event["solve"]):
                groups.setdefault(key, []).append(event)
    result = {}
    for key, events in sorted(groups.items()):
        ratios = [float(e["pivot-ratio"]) for e in events if float(e["pivot-ratio"]) > 0]
        result[key] = {
            "event_count": len(events),
            "guard_count": sum(e["guard"] == "1" for e in events),
            "sentinel_count": sum(float(e["pivot-ratio"]) == 0 for e in events),
            "pivot_ratio": {"count": len(ratios), "minimum": min(ratios) if ratios else None,
                            "p01": _percentile(ratios, .01), **_distribution(ratios)},
            "ridge_multiplier_max": _distribution(float(e["ridge-max"]) for e in events),
        }
    solve_groups: dict[str, Counter] = {}
    availability_groups: dict[str, Counter] = {}
    populations: dict[str, list[list[int]]] = {}
    for case in cases:
        diagnostics = case.get("diagnostics", {})
        for event in diagnostics.get("populations", []):
            for prefix in ("all", case["family"] + "/" + case["topology"]):
                for field in ("accepted-active-population", "operator-nominal-population"):
                    populations.setdefault(prefix + "/" + field, []).append([int(v) for v in event[field].split("/")])
        for event in diagnostics.get("solves", []):
            for prefix in ("all", case["family"] + "/" + case["topology"]):
                key = prefix + "/" + event["phase"] + "/" + event["solve"]
                solve_groups.setdefault(key, Counter())[event["status"]] += 1
        for event in diagnostics.get("availability", []):
            for prefix in ("all", case["family"] + "/" + case["topology"]):
                key = prefix + "/" + event["phase"]
                counts = availability_groups.setdefault(key, Counter())
                counts["event_count"] += 1
                for field in ("nominal-atoms", "shape-unavailable", "offset-unavailable"):
                    counts[field] += int(event[field])
    return {"conditioning_case_count": sum(bool(c.get("diagnostics", {}).get("conditioning")) for c in cases),
            "case_count": len(cases),
            "populations": {k: {"event_count": len(v), "coordinates": [_distribution(row[i] for row in v) for i in range(3)]} for k, v in sorted(populations.items())},
            "conditioning": result,
            "solves": {k: dict(v) for k, v in sorted(solve_groups.items())},
            "availability": {k: dict(v) for k, v in sorted(availability_groups.items())}}


def compare(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    if before.get("schema_version") != SCHEMA_VERSION or after.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported corpus aggregate schema")
    pairs = _paired_rows(before, after)
    complete = (len(pairs) == before.get("case_count") == after.get("case_count") == 600
                and before.get("failed_case_count", 0) == after.get("failed_case_count", 0) == 0)
    identity = all(all(a.get(k) == b.get(k) for k in ("family", "topology", "seed"))
                   and digest(a.get("frozen_truth_sha256"))
                   and a["frozen_truth_sha256"] == b.get("frozen_truth_sha256") for a, b in pairs)
    def matches(field: str) -> int:
        return sum(digest(a.get(field)) and a[field] == b.get(field) for a, b in pairs)
    semantic_matches = matches("semantic_trajectory_sha256")
    terminal_matches = matches("terminal_state_sha256")
    safety_count = sum(a.get("safety_regression") is not False or b.get("safety_regression") is not False for a, b in pairs)
    deltas = {field: _numeric_delta(pairs, field) for field in ("objective", "truth_rmse", "accepted_iteration", "elapsed_seconds")}
    exact = all(all(finite(a.get(k)) and finite(b.get(k)) and a[k] == b[k]
                    for k in ("objective", "truth_rmse", "accepted_iteration"))
                and isinstance(a.get("stop_reason"), str) and a["stop_reason"] == b.get("stop_reason") for a, b in pairs)
    timing_valid = all(finite(a.get("elapsed_seconds")) and a["elapsed_seconds"] > 0
                       and finite(b.get("elapsed_seconds")) and b["elapsed_seconds"] > 0 for a, b in pairs)
    elapsed_before = _distribution(a["elapsed_seconds"] for a, _ in pairs if finite(a.get("elapsed_seconds")))
    elapsed_after = _distribution(b["elapsed_seconds"] for _, b in pairs if finite(b.get("elapsed_seconds")))
    faster = bool(pairs) and timing_valid and all(elapsed_after[k] < elapsed_before[k] for k in ("median", "p90"))
    def gate(conditions: dict[str, bool]) -> dict[str, Any]:
        return {"passed": all(conditions.values()), "conditions": conditions}
    strata = {}
    for key in sorted({a["family"] + "/" + a["topology"] for a, _ in pairs}):
        subset = [(a, b) for a, b in pairs if a["family"] + "/" + a["topology"] == key]
        strata[key] = {"case_count": len(subset), "elapsed_delta": _numeric_delta(subset, "elapsed_seconds")}
    return {
        "schema_version": COMPARISON_SCHEMA_VERSION, "case_count": len(pairs),
        "production_semantic_match_count": semantic_matches,
        "terminal_state_match_count": terminal_matches,
        "before_termination_counts": before.get("termination_counts", {}),
        "after_termination_counts": after.get("termination_counts", {}),
        "stop_distribution_matches": before.get("termination_counts") == after.get("termination_counts"),
        "objective_delta": deltas["objective"], "truth_rmse_delta": deltas["truth_rmse"],
        "accepted_iteration_delta": deltas["accepted_iteration"],
        "safety_regression_count": safety_count,
        "elapsed_seconds": {"before": elapsed_before, "after": elapsed_after, "strictly_lower": faster,
                            "delta": deltas["elapsed_seconds"],
                            "improved_fraction": sum(finite(a.get("elapsed_seconds")) and finite(b.get("elapsed_seconds")) and b["elapsed_seconds"] < a["elapsed_seconds"] for a, b in pairs) / len(pairs) if pairs else None},
        "per_case": [{"case_id": a["case_id"], "elapsed_delta": b["elapsed_seconds"] - a["elapsed_seconds"] if finite(a.get("elapsed_seconds")) and finite(b.get("elapsed_seconds")) else None} for a, b in pairs],
        "strata": strata,
        "safety_gate": gate({"complete-pair": complete, "frozen-identities": identity, "zero-safety-regression": safety_count == 0}),
        "quality_gate": gate({"complete-pair": complete, "frozen-identities": identity,
                              "production-semantic-digest": semantic_matches == len(pairs),
                              "terminal-state-digest": terminal_matches == len(pairs), "per-case-outcomes-identical": exact}),
        "efficiency_gate": gate({"complete-pair": complete, "frozen-identities": identity, "complete-timing": timing_valid, "audit-cost-strictly-lower": faster}),
    }


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summaries", nargs="*", type=Path)
    parser.add_argument("--before", type=Path)
    parser.add_argument("--after", type=Path)
    parser.add_argument("--json", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_arguments(argv)
    if args.before is not None or args.after is not None:
        if args.before is None or args.after is None or args.summaries:
            raise ValueError("Use both --before and --after without summaries")
        report = compare(
            json.loads(args.before.read_text(encoding="utf-8")),
            json.loads(args.after.read_text(encoding="utf-8")))
    else:
        report = analyze(
            json.loads(path.read_text(encoding="utf-8"))
            for path in args.summaries)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json is None:
        print(payload, end="")
    else:
        args.json.write_text(payload, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
