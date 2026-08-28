#!/usr/bin/env python3
"""Aggregate and compare production-only convergence corpus summaries."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import statistics
from typing import Any, Iterable, Sequence


SCHEMA_VERSION = 8
COMPARISON_SCHEMA_VERSION = 4
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
    rows = list(values)
    return {
        "median": statistics.median(rows) if rows else None,
        "p90": _percentile(rows, 0.90),
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
        "safety_regression": bool(summary.get("safety_regression", False)),
    }


def analyze(summaries: Iterable[dict[str, Any]]) -> dict[str, Any]:
    complete = [
        summary for summary in summaries if summary.get("status") == "complete"
    ]
    cases = [_case_row(summary) for summary in complete]
    return {
        "schema_version": SCHEMA_VERSION,
        "case_count": len(cases),
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
    }


def _paired_rows(
    before: dict[str, Any],
    after: dict[str, Any],
) -> list[tuple[dict[str, Any], dict[str, Any]]]:
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
        if before.get(field) is not None and after.get(field) is not None
    ]
    return {"count": len(deltas), **_distribution(deltas)}


def compare(
    before: dict[str, Any],
    after: dict[str, Any],
) -> dict[str, Any]:
    pairs = _paired_rows(before, after)
    elapsed_before = _distribution(
        before_row["elapsed_seconds"] for before_row, _ in pairs
        if before_row.get("elapsed_seconds") is not None)
    elapsed_after = _distribution(
        after_row["elapsed_seconds"] for _, after_row in pairs
        if after_row.get("elapsed_seconds") is not None)
    elapsed_strictly_lower = (
        elapsed_before["median"] is not None and
        elapsed_before["p90"] is not None and
        elapsed_after["median"] is not None and
        elapsed_after["p90"] is not None and
        elapsed_after["median"] < elapsed_before["median"] and
        elapsed_after["p90"] < elapsed_before["p90"]
    )
    semantic_matches = sum(
        before_row.get("semantic_trajectory_sha256") ==
        after_row.get("semantic_trajectory_sha256")
        for before_row, after_row in pairs)
    terminal_matches = sum(
        before_row.get("terminal_state_sha256") ==
        after_row.get("terminal_state_sha256")
        for before_row, after_row in pairs)
    safety_regressions = sum(
        after_row.get("safety_regression", False) for _, after_row in pairs)
    objective_delta = _numeric_delta(pairs, "objective")
    truth_delta = _numeric_delta(pairs, "truth_rmse")
    accepted_delta = _numeric_delta(pairs, "accepted_iteration")
    stop_distribution_matches = (
        before.get("termination_counts") == after.get("termination_counts"))
    zero_outcome_deltas = all(
        report[statistic] == 0.0
        for report in (objective_delta, truth_delta, accepted_delta)
        for statistic in ("median", "p90")
    )
    blocking_conditions = {
        "complete-pair": (
            len(pairs) == before.get("case_count") == after.get("case_count")),
        "production-semantic-digest": semantic_matches == len(pairs),
        "terminal-state-digest": terminal_matches == len(pairs),
        "stop-distribution": stop_distribution_matches,
        "zero-outcome-deltas": zero_outcome_deltas,
        "zero-safety-regression": safety_regressions == 0,
        "audit-cost-strictly-lower": elapsed_strictly_lower,
    }
    return {
        "schema_version": COMPARISON_SCHEMA_VERSION,
        "case_count": len(pairs),
        "production_semantic_match_count": semantic_matches,
        "terminal_state_match_count": terminal_matches,
        "before_termination_counts": before.get("termination_counts", {}),
        "after_termination_counts": after.get("termination_counts", {}),
        "stop_distribution_matches": stop_distribution_matches,
        "objective_delta": objective_delta,
        "truth_rmse_delta": truth_delta,
        "accepted_iteration_delta": accepted_delta,
        "safety_regression_count": safety_regressions,
        "elapsed_seconds": {
            "before": elapsed_before,
            "after": elapsed_after,
            "strictly_lower": elapsed_strictly_lower,
        },
        "blocking_gate": {
            "passed": all(blocking_conditions.values()),
            "conditions": blocking_conditions,
        },
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
