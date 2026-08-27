#!/usr/bin/env python3
"""Aggregate case summaries from the convergence-exposure corpus."""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import statistics
from typing import Any, Iterable


ABSOLUTE_TOLERANCE = 1.0e-8
RELATIVE_TOLERANCE = 1.0e-3
POLICY_EXPOSURES = {
    "legacy-population": "legacy_population",
    "legacy-maximum": "maximum_gate",
    "solver-qualified": "solver_qualification",
}
MINIMUM_TOTAL_EXPOSURES = 15
MINIMUM_EXPOSURES_PER_FAMILY = 5
MINIMUM_BENEFIT_RATIO = 0.70
MAXIMUM_HARM_RATIO = 0.10
REQUIRED_FAMILIES = ("natural", "stationarity", "population")


def _materially_lower(candidate: float | None, reference: float | None) -> bool:
    if candidate is None or reference is None:
        return False
    tolerance = ABSOLUTE_TOLERANCE + RELATIVE_TOLERANCE * abs(reference)
    return candidate < reference - tolerance


def _materially_higher(candidate: float | None, reference: float | None) -> bool:
    if candidate is None or reference is None:
        return False
    tolerance = ABSOLUTE_TOLERANCE + RELATIVE_TOLERANCE * abs(reference)
    return candidate > reference + tolerance


def classify_exposure(case_summary: dict[str, Any]) -> str:
    audit = case_summary["audit"]
    if audit["status"] == "no_convergence_trigger":
        return "no-trigger"
    experiment = audit["experiments"][0]
    exposed = [
        name for name in POLICY_EXPOSURES.values()
        if experiment["exposures"][name]
    ]
    return "+".join(exposed) if exposed else "policy-agreement"


def classify_policy_outcome(
    case_summary: dict[str, Any],
    policy: str,
) -> dict[str, Any]:
    experiment = case_summary["audit"]["experiments"][0]
    production = experiment["policies"]["production"]
    candidate = experiment["policies"][policy]
    if not candidate["reached"]:
        termination_reason = experiment.get("termination_reason", "unrecorded")
        if termination_reason == "budget-exhausted":
            category = "unresolved-budget-exhausted"
        elif termination_reason == "unrecorded":
            category = "unresolved"
        else:
            category = "terminated-existing-safeguard"
        return {
            "category": category,
            "truth_improvement": False,
            "termination_reason": termination_reason,
        }

    production_truth = production.get("truth_metrics") or {}
    candidate_truth = candidate.get("truth_metrics") or {}
    production_truth_error = production_truth.get("transformed_aggregate_rmse")
    candidate_truth_error = candidate_truth.get("transformed_aggregate_rmse")
    truth_improvement = _materially_lower(
        candidate_truth_error, production_truth_error)
    truth_harm = _materially_higher(candidate_truth_error, production_truth_error)
    objective_improvement = bool(candidate["material_objective_improvement"])
    objective_harm = _materially_higher(
        candidate.get("objective"), production.get("objective"))
    safety_regression = bool(case_summary.get("safety_regression", False))
    raw_median = max(candidate.get("raw_median", [math.inf]))
    production_raw_median = max(production.get("raw_median", [math.inf]))
    raw_median_regression = _materially_higher(raw_median, production_raw_median)
    raw_median_ratio = (
        raw_median / production_raw_median
        if math.isfinite(raw_median) and
        math.isfinite(production_raw_median) and production_raw_median > 0.0
        else None)

    if objective_harm or truth_harm or safety_regression:
        category = "material-harm"
    elif objective_improvement or truth_improvement:
        category = "material-benefit"
    else:
        category = "neutral-churn"
    return {
        "category": category,
        "objective_improvement": objective_improvement,
        "truth_improvement": truth_improvement,
        "raw_median_regression": raw_median_regression,
        "raw_median_ratio": raw_median_ratio,
        "safety_regression": safety_regression,
        "extra_attempts": candidate.get("extra_attempts"),
        "extra_accepted_iterations": candidate.get("extra_accepted_iterations"),
        "objective_delta": candidate.get("objective_delta"),
        "objective_relative_delta": candidate.get("objective_relative_delta"),
        "production_truth_error": production_truth_error,
        "candidate_truth_error": candidate_truth_error,
    }


def select_replay_cases(case_rows: Iterable[dict[str, Any]]) -> list[str]:
    by_exposure: dict[str, list[str]] = {
        name: [] for name in POLICY_EXPOSURES.values()
    }
    for row in case_rows:
        for name, exposed in row["exposures"].items():
            if exposed:
                by_exposure[name].append(row["case_id"])
    for values in by_exposure.values():
        values.sort()
    selected: list[str] = []
    for exposure_name in by_exposure:
        selected.extend(by_exposure[exposure_name][:10])
    selected = sorted(set(selected))
    if len(selected) < 30:
        remaining = sorted(
            case_id for values in by_exposure.values() for case_id in values
            if case_id not in selected)
        selected.extend(remaining[:30 - len(selected)])
    return selected[:30]


def analyze(case_summaries: Iterable[dict[str, Any]]) -> dict[str, Any]:
    summaries = sorted(case_summaries, key=lambda value: value["case"]["case_id"])
    rows: list[dict[str, Any]] = []
    exposure_counts: Counter[str] = Counter()
    comparator_exposure_counts: Counter[str] = Counter()
    family_exposure_counts: dict[str, Counter[str]] = {}
    topology_exposure_counts: dict[str, Counter[str]] = {}
    termination_counts: Counter[str] = Counter()
    production_convergence_count = 0
    accepted_only_shadow_count = 0
    maximum_evidence: dict[str, Counter[str]] = {
        population: Counter() for population in (
            "production", "legacy_population", "solver_qualified")
    }
    for summary in summaries:
        exposure_class = classify_exposure(summary)
        exposure_counts[exposure_class] += 1
        family = summary["case"]["family"]
        topology = summary["case"]["topology"]
        family_exposure_counts.setdefault(family, Counter())[exposure_class] += 1
        topology_exposure_counts.setdefault(topology, Counter())[exposure_class] += 1
        audit = summary["audit"]
        production_convergence_count += int(
            audit["status"] == "counterfactual_records")
        shadow = audit.get("accepted_only_shadow", {"reached": False})
        accepted_only_shadow_count += int(shadow.get("reached", False))
        terminal = audit.get("terminal")
        exposures = {
            name: False for name in POLICY_EXPOSURES.values()
        }
        if terminal is not None:
            termination_counts[terminal["reason"]] += 1
        elif audit["status"] == "no_convergence_trigger":
            termination_counts["no-convergence-trigger"] += 1
        else:
            experiment = audit["experiments"][0]
            termination_counts[
                experiment.get("termination_reason", "unrecorded")
            ] += 1
            exposures = {
                name: bool(experiment["exposures"][name])
                for name in POLICY_EXPOSURES.values()
            }
            for name, exposed in exposures.items():
                comparator_exposure_counts[name] += int(exposed)
        trajectory = audit.get("trajectory_audit", {})
        p99_implies_max = trajectory.get("p99_implies_max", {})
        unique_blockers = trajectory.get("unique_blocker_count", {})
        for population in maximum_evidence:
            for name, count in p99_implies_max.get(population, {}).items():
                maximum_evidence[population][name] += int(count)
        legacy_maximum_unique = unique_blockers.get("legacy_maximum", {})
        maximum_evidence["production"]["accepted_max_unique_catches"] += int(
            legacy_maximum_unique.get("accepted_max", 0))
        maximum_evidence["production"]["raw_max_unique_catches"] += int(
            legacy_maximum_unique.get("raw_max", 0))
        outcomes = {}
        if exposure_class not in ("no-trigger", "policy-agreement"):
            outcomes = {
                policy: classify_policy_outcome(summary, policy)
                for policy, exposure_name in POLICY_EXPOSURES.items()
                if exposures[exposure_name]
            }
        rows.append({
            "case_id": summary["case"]["case_id"],
            "family": family,
            "topology": topology,
            "level": summary["case"]["level"],
            "exposure_class": exposure_class,
            "exposures": exposures,
            "outcomes": outcomes,
            "accepted_only_shadow": shadow,
            "terminal": terminal,
            "production_converged": audit["status"] == "counterfactual_records",
            "safety_regression": bool(summary.get("safety_regression", False)),
        })

    decisions: dict[str, Any] = {}
    for policy, exposure_name in POLICY_EXPOSURES.items():
        applicable = [
            row for row in rows if row["exposures"][exposure_name]
        ]
        outcome_counts = Counter(
            row["outcomes"][policy]["category"] for row in applicable)
        objective_improvement_count = sum(
            row["outcomes"][policy].get("objective_improvement", False)
            for row in applicable)
        truth_improvement_count = sum(
            row["outcomes"][policy].get("truth_improvement", False)
            for row in applicable)
        family_counts = Counter(row["family"] for row in applicable)
        total = len(applicable)
        benefit_ratio = outcome_counts["material-benefit"] / total if total else 0.0
        harm_ratio = outcome_counts["material-harm"] / total if total else 0.0
        raw_median_regressions = sum(
            row["outcomes"][policy].get("raw_median_regression", False)
            for row in applicable)
        safety_regressions = sum(
            row["outcomes"][policy].get("safety_regression", False)
            for row in applicable)
        raw_median_ratios = [
            row["outcomes"][policy].get("raw_median_ratio")
            for row in applicable
            if row["outcomes"][policy].get("raw_median_ratio") is not None
        ]
        aggregate_raw_median_ratio = (
            statistics.median(raw_median_ratios) if raw_median_ratios else None)
        aggregate_raw_median_regression = (
            aggregate_raw_median_ratio is None or
            _materially_higher(aggregate_raw_median_ratio, 1.0))
        family_quota_met = all(
            family_counts[name] >= MINIMUM_EXPOSURES_PER_FAMILY
            for name in REQUIRED_FAMILIES)
        candidate = (
            total >= MINIMUM_TOTAL_EXPOSURES and
            family_quota_met and
            benefit_ratio >= MINIMUM_BENEFIT_RATIO and
            harm_ratio <= MAXIMUM_HARM_RATIO and
            not aggregate_raw_median_regression and
            safety_regressions == 0
        )
        decisions[policy] = {
            "applicable_exposure_count": total,
            "family_counts": dict(sorted(family_counts.items())),
            "outcome_counts": dict(sorted(outcome_counts.items())),
            "objective_improvement_count": objective_improvement_count,
            "truth_improvement_count": truth_improvement_count,
            "benefit_ratio": benefit_ratio,
            "harm_ratio": harm_ratio,
            "raw_median_regression_count": raw_median_regressions,
            "aggregate_raw_median_ratio": aggregate_raw_median_ratio,
            "aggregate_raw_median_regression": aggregate_raw_median_regression,
            "safety_regression_count": safety_regressions,
            (
                "redesign_candidate" if policy == "solver-qualified" else
                "rollback_candidate"
            ): candidate,
        }

    genuine_exposure_count = sum(
        1 for row in rows if any(row["exposures"].values()))
    shortfall = {}
    for exposure_name in POLICY_EXPOSURES.values():
        exposed_rows = [row for row in rows if row["exposures"][exposure_name]]
        family_counts = Counter(row["family"] for row in exposed_rows)
        shortfall[exposure_name] = {
            "total": max(
                0,
                MINIMUM_TOTAL_EXPOSURES - comparator_exposure_counts[exposure_name]),
            "families": {
                family: max(
                    0,
                    MINIMUM_EXPOSURES_PER_FAMILY - family_counts[family])
                for family in REQUIRED_FAMILIES
            },
        }
    corpus_target_met = all(
        values["total"] == 0 and not any(values["families"].values())
        for values in shortfall.values()
    )
    accepted_only_saved_attempts = [
        row["terminal"]["attempts_after_accepted_only"]
        for row in rows
        if row["terminal"] is not None and
        "attempts_after_accepted_only" in row["terminal"]
    ]
    return {
        "schema_version": 4,
        "case_count": len(summaries),
        "production_convergence_count": production_convergence_count,
        "accepted_only_shadow_count": accepted_only_shadow_count,
        "accepted_only_saved_attempts": {
            "median": (
                statistics.median(accepted_only_saved_attempts)
                if accepted_only_saved_attempts else None),
            "p90": _percentile(accepted_only_saved_attempts, 0.90),
            "maximum": (
                max(accepted_only_saved_attempts)
                if accepted_only_saved_attempts else None),
        },
        "genuine_exposure_count": genuine_exposure_count,
        "exposure_counts": dict(sorted(exposure_counts.items())),
        "comparator_exposure_counts": {
            name: comparator_exposure_counts.get(name, 0)
            for name in POLICY_EXPOSURES.values()
        },
        "exposure_counts_by_family": {
            name: dict(sorted(counts.items()))
            for name, counts in sorted(family_exposure_counts.items())
        },
        "exposure_counts_by_topology": {
            name: dict(sorted(counts.items()))
            for name, counts in sorted(topology_exposure_counts.items())
        },
        "termination_counts": dict(sorted(termination_counts.items())),
        "maximum_evidence": {
            name: dict(sorted(counts.items()))
            for name, counts in maximum_evidence.items()
        },
        "corpus_target_met": corpus_target_met,
        "corpus_shortfall": shortfall,
        "replay_case_ids": select_replay_cases(rows),
        "policy_decisions": decisions,
        "cases": rows,
    }


def _percentile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def compare(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    before_by_id = {row["case_id"]: row for row in before["cases"]}
    after_by_id = {row["case_id"]: row for row in after["cases"]}
    if set(before_by_id) != set(after_by_id):
        raise ValueError("Before and after corpus case IDs do not match")

    rows = []
    objective_deltas: list[float] = []
    truth_deltas: list[float] = []
    objective_harm_count = 0
    objective_benefit_count = 0
    truth_harm_count = 0
    truth_benefit_count = 0
    for case_id in sorted(before_by_id):
        baseline = before_by_id[case_id]
        candidate = after_by_id[case_id]
        baseline_terminal = baseline.get("terminal") or {}
        candidate_terminal = candidate.get("terminal") or {}
        baseline_objective = baseline_terminal.get("objective")
        candidate_objective = candidate_terminal.get("objective")
        objective_delta = None
        if baseline_objective is not None and candidate_objective is not None:
            objective_delta = candidate_objective - baseline_objective
            objective_deltas.append(objective_delta)
            objective_harm_count += int(
                _materially_higher(candidate_objective, baseline_objective))
            objective_benefit_count += int(
                _materially_lower(candidate_objective, baseline_objective))
        baseline_truth = (
            baseline_terminal.get("truth_metrics") or {}).get(
                "transformed_aggregate_rmse")
        candidate_truth = (
            candidate_terminal.get("truth_metrics") or {}).get(
                "transformed_aggregate_rmse")
        truth_delta = None
        if baseline_truth is not None and candidate_truth is not None:
            truth_delta = candidate_truth - baseline_truth
            truth_deltas.append(truth_delta)
            truth_harm_count += int(
                _materially_higher(candidate_truth, baseline_truth))
            truth_benefit_count += int(
                _materially_lower(candidate_truth, baseline_truth))
        rows.append({
            "case_id": case_id,
            "family": candidate["family"],
            "topology": candidate["topology"],
            "before_converged": baseline.get("production_converged", False),
            "after_converged": candidate.get("production_converged", False),
            "before_stop_reason": baseline_terminal.get("reason"),
            "after_stop_reason": candidate_terminal.get("reason"),
            "objective_delta": objective_delta,
            "truth_rmse_delta": truth_delta,
            "accepted_only_shadow": candidate.get("accepted_only_shadow"),
            "safety_regression": candidate.get("safety_regression", False),
        })

    def summarize_group(group_rows: list[dict[str, Any]]) -> dict[str, Any]:
        group_objective_deltas = [
            row["objective_delta"] for row in group_rows
            if row["objective_delta"] is not None]
        group_truth_deltas = [
            row["truth_rmse_delta"] for row in group_rows
            if row["truth_rmse_delta"] is not None]
        return {
            "case_count": len(group_rows),
            "before_convergence_count": sum(
                row["before_converged"] for row in group_rows),
            "after_convergence_count": sum(
                row["after_converged"] for row in group_rows),
            "accepted_only_shadow_count": sum(
                bool((row.get("accepted_only_shadow") or {}).get("reached"))
                for row in group_rows),
            "objective_delta_median": (
                statistics.median(group_objective_deltas)
                if group_objective_deltas else None),
            "objective_delta_p90": _percentile(group_objective_deltas, 0.90),
            "truth_rmse_delta_median": (
                statistics.median(group_truth_deltas)
                if group_truth_deltas else None),
            "truth_rmse_delta_p90": _percentile(group_truth_deltas, 0.90),
            "safety_regression_count": sum(
                row["safety_regression"] for row in group_rows),
        }

    family_names = sorted({row["family"] for row in rows})
    topology_names = sorted({
        (row["family"], row["topology"]) for row in rows})
    by_family = {
        family: summarize_group([
            row for row in rows if row["family"] == family])
        for family in family_names
    }
    by_topology = {
        f"{family}/{topology}": summarize_group([
            row for row in rows
            if row["family"] == family and row["topology"] == topology])
        for family, topology in topology_names
    }

    return {
        "schema_version": 1,
        "case_count": len(rows),
        "before_convergence_count": before.get("production_convergence_count", 0),
        "after_convergence_count": after.get("production_convergence_count", 0),
        "accepted_only_shadow_count": after.get("accepted_only_shadow_count", 0),
        "before_stop_reasons": before.get("termination_counts", {}),
        "after_stop_reasons": after.get("termination_counts", {}),
        "objective_delta": {
            "median": statistics.median(objective_deltas) if objective_deltas else None,
            "p90": _percentile(objective_deltas, 0.90),
            "material_benefit_count": objective_benefit_count,
            "material_harm_count": objective_harm_count,
        },
        "truth_rmse_delta": {
            "median": statistics.median(truth_deltas) if truth_deltas else None,
            "p90": _percentile(truth_deltas, 0.90),
            "material_benefit_count": truth_benefit_count,
            "material_harm_count": truth_harm_count,
        },
        "safety_regression_count": sum(
            row["safety_regression"] for row in rows),
        "by_family": by_family,
        "by_topology": by_topology,
        "cases": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path, nargs="+")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    summaries = [
        json.loads(path.read_text(encoding="utf-8")) for path in args.summary
    ]
    output = json.dumps(analyze(summaries), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    else:
        print(output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
