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
    "fixed-point-operator": "fixed_point_operator",
    "fixed-point-operator-maximum": "fixed_point_operator_maximum",
}
MINIMUM_TOTAL_EXPOSURES = 15
MINIMUM_EXPOSURES_PER_FAMILY = 5
MINIMUM_BENEFIT_RATIO = 0.70
MAXIMUM_HARM_RATIO = 0.10
REQUIRED_FAMILIES = ("natural", "stationarity", "population")
SHADOW_REQUIRED_FAMILIES = ("natural", "stationarity")
SHADOW_POLICIES = (
    "accepted-only-k2",
    "accepted-only-k3",
    "accepted-only-k5",
    "dynamic-raw",
)
TRUTH_RMSE_METRICS = (
    "amplitude_rmse",
    "width_rmse",
    "offset_rmse",
    "transformed_aggregate_rmse",
)
GUARD_TRUST_DECOUPLING_CASE_COUNT = 600


def summarize_trust_model_shadow(rows: list[dict[str, Any]]) -> dict[str, Any]:
    records = [
        {
            **record,
            "case-id": row["case_id"],
            "family": row["family"],
            "topology": f'{row["family"]}/{row["topology"]}',
        }
        for row in rows
        for record in row.get("trust_model_shadow", {}).get("records", [])
    ]

    def rho_bin(record: dict[str, Any]) -> str:
        rho = record.get("rho")
        if rho is None:
            return "unavailable"
        if rho < 0.25:
            return "low"
        if rho <= 0.75:
            return "mid"
        return (
            "high-boundary"
            if float(record["boundary-utilization"]) >= 0.8 else
            "high-interior")

    def summarize(group: list[dict[str, Any]]) -> dict[str, Any]:
        statuses = Counter(str(record["status"]) for record in group)
        rho_bins = Counter(rho_bin(record) for record in group)
        eligible_count = sum(
            record["status"] != "nonmaterial-step" for record in group)
        action_confusion = Counter(
            f'{record["current-action"]}->{record["shadow-action"]}'
            for record in group
            if record["readiness-eligible"] and
            record["shadow-action"] != "suppressed")
        rho_values = [
            float(record["rho"]) for record in group
            if record.get("rho") is not None]
        calibration_errors = [abs(value - 1.0) for value in rho_values]
        return {
            "record_count": len(group),
            "eligible_count": eligible_count,
            "available_count": statuses["available"],
            "availability_ratio": (
                statuses["available"] / eligible_count
                if eligible_count else None),
            "status_counts": dict(sorted(statuses.items())),
            "rho_bins": dict(sorted(rho_bins.items())),
            "rho": {
                "median": statistics.median(rho_values) if rho_values else None,
                "p90": _percentile(rho_values, 0.90),
                "absolute_calibration_error_median": (
                    statistics.median(calibration_errors)
                    if calibration_errors else None),
                "absolute_calibration_error_p90": _percentile(
                    calibration_errors, 0.90),
            },
            "action_confusion": dict(sorted(action_confusion.items())),
            "elapsed_milliseconds": sum(
                float(record["elapsed-ms"]) for record in group),
        }

    stratum_value = {
        "family": lambda record: str(record["family"]),
        "topology": lambda record: str(record["topology"]),
        "source": lambda record: str(record["source"]),
        "trial-disposition": lambda record: str(record["trial-disposition"]),
        "factor": lambda record: f'{float(record["factor"]):.8g}',
        "prediction-status": lambda record: str(record["status"]),
        "boundary": lambda record: (
            "rescued" if record["boundary-rescued"] else
            "touched" if record["boundary-touched"] else "none"),
        "cluster-size": lambda record: str(record["atoms"]),
        "unselected-dependencies": lambda record: (
            "0" if int(record["unselected-dependencies"]) == 0 else
            "1" if int(record["unselected-dependencies"]) == 1 else "2+"),
    }
    overall = summarize(records)
    strata = {
        name: {
            value: summarize([
                record for record in records if classifier(record) == value])
            for value in sorted({classifier(record) for record in records})
        }
        for name, classifier in stratum_value.items()
    }
    action_records = [
        record for record in records
        if record["readiness-eligible"] and
        record["shadow-action"] != "suppressed"
    ]
    action_differences = [
        record for record in action_records
        if record["current-action"] != record["shadow-action"]
    ]
    shrink_opportunities = sum(
        record["current-action"] != "shrink" and
        record["shadow-action"] == "shrink"
        for record in action_differences)
    growth_related_opportunities = sum(
        "grow" in (record["current-action"], record["shadow-action"])
        for record in action_differences)
    performance_ratios = [
        float(value)
        for row in rows
        if (value := row.get("trust_model_shadow", {}).get(
            "candidate_phase_ratio")) is not None
    ]
    performance = {
        "case_count": len(performance_ratios),
        "candidate_phase_ratio_median": (
            statistics.median(performance_ratios)
            if performance_ratios else None),
        "candidate_phase_ratio_p90": _percentile(performance_ratios, 0.90),
    }
    family_summaries = strata["family"]
    topology_summaries = strata["topology"]
    high_count = (
        overall["rho_bins"].get("high-interior", 0) +
        overall["rho_bins"].get("high-boundary", 0))
    diversity_counts = {
        "low": overall["rho_bins"].get("low", 0),
        "mid": overall["rho_bins"].get("mid", 0),
        "high": high_count,
        "high-boundary": overall["rho_bins"].get("high-boundary", 0),
    }
    family_diversity_passed = all(
        summary["rho_bins"].get("low", 0) >= 20 and
        summary["rho_bins"].get("mid", 0) >= 20 and
        summary["rho_bins"].get("high-interior", 0) +
        summary["rho_bins"].get("high-boundary", 0) >= 20
        for family, summary in family_summaries.items()
        if family in REQUIRED_FAMILIES
    ) and all(family in family_summaries for family in REQUIRED_FAMILIES)
    family_calibration_passed = all(
        summary["rho"]["absolute_calibration_error_median"] is not None and
        summary["rho"]["absolute_calibration_error_median"] <= 0.75
        for family, summary in family_summaries.items()
        if family in REQUIRED_FAMILIES
    ) and all(family in family_summaries for family in REQUIRED_FAMILIES)
    integrity_passed = (
        len(rows) == GUARD_TRUST_DECOUPLING_CASE_COUNT and
        not any(row.get("safety_regression", False) for row in rows) and
        all(row.get("production_artifacts_identical") is True for row in rows)
    )
    conditions = {
        "paired-production-integrity": integrity_passed,
        "overall-coverage": (
            overall["availability_ratio"] is not None and
            overall["availability_ratio"] >= 0.70),
        "family-coverage": (
            all(family in family_summaries for family in REQUIRED_FAMILIES) and
            all(
                family_summaries[family]["availability_ratio"] is not None and
                family_summaries[family]["availability_ratio"] >= 0.60
                for family in REQUIRED_FAMILIES)),
        "topology-coverage": (
            bool(topology_summaries) and all(
                summary["availability_ratio"] is not None and
                summary["availability_ratio"] >= 0.50
                for summary in topology_summaries.values())),
        "rho-diversity": (
            diversity_counts["low"] >= 100 and
            diversity_counts["mid"] >= 100 and
            diversity_counts["high"] >= 100 and
            diversity_counts["high-boundary"] >= 50 and
            family_diversity_passed),
        "rho-calibration": (
            overall["rho"]["absolute_calibration_error_median"] is not None and
            overall["rho"]["absolute_calibration_error_median"] <= 0.50 and
            overall["rho"]["absolute_calibration_error_p90"] is not None and
            overall["rho"]["absolute_calibration_error_p90"] <= 2.0 and
            family_calibration_passed),
        "action-opportunity": (
            len(action_records) > 0 and
            len(action_differences) >= 100 and
            len(action_differences) / len(action_records) >= 0.01 and
            shrink_opportunities >= 25 and
            growth_related_opportunities >= 25),
        "instrumentation-cost": (
            len(performance_ratios) == len(rows) and bool(rows) and
            performance["candidate_phase_ratio_median"] is not None and
            performance["candidate_phase_ratio_median"] <= 0.25 and
            performance["candidate_phase_ratio_p90"] is not None and
            performance["candidate_phase_ratio_p90"] <= 0.40),
    }
    failed_conditions = [
        name for name, passed in conditions.items() if not passed]
    replay_reason_by_case: dict[str, set[str]] = {}
    for record in records:
        reasons = replay_reason_by_case.setdefault(record["case-id"], set())
        bin_name = rho_bin(record)
        if bin_name in ("low", "mid"):
            reasons.add(f"rho-{bin_name}")
        if bin_name == "unavailable":
            reasons.add(f'status-{record["status"]}')
        if (record["readiness-eligible"] and
                record["shadow-action"] != "suppressed" and
                record["current-action"] != record["shadow-action"]):
            reasons.add("action-divergence")
    reason_priority = {
        "action-divergence": 0,
        "rho-low": 1,
        "rho-mid": 2,
    }
    prioritized_replays = sorted(
        (
            min((reason_priority.get(reason, 3) for reason in reasons), default=3),
            case_id,
            reasons,
        )
        for case_id, reasons in replay_reason_by_case.items()
        if reasons)
    replay_cases = [
        {"case_id": case_id, "reasons": sorted(reasons)}
        for _, case_id, reasons in prioritized_replays[:30]
    ]
    recommendation = not failed_conditions
    return {
        "diagnostic_only": True,
        **overall,
        "strata": strata,
        "candidate_funnel": {
            name: sum(
                int(row.get("trust_model_shadow", {}).get("funnel", {}).get(
                    name, 0))
                for row in rows)
            for name in (
                "generated", "invalid", "trust-skipped", "guard-rejected",
                "nonmaterial", "objective-evaluated",
                "polish-objective-evaluated")
        },
        "diversity_counts": diversity_counts,
        "action_opportunity": {
            "eligible_count": len(action_records),
            "difference_count": len(action_differences),
            "difference_ratio": (
                len(action_differences) / len(action_records)
                if action_records else None),
            "shrink_count": shrink_opportunities,
            "growth_related_count": growth_related_opportunities,
        },
        "performance": performance,
        "evidence_gate": {
            "passed": recommendation,
            "conditions": conditions,
            "failed_conditions": failed_conditions,
        },
        "priority_replay_cases": replay_cases,
        "model_based_controller_experiment_recommended": recommendation,
        "production_promotion_recommended": False,
    }


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
        if experiment["exposures"].get(name, False)
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


def classify_shadow_policy_outcome(
    policy_report: dict[str, Any],
    terminal: dict[str, Any] | None,
) -> dict[str, Any]:
    if not policy_report.get("reached", False) or terminal is None:
        return {"reached": False, "effective_exposure": False}
    attempts_saved = int(policy_report.get("attempts_saved", 0))
    effective_exposure = (
        policy_report.get("source") == "shadow-policy" and attempts_saved > 0)
    candidate_objective = policy_report.get("objective")
    terminal_objective = terminal.get("objective")
    objective_delta = (
        candidate_objective - terminal_objective
        if candidate_objective is not None and terminal_objective is not None
        else None)
    objective_harm = _materially_higher(
        candidate_objective, terminal_objective)
    objective_benefit = _materially_lower(
        candidate_objective, terminal_objective)
    candidate_truth = policy_report.get("truth_metrics") or {}
    terminal_truth = terminal.get("truth_metrics") or {}
    truth_deltas = {}
    truth_harm = {}
    truth_benefit = {}
    for metric in TRUTH_RMSE_METRICS:
        candidate_value = candidate_truth.get(metric)
        terminal_value = terminal_truth.get(metric)
        truth_deltas[metric] = (
            candidate_value - terminal_value
            if candidate_value is not None and terminal_value is not None
            else None)
        truth_harm[metric] = _materially_higher(candidate_value, terminal_value)
        truth_benefit[metric] = _materially_lower(candidate_value, terminal_value)
    comparison_complete = (
        candidate_objective is not None and terminal_objective is not None and
        all(
            candidate_truth.get(metric) is not None and
            terminal_truth.get(metric) is not None
            for metric in TRUTH_RMSE_METRICS)
    )
    return {
        "reached": True,
        "source": policy_report.get("source"),
        "effective_exposure": effective_exposure,
        "attempts_saved": attempts_saved,
        "accepted_iterations_saved": int(
            policy_report.get("accepted_iterations_saved", 0)),
        "objective_delta": objective_delta,
        "objective_harm": objective_harm,
        "objective_benefit": objective_benefit,
        "truth_rmse_delta": truth_deltas,
        "truth_rmse_harm": truth_harm,
        "truth_rmse_benefit": truth_benefit,
        "endpoint_safe": bool(policy_report.get("endpoint_safe", False)),
        "comparison_complete": comparison_complete,
        "continuation_safety_events": policy_report.get(
            "continuation_safety_events", {}),
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
            "production", "legacy_population", "solver_qualified",
            "fixed_point_operator")
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
        trust_model_shadow = audit.get(
            "trust_model_shadow", {"diagnostic_only": True, "records": []})
        terminal = audit.get("terminal")
        shadow_policy_outcomes = {
            policy: classify_shadow_policy_outcome(
                audit.get("shadow_policies", {}).get(policy, {}), terminal)
            for policy in SHADOW_POLICIES
        }
        accepted_only_shadow_count += int(shadow.get("reached", False))
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
        if audit["status"] != "no_convergence_trigger":
            experiment = audit["experiments"][0]
            exposures = {
                name: bool(experiment["exposures"].get(name, False))
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
            legacy_maximum_unique.get(
                "guarded_max", legacy_maximum_unique.get("raw_max", 0)))
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
            "trust_model_shadow": trust_model_shadow,
            "shadow_policy_outcomes": shadow_policy_outcomes,
            "terminal": terminal,
            "production_converged": audit["status"] == "counterfactual_records",
            "safety_regression": bool(summary.get("safety_regression", False)),
            "production_artifacts_identical": summary.get(
                "production_artifacts_identical"),
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
                "promotion_candidate" if policy.startswith("fixed-point-") else
                "rollback_candidate"
            ): candidate,
        }

    shadow_decisions: dict[str, Any] = {}
    for policy in SHADOW_POLICIES:
        reached_rows = [
            row for row in rows
            if row["shadow_policy_outcomes"][policy].get("reached", False)
        ]
        exposed_rows = [
            row for row in reached_rows
            if row["shadow_policy_outcomes"][policy]["effective_exposure"]
        ]
        family_counts = Counter(row["family"] for row in exposed_rows)
        topology_counts = Counter(row["topology"] for row in exposed_rows)
        attempts_saved = [
            row["shadow_policy_outcomes"][policy]["attempts_saved"]
            for row in exposed_rows
        ]
        accepted_iterations_saved = [
            row["shadow_policy_outcomes"][policy]["accepted_iterations_saved"]
            for row in exposed_rows
        ]
        objective_deltas = [
            row["shadow_policy_outcomes"][policy]["objective_delta"]
            for row in exposed_rows
            if row["shadow_policy_outcomes"][policy]["objective_delta"] is not None
        ]
        objective_harm_count = sum(
            row["shadow_policy_outcomes"][policy]["objective_harm"]
            for row in exposed_rows)
        objective_benefit_count = sum(
            row["shadow_policy_outcomes"][policy]["objective_benefit"]
            for row in exposed_rows)
        truth_harm_counts = {
            metric: sum(
                row["shadow_policy_outcomes"][policy]["truth_rmse_harm"][metric]
                for row in exposed_rows)
            for metric in TRUTH_RMSE_METRICS
        }
        truth_benefit_counts = {
            metric: sum(
                row["shadow_policy_outcomes"][policy]["truth_rmse_benefit"][metric]
                for row in exposed_rows)
            for metric in TRUTH_RMSE_METRICS
        }
        truth_harm_case_count = sum(
            any(row["shadow_policy_outcomes"][policy]["truth_rmse_harm"].values())
            for row in exposed_rows)
        endpoint_safety_violation_count = sum(
            not row["shadow_policy_outcomes"][policy]["endpoint_safe"]
            for row in exposed_rows)
        incomplete_comparison_count = sum(
            not row["shadow_policy_outcomes"][policy]["comparison_complete"]
            for row in exposed_rows)
        continuation_safety_event_counts: Counter[str] = Counter()
        for row in exposed_rows:
            continuation_safety_event_counts.update(
                row["shadow_policy_outcomes"][policy]
                ["continuation_safety_events"])

        def summarize_shadow_group(
            group_rows: list[dict[str, Any]],
        ) -> dict[str, Any]:
            group_attempts = [
                row["shadow_policy_outcomes"][policy]["attempts_saved"]
                for row in group_rows]
            group_objective_deltas = [
                row["shadow_policy_outcomes"][policy]["objective_delta"]
                for row in group_rows
                if row["shadow_policy_outcomes"][policy]["objective_delta"]
                is not None]
            group_truth_deltas = {
                metric: [
                    row["shadow_policy_outcomes"][policy]
                    ["truth_rmse_delta"][metric]
                    for row in group_rows
                    if row["shadow_policy_outcomes"][policy]
                    ["truth_rmse_delta"][metric] is not None]
                for metric in TRUTH_RMSE_METRICS
            }
            return {
                "effective_exposure_count": len(group_rows),
                "attempts_saved_total": sum(group_attempts),
                "attempts_saved_median": (
                    statistics.median(group_attempts) if group_attempts else None),
                "objective_delta_median": (
                    statistics.median(group_objective_deltas)
                    if group_objective_deltas else None),
                "objective_delta_p90": _percentile(
                    group_objective_deltas, 0.90),
                "truth_rmse_delta": {
                    metric: {
                        "median": statistics.median(values) if values else None,
                        "p90": _percentile(values, 0.90),
                    }
                    for metric, values in group_truth_deltas.items()
                },
                "objective_material_harm_count": sum(
                    row["shadow_policy_outcomes"][policy]["objective_harm"]
                    for row in group_rows),
                "truth_material_harm_counts": {
                    metric: sum(
                        row["shadow_policy_outcomes"][policy]
                        ["truth_rmse_harm"][metric]
                        for row in group_rows)
                    for metric in TRUTH_RMSE_METRICS
                },
                "endpoint_safety_violation_count": sum(
                    not row["shadow_policy_outcomes"][policy]["endpoint_safe"]
                    for row in group_rows),
            }
        exposure_quota_met = (
            len(exposed_rows) >= MINIMUM_TOTAL_EXPOSURES and
            all(
                family_counts[family] >= MINIMUM_EXPOSURES_PER_FAMILY
                for family in SHADOW_REQUIRED_FAMILIES)
        )
        promotion_candidate = (
            exposure_quota_met and
            endpoint_safety_violation_count == 0 and
            incomplete_comparison_count == 0 and
            objective_harm_count == 0 and
            not any(truth_harm_counts.values())
        )
        shadow_decisions[policy] = {
            "reached_count": len(reached_rows),
            "production_agreement_count": sum(
                row["shadow_policy_outcomes"][policy].get("source") ==
                "production-agreement" for row in reached_rows),
            "effective_exposure_count": len(exposed_rows),
            "unreached_count": sum(
                row["accepted_only_shadow"].get("reached", False)
                for row in rows) - len(reached_rows),
            "family_counts": dict(sorted(family_counts.items())),
            "topology_counts": dict(sorted(topology_counts.items())),
            "attempts_saved": {
                "median": statistics.median(attempts_saved)
                if attempts_saved else None,
                "p90": _percentile(attempts_saved, 0.90),
                "maximum": max(attempts_saved) if attempts_saved else None,
                "total": sum(attempts_saved),
            },
            "accepted_iterations_saved": {
                "median": statistics.median(accepted_iterations_saved)
                if accepted_iterations_saved else None,
                "p90": _percentile(accepted_iterations_saved, 0.90),
                "maximum": max(accepted_iterations_saved)
                if accepted_iterations_saved else None,
                "total": sum(accepted_iterations_saved),
            },
            "objective_delta": {
                "median": statistics.median(objective_deltas)
                if objective_deltas else None,
                "p90": _percentile(objective_deltas, 0.90),
            },
            "truth_rmse_delta": {
                metric: {
                    "median": statistics.median(values) if values else None,
                    "p90": _percentile(values, 0.90),
                }
                for metric, values in {
                    metric: [
                        row["shadow_policy_outcomes"][policy]
                        ["truth_rmse_delta"][metric]
                        for row in exposed_rows
                        if row["shadow_policy_outcomes"][policy]
                        ["truth_rmse_delta"][metric] is not None]
                    for metric in TRUTH_RMSE_METRICS
                }.items()
            },
            "objective_material_benefit_count": objective_benefit_count,
            "objective_material_harm_count": objective_harm_count,
            "truth_material_benefit_counts": truth_benefit_counts,
            "truth_material_harm_counts": truth_harm_counts,
            "truth_material_harm_case_count": truth_harm_case_count,
            "endpoint_safety_violation_count": endpoint_safety_violation_count,
            "incomplete_comparison_count": incomplete_comparison_count,
            "continuation_safety_event_counts": dict(
                sorted(continuation_safety_event_counts.items())),
            "exposure_quota_met": exposure_quota_met,
            "promotion_candidate": promotion_candidate,
            "by_family": {
                family: summarize_shadow_group([
                    row for row in exposed_rows if row["family"] == family])
                for family in sorted({row["family"] for row in exposed_rows})
            },
            "by_topology": {
                topology: summarize_shadow_group([
                    row for row in exposed_rows if row["topology"] == topology])
                for topology in sorted({row["topology"] for row in exposed_rows})
            },
        }

    promotable = [
        policy for policy in SHADOW_POLICIES
        if shadow_decisions[policy]["promotion_candidate"]
    ]
    recommended_shadow_policy = None
    if promotable:
        maximum_total = max(
            shadow_decisions[policy]["attempts_saved"]["total"]
            for policy in promotable)
        near_maximum = [
            policy for policy in promotable
            if shadow_decisions[policy]["attempts_saved"]["total"] >=
            0.95 * maximum_total
        ]
        if "dynamic-raw" in near_maximum:
            recommended_shadow_policy = "dynamic-raw"
        else:
            recommended_shadow_policy = next(
                policy for policy in SHADOW_POLICIES if policy in near_maximum)

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
        "schema_version": 6,
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
        "shadow_policy_decisions": shadow_decisions,
        "recommended_shadow_policy": recommended_shadow_policy,
        "shadow_policy_production_change_recommended": (
            recommended_shadow_policy is not None),
        "trust_model_shadow": summarize_trust_model_shadow(rows),
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
    before_objectives: list[float] = []
    after_objectives: list[float] = []
    before_accepted_iterations: list[int] = []
    after_accepted_iterations: list[int] = []
    accepted_iteration_deltas: list[int] = []
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
            before_objectives.append(baseline_objective)
            after_objectives.append(candidate_objective)
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
        baseline_accepted_iteration = baseline_terminal.get("accepted_iteration")
        candidate_accepted_iteration = candidate_terminal.get("accepted_iteration")
        accepted_iteration_delta = None
        if (baseline_accepted_iteration is not None and
                candidate_accepted_iteration is not None):
            baseline_accepted_iteration = int(baseline_accepted_iteration)
            candidate_accepted_iteration = int(candidate_accepted_iteration)
            accepted_iteration_delta = (
                candidate_accepted_iteration - baseline_accepted_iteration)
            before_accepted_iterations.append(baseline_accepted_iteration)
            after_accepted_iterations.append(candidate_accepted_iteration)
            accepted_iteration_deltas.append(accepted_iteration_delta)
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
            "before_accepted_iteration": baseline_accepted_iteration,
            "after_accepted_iteration": candidate_accepted_iteration,
            "accepted_iteration_delta": accepted_iteration_delta,
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
        group_accepted_iteration_deltas = [
            row["accepted_iteration_delta"] for row in group_rows
            if row["accepted_iteration_delta"] is not None]
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
            "accepted_iteration_delta_median": (
                statistics.median(group_accepted_iteration_deltas)
                if group_accepted_iteration_deltas else None),
            "accepted_iteration_delta_p90": _percentile(
                group_accepted_iteration_deltas, 0.90),
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

    before_objective_median = (
        statistics.median(before_objectives) if before_objectives else None)
    after_objective_median = (
        statistics.median(after_objectives) if after_objectives else None)
    before_accepted_iteration_median = (
        statistics.median(before_accepted_iterations)
        if before_accepted_iterations else None)
    after_accepted_iteration_median = (
        statistics.median(after_accepted_iterations)
        if after_accepted_iterations else None)
    safety_regression_count = sum(row["safety_regression"] for row in rows)
    objective_harm_ratio = (
        objective_harm_count / len(objective_deltas)
        if objective_deltas else None)
    complete_pair_count = min(
        len(objective_deltas), len(accepted_iteration_deltas))
    blocking_gate = {
        "expected_case_count": GUARD_TRUST_DECOUPLING_CASE_COUNT,
        "complete_pair_count": complete_pair_count,
        "complete_corpus": (
            len(rows) == GUARD_TRUST_DECOUPLING_CASE_COUNT and
            complete_pair_count == GUARD_TRUST_DECOUPLING_CASE_COUNT),
        "safety_regression_free": safety_regression_count == 0,
        "objective_median_regression": _materially_higher(
            after_objective_median, before_objective_median),
        "objective_harm_ratio": objective_harm_ratio,
        "objective_harm_ratio_passed": (
            objective_harm_ratio is not None and
            objective_harm_ratio <= MAXIMUM_HARM_RATIO),
        "accepted_iteration_median_regression": (
            before_accepted_iteration_median is None or
            after_accepted_iteration_median is None or
            after_accepted_iteration_median > before_accepted_iteration_median),
    }
    blocking_gate["passed"] = (
        blocking_gate["complete_corpus"] and
        blocking_gate["safety_regression_free"] and
        not blocking_gate["objective_median_regression"] and
        blocking_gate["objective_harm_ratio_passed"] and
        not blocking_gate["accepted_iteration_median_regression"])

    return {
        "schema_version": 2,
        "case_count": len(rows),
        "before_convergence_count": before.get("production_convergence_count", 0),
        "after_convergence_count": after.get("production_convergence_count", 0),
        "accepted_only_shadow_count": after.get("accepted_only_shadow_count", 0),
        "before_stop_reasons": before.get("termination_counts", {}),
        "after_stop_reasons": after.get("termination_counts", {}),
        "objective_delta": {
            "before_median": before_objective_median,
            "after_median": after_objective_median,
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
        "accepted_iteration_delta": {
            "before_median": before_accepted_iteration_median,
            "after_median": after_accepted_iteration_median,
            "median": (
                statistics.median(accepted_iteration_deltas)
                if accepted_iteration_deltas else None),
            "p90": _percentile(accepted_iteration_deltas, 0.90),
        },
        "safety_regression_count": safety_regression_count,
        "guard_trust_decoupling_blocking_gate": blocking_gate,
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
