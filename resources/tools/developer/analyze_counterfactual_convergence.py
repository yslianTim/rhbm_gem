#!/usr/bin/env python3
"""Aggregate developer-only counterfactual convergence continuation records."""

from __future__ import annotations

import argparse
from collections import defaultdict
import importlib.util
import json
import math
from pathlib import Path
import re
import statistics
from typing import Any


CHECKPOINT_MARKER = "Counterfactual convergence checkpoint:"
ATOM_MARKER = "Counterfactual convergence atom:"
TERMINATION_MARKER = "Counterfactual convergence termination:"
SHADOW_CHECKPOINT_MARKER = "Accepted-only shadow checkpoint:"
SHADOW_ATOM_MARKER = "Accepted-only shadow atom:"
AUDIT_TERMINAL_MARKER = "Second-stage audit terminal:"
AUDIT_TERMINAL_ATOM_MARKER = "Second-stage audit terminal atom:"
TRUST_MODEL_SHADOW_MARKER = "Trust-model shadow:"
TRUST_MODEL_FUNNEL_MARKER = "Trust-model funnel:"
TRUST_MODEL_PERFORMANCE_MARKER = "Trust-model performance:"
FIELD_PATTERN = re.compile(
    r"(?:^|, )(?P<name>[a-z][a-z0-9-]*)=(?P<value>[^,]+)"
)
POLICIES = (
    "production",
    "historical-all-selected",
    "historical-cluster-active-proposal-maximum",
    "historical-active-proposal",
    "production-maximum",
)
LEGACY_POLICY_MAP = {
    "production": "production",
    "legacy-population": "historical-all-selected",
    "legacy-maximum": "historical-cluster-active-proposal-maximum",
    "solver-qualified": "historical-active-proposal",
    "fixed-point-operator-maximum": "production-maximum",
}
SHADOW_POLICIES = (
    "accepted-only-k2",
    "accepted-only-k3",
    "accepted-only-k5",
    "dynamic-raw",
)
OBJECTIVE_ABSOLUTE_TOLERANCE = 1.0e-8
OBJECTIVE_RELATIVE_TOLERANCE = 1.0e-3

CONVERGENCE_ANALYZER_PATH = Path(__file__).with_name("analyze_convergence_audit.py")
CONVERGENCE_SPEC = importlib.util.spec_from_file_location(
    "convergence_safeguard_analyzer", CONVERGENCE_ANALYZER_PATH)
assert CONVERGENCE_SPEC is not None and CONVERGENCE_SPEC.loader is not None
CONVERGENCE_ANALYZER = importlib.util.module_from_spec(CONVERGENCE_SPEC)
CONVERGENCE_SPEC.loader.exec_module(CONVERGENCE_ANALYZER)


def _fields(
    line: str,
    marker: str,
    schema: str = "3",
) -> dict[str, str] | None:
    position = line.find(marker)
    if position < 0:
        return None
    payload = line[position + len(marker):].strip()
    fields = {
        match.group("name"): match.group("value").strip().rstrip(".")
        for match in FIELD_PATTERN.finditer(payload)
    }
    if fields.get("schema") != schema:
        return None
    return fields


def _counterfactual_fields(line: str, marker: str) -> dict[str, str] | None:
    if fields := _fields(line, marker, "4"):
        if fields.get("comparator-set") != "1":
            raise ValueError("Unsupported counterfactual comparator set")
        required_by_marker = {
            CHECKPOINT_MARKER: {
                "experiment", "policy", "try", "acc", "extra-try",
                "extra-acc", "extra-ms", "final-polish",
                "solver-qualified", "restricted", "all-fixed", "population",
                "objective", "accepted-median", "accepted-p99",
                "accepted-max", "raw-median", "raw-p99", "raw-max",
            },
            ATOM_MARKER: {
                "experiment", "policy", "serial", "group", "shape-active",
                "offset-active", "quarantined", "amplitude", "width", "offset",
            },
            TERMINATION_MARKER: {
                "experiment", "reason", "try", "acc", "extra-try",
                "extra-acc", "extra-ms", "checkpoints",
            },
        }
        required = required_by_marker[marker]
        if not required.issubset(fields):
            missing = ", ".join(sorted(required - fields.keys()))
            raise ValueError(
                f"Counterfactual schema-4 record is missing: {missing}")
        policy = fields.get("policy")
        if policy is not None and policy not in POLICIES:
            raise ValueError(f"Unknown counterfactual policy: {policy}")
        if "checkpoints" in fields and len(fields["checkpoints"].split("/")) != 5:
            raise ValueError("Checkpoint vector must contain five values")
        return fields
    fields = _fields(line, marker, "3")
    if fields is None:
        return None
    policy = fields.get("policy")
    if policy == "fixed-point-operator":
        return None
    if policy is not None:
        if policy not in LEGACY_POLICY_MAP:
            raise ValueError(f"Unknown legacy counterfactual policy: {policy}")
        fields["policy"] = LEGACY_POLICY_MAP[policy]
    if "checkpoints" in fields:
        checkpoints = fields["checkpoints"].split("/")
        if len(checkpoints) != 6:
            raise ValueError("Legacy checkpoint vector must contain six values")
        fields["checkpoints"] = "/".join(
            checkpoints[index] for index in (0, 1, 2, 3, 5))
    fields["source-schema"] = "3"
    fields["comparator-set"] = "1"
    return fields


def parse_log(text: str) -> dict[str, Any]:
    checkpoints: list[dict[str, str]] = []
    atoms: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    terminations: list[dict[str, str]] = []
    shadow_checkpoint: dict[str, str] | None = None
    shadow_atoms: list[dict[str, str]] = []
    shadow_policy_checkpoints: list[dict[str, str]] = []
    shadow_policy_atoms: dict[str, list[dict[str, str]]] = defaultdict(list)
    audit_terminal: dict[str, str] | None = None
    audit_terminal_atoms: list[dict[str, str]] = []
    trajectory_records: list[dict[str, str]] = []
    trust_model_shadow_records: list[dict[str, Any]] = []
    trust_model_funnel_records: list[dict[str, Any]] = []
    trust_model_performance: dict[str, float] | None = None
    for line in text.splitlines():
        trajectory_record = CONVERGENCE_ANALYZER.parse_record(line)
        if trajectory_record is not None:
            trajectory_records.append(trajectory_record)
        if trust_model := _fields(line, TRUST_MODEL_SHADOW_MARKER, "2"):
            optional_float_fields = (
                "actual-reduction", "polish-reduction",
                "predicted-residual-reduction",
                "predicted-penalty-reduction", "predicted-reduction", "rho")
            trust_model_shadow_records.append({
                **trust_model,
                **{
                    name: (
                        None if trust_model[name] == "-" else
                        float(trust_model[name]))
                    for name in optional_float_fields
                },
                **{
                    name: int(trust_model[name])
                    for name in (
                        "try", "acc", "atoms", "key-first", "key-last",
                        "boundary-touched", "boundary-rescued",
                        "readiness-eligible", "final-local-candidate",
                        "search-pass", "trial", "rejected-by-previous",
                        "rejected-by-best", "rejected-by-strict-polish",
                        "objective-backtracked",
                        "unselected-dependencies")
                },
                "factor": float(trust_model["factor"]),
                "step-norm": float(trust_model["step-norm"]),
                "boundary-utilization": float(
                    trust_model["boundary-utilization"]),
                "elapsed-ms": float(trust_model["elapsed-ms"]),
            })
        if funnel := _fields(line, TRUST_MODEL_FUNNEL_MARKER, "1"):
            trust_model_funnel_records.append({
                **funnel,
                **{
                    name: int(funnel[name])
                    for name in (
                        "try", "acc", "atoms", "key-first", "key-last",
                        "generated", "invalid", "trust-skipped",
                        "guard-rejected", "nonmaterial", "objective-evaluated",
                        "polish-objective-evaluated")
                },
            })
        if performance := _fields(
                line, TRUST_MODEL_PERFORMANCE_MARKER, "1"):
            trust_model_performance = {
                "candidate-ms": float(performance["candidate-ms"]),
                "total-ms": float(performance["total-ms"]),
            }
        if checkpoint := _counterfactual_fields(line, CHECKPOINT_MARKER):
            checkpoints.append(checkpoint)
        elif atom := _counterfactual_fields(line, ATOM_MARKER):
            atoms[(atom["experiment"], atom["policy"])].append(atom)
        elif termination := _counterfactual_fields(line, TERMINATION_MARKER):
            terminations.append(termination)
        elif shadow := _fields(line, SHADOW_CHECKPOINT_MARKER, "1"):
            shadow_checkpoint = shadow
        elif shadow_policy := _fields(line, SHADOW_CHECKPOINT_MARKER, "2"):
            shadow_policy_checkpoints.append(shadow_policy)
        elif shadow_atom := _fields(line, SHADOW_ATOM_MARKER, "1"):
            shadow_atoms.append(shadow_atom)
        elif shadow_policy_atom := _fields(line, SHADOW_ATOM_MARKER, "2"):
            shadow_policy_atoms[shadow_policy_atom["policy"]].append(
                shadow_policy_atom)
        elif terminal := _fields(line, AUDIT_TERMINAL_MARKER, "1"):
            audit_terminal = terminal
        elif terminal_atom := _fields(line, AUDIT_TERMINAL_ATOM_MARKER, "1"):
            audit_terminal_atoms.append(terminal_atom)
    return {
        "checkpoints": checkpoints,
        "atoms": atoms,
        "terminations": terminations,
        "shadow_checkpoint": shadow_checkpoint,
        "shadow_atoms": shadow_atoms,
        "shadow_policy_checkpoints": shadow_policy_checkpoints,
        "shadow_policy_atoms": shadow_policy_atoms,
        "audit_terminal": audit_terminal,
        "audit_terminal_atoms": audit_terminal_atoms,
        "trajectory_records": trajectory_records,
        "trust_model_shadow_records": trust_model_shadow_records,
        "trust_model_funnel_records": trust_model_funnel_records,
        "trust_model_performance": trust_model_performance,
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


def analyze_trust_model_shadow(
    records: list[dict[str, Any]],
    funnel_records: list[dict[str, Any]],
    performance: dict[str, float] | None,
) -> dict[str, Any]:
    status_counts: dict[str, int] = defaultdict(int)
    rho_bins: dict[str, int] = defaultdict(int)
    action_confusion: dict[str, int] = defaultdict(int)
    strata: dict[str, dict[str, dict[str, int]]] = {
        name: defaultdict(lambda: defaultdict(int))
        for name in (
            "source", "trial-disposition", "factor", "prediction-status",
            "boundary", "cluster-size", "unselected-dependencies")
    }
    rho_values: list[float] = []
    for record in records:
        status = str(record["status"])
        status_counts[status] += 1
        rho = record.get("rho")
        if rho is not None:
            rho_values.append(float(rho))
            if rho < 0.25:
                rho_bin = "low"
            elif rho <= 0.75:
                rho_bin = "mid"
            elif float(record["boundary-utilization"]) >= 0.8:
                rho_bin = "high-boundary"
            else:
                rho_bin = "high-interior"
            rho_bins[rho_bin] += 1
        else:
            rho_bin = "unavailable"
            rho_bins[rho_bin] += 1
        if record["readiness-eligible"] and record["shadow-action"] != "suppressed":
            action_confusion[
                f'{record["current-action"]}->{record["shadow-action"]}'] += 1
        boundary = (
            "rescued" if record["boundary-rescued"] else
            "touched" if record["boundary-touched"] else "none")
        unselected_count = int(record["unselected-dependencies"])
        unselected_stratum = (
            "0" if unselected_count == 0 else
            "1" if unselected_count == 1 else "2+")
        values = {
            "source": str(record["source"]),
            "trial-disposition": str(record["trial-disposition"]),
            "factor": f'{float(record["factor"]):.8g}',
            "prediction-status": status,
            "boundary": boundary,
            "cluster-size": str(record["atoms"]),
            "unselected-dependencies": unselected_stratum,
        }
        for stratum_name, stratum_value in values.items():
            group = strata[stratum_name][stratum_value]
            group["records"] += 1
            group[f"status:{status}"] += 1
            group[f"rho:{rho_bin}"] += 1
    eligible_count = sum(
        record["status"] != "nonmaterial-step" for record in records)
    elapsed_milliseconds = sum(
        float(record["elapsed-ms"]) for record in records)
    candidate_milliseconds = (
        float(performance["candidate-ms"]) if performance is not None else None)
    funnel_totals = {
        name: sum(int(record[name]) for record in funnel_records)
        for name in (
            "generated", "invalid", "trust-skipped", "guard-rejected",
            "nonmaterial", "objective-evaluated",
            "polish-objective-evaluated")
    }
    calibration_errors = [abs(value - 1.0) for value in rho_values]
    return {
        "diagnostic_only": True,
        "record_count": len(records),
        "eligible_count": eligible_count,
        "available_count": status_counts.get("available", 0),
        "availability_ratio": (
            status_counts.get("available", 0) / eligible_count
            if eligible_count else None),
        "status_counts": dict(sorted(status_counts.items())),
        "rho_bins": dict(sorted(rho_bins.items())),
        "action_confusion": dict(sorted(action_confusion.items())),
        "rho": {
            "median": statistics.median(rho_values) if rho_values else None,
            "p90": _percentile(rho_values, 0.90),
            "absolute_calibration_error_median": (
                statistics.median(calibration_errors)
                if calibration_errors else None),
            "absolute_calibration_error_p90": _percentile(
                calibration_errors, 0.90),
        },
        "elapsed_milliseconds": elapsed_milliseconds,
        "candidate_milliseconds": candidate_milliseconds,
        "candidate_phase_ratio": (
            elapsed_milliseconds / candidate_milliseconds
            if candidate_milliseconds is not None and candidate_milliseconds > 0.0
            else None),
        "funnel": funnel_totals,
        "strata": {
            name: {
                value: dict(sorted(counts.items()))
                for value, counts in sorted(groups.items())
            }
            for name, groups in strata.items()
        },
        "funnels": funnel_records,
        "records": records,
    }


def _numbers(value: str) -> list[float]:
    return [float(item) for item in value.split("/")]


def _objective_total(value: str) -> float | None:
    if value == "-/-/-/-":
        return None
    numbers = _numbers(value)
    return numbers[-1] if numbers and math.isfinite(numbers[-1]) else None


def _load_truth(path: Path | None) -> dict[int, dict[str, float]]:
    if path is None:
        return {}
    value = json.loads(path.read_text(encoding="utf-8"))
    rows = value.get("atoms", value) if isinstance(value, dict) else value
    if isinstance(rows, dict):
        rows = [dict(parameters, serial=int(serial)) for serial, parameters in rows.items()]
    truth: dict[int, dict[str, float]] = {}
    for row in rows:
        serial = int(row["serial"])
        truth[serial] = {
            name: float(row[name]) for name in ("amplitude", "width", "offset")
        }
    return truth


def _truth_metrics(
    atom_rows: list[dict[str, str]],
    truth: dict[int, dict[str, float]],
) -> dict[str, float] | None:
    if not truth:
        return None
    squared = {name: 0.0 for name in ("amplitude", "width", "offset")}
    transformed_squared = [0.0, 0.0, 0.0]
    maximum_absolute_offset_error = 0.0
    seen: set[int] = set()
    for row in atom_rows:
        serial = int(row["serial"])
        if serial not in truth:
            continue
        seen.add(serial)
        for name in squared:
            error = float(row[name]) - truth[serial][name]
            squared[name] += error * error
            if name == "offset":
                maximum_absolute_offset_error = max(
                    maximum_absolute_offset_error, abs(error))
        fitted_amplitude = float(row["amplitude"])
        fitted_width = float(row["width"])
        truth_amplitude = truth[serial]["amplitude"]
        truth_width = truth[serial]["width"]
        if min(fitted_amplitude, fitted_width, truth_amplitude, truth_width) <= 0.0:
            return None
        fitted_height = fitted_amplitude / (
            2.0 * math.pi * fitted_width * fitted_width) ** 1.5
        truth_height = truth_amplitude / (
            2.0 * math.pi * truth_width * truth_width) ** 1.5
        transformed_error = (
            math.log(fitted_height) - math.log(truth_height),
            math.log(fitted_width) - math.log(truth_width),
            float(row["offset"]) / fitted_height -
            truth[serial]["offset"] / truth_height,
        )
        if not all(math.isfinite(value) for value in transformed_error):
            return None
        for index, error in enumerate(transformed_error):
            transformed_squared[index] += error * error
    if seen != set(truth):
        return None
    count = len(seen)
    result = {
        f"{name}_rmse": math.sqrt(total / count)
        for name, total in squared.items()
    }
    transformed_names = ("log_peak", "log_width", "offset_over_peak")
    transformed_rmse = [
        math.sqrt(total / count) for total in transformed_squared
    ]
    result.update({
        f"transformed_{name}_rmse": value
        for name, value in zip(transformed_names, transformed_rmse)
    })
    result["transformed_aggregate_rmse"] = math.sqrt(
        sum(transformed_squared) / (3 * count))
    result["maximum_absolute_offset_error"] = maximum_absolute_offset_error
    return result


def _continuation_safety_events(
    records: list[dict[str, str]],
    checkpoint_try: int,
) -> dict[str, int]:
    counts = defaultdict(int)
    for record in records:
        if int(record["try"]) <= checkpoint_try:
            continue
        qualification = [
            int(float(value)) for value in record.get("qualification", "").split("/")
            if value]
        blockers = [
            int(float(value.rstrip(".")))
            for value in record.get("blockers", "").split("/") if value]
        if len(qualification) >= 17:
            counts["hard_failure"] += qualification[7] + qualification[13]
            counts["mixed_activity"] += qualification[16]
        if len(blockers) >= 4:
            counts["suspicious"] += blockers[0]
            counts["rejection"] += blockers[1]
            counts["quarantine_transition"] += blockers[2]
            counts["domain_change"] += blockers[3]
        summary_fields = (
            "accepted-active-p99", "accepted-active-max",
            "operator-nominal-residual-p99", "operator-nominal-residual-max")
        if any(
                not math.isfinite(value)
                for field in summary_fields
                for value in _numbers(record.get(field, "inf"))):
            counts["nonfinite"] += 1
    return {
        name: counts.get(name, 0) for name in (
            "nonfinite", "hard_failure", "mixed_activity", "suspicious",
            "rejection", "quarantine_transition", "domain_change")
    }


def analyze(parsed: dict[str, Any], truth: dict[int, dict[str, float]]) -> dict[str, Any]:
    checkpoints = parsed["checkpoints"]
    trajectory_audit = CONVERGENCE_ANALYZER.analyze_records(
        parsed.get("trajectory_records", []))
    trust_model_shadow = analyze_trust_model_shadow(
        parsed.get("trust_model_shadow_records", []),
        parsed.get("trust_model_funnel_records", []),
        parsed.get("trust_model_performance"))
    shadow_record = parsed.get("shadow_checkpoint")
    terminal_record = parsed.get("audit_terminal")
    shadow_report = {
        "reached": shadow_record is not None,
    }
    if shadow_record is not None:
        shadow_report.update({
            "try": int(shadow_record["try"]),
            "accepted_iteration": int(shadow_record["acc"]),
            "objective": _objective_total(shadow_record["objective"]),
            "accepted_p99": _numbers(shadow_record["accepted-p99"]),
            "raw_p99": _numbers(shadow_record["raw-p99"]),
            "truth_metrics": _truth_metrics(
                parsed.get("shadow_atoms", []), truth),
        })
    terminal_report = None
    if terminal_record is not None:
        terminal_report = {
            "reason": terminal_record["reason"],
            "try": int(terminal_record["try"]),
            "accepted_iteration": int(terminal_record["acc"]),
            "objective": _objective_total(terminal_record["objective"]),
            "truth_metrics": _truth_metrics(
                parsed.get("audit_terminal_atoms", []), truth),
        }
        if shadow_record is not None:
            terminal_report["attempts_after_accepted_only"] = (
                terminal_report["try"] - int(shadow_record["try"]))
            terminal_report["accepted_iterations_after_accepted_only"] = (
                terminal_report["accepted_iteration"] -
                int(shadow_record["acc"]))
    shadow_policy_records = {
        record["policy"]: record
        for record in parsed.get("shadow_policy_checkpoints", [])
    }
    shadow_policy_reports: dict[str, dict[str, Any]] = {}
    for policy in SHADOW_POLICIES:
        record = shadow_policy_records.get(policy)
        if record is not None:
            checkpoint_try = int(record["try"])
            checkpoint_acc = int(record["acc"])
            checkpoint_objective = _objective_total(record["objective"])
            endpoint_safe = record.get("checkpoint-safe") == "1"
            blockers = [
                int(float(value)) for value in record.get(
                    "blockers", "1/1/1/1").split("/")]
            endpoint_safe = endpoint_safe and not any(blockers)
            report = {
                "reached": True,
                "source": "shadow-policy",
                "try": checkpoint_try,
                "accepted_iteration": checkpoint_acc,
                "streak": int(record["streak"]),
                "effective_raw_threshold": float(record["raw-threshold"]),
                "objective": checkpoint_objective,
                "accepted_p99": _numbers(record["accepted-p99"]),
                "raw_p99": _numbers(record["raw-p99"]),
                "endpoint_safe": endpoint_safe,
                "truth_metrics": _truth_metrics(
                    parsed.get("shadow_policy_atoms", {}).get(policy, []), truth),
                "continuation_safety_events": _continuation_safety_events(
                    parsed.get("trajectory_records", []), checkpoint_try),
            }
            if terminal_report is not None:
                report["attempts_saved"] = terminal_report["try"] - checkpoint_try
                report["accepted_iterations_saved"] = (
                    terminal_report["accepted_iteration"] - checkpoint_acc)
            shadow_policy_reports[policy] = report
        elif terminal_report is not None and terminal_report["reason"] == "converged":
            shadow_policy_reports[policy] = {
                "reached": True,
                "source": "production-agreement",
                "try": terminal_report["try"],
                "accepted_iteration": terminal_report["accepted_iteration"],
                "objective": terminal_report["objective"],
                "endpoint_safe": terminal_report["objective"] is not None,
                "truth_metrics": terminal_report["truth_metrics"],
                "attempts_saved": 0,
                "accepted_iterations_saved": 0,
                "continuation_safety_events": _continuation_safety_events(
                    parsed.get("trajectory_records", []), terminal_report["try"]),
            }
        else:
            shadow_policy_reports[policy] = {"reached": False}
    if not checkpoints:
        return {
            "status": "no_convergence_trigger",
            "experiment_count": 0,
            "experiments": [],
            "accepted_only_shadow": shadow_report,
            "shadow_policies": shadow_policy_reports,
            "terminal": terminal_report,
            "trajectory_audit": trajectory_audit,
            "trust_model_shadow": trust_model_shadow,
        }

    by_experiment: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for checkpoint in checkpoints:
        by_experiment[checkpoint["experiment"]][checkpoint["policy"]] = checkpoint
    termination_by_experiment = {
        record["experiment"]: record for record in parsed["terminations"]
    }
    reports: list[dict[str, Any]] = []
    for experiment, policy_records in sorted(by_experiment.items()):
        production = policy_records.get("production")
        if production is None:
            raise ValueError(f"Experiment {experiment} has no production checkpoint")
        trigger_try = int(production["try"])
        production_objective = _objective_total(production["objective"])
        policy_reports: dict[str, Any] = {}
        for policy in POLICIES:
            record = policy_records.get(policy)
            if record is None:
                policy_reports[policy] = {"reached": False}
                continue
            objective = _objective_total(record["objective"])
            delta = (
                objective - production_objective
                if objective is not None and production_objective is not None else None)
            tolerance = (
                OBJECTIVE_ABSOLUTE_TOLERANCE +
                OBJECTIVE_RELATIVE_TOLERANCE * abs(production_objective)
                if production_objective is not None else None)
            policy_reports[policy] = {
                "reached": True,
                "try": int(record["try"]),
                "accepted_iteration": int(record["acc"]),
                "extra_attempts": int(record["extra-try"]),
                "extra_accepted_iterations": int(record["extra-acc"]),
                "elapsed_milliseconds": (
                    float(record["extra-ms"]) if "extra-ms" in record else None),
                "objective": objective,
                "objective_delta": delta,
                "objective_relative_delta": (
                    delta / abs(production_objective)
                    if delta is not None and production_objective not in (None, 0.0)
                    else None),
                "material_objective_improvement": (
                    objective < production_objective - tolerance
                    if objective is not None and production_objective is not None and
                    tolerance is not None else False),
                "accepted_median": _numbers(record["accepted-median"]),
                "accepted_p99": _numbers(record["accepted-p99"]),
                "accepted_max": _numbers(record["accepted-max"]),
                "raw_median": _numbers(record["raw-median"]),
                "raw_p99": _numbers(record["raw-p99"]),
                "raw_max": _numbers(record["raw-max"]),
                "truth_metrics": _truth_metrics(
                    parsed["atoms"].get((experiment, policy), []), truth),
            }

        def exposed(policy: str) -> bool:
            record = policy_records.get(policy)
            return record is None or int(record["try"]) > trigger_try

        termination = termination_by_experiment.get(experiment)
        reports.append({
            "experiment": experiment,
            "trigger_try": trigger_try,
            "trigger_accepted_iteration": int(production["acc"]),
            "exposures": {
                "historical-all-selected": exposed(
                    "historical-all-selected"),
                "historical-cluster-active-proposal-maximum": exposed(
                    "historical-cluster-active-proposal-maximum"),
                "historical-active-proposal": exposed(
                    "historical-active-proposal"),
                "production-maximum": exposed("production-maximum"),
            },
            "actual_continuation": any(
                report.get("extra_attempts", 0) > 0
                for report in policy_reports.values()),
            "termination_reason": (
                termination.get("reason") if termination else "unrecorded"),
            "termination_elapsed_milliseconds": (
                float(termination["extra-ms"])
                if termination and "extra-ms" in termination else None),
            "policies": policy_reports,
        })
    exposure_names = (
        "historical-all-selected",
        "historical-cluster-active-proposal-maximum",
        "historical-active-proposal",
        "production-maximum",
    )
    exposure_counts = {
        name: sum(report["exposures"][name] for report in reports)
        for name in exposure_names
    }
    exposure_overlap: dict[str, int] = defaultdict(int)
    termination_counts: dict[str, int] = defaultdict(int)
    termination_category_counts: dict[str, int] = defaultdict(int)
    material_improvement_counts = {policy: 0 for policy in POLICIES}
    unresolved_policy_counts = {policy: 0 for policy in POLICIES}
    for report in reports:
        signature = "+".join(
            name for name in exposure_names if report["exposures"][name]
        ) or "none"
        exposure_overlap[signature] += 1
        termination_counts[report["termination_reason"]] += 1
        termination_reason = report["termination_reason"]
        if termination_reason == "policy-agreement":
            termination_category = "policy_agreement"
        elif termination_reason == "all-policies-reached":
            termination_category = "all_candidate_policies_reached"
        elif termination_reason == "budget-exhausted":
            termination_category = "budget_exhaustion"
        else:
            termination_category = "other_safeguard"
        termination_category_counts[termination_category] += 1
        for policy, policy_report in report["policies"].items():
            if not policy_report["reached"]:
                unresolved_policy_counts[policy] += 1
            elif policy_report["material_objective_improvement"]:
                material_improvement_counts[policy] += 1
    return {
        "status": "counterfactual_records",
        "experiment_count": len(reports),
        "exposure_counts": exposure_counts,
        "exposure_overlap": dict(sorted(exposure_overlap.items())),
        "actual_continuation_count": sum(
            report["actual_continuation"] for report in reports),
        "termination_counts": dict(sorted(termination_counts.items())),
        "termination_category_counts": dict(
            sorted(termination_category_counts.items())),
        "material_objective_improvement_counts": material_improvement_counts,
        "unresolved_policy_counts": unresolved_policy_counts,
        "experiments": reports,
        "accepted_only_shadow": shadow_report,
        "shadow_policies": shadow_policy_reports,
        "terminal": terminal_report,
        "trajectory_audit": trajectory_audit,
        "trust_model_shadow": trust_model_shadow,
    }


def format_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Counterfactual convergence continuation audit",
        "",
        f"- Status: `{report['status']}`",
        f"- Experiments: {report['experiment_count']}",
    ]
    if report["status"] == "counterfactual_records":
        lines.extend((
            f"- Exposure counts: {json.dumps(report['exposure_counts'], sort_keys=True)}",
            f"- Exposure overlap: {json.dumps(report['exposure_overlap'], sort_keys=True)}",
            f"- Actual continuations: {report['actual_continuation_count']}",
            f"- Terminations: {json.dumps(report['termination_counts'], sort_keys=True)}",
        ))
    for experiment in report.get("experiments", []):
        lines.extend((
            "",
            f"## Experiment {experiment['experiment']}",
            "",
            f"- Termination: `{experiment['termination_reason']}`",
            f"- Actual continuation: {experiment['actual_continuation']}",
            f"- Exposures: {json.dumps(experiment['exposures'], sort_keys=True)}",
            "",
            "| Policy | Reached | Extra attempts | Objective delta | Material improvement |",
            "| --- | --- | ---: | ---: | --- |",
        ))
        for policy in POLICIES:
            value = experiment["policies"][policy]
            lines.append(
                f"| {policy} | {value['reached']} | "
                f"{value.get('extra_attempts', '-')} | "
                f"{value.get('objective_delta', '-')} | "
                f"{value.get('material_objective_improvement', '-')} |")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--format", choices=("json", "markdown"), default="markdown")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    parsed = parse_log(args.log.read_text(encoding="utf-8", errors="replace"))
    report = analyze(parsed, _load_truth(args.truth))
    output = (
        json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.format == "json" else format_markdown(report))
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    else:
        print(output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
