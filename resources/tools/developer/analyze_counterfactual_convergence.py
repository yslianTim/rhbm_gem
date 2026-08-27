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
from typing import Any


CHECKPOINT_MARKER = "Counterfactual convergence checkpoint:"
ATOM_MARKER = "Counterfactual convergence atom:"
TERMINATION_MARKER = "Counterfactual convergence termination:"
FIELD_PATTERN = re.compile(
    r"(?:^|, )(?P<name>[a-z][a-z0-9-]*)=(?P<value>[^,]+)"
)
POLICIES = (
    "production",
    "strict-current",
    "current-dof",
    "strict-dof",
    "strict-member",
)
OBJECTIVE_ABSOLUTE_TOLERANCE = 1.0e-8
OBJECTIVE_RELATIVE_TOLERANCE = 1.0e-3

CONVERGENCE_ANALYZER_PATH = Path(__file__).with_name("analyze_convergence_audit.py")
CONVERGENCE_SPEC = importlib.util.spec_from_file_location(
    "convergence_safeguard_analyzer", CONVERGENCE_ANALYZER_PATH)
assert CONVERGENCE_SPEC is not None and CONVERGENCE_SPEC.loader is not None
CONVERGENCE_ANALYZER = importlib.util.module_from_spec(CONVERGENCE_SPEC)
CONVERGENCE_SPEC.loader.exec_module(CONVERGENCE_ANALYZER)


def _fields(line: str, marker: str) -> dict[str, str] | None:
    position = line.find(marker)
    if position < 0:
        return None
    payload = line[position + len(marker):].strip()
    fields = {
        match.group("name"): match.group("value").strip().rstrip(".")
        for match in FIELD_PATTERN.finditer(payload)
    }
    if fields.get("schema") != "1":
        return None
    return fields


def parse_log(text: str) -> dict[str, Any]:
    checkpoints: list[dict[str, str]] = []
    atoms: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    terminations: list[dict[str, str]] = []
    trajectory_records: list[dict[str, str]] = []
    for line in text.splitlines():
        trajectory_record = CONVERGENCE_ANALYZER.parse_record(line)
        if trajectory_record is not None:
            trajectory_records.append(trajectory_record)
        if checkpoint := _fields(line, CHECKPOINT_MARKER):
            checkpoints.append(checkpoint)
        elif atom := _fields(line, ATOM_MARKER):
            atoms[(atom["experiment"], atom["policy"])].append(atom)
        elif termination := _fields(line, TERMINATION_MARKER):
            terminations.append(termination)
    return {
        "checkpoints": checkpoints,
        "atoms": atoms,
        "terminations": terminations,
        "trajectory_records": trajectory_records,
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


def analyze(parsed: dict[str, Any], truth: dict[int, dict[str, float]]) -> dict[str, Any]:
    checkpoints = parsed["checkpoints"]
    trajectory_audit = CONVERGENCE_ANALYZER.analyze_records(
        parsed.get("trajectory_records", []))
    if not checkpoints:
        return {
            "status": "no_convergence_trigger",
            "experiment_count": 0,
            "experiments": [],
            "trajectory_audit": trajectory_audit,
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
                "stationarity": exposed("strict-current"),
                "active_dof_population": exposed("current-dof"),
                "combined": exposed("strict-dof"),
                "active_member_diagnostic": exposed("strict-member"),
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
        "stationarity", "active_dof_population", "combined",
        "active_member_diagnostic")
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
            name for name in exposure_names[:3] if report["exposures"][name]
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
        "trajectory_audit": trajectory_audit,
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
