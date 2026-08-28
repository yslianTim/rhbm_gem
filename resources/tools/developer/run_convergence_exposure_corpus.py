#!/usr/bin/env python3
"""Run and aggregate the build-gated convergence-exposure corpus."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import importlib.util
import json
import math
from pathlib import Path
import re
import subprocess
from typing import Any, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[3]
COUNTERFACTUAL_ANALYZER_PATH = Path(__file__).with_name(
    "analyze_counterfactual_convergence.py")
CORPUS_ANALYZER_PATH = Path(__file__).with_name(
    "analyze_convergence_exposure_corpus.py")
TRUTH_MARKER = "Convergence exposure truth:"
CASE_SUMMARY_SCHEMA_VERSION = 9
FIELD_PATTERN = re.compile(
    r"(?:^|, )(?P<name>[a-z][a-z0-9-]*)=(?P<value>[^,]+)")


def _load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


COUNTERFACTUAL_ANALYZER = _load_module(
    "counterfactual_convergence_analyzer", COUNTERFACTUAL_ANALYZER_PATH)
CORPUS_ANALYZER = _load_module(
    "convergence_exposure_corpus_analyzer", CORPUS_ANALYZER_PATH)


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")


def expand_manifest(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    if manifest.get("schema_version") != 1:
        raise ValueError("Exposure manifest schema_version must be 1")
    replicas = int(manifest["replicas_per_variant"])
    seed_base = int(manifest["seed_base"])
    cases: list[dict[str, Any]] = []
    for family_index, family in enumerate(manifest["families"]):
        topologies = family["topologies"]
        levels = family["levels"]
        if len(topologies) != 5 or len(levels) != 5:
            raise ValueError("Each exposure family must define five topologies and levels")
        for topology_index, topology in enumerate(topologies):
            for level_index, level in enumerate(levels):
                variant = topology_index * len(levels) + level_index
                for replica in range(replicas):
                    family_name = str(family["name"])
                    cases.append({
                        "case_id": f"{family_name}-v{variant:02d}-r{replica}",
                        "family": family_name,
                        "family_index": family_index,
                        "variant": variant,
                        "topology": topology,
                        "topology_index": topology_index,
                        "level": level_index,
                        "replica": replica,
                        "seed": seed_base + 10000 * family_index +
                        100 * variant + replica,
                        **level,
                    })
    case_ids = [case["case_id"] for case in cases]
    if len(cases) != 600 or len(set(case_ids)) != len(case_ids):
        raise ValueError("Exposure manifest must expand to 600 unique cases")
    return cases


def build_command(
    executable: Path,
    case: dict[str, Any],
    thread_count: int,
) -> list[str]:
    command = [
        str(executable),
        "--case-id", case["case_id"],
        "--family", case["family"],
        "--topology", case["topology"],
        "--level", str(case["level"]),
        "--replica", str(case["replica"]),
        "--seed", str(case["seed"]),
        "--threads", str(thread_count),
    ]
    for name in (
        "noise_sigma", "separation", "initial_perturbation",
        "active_target", "largest_group_ratio"):
        if name in case:
            command.extend(("--" + name.replace("_", "-"), str(case[name])))
    return command


def parse_truth(log_text: str) -> dict[int, dict[str, float]]:
    truth: dict[int, dict[str, float]] = {}
    for line in log_text.splitlines():
        marker_position = line.find(TRUTH_MARKER)
        if marker_position < 0:
            continue
        payload = line[marker_position + len(TRUTH_MARKER):].strip()
        fields = {
            match.group("name"): match.group("value").strip().rstrip(".")
            for match in FIELD_PATTERN.finditer(payload)
        }
        if fields.get("schema") != "1":
            continue
        truth[int(fields["serial"])] = {
            name: float(fields[name]) for name in ("amplitude", "width", "offset")
        }
    return truth


def load_reference_truth(
    reference_truth_directory: Path | None,
    case_id: str,
) -> dict[int, dict[str, float]] | None:
    if reference_truth_directory is None:
        return None
    path = reference_truth_directory / "cases" / case_id / "scenario-truth.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    return {
        int(row["serial_id"]): {
            name: float(row[name]) for name in ("amplitude", "width", "offset")
        }
        for row in value["atoms"]
    }


def _numbers(value: str) -> list[float]:
    return [float(item) for item in value.split("/")]


def detect_safety_regression(parsed: dict[str, Any], audit: dict[str, Any]) -> bool:
    terminal = parsed.get("audit_terminal")
    if terminal is not None:
        objective = terminal.get("objective", "-/-/-/-")
        if objective == "-/-/-/-" or any(
                not math.isfinite(value) for value in _numbers(objective)):
            return True
        for atom in parsed.get("audit_terminal_atoms", []):
            amplitude = float(atom["amplitude"])
            width = float(atom["width"])
            offset = float(atom["offset"])
            if (not all(math.isfinite(value) for value in
                        (amplitude, width, offset)) or
                    amplitude <= 0.0 or width <= 0.0):
                return True
    if audit["status"] == "no_convergence_trigger":
        return False
    trigger_try = int(audit["experiments"][0]["trigger_try"])
    records = parsed.get("trajectory_records", [])
    baseline = next(
        (record for record in records if int(record["try"]) == trigger_try), None)
    if baseline is None:
        return True

    def safety_signature(record: dict[str, str]) -> tuple[int, int, int, bool]:
        qualification = [
            int(float(value)) for value in record["qualification"].split("/")]
        blockers = [
            int(float(value)) for value in record["blockers"].split("/")]
        summary_fields = (
            "production-accepted-p99", "production-accepted-max",
            "legacy-accepted-p99", "legacy-accepted-max",
            "fixed-point-residual-p99", "fixed-point-residual-max")
        nonfinite = any(
            not math.isfinite(value)
            for field in summary_fields for value in _numbers(record[field]))
        return (
            qualification[7] + qualification[13],
            qualification[16],
            blockers[2],
            nonfinite,
        )

    baseline_signature = safety_signature(baseline)
    for record in records:
        if int(record["try"]) <= trigger_try:
            continue
        signature = safety_signature(record)
        if any(signature[index] > baseline_signature[index] for index in range(3)):
            return True
        if signature[3] and not baseline_signature[3]:
            return True
    return False


def run_case(
    executable: Path,
    case: dict[str, Any],
    output_directory: Path,
    thread_count: int,
    reference_truth_directory: Path | None,
) -> dict[str, Any]:
    case_directory = output_directory / "cases" / case["case_id"]
    case_directory.mkdir(parents=True, exist_ok=True)
    summary_path = case_directory / "case-summary.json"
    if summary_path.is_file():
        value = json.loads(summary_path.read_text(encoding="utf-8"))
        trajectory_path = case_directory / "trajectory-schema-7.json"
        trajectory = (
            json.loads(trajectory_path.read_text(encoding="utf-8"))
            if trajectory_path.is_file() else {})
        trust_model_path = case_directory / "trust-model-shadow-schema-1.json"
        trust_model = (
            json.loads(trust_model_path.read_text(encoding="utf-8"))
            if trust_model_path.is_file() else {})
        if (value.get("schema_version") == CASE_SUMMARY_SCHEMA_VERSION and
                value.get("case") == case and
                value.get("thread_count") == thread_count and
                value.get("reference_truth_directory") == (
                    str(reference_truth_directory)
                    if reference_truth_directory is not None else None) and
                value.get("status") == "complete" and
                trajectory.get("schema_version") == 7 and
                trust_model.get("schema_version") == 1):
            return value

    command = build_command(executable, case, thread_count)
    try:
        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
    except OSError as error:
        summary = {
            "schema_version": CASE_SUMMARY_SCHEMA_VERSION,
            "status": "failed",
            "case": case,
            "thread_count": thread_count,
            "return_code": None,
            "command": command,
            "error": str(error),
        }
        write_json(summary_path, summary)
        return summary
    log_text = completed.stdout.decode("utf-8", errors="replace")
    (case_directory / "run.log").write_text(log_text, encoding="utf-8")
    if completed.returncode != 0:
        summary = {
            "schema_version": CASE_SUMMARY_SCHEMA_VERSION,
            "status": "failed",
            "case": case,
            "thread_count": thread_count,
            "return_code": completed.returncode,
            "command": command,
        }
        write_json(summary_path, summary)
        return summary
    try:
        parsed = COUNTERFACTUAL_ANALYZER.parse_log(log_text)
        truth = load_reference_truth(
            reference_truth_directory, case["case_id"])
        if truth is None:
            truth = parse_truth(log_text)
    except (KeyError, TypeError, ValueError) as error:
        summary = {
            "schema_version": CASE_SUMMARY_SCHEMA_VERSION,
            "status": "failed",
            "case": case,
            "thread_count": thread_count,
            "return_code": completed.returncode,
            "command": command,
            "error": f"audit record parsing failed: {error}",
        }
        write_json(summary_path, summary)
        return summary
    write_json(case_directory / "scenario-truth.json", {
        "schema_version": 1,
        "case": case,
        "atoms": [
            {"serial_id": serial_id, **parameters}
            for serial_id, parameters in sorted(truth.items())
        ],
    })
    write_json(
        case_directory / "trajectory-schema-7.json",
        {"schema_version": 7, "records": parsed["trajectory_records"]})
    write_json(
        case_directory / "trust-model-shadow-schema-1.json",
        {
            "schema_version": 1,
            "diagnostic_only": True,
            "records": [
                {
                    name: value for name, value in record.items()
                    if name != "elapsed-ms"
                }
                for record in parsed["trust_model_shadow_records"]
            ],
        })
    write_json(case_directory / "counterfactual-schema-3.json", {
        "schema_version": 3,
        "checkpoints": parsed["checkpoints"],
        "terminations": parsed["terminations"],
        "atoms": [
            atom
            for records in parsed["atoms"].values()
            for atom in records
        ],
    })
    write_json(case_directory / "shadow-terminal-schema-1.json", {
        "schema_version": 1,
        "accepted_only_checkpoint": parsed["shadow_checkpoint"],
        "accepted_only_atoms": parsed["shadow_atoms"],
        "terminal": parsed["audit_terminal"],
        "terminal_atoms": parsed["audit_terminal_atoms"],
    })
    write_json(case_directory / "shadow-continuation-schema-2.json", {
        "schema_version": 2,
        "checkpoints": parsed["shadow_policy_checkpoints"],
        "atoms": [
            atom
            for records in parsed["shadow_policy_atoms"].values()
            for atom in records
        ],
        "terminal": parsed["audit_terminal"],
        "terminal_atoms": parsed["audit_terminal_atoms"],
    })
    audit = COUNTERFACTUAL_ANALYZER.analyze(parsed, truth)
    summary = {
        "schema_version": CASE_SUMMARY_SCHEMA_VERSION,
        "status": "complete",
        "case": case,
        "thread_count": thread_count,
        "reference_truth_directory": (
            str(reference_truth_directory)
            if reference_truth_directory is not None else None),
        "command": command,
        "truth_atom_count": len(truth),
        "safety_regression": detect_safety_regression(parsed, audit),
        "audit": audit,
    }
    write_json(summary_path, summary)
    return summary


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument(
        "--manifest", type=Path,
        default=PROJECT_ROOT / "tests" / "benchmarks" /
        "convergence_exposure_manifest.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--case-id")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--reference-truth-dir", type=Path)
    parser.add_argument(
        "--case-list", type=Path,
        help="JSON object containing the sorted case_ids to replay")
    return parser.parse_args(argv)


def run(argv: Sequence[str] | None = None) -> int:
    args = parse_arguments(argv)
    executable = args.executable.resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"Exposure runner does not exist: {executable}")
    reference_truth_directory = (
        args.reference_truth_dir.resolve()
        if args.reference_truth_dir is not None else None)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    cases = expand_manifest(manifest)
    if args.case_list is not None:
        selection = json.loads(args.case_list.read_text(encoding="utf-8"))
        if selection.get("schema_version") != 1:
            raise ValueError("Case-list schema_version must be 1")
        requested_case_ids = list(selection.get("case_ids", []))
        if (requested_case_ids != sorted(requested_case_ids) or
                len(requested_case_ids) != len(set(requested_case_ids))):
            raise ValueError("Case-list case_ids must be sorted and unique")
        case_by_id = {case["case_id"]: case for case in cases}
        unknown = [case_id for case_id in requested_case_ids if case_id not in case_by_id]
        if unknown:
            raise ValueError(f"Unknown case-list case_id: {unknown[0]}")
        cases = [case_by_id[case_id] for case_id in requested_case_ids]
    if args.case_id:
        cases = [case for case in cases if case["case_id"] == args.case_id]
        if not cases:
            raise ValueError(f"Unknown exposure case_id: {args.case_id}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.jobs <= 0:
        raise ValueError("--jobs must be positive")
    summaries_by_id: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                run_case,
                executable,
                case,
                args.output_dir,
                args.threads,
                reference_truth_directory): case
            for case in cases
        }
        for completed_count, future in enumerate(as_completed(futures), start=1):
            case = futures[future]
            summaries_by_id[case["case_id"]] = future.result()
            if (completed_count == 1 or completed_count % 25 == 0 or
                    completed_count == len(cases)):
                print(
                    f"Convergence exposure corpus progress: "
                    f"{completed_count}/{len(cases)} ({case['case_id']})",
                    flush=True)
    summaries = [summaries_by_id[case["case_id"]] for case in cases]
    complete = [summary for summary in summaries if summary["status"] == "complete"]
    failed = [summary for summary in summaries if summary["status"] == "failed"]
    aggregate = CORPUS_ANALYZER.analyze(complete)
    aggregate["failed_case_count"] = len(failed)
    aggregate["failed_case_ids"] = [summary["case"]["case_id"] for summary in failed]
    aggregate["incomplete_corpus"] = bool(failed)
    if failed:
        aggregate["corpus_target_met"] = False
        for policy, decision in aggregate["policy_decisions"].items():
            decision[
                "redesign_candidate" if policy == "solver-qualified" else
                "rollback_candidate"
            ] = False
    write_json(args.output_dir / "aggregate.json", aggregate)
    if reference_truth_directory is not None:
        baseline_aggregate_path = reference_truth_directory / "aggregate.json"
        if baseline_aggregate_path.is_file():
            baseline_aggregate = json.loads(
                baseline_aggregate_path.read_text(encoding="utf-8"))
            aggregate_case_ids = {
                row["case_id"] for row in aggregate["cases"]
            }
            baseline_case_ids = {
                row["case_id"] for row in baseline_aggregate["cases"]
            }
            if aggregate_case_ids != baseline_case_ids:
                if not aggregate_case_ids.issubset(baseline_case_ids):
                    raise ValueError(
                        "Reference truth aggregate does not contain every replay case")
                baseline_aggregate = dict(baseline_aggregate)
                baseline_aggregate["cases"] = [
                    row for row in baseline_aggregate["cases"]
                    if row["case_id"] in aggregate_case_ids
                ]
                baseline_aggregate["production_convergence_count"] = sum(
                    row.get("production_converged", False)
                    for row in baseline_aggregate["cases"])
                baseline_aggregate["accepted_only_shadow_count"] = sum(
                    bool(row.get("accepted_only_shadow", {}).get("reached"))
                    for row in baseline_aggregate["cases"])
            comparison = CORPUS_ANALYZER.compare(baseline_aggregate, aggregate)
            write_json(args.output_dir / "comparison.json", comparison)
            objective = comparison["objective_delta"]
            truth_rmse = comparison["truth_rmse_delta"]
            accepted_iteration = comparison["accepted_iteration_delta"]
            blocking_gate = comparison[
                "guard_trust_decoupling_blocking_gate"]
            markdown_lines = [
                "# Convergence exposure before/after comparison",
                "",
                f"- Cases: {comparison['case_count']}",
                f"- Production convergence: "
                f"{comparison['before_convergence_count']} → "
                f"{comparison['after_convergence_count']}",
                f"- Accepted-only shadow checkpoints: "
                f"{comparison['accepted_only_shadow_count']}",
                f"- Objective delta median/p90: "
                f"{objective['median']} / {objective['p90']}",
                f"- Objective material benefit/harm: "
                f"{objective['material_benefit_count']} / "
                f"{objective['material_harm_count']}",
                f"- Accepted iterations before/after median: "
                f"{accepted_iteration['before_median']} / "
                f"{accepted_iteration['after_median']}",
                f"- Accepted-iteration delta median/p90: "
                f"{accepted_iteration['median']} / "
                f"{accepted_iteration['p90']}",
                f"- Truth RMSE delta median/p90: "
                f"{truth_rmse['median']} / {truth_rmse['p90']}",
                f"- Truth material benefit/harm: "
                f"{truth_rmse['material_benefit_count']} / "
                f"{truth_rmse['material_harm_count']}",
                f"- Safety regressions: "
                f"{comparison['safety_regression_count']}",
                f"- Guard/trust decoupling blocking gate: "
                f"{'PASS' if blocking_gate['passed'] else 'FAIL'}",
                "",
                "## By family",
                "",
                "| Family | Cases | Convergence before→after | Shadow | "
                "Objective Δ median/p90 | Truth RMSE Δ median/p90 | Safety |",
                "|---|---:|---:|---:|---:|---:|---:|",
            ]
            for family, values in comparison["by_family"].items():
                markdown_lines.append(
                    f"| {family} | {values['case_count']} | "
                    f"{values['before_convergence_count']}→"
                    f"{values['after_convergence_count']} | "
                    f"{values['accepted_only_shadow_count']} | "
                    f"{values['objective_delta_median']} / "
                    f"{values['objective_delta_p90']} | "
                    f"{values['truth_rmse_delta_median']} / "
                    f"{values['truth_rmse_delta_p90']} | "
                    f"{values['safety_regression_count']} |"
                )
            markdown_lines.extend((
                "",
                "## By topology",
                "",
                "| Topology | Cases | Convergence before→after | Shadow | Safety |",
                "|---|---:|---:|---:|---:|",
            ))
            for topology, values in comparison["by_topology"].items():
                markdown_lines.append(
                    f"| {topology} | {values['case_count']} | "
                    f"{values['before_convergence_count']}→"
                    f"{values['after_convergence_count']} | "
                    f"{values['accepted_only_shadow_count']} | "
                    f"{values['safety_regression_count']} |"
                )
            markdown_lines.append("")
            markdown = "\n".join(markdown_lines)
            (args.output_dir / "comparison.md").write_text(
                markdown, encoding="utf-8")
    write_json(args.output_dir / "replay-manifest.json", {
        "schema_version": 1,
        "case_ids": aggregate["replay_case_ids"],
    })
    shadow_case_ids = sorted(
        row["case_id"] for row in aggregate["cases"]
        if row.get("accepted_only_shadow", {}).get("reached", False))
    write_json(args.output_dir / "shadow-replay-manifest.json", {
        "schema_version": 1,
        "case_ids": shadow_case_ids,
    })
    decision = {
        "schema_version": 1,
        "case_count": aggregate["case_count"],
        "failed_case_count": len(failed),
        "accepted_only_shadow_count": aggregate["accepted_only_shadow_count"],
        "policy_decisions": aggregate["shadow_policy_decisions"],
        "recommended_policy": aggregate["recommended_shadow_policy"],
        "production_change_recommended": aggregate[
            "shadow_policy_production_change_recommended"],
    }
    write_json(args.output_dir / "shadow-continuation-decision.json", decision)
    decision_lines = [
        "# Accepted-only shadow continuation decision",
        "",
        f"- Cases: {decision['case_count']}",
        f"- Accepted-only shadow cases: {decision['accepted_only_shadow_count']}",
        f"- Failed cases: {decision['failed_case_count']}",
        f"- Recommended policy: {decision['recommended_policy'] or 'none'}",
        f"- Production change recommended: "
        f"{str(decision['production_change_recommended']).lower()}",
        "",
        "| Policy | Exposure | Natural | Stationarity | Attempts saved total | "
        "Objective Δ median/p90 | Truth aggregate Δ median/p90 | "
        "Objective harm | Truth harm | Safety | Incomplete | Promote |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for policy, values in decision["policy_decisions"].items():
        decision_lines.append(
            f"| {policy} | {values['effective_exposure_count']} | "
            f"{values['family_counts'].get('natural', 0)} | "
            f"{values['family_counts'].get('stationarity', 0)} | "
            f"{values['attempts_saved']['total']} | "
            f"{values['objective_delta']['median']} / "
            f"{values['objective_delta']['p90']} | "
            f"{values['truth_rmse_delta']['transformed_aggregate_rmse']['median']} / "
            f"{values['truth_rmse_delta']['transformed_aggregate_rmse']['p90']} | "
            f"{values['objective_material_harm_count']} | "
            f"{values['truth_material_harm_case_count']} | "
            f"{values['endpoint_safety_violation_count']} | "
            f"{values['incomplete_comparison_count']} | "
            f"{str(values['promotion_candidate']).lower()} |"
        )
    decision_lines.append("")
    (args.output_dir / "shadow-continuation-decision.md").write_text(
        "\n".join(decision_lines), encoding="utf-8")
    print(
        f"Convergence exposure corpus: cases={len(summaries)}, "
        f"exposures={aggregate['genuine_exposure_count']}, "
        f"failed={len(failed)}, target_met={aggregate['corpus_target_met']}.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(run())
