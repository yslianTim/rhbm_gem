#!/usr/bin/env python3
"""Run and aggregate the production-only convergence corpus."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import subprocess
import time
from typing import Any, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[3]
CONVERGENCE_ANALYZER_PATH = Path(__file__).with_name(
    "analyze_convergence_audit.py")
CORPUS_ANALYZER_PATH = Path(__file__).with_name(
    "analyze_convergence_exposure_corpus.py")
TRUTH_MARKER = "Convergence exposure truth:"
CASE_SUMMARY_SCHEMA_VERSION = 12
FIELD_PATTERN = re.compile(
    r"(?:^|, )(?P<name>[a-z][a-z0-9-]*)=(?P<value>[^,]+)")


def _load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CONVERGENCE_ANALYZER = _load_module(
    "convergence_audit_analyzer", CONVERGENCE_ANALYZER_PATH)
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
            raise ValueError(
                "Each exposure family must define five topologies and levels")
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
        "active_target", "largest_group_ratio",
    ):
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
        if fields.get("schema") == "1":
            truth[int(fields["serial"])] = {
                name: float(fields[name])
                for name in ("amplitude", "width", "offset")
            }
    return truth


def load_reference_truth(
    reference_directory: Path | None,
    case_id: str,
) -> dict[int, dict[str, float]] | None:
    if reference_directory is None:
        return None
    path = reference_directory / "cases" / case_id / "scenario-truth.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    return {
        int(row["serial_id"]): {
            name: float(row[name]) for name in ("amplitude", "width", "offset")
        }
        for row in value["atoms"]
    }


def semantic_digest(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def semantic_trajectory(records: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        {
            name: normalized[name]
            for name in CONVERGENCE_ANALYZER.PRODUCTION_FIELDS
        }
        for record in records
        for normalized in [CONVERGENCE_ANALYZER.normalize_record(record)]
    ]


def normalized_terminal(parsed: dict[str, Any]) -> dict[str, Any]:
    terminal = parsed.get("terminal")
    if terminal is None:
        raise ValueError("Missing terminal record")
    atoms = [
        {
            name: atom[name]
            for name in ("serial", "group", "amplitude", "width", "offset")
        }
        for atom in sorted(
            parsed["terminal_atoms"], key=lambda row: int(row["serial"]))
    ]
    return {
        "reason": terminal["reason"],
        "try": terminal["try"],
        "acc": terminal["acc"],
        "objective": terminal["objective"],
        "atoms": atoms,
    }


def _objective_total(value: str) -> float:
    values = [float(item) for item in value.split("/")]
    return values[-1] if len(values) == 4 else sum(values)


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
        amplitude = float(row["amplitude"])
        width = float(row["width"])
        truth_amplitude = truth[serial]["amplitude"]
        truth_width = truth[serial]["width"]
        if min(amplitude, width, truth_amplitude, truth_width) <= 0.0:
            return None
        height = amplitude / (2.0 * math.pi * width * width) ** 1.5
        truth_height = truth_amplitude / (
            2.0 * math.pi * truth_width * truth_width) ** 1.5
        errors = (
            math.log(height) - math.log(truth_height),
            math.log(width) - math.log(truth_width),
            float(row["offset"]) / height -
            truth[serial]["offset"] / truth_height,
        )
        if not all(math.isfinite(value) for value in errors):
            return None
        for index, error in enumerate(errors):
            transformed_squared[index] += error * error
    if seen != set(truth):
        return None
    count = len(seen)
    result = {
        f"{name}_rmse": math.sqrt(total / count)
        for name, total in squared.items()
    }
    for name, total in zip(
        ("log_peak", "log_width", "offset_over_peak"), transformed_squared
    ):
        result[f"transformed_{name}_rmse"] = math.sqrt(total / count)
    result["transformed_aggregate_rmse"] = math.sqrt(
        sum(transformed_squared) / (3 * count))
    result["maximum_absolute_offset_error"] = maximum_absolute_offset_error
    return result


def _safety_regression(parsed: dict[str, Any]) -> bool:
    terminal = parsed.get("terminal")
    if terminal is None:
        return True
    try:
        if not all(math.isfinite(float(value)) for value in
                   terminal["objective"].split("/")):
            return True
        return any(
            not all(math.isfinite(float(atom[name])) for name in
                    ("amplitude", "width", "offset")) or
            float(atom["amplitude"]) <= 0.0 or float(atom["width"]) <= 0.0
            for atom in parsed["terminal_atoms"])
    except (KeyError, TypeError, ValueError):
        return True


def _failed_summary(
    case: dict[str, Any],
    thread_count: int,
    command: list[str],
    elapsed_seconds: float,
    return_code: int | None,
    error: str | None = None,
) -> dict[str, Any]:
    return {
        "schema_version": CASE_SUMMARY_SCHEMA_VERSION,
        "status": "failed",
        "case": case,
        "thread_count": thread_count,
        "return_code": return_code,
        "command": command,
        "elapsed_seconds": elapsed_seconds,
        **({"error": error} if error is not None else {}),
    }


def run_case(
    executable: Path,
    case: dict[str, Any],
    output_directory: Path,
    thread_count: int,
    reference_directory: Path | None,
) -> dict[str, Any]:
    case_directory = output_directory / "cases" / case["case_id"]
    case_directory.mkdir(parents=True, exist_ok=True)
    summary_path = case_directory / "case-summary.json"
    expected_reference = (
        str(reference_directory) if reference_directory is not None else None)
    required_artifacts = (
        "run.log", "scenario-truth.json", "trajectory-schema-9.json",
        "terminal-schema-2.json", "case-summary.json")
    if summary_path.is_file() and all(
        (case_directory / name).is_file() for name in required_artifacts
    ):
        cached = json.loads(summary_path.read_text(encoding="utf-8"))
        if (cached.get("schema_version") == CASE_SUMMARY_SCHEMA_VERSION and
                cached.get("status") == "complete" and
                cached.get("case") == case and
                cached.get("thread_count") == thread_count and
                cached.get("reference_truth_directory") == expected_reference):
            return cached
    command = build_command(executable, case, thread_count)
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=False)
    except OSError as error:
        summary = _failed_summary(
            case, thread_count, command, time.perf_counter() - started,
            None, str(error))
        write_json(summary_path, summary)
        return summary
    elapsed_seconds = time.perf_counter() - started
    log_text = completed.stdout.decode("utf-8", errors="replace")
    (case_directory / "run.log").write_text(log_text, encoding="utf-8")
    if completed.returncode != 0:
        summary = _failed_summary(
            case, thread_count, command, elapsed_seconds, completed.returncode)
        write_json(summary_path, summary)
        return summary
    try:
        parsed = CONVERGENCE_ANALYZER.parse_log(log_text)
        terminal_digest_value = normalized_terminal(parsed)
        truth = load_reference_truth(reference_directory, case["case_id"])
        if truth is None:
            truth = parse_truth(log_text)
        truth_metrics = _truth_metrics(parsed["terminal_atoms"], truth)
        if truth_metrics is None:
            raise ValueError("Terminal atoms do not cover frozen truth")
    except (KeyError, TypeError, ValueError) as error:
        summary = _failed_summary(
            case, thread_count, command, elapsed_seconds, completed.returncode,
            f"audit record parsing failed: {error}")
        write_json(summary_path, summary)
        return summary
    scenario_truth = {
        "schema_version": 1,
        "case": case,
        "atoms": [
            {"serial_id": serial_id, **parameters}
            for serial_id, parameters in sorted(truth.items())
        ],
    }
    trajectory = semantic_trajectory(parsed["trajectory_records"])
    terminal = parsed["terminal"]
    terminal_report = {
        "reason": terminal["reason"],
        "try": int(terminal["try"]),
        "accepted_iteration": int(terminal["acc"]),
        "objective": _objective_total(terminal["objective"]),
        "truth_metrics": truth_metrics,
    }
    trajectory_digest = semantic_digest(trajectory)
    terminal_digest = semantic_digest(terminal_digest_value)
    reference_matches: bool | None = None
    if reference_directory is not None:
        baseline_summary = json.loads((
            reference_directory / "cases" / case["case_id"] /
            "case-summary.json").read_text(encoding="utf-8"))
        reference_matches = (
            baseline_summary.get("semantic_trajectory_sha256") ==
            trajectory_digest and
            baseline_summary.get("terminal_state_sha256") == terminal_digest)
    write_json(case_directory / "scenario-truth.json", scenario_truth)
    write_json(case_directory / "trajectory-schema-9.json", {
        "schema_version": 9,
        "certificate_definition": 1,
        "records": trajectory,
    })
    write_json(case_directory / "terminal-schema-2.json", {
        "schema_version": 2,
        **terminal_digest_value,
    })
    summary = {
        "schema_version": CASE_SUMMARY_SCHEMA_VERSION,
        "certificate_definition": 1,
        "status": "complete",
        "case": case,
        "thread_count": thread_count,
        "reference_truth_directory": (
            expected_reference),
        "command": command,
        "elapsed_seconds": elapsed_seconds,
        "truth_atom_count": len(truth),
        "production_converged": terminal["reason"] == "converged",
        "safety_regression": _safety_regression(parsed),
        "production_artifacts_identical": reference_matches,
        "semantic_trajectory_sha256": trajectory_digest,
        "terminal_state_sha256": terminal_digest,
        "frozen_truth_sha256": semantic_digest(scenario_truth),
        "terminal": terminal_report,
    }
    write_json(summary_path, summary)
    return summary


def build_compact_baseline(
    manifest_bytes: bytes,
    summaries: list[dict[str, Any]],
    aggregate: dict[str, Any],
) -> dict[str, Any]:
    complete = [
        summary for summary in summaries if summary.get("status") == "complete"
    ]
    case_identity = [
        {"case_id": row["case"]["case_id"], "seed": row["case"]["seed"]}
        for row in complete
    ]
    truth_identity = [
        {
            "case_id": row["case"]["case_id"],
            "sha256": row["frozen_truth_sha256"],
        }
        for row in complete
    ]
    return {
        "schema_version": 2,
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "case_identity_sha256": semantic_digest(case_identity),
        "frozen_truth_sha256": semantic_digest(truth_identity),
        "schema_contract": {
            "certificate_definition": 1,
            "trajectory": 9,
            "terminal": 2,
            "case_summary": CASE_SUMMARY_SCHEMA_VERSION,
            "aggregate": CORPUS_ANALYZER.SCHEMA_VERSION,
            "comparison": CORPUS_ANALYZER.COMPARISON_SCHEMA_VERSION,
        },
        "production_definition": (
            "solver-qualified && accepted-active-p99 && operator-complete && "
            "operator-nominal-p99 && invariants-clear && orthogonal-clear"),
        "case_count": aggregate["case_count"],
        "failed_case_count": aggregate.get("failed_case_count", 0),
        "production_convergence_count": aggregate[
            "production_convergence_count"],
        "termination_counts": aggregate["termination_counts"],
        "safety_regression_count": aggregate["safety_regression_count"],
        "case_fields": [
            "case_id", "seed", "semantic_trajectory_sha256",
            "terminal_state_sha256", "stop_reason"],
        "cases": [
            [
                summary["case"]["case_id"],
                summary["case"]["seed"],
                summary["semantic_trajectory_sha256"],
                summary["terminal_state_sha256"],
                summary["terminal"]["reason"],
            ]
            for summary in complete
        ],
    }


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
    parser.add_argument("--case-list", type=Path)
    return parser.parse_args(argv)


def run(argv: Sequence[str] | None = None) -> int:
    args = parse_arguments(argv)
    executable = args.executable.resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"Exposure runner does not exist: {executable}")
    if args.jobs <= 0:
        raise ValueError("--jobs must be positive")
    reference_directory = (
        args.reference_truth_dir.resolve()
        if args.reference_truth_dir is not None else None)
    manifest_bytes = args.manifest.read_bytes()
    cases = expand_manifest(json.loads(manifest_bytes.decode("utf-8")))
    if args.case_list is not None:
        selection = json.loads(args.case_list.read_text(encoding="utf-8"))
        if selection.get("schema_version") != 1:
            raise ValueError("Case-list schema_version must be 1")
        requested = list(selection.get("case_ids", []))
        if requested != sorted(set(requested)):
            raise ValueError("Case-list case_ids must be sorted and unique")
        case_by_id = {case["case_id"]: case for case in cases}
        unknown = [case_id for case_id in requested if case_id not in case_by_id]
        if unknown:
            raise ValueError(f"Unknown case-list case_id: {unknown[0]}")
        cases = [case_by_id[case_id] for case_id in requested]
    if args.case_id is not None:
        cases = [case for case in cases if case["case_id"] == args.case_id]
        if not cases:
            raise ValueError(f"Unknown exposure case_id: {args.case_id}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries_by_id: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                run_case, executable, case, args.output_dir, args.threads,
                reference_directory): case
            for case in cases
        }
        for completed_count, future in enumerate(as_completed(futures), start=1):
            case = futures[future]
            summaries_by_id[case["case_id"]] = future.result()
            if (completed_count == 1 or completed_count % 25 == 0 or
                    completed_count == len(cases)):
                print(
                    f"Convergence corpus progress: {completed_count}/"
                    f"{len(cases)} ({case['case_id']})", flush=True)
    summaries = [summaries_by_id[case["case_id"]] for case in cases]
    complete = [row for row in summaries if row["status"] == "complete"]
    failed = [row for row in summaries if row["status"] == "failed"]
    aggregate = CORPUS_ANALYZER.analyze(complete)
    aggregate.update({
        "failed_case_count": len(failed),
        "failed_case_ids": [row["case"]["case_id"] for row in failed],
        "incomplete_corpus": bool(failed),
    })
    write_json(args.output_dir / "aggregate.json", aggregate)
    if reference_directory is not None:
        baseline = json.loads((
            reference_directory / "aggregate.json").read_text(encoding="utf-8"))
        comparison = CORPUS_ANALYZER.compare(baseline, aggregate)
        comparison["baseline_complete_count"] = baseline.get("case_count", 0)
        comparison["candidate_complete_count"] = aggregate["case_count"]
        comparison["baseline_failed_case_count"] = baseline.get(
            "failed_case_count", 0)
        comparison["candidate_failed_case_count"] = len(failed)
        if (comparison["baseline_complete_count"] != 600 or
                comparison["candidate_complete_count"] != 600 or
                comparison["baseline_failed_case_count"] != 0 or
                comparison["candidate_failed_case_count"] != 0):
            comparison["blocking_gate"]["conditions"]["complete-pair"] = False
            comparison["blocking_gate"]["passed"] = False
        write_json(args.output_dir / "comparison.json", comparison)
    write_json(
        args.output_dir / "compact-baseline.json",
        build_compact_baseline(manifest_bytes, summaries, aggregate))
    print(
        f"Convergence corpus: cases={len(summaries)}, complete={len(complete)}, "
        f"failed={len(failed)}.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(run())
