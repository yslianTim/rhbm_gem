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
CASE_SUMMARY_SCHEMA_VERSION = 4
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


def _numbers(value: str) -> list[float]:
    return [float(item) for item in value.split("/")]


def detect_safety_regression(parsed: dict[str, Any], audit: dict[str, Any]) -> bool:
    if audit["status"] == "no_convergence_trigger":
        return False
    trigger_try = int(audit["experiments"][0]["trigger_try"])
    records = parsed.get("trajectory_records", [])
    baseline = next(
        (record for record in records if int(record["try"]) == trigger_try), None)
    if baseline is None:
        return True

    def safety_signature(record: dict[str, str]) -> tuple[int, int, int, bool]:
        strict = [int(value) for value in record["strict-stationarity"].split("/")]
        blockers = [int(value) for value in record["blockers"].split("/")]
        summary_fields = (
            "accepted-p99", "accepted-max", "raw-p99", "raw-max",
            "shadow-accepted-p99", "shadow-accepted-max",
            "shadow-raw-p99", "shadow-raw-max",
            "dof-accepted-p99", "dof-accepted-max",
            "dof-raw-p99", "dof-raw-max")
        nonfinite = any(
            not math.isfinite(value)
            for field in summary_fields for value in _numbers(record[field]))
        return strict[7] + strict[13], strict[16], blockers[2], nonfinite

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
) -> dict[str, Any]:
    case_directory = output_directory / "cases" / case["case_id"]
    case_directory.mkdir(parents=True, exist_ok=True)
    summary_path = case_directory / "case-summary.json"
    if summary_path.is_file():
        value = json.loads(summary_path.read_text(encoding="utf-8"))
        if (value.get("schema_version") == CASE_SUMMARY_SCHEMA_VERSION and
                value.get("case") == case and
                value.get("thread_count") == thread_count and
                value.get("status") == "complete"):
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
        case_directory / "trajectory-schema-2.json",
        {"schema_version": 2, "records": parsed["trajectory_records"]})
    write_json(case_directory / "counterfactual-schema-1.json", {
        "schema_version": 1,
        "checkpoints": parsed["checkpoints"],
        "terminations": parsed["terminations"],
        "atoms": [
            atom
            for records in parsed["atoms"].values()
            for atom in records
        ],
    })
    audit = COUNTERFACTUAL_ANALYZER.analyze(parsed, truth)
    summary = {
        "schema_version": CASE_SUMMARY_SCHEMA_VERSION,
        "status": "complete",
        "case": case,
        "thread_count": thread_count,
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
    return parser.parse_args(argv)


def run(argv: Sequence[str] | None = None) -> int:
    args = parse_arguments(argv)
    executable = args.executable.resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"Exposure runner does not exist: {executable}")
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    cases = expand_manifest(manifest)
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
                run_case, executable, case, args.output_dir, args.threads): case
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
        for decision in aggregate["policy_decisions"].values():
            decision["production_redesign_candidate"] = False
    write_json(args.output_dir / "aggregate.json", aggregate)
    write_json(args.output_dir / "replay-manifest.json", {
        "schema_version": 1,
        "case_ids": aggregate["replay_case_ids"],
    })
    print(
        f"Convergence exposure corpus: cases={len(summaries)}, "
        f"exposures={aggregate['genuine_exposure_count']}, "
        f"failed={len(failed)}, target_met={aggregate['corpus_target_met']}.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(run())
