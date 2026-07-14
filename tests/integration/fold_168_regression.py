#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sqlite3
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any, Sequence


SCHEMA_VERSION = 2
SAVED_KEY = "fold_gaus_charge1"
EXPECTED_INPUT_HASHES = {
    "model": "156d35aa326f0d4408d726a999329d2ffede775489aeaa5d99a2cc9b9f663cab",
    "map": "5e0dbb13fc3a76f8a944e6e2b18393d1896fafc2ec9020457cca8e8a421f120e",
}
COMMAND_ARGUMENT_TEMPLATE = [
    "potential_analysis",
    "-d", "{database}",
    "-a", "{model}",
    "-m", "{map}",
    "-k", SAVED_KEY,
    "-j", "4",
    "-v", "3",
    "--simulation", "true",
    "-r", "0.50",
    "--fit-max", "1.0",
    "--exclude-hydrogen", "true",
]
ATOM_RELATIVE_TOLERANCE = 1.0e-6
ATOM_ABSOLUTE_TOLERANCE = 1.0e-8
STRICT_FLOAT_TOLERANCE = 1.0e-12
ANSI_ESCAPE_PATTERN = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


class RegressionError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_log(text: str) -> str:
    return ANSI_ESCAPE_PATTERN.sub("", text).replace("\r", "\n")


def require_unique_match(pattern: str, text: str, description: str) -> re.Match[str]:
    matches = list(re.finditer(pattern, text, flags=re.MULTILINE))
    if len(matches) != 1:
        raise RegressionError(
            f"Expected exactly one {description}; found {len(matches)}.")
    return matches[0]


def parse_joint_statuses(text: str) -> dict[str, dict[str, int]]:
    statuses: dict[str, dict[str, int]] = {}
    for entry in text.split(","):
        match = re.fullmatch(r"\s*([a-z0-9-]+):(\d+)/(\d+)\s*", entry)
        if match is None:
            raise RegressionError(f"Invalid joint-offset status entry: {entry!r}.")
        name = match.group(1)
        if name in statuses:
            raise RegressionError(f"Duplicate joint-offset status: {name}.")
        statuses[name] = {
            "clusters": int(match.group(2)),
            "atoms": int(match.group(3)),
        }
    return statuses


def parse_iteration_line(line: str) -> dict[str, Any]:
    prefix = require_unique_match(
        rf"Iter\. (\d+)/(\d+), active/frozen atoms = (\d+)/(\d+)",
        line,
        "iteration prefix")
    component = require_unique_match(
        rf"iteration components/max-atoms/active-ratio = "
        rf"(\d+)/(\d+)/({FLOAT_PATTERN})",
        line,
        "iteration component summary")
    freeze = require_unique_match(
        r"freeze outcomes ineligible/above-threshold/stabilizing/newly-frozen = "
        r"(\d+)/(\d+)/(\d+)/(\d+)",
        line,
        "freeze outcome summary")
    thaw = require_unique_match(
        r"thaw events dependency/suspicious = (\d+)/(\d+)",
        line,
        "thaw event summary")

    status_marker = "joint-offset statuses clusters/atoms = "
    status_start = line.find(status_marker)
    if status_start < 0:
        raise RegressionError("Iteration line is missing joint-offset statuses.")
    status_start += len(status_marker)
    status_end_candidates = [
        index for marker in (", health-unhealthy", ", freeze outcomes")
        if (index := line.find(marker, status_start)) >= 0
    ]
    if not status_end_candidates:
        raise RegressionError("Iteration line has no joint-offset status terminator.")
    statuses = parse_joint_statuses(line[status_start:min(status_end_candidates)])

    return {
        "iteration": int(prefix.group(1)),
        "iteration_limit": int(prefix.group(2)),
        "active_atoms_after_update": int(prefix.group(3)),
        "frozen_atoms_after_update": int(prefix.group(4)),
        "components_at_start": int(component.group(1)),
        "maximum_component_atoms_at_start": int(component.group(2)),
        "maximum_component_active_ratio": float(component.group(3)),
        "joint_statuses": statuses,
        "freeze_outcomes": {
            "ineligible": int(freeze.group(1)),
            "above_threshold": int(freeze.group(2)),
            "stabilizing": int(freeze.group(3)),
            "newly_frozen": int(freeze.group(4)),
        },
        "thaw_events": {
            "dependency": int(thaw.group(1)),
            "suspicious": int(thaw.group(2)),
        },
    }


def parse_log(text: str) -> dict[str, Any]:
    normalized = normalize_log(text)
    topology_match = require_unique_match(
        rf"^Local-fitting coupling graph mode = ([a-z-]+), minimum weight = "
        rf"({FLOAT_PATTERN}), candidate/retained/cut edges = (\d+)/(\d+)/(\d+),"
        rf".*initial components/max atoms/ratio = (\d+)/(\d+)/({FLOAT_PATTERN})\.$",
        normalized,
        "coupling topology summary")
    topology = {
        "mode": topology_match.group(1),
        "minimum_weight": float(topology_match.group(2)),
        "candidate_edges": int(topology_match.group(3)),
        "retained_edges": int(topology_match.group(4)),
        "cut_edges": int(topology_match.group(5)),
        "initial_components": int(topology_match.group(6)),
        "initial_maximum_component_atoms": int(topology_match.group(7)),
        "initial_maximum_component_ratio": float(topology_match.group(8)),
    }

    sensitivity_pattern = re.compile(
        rf"^Coupling sensitivity: threshold=({FLOAT_PATTERN}), retained/cut=(\d+)/(\d+), "
        rf"components/max-atoms/ratio=(\d+)/(\d+)/({FLOAT_PATTERN})\.$",
        flags=re.MULTILINE)
    sensitivity: list[dict[str, Any]] = []
    seen_thresholds: set[float] = set()
    for match in sensitivity_pattern.finditer(normalized):
        threshold = float(match.group(1))
        if threshold in seen_thresholds:
            raise RegressionError(f"Duplicate coupling sensitivity threshold: {threshold}.")
        seen_thresholds.add(threshold)
        sensitivity.append({
            "threshold": threshold,
            "retained_edges": int(match.group(2)),
            "cut_edges": int(match.group(3)),
            "components": int(match.group(4)),
            "maximum_component_atoms": int(match.group(5)),
            "maximum_component_ratio": float(match.group(6)),
        })

    iteration_by_number: dict[int, dict[str, Any]] = {}
    for line in normalized.splitlines():
        if not line.startswith("Iter. "):
            continue
        iteration = parse_iteration_line(line)
        number = iteration["iteration"]
        if number in iteration_by_number:
            raise RegressionError(f"Duplicate local-fitting iteration: {number}.")
        iteration_by_number[number] = iteration

    maximum_iteration_matches = list(re.finditer(
        rf"^\[Warning\] Reached maximum iteration size; applying best validated audit state; "
        rf"audit best source = accepted iteration (\d+), fixed audit objective "
        rf"residual/width/offset/total = ({FLOAT_PATTERN})/({FLOAT_PATTERN})/"
        rf"({FLOAT_PATTERN})/({FLOAT_PATTERN});",
        normalized,
        flags=re.MULTILINE))
    stagnation_matches = list(re.finditer(
        rf"^\[Warning\] Stopped local fitting because forced fixed-point produced no "
        rf"objective improvement after (\d+) accepted iterations; stagnation "
        rf"forced/stalled clusters = (\d+)/(\d+), atoms = (\d+)/(\d+); "
        rf"applying best validated audit state; audit best source = accepted iteration "
        rf"(\d+), fixed audit objective residual/width/offset/total = "
        rf"({FLOAT_PATTERN})/({FLOAT_PATTERN})/({FLOAT_PATTERN})/({FLOAT_PATTERN});",
        normalized,
        flags=re.MULTILINE))
    if len(maximum_iteration_matches) + len(stagnation_matches) != 1:
        raise RegressionError(
            "Expected exactly one terminal audit summary; found "
            f"{len(maximum_iteration_matches) + len(stagnation_matches)}.")

    if stagnation_matches:
        terminal_match = stagnation_matches[0]
        terminal = {
            "outcome": "forced-fixed-point-stagnation",
            "accepted_iterations": int(terminal_match.group(1)),
            "forced_clusters": int(terminal_match.group(2)),
            "stalled_clusters": int(terminal_match.group(3)),
            "forced_atoms": int(terminal_match.group(4)),
            "stalled_atoms": int(terminal_match.group(5)),
            "applied_state": "best-validated-audit",
            "audit_source": "accepted-iteration",
            "audit_iteration": int(terminal_match.group(6)),
            "audit_objective": {
                "residual": float(terminal_match.group(7)),
                "width": float(terminal_match.group(8)),
                "offset": float(terminal_match.group(9)),
                "total": float(terminal_match.group(10)),
            },
        }
    else:
        terminal_match = maximum_iteration_matches[0]
        terminal = {
            "outcome": "maximum-iteration",
            "applied_state": "best-validated-audit",
            "audit_source": "accepted-iteration",
            "audit_iteration": int(terminal_match.group(1)),
            "audit_objective": {
                "residual": float(terminal_match.group(2)),
                "width": float(terminal_match.group(3)),
                "offset": float(terminal_match.group(4)),
                "total": float(terminal_match.group(5)),
            },
        }
    return {
        "topology": topology,
        "sensitivity": sensitivity,
        "iterations": [iteration_by_number[key] for key in sorted(iteration_by_number)],
        "terminal": terminal,
    }


def read_atom_results(database_path: Path) -> list[dict[str, Any]]:
    query = """
        SELECT serial_id,
               amplitude_estimate_mdpde,
               width_estimate_mdpde,
               intercept_estimate_mdpde,
               alpha_r
          FROM model_atom_local_potential
         WHERE key_tag = ?
         ORDER BY serial_id
    """
    with sqlite3.connect(database_path) as connection:
        rows = connection.execute(query, (SAVED_KEY,)).fetchall()

    atoms: list[dict[str, Any]] = []
    seen_serial_ids: set[int] = set()
    for row in rows:
        serial_id = int(row[0])
        if serial_id in seen_serial_ids:
            raise RegressionError(f"Duplicate atom serial ID in output database: {serial_id}.")
        seen_serial_ids.add(serial_id)
        try:
            values = [float(value) for value in row[1:]]
        except (TypeError, ValueError) as error:
            raise RegressionError(
                f"Invalid atom result for serial ID {serial_id}.") from error
        if not all(math.isfinite(value) for value in values):
            raise RegressionError(f"Non-finite atom result for serial ID {serial_id}.")
        atoms.append({
            "serial_id": serial_id,
            "amplitude_mdpde": values[0],
            "width_mdpde": values[1],
            "intercept_mdpde": values[2],
            "alpha_r": values[3],
        })
    return atoms


def validate_input_hashes(paths: dict[str, Path]) -> dict[str, str]:
    actual_hashes: dict[str, str] = {}
    failures: list[str] = []
    for name, expected_hash in EXPECTED_INPUT_HASHES.items():
        path = paths[name]
        if not path.is_file():
            failures.append(f"{name} input does not exist: {path}")
            continue
        actual_hash = sha256_file(path)
        actual_hashes[name] = actual_hash
        if actual_hash != expected_hash:
            failures.append(
                f"{name} SHA-256 mismatch: expected {expected_hash}, got {actual_hash}")
    if failures:
        raise RegressionError("\n".join(failures))
    return actual_hashes


def build_command(
    executable: Path,
    database: Path,
    model: Path,
    map_path: Path,
) -> list[str]:
    substitutions = {
        "database": str(database),
        "model": str(model),
        "map": str(map_path),
    }
    return [
        str(executable),
        *(argument.format(**substitutions) for argument in COMMAND_ARGUMENT_TEMPLATE),
    ]


def make_empty_actual(input_hashes: dict[str, str]) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "input_hashes": input_hashes,
        "command_arguments": COMMAND_ARGUMENT_TEMPLATE,
        "topology": None,
        "sensitivity": [],
        "iterations": [],
        "terminal": None,
        "atoms": [],
    }


def compare_values(
    expected: Any,
    actual: Any,
    path: str = "root",
) -> list[str]:
    differences: list[str] = []
    if isinstance(expected, dict):
        if not isinstance(actual, dict):
            return [f"{path}: expected object, got {type(actual).__name__}"]
        expected_keys = set(expected)
        actual_keys = set(actual)
        for key in sorted(expected_keys - actual_keys):
            differences.append(f"{path}: missing key {key!r}")
        for key in sorted(actual_keys - expected_keys):
            differences.append(f"{path}: unexpected key {key!r}")
        for key in sorted(expected_keys & actual_keys):
            differences.extend(compare_values(expected[key], actual[key], f"{path}.{key}"))
        return differences
    if isinstance(expected, list):
        if not isinstance(actual, list):
            return [f"{path}: expected list, got {type(actual).__name__}"]
        if len(expected) != len(actual):
            differences.append(
                f"{path}: expected {len(expected)} entries, got {len(actual)}")
        for index, (expected_item, actual_item) in enumerate(zip(expected, actual)):
            differences.extend(
                compare_values(expected_item, actual_item, f"{path}[{index}]"))
        return differences
    if isinstance(expected, bool) or isinstance(actual, bool):
        if expected != actual:
            differences.append(f"{path}: expected {expected!r}, got {actual!r}")
        return differences
    if isinstance(expected, float):
        if not isinstance(actual, (int, float)) or not math.isfinite(float(actual)):
            return [f"{path}: expected finite number, got {actual!r}"]
        is_atom_value = path.startswith("root.atoms[") and not path.endswith(".serial_id")
        relative_tolerance = (
            ATOM_RELATIVE_TOLERANCE if is_atom_value else STRICT_FLOAT_TOLERANCE)
        absolute_tolerance = (
            ATOM_ABSOLUTE_TOLERANCE if is_atom_value else STRICT_FLOAT_TOLERANCE)
        if not math.isclose(
            expected,
            float(actual),
            rel_tol=relative_tolerance,
            abs_tol=absolute_tolerance,
        ):
            differences.append(f"{path}: expected {expected:.17g}, got {float(actual):.17g}")
        return differences
    if expected != actual:
        differences.append(f"{path}: expected {expected!r}, got {actual!r}")
    return differences


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=False, allow_nan=False) + "\n",
        encoding="utf-8")


def load_baseline(path: Path) -> dict[str, Any]:
    try:
        baseline = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RegressionError(f"Failed to read benchmark baseline {path}: {error}") from error
    if baseline.get("schema_version") != SCHEMA_VERSION:
        raise RegressionError(
            f"Unsupported baseline schema version: {baseline.get('schema_version')!r}.")
    if baseline.get("input_hashes") != EXPECTED_INPUT_HASHES:
        raise RegressionError("Baseline input hashes do not match the fixed fold-168 fixture.")
    if baseline.get("command_arguments") != COMMAND_ARGUMENT_TEMPLATE:
        raise RegressionError("Baseline command arguments do not match the fixed benchmark command.")
    return baseline


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the external fold-168 regression benchmark.")
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--map", dest="map_path", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def run(argv: Sequence[str] | None = None) -> int:
    args = parse_arguments(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.output_dir / "run.log"
    actual_path = args.output_dir / "actual.json"
    report_path = args.output_dir / "report.json"
    input_paths = {
        "model": args.model.resolve(),
        "map": args.map_path.resolve(),
    }
    actual_hashes: dict[str, str] = {}
    actual = make_empty_actual(actual_hashes)
    command: list[str] = []
    errors: list[str] = []
    differences: list[str] = []
    wall_time_seconds: float | None = None
    log_text = ""

    try:
        baseline = load_baseline(args.baseline.resolve())
        executable = args.executable.resolve()
        if not executable.is_file():
            raise RegressionError(f"Benchmark executable does not exist: {executable}")
        actual_hashes = validate_input_hashes(input_paths)
        actual = make_empty_actual(actual_hashes)

        with tempfile.TemporaryDirectory(prefix="rhbm_fold_168_regression_") as temp_dir:
            temporary_database = Path(temp_dir) / "database.sqlite"
            if temporary_database.exists():
                raise RegressionError(
                    f"Temporary output database already exists: {temporary_database}")
            command = build_command(
                executable,
                temporary_database,
                input_paths["model"],
                input_paths["map"])
            start_time = time.perf_counter()
            completed = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False)
            wall_time_seconds = time.perf_counter() - start_time
            log_text = completed.stdout.decode("utf-8", errors="replace")
            log_path.write_text(log_text, encoding="utf-8")
            if completed.returncode != 0:
                raise RegressionError(
                    f"Benchmark command exited with status {completed.returncode}.")
            if not temporary_database.is_file():
                raise RegressionError(
                    "Benchmark command did not create the temporary output database.")

            parsed_log = parse_log(log_text)
            actual.update(parsed_log)
            actual["atoms"] = read_atom_results(temporary_database)

        differences = compare_values(baseline, actual)
    except Exception as error:  # Preserve artifacts for all benchmark failures.
        errors.append(str(error))

    log_path.write_text(log_text, encoding="utf-8")
    write_json(actual_path, actual)
    passed = not errors and not differences
    report = {
        "schema_version": SCHEMA_VERSION,
        "passed": passed,
        "performance_gate": False,
        "wall_time_seconds": wall_time_seconds,
        "command": command,
        "errors": errors,
        "differences": differences,
        "artifacts": {
            "log": str(log_path),
            "actual": str(actual_path),
        },
    }
    write_json(report_path, report)

    if passed:
        print(
            "fold-168 regression passed; "
            f"wall time = {wall_time_seconds:.3f} s (observation only).")
        return 0
    print(f"fold-168 regression failed; see {report_path}.")
    for message in [*errors, *differences[:20]]:
        print(f"- {message}")
    if len(differences) > 20:
        print(f"- ... {len(differences) - 20} additional differences")
    return 1


def main() -> int:
    return run()


if __name__ == "__main__":
    raise SystemExit(main())
