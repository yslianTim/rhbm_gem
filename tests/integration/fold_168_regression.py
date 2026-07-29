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


SCHEMA_VERSION = 4
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
EXPECTED_ATOM_COUNT = 168
QUALITY_TOLERANCE_RATIO = 1.05
MAXIMUM_ACCEPTED_ITERATIONS = 10
TRUTH_WIDTH = 0.50
TRUTH_OFFSET = 1.0
SECOND_STAGE_SUMMARY_PATTERN = re.compile(
    r"Second-stage local fitting summary: "
    r"accepted_iterations=(?P<accepted_iterations>\d+), "
    r"best_iteration=(?P<best_iteration>initial|unavailable|\d+), "
    r"stop_reason=(?P<stop_reason>[a-z-]+), "
    r"best_audit_objective=(?P<best_audit_objective>unavailable|"
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?), "
    r"final_uses_polish=(?P<final_uses_polish>yes|no|unavailable), "
    r"final_state_source=(?P<final_state_source>"
    r"best-audit|latest-validated|unavailable)\."
)


class RegressionError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_atom_results(database_path: Path) -> list[dict[str, Any]]:
    query = """
        SELECT local.serial_id,
               atom.element,
               local.amplitude_estimate_mdpde,
               local.width_estimate_mdpde,
               local.intercept_estimate_mdpde,
               local.alpha_r
          FROM model_atom_local_potential AS local
          JOIN model_atom AS atom
            ON atom.key_tag = local.key_tag
           AND atom.serial_id = local.serial_id
         WHERE local.key_tag = ?
         ORDER BY local.serial_id
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
            element = int(row[1])
            values = [float(value) for value in row[2:]]
        except (TypeError, ValueError) as error:
            raise RegressionError(
                f"Invalid atom result for serial ID {serial_id}.") from error
        if not all(math.isfinite(value) for value in values):
            raise RegressionError(f"Non-finite atom result for serial ID {serial_id}.")
        atoms.append({
            "serial_id": serial_id,
            "element": element,
            "amplitude_mdpde": values[0],
            "width_mdpde": values[1],
            "intercept_mdpde": values[2],
            "alpha_r": values[3],
        })
    return atoms


def calculate_quality_metrics(atoms: Sequence[dict[str, Any]]) -> dict[str, float]:
    if not atoms:
        raise RegressionError("Cannot calculate fold-168 quality metrics without atoms.")

    amplitude_error_square_sum = 0.0
    width_error_square_sum = 0.0
    offset_error_square_sum = 0.0
    maximum_absolute_offset = 0.0
    for atom in atoms:
        amplitude = float(atom["amplitude_mdpde"])
        width = float(atom["width_mdpde"])
        offset = float(atom["intercept_mdpde"])
        atomic_number = int(atom["element"])
        amplitude_error_square_sum += (amplitude - atomic_number) ** 2
        width_error_square_sum += (width - TRUTH_WIDTH) ** 2
        offset_error_square_sum += (offset - TRUTH_OFFSET) ** 2
        maximum_absolute_offset = max(maximum_absolute_offset, abs(offset))

    atom_count = len(atoms)
    return {
        "amplitude_rmse": math.sqrt(amplitude_error_square_sum / atom_count),
        "width_rmse": math.sqrt(width_error_square_sum / atom_count),
        "offset_rmse": math.sqrt(offset_error_square_sum / atom_count),
        "maximum_absolute_offset": maximum_absolute_offset,
    }


def parse_second_stage_summary(log_text: str) -> dict[str, Any]:
    matches = list(SECOND_STAGE_SUMMARY_PATTERN.finditer(log_text))
    if len(matches) != 1:
        raise RegressionError(
            "Expected exactly one parseable second-stage final summary, "
            f"found {len(matches)}.")
    values = matches[0].groupdict()
    best_objective_text = values["best_audit_objective"]
    best_objective = (
        None if best_objective_text == "unavailable" else float(best_objective_text))
    if best_objective is not None and not math.isfinite(best_objective):
        raise RegressionError("Second-stage best audit objective is not finite.")
    final_uses_polish_text = values["final_uses_polish"]
    final_uses_polish = (
        None if final_uses_polish_text == "unavailable" else
        final_uses_polish_text == "yes")
    return {
        "accepted_iterations": int(values["accepted_iterations"]),
        "best_iteration": values["best_iteration"],
        "stop_reason": values["stop_reason"],
        "best_audit_objective": best_objective,
        "final_uses_polish": final_uses_polish,
        "final_state_source": values["final_state_source"],
    }


def validate_quality_gate(
    baseline: dict[str, Any],
    actual: dict[str, Any],
) -> list[str]:
    differences: list[str] = []
    atoms = actual.get("atoms", [])
    expected_serial_ids = baseline["serial_ids"]
    actual_serial_ids = [atom["serial_id"] for atom in atoms]
    if len(atoms) != EXPECTED_ATOM_COUNT:
        differences.append(
            f"atoms: expected {EXPECTED_ATOM_COUNT} entries, got {len(atoms)}")
    if actual_serial_ids != expected_serial_ids:
        differences.append("atoms: serial IDs differ from the fixed fold-168 baseline")

    for atom in atoms:
        serial_id = atom["serial_id"]
        if int(atom["element"]) <= 0:
            differences.append(f"atom {serial_id}: invalid atomic number")
        if float(atom["amplitude_mdpde"]) <= 0.0:
            differences.append(f"atom {serial_id}: amplitude is not positive")
        if float(atom["width_mdpde"]) <= 0.0:
            differences.append(f"atom {serial_id}: width is not positive")
        alpha_r = float(atom["alpha_r"])
        if not 0.0 <= alpha_r <= 1.0:
            differences.append(f"atom {serial_id}: alpha_r is outside [0, 1]")

    actual_metrics = actual.get("quality_metrics")
    if not isinstance(actual_metrics, dict):
        differences.append("quality_metrics: missing")
    else:
        for name, reference_value in baseline["reference_quality_metrics"].items():
            actual_value = actual_metrics.get(name)
            if not isinstance(actual_value, (int, float)) or not math.isfinite(actual_value):
                differences.append(f"quality_metrics.{name}: missing or non-finite")
                continue
            limit = float(reference_value) * QUALITY_TOLERANCE_RATIO
            if float(actual_value) > limit:
                differences.append(
                    f"quality_metrics.{name}: {float(actual_value):.17g} exceeds "
                    f"105% reference limit {limit:.17g}")

    summary = actual.get("second_stage_summary")
    if not isinstance(summary, dict):
        differences.append("second_stage_summary: missing")
    elif int(summary.get("accepted_iterations", MAXIMUM_ACCEPTED_ITERATIONS + 1)) > \
            MAXIMUM_ACCEPTED_ITERATIONS:
        differences.append(
            "second_stage_summary.accepted_iterations: "
            f"expected <= {MAXIMUM_ACCEPTED_ITERATIONS}, "
            f"got {summary.get('accepted_iterations')}")
    return differences


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
        "atoms": [],
        "quality_metrics": None,
        "second_stage_summary": None,
    }


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
    serial_ids = baseline.get("serial_ids")
    if not isinstance(serial_ids, list) or len(serial_ids) != EXPECTED_ATOM_COUNT:
        raise RegressionError("Baseline does not contain 168 serial IDs.")
    reference_metrics = baseline.get("reference_quality_metrics")
    expected_metric_names = {
        "amplitude_rmse",
        "width_rmse",
        "offset_rmse",
        "maximum_absolute_offset",
    }
    if not isinstance(reference_metrics, dict) or set(reference_metrics) != expected_metric_names:
        raise RegressionError("Baseline quality metric schema is invalid.")
    if not all(
        isinstance(value, (int, float)) and math.isfinite(value) and value >= 0.0
        for value in reference_metrics.values()
    ):
        raise RegressionError("Baseline quality metrics must be finite and non-negative.")
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

            actual["atoms"] = read_atom_results(temporary_database)
            actual["quality_metrics"] = calculate_quality_metrics(actual["atoms"])
            actual["second_stage_summary"] = parse_second_stage_summary(log_text)

        differences = validate_quality_gate(baseline, actual)
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
