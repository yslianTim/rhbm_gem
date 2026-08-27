#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import time
from typing import Any, Sequence

import fold_168_regression as fold


PROJECT_ROOT = Path(__file__).resolve().parents[2]
ANALYZER_PATH = (
    PROJECT_ROOT / "resources" / "tools" / "developer"
    / "analyze_counterfactual_convergence.py"
)
SPEC = importlib.util.spec_from_file_location("counterfactual_analyzer", ANALYZER_PATH)
assert SPEC is not None and SPEC.loader is not None
ANALYZER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYZER)


def build_debug_command(
    executable: Path,
    database: Path,
    model: Path,
    map_path: Path,
) -> list[str]:
    command = fold.build_command(executable, database, model, map_path)
    verbosity_position = command.index("-v") + 1
    command[verbosity_position] = "4"
    return command


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the build-gated fold-168 counterfactual convergence audit.")
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
    audit_path = args.output_dir / "counterfactual-audit.json"
    report_path = args.output_dir / "report.json"
    errors: list[str] = []
    differences: list[str] = []
    command: list[str] = []
    wall_time_seconds: float | None = None
    summary: dict[str, Any] | None = None
    audit_report: dict[str, Any] = {
        "status": "not_run", "experiment_count": 0, "experiments": []}
    quality_metrics: dict[str, float] | None = None
    log_text = ""
    try:
        baseline = fold.load_baseline(args.baseline.resolve())
        input_paths = {
            "model": args.model.resolve(),
            "map": args.map_path.resolve(),
        }
        fold.validate_input_hashes(input_paths)
        executable = args.executable.resolve()
        if not executable.is_file():
            raise fold.RegressionError(f"Audit executable does not exist: {executable}")

        with tempfile.TemporaryDirectory(
                prefix="rhbm_fold_168_counterfactual_") as temp_dir:
            database = Path(temp_dir) / "database.sqlite"
            command = build_debug_command(
                executable, database, input_paths["model"], input_paths["map"])
            started = time.perf_counter()
            completed = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
                cwd=temp_dir,
            )
            wall_time_seconds = time.perf_counter() - started
            log_text = completed.stdout.decode("utf-8", errors="replace")
            if completed.returncode != 0:
                raise fold.RegressionError(
                    f"Audit command exited with status {completed.returncode}.")
            atoms = fold.read_atom_results(database)
            quality_metrics = fold.calculate_quality_metrics(atoms)
            summary = fold.parse_second_stage_summary(log_text)
            actual = {
                "schema_version": fold.SCHEMA_VERSION,
                "input_hashes": fold.EXPECTED_INPUT_HASHES,
                "command_arguments": fold.COMMAND_ARGUMENT_TEMPLATE,
                "atoms": atoms,
                "quality_metrics": quality_metrics,
                "second_stage_summary": summary,
                "residue_cutoff_summary": fold.parse_residue_cutoff_summary(log_text),
            }
            write_json(actual_path, actual)
            truth = {
                int(atom["serial_id"]): {
                    "amplitude": float(atom["element"]),
                    "width": fold.TRUTH_WIDTH,
                    "offset": fold.TRUTH_OFFSET,
                }
                for atom in atoms
            }
            audit_report = ANALYZER.analyze(ANALYZER.parse_log(log_text), truth)
            differences = fold.validate_quality_gate(baseline, actual)
            if audit_report["status"] == "counterfactual_records":
                differences = [
                    difference for difference in differences
                    if not difference.startswith(
                        "second_stage_summary.accepted_iterations:")
                ]
            if summary["stop_reason"] == "converged" and \
                    audit_report["status"] == "no_convergence_trigger":
                errors.append(
                    "The audit build reached convergence without a counterfactual checkpoint.")
    except Exception as error:  # Preserve diagnostic artifacts on every failure.
        errors.append(str(error))

    log_path.write_text(log_text, encoding="utf-8")
    write_json(audit_path, audit_report)
    passed = not errors and not differences
    report = {
        "schema_version": 1,
        "passed": passed,
        "wall_time_seconds": wall_time_seconds,
        "command": command,
        "second_stage_summary": summary,
        "quality_metrics": quality_metrics,
        "counterfactual_status": audit_report["status"],
        "errors": errors,
        "differences": differences,
        "artifacts": {
            "log": str(log_path),
            "actual": str(actual_path),
            "audit": str(audit_path),
        },
    }
    write_json(report_path, report)
    if passed:
        print(
            "fold-168 counterfactual audit passed; "
            f"status={audit_report['status']}, wall time={wall_time_seconds:.3f}s.")
        return 0
    print(f"fold-168 counterfactual audit failed; see {report_path}.")
    for message in [*errors, *differences[:20]]:
        print(f"- {message}")
    return 1


if __name__ == "__main__":
    raise SystemExit(run())
