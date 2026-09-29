#!/usr/bin/env python3

from __future__ import annotations

import argparse
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Sequence

from fold_168_support import *


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the external fold-168 regression benchmark.")
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--map", dest="map_path", type=Path, required=True)
    parser.add_argument("--simulation-manifest", type=Path,
                        help="Generation record; defaults to <map>.simulation.json")
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
        "manifest": (args.simulation_manifest or Path(str(args.map_path) + ".simulation.json")).resolve(),
    }
    actual_hashes: dict[str, str] = {}
    actual = make_empty_actual(actual_hashes)
    command: list[str] = []
    errors: list[str] = []
    differences: list[str] = []
    wall_time_seconds: float | None = None
    log_text = ""
    scoring_status = "failed"
    gates = {
        "quality_gate": {"status": "uncalibrated", "passed": False},
        "iteration_gate": {"passed": None},
        "atom_cutoff_gate": {"passed": None},
    }

    try:
        baseline = load_baseline(args.baseline.resolve())
        executable = args.executable.resolve()
        if not executable.is_file():
            raise RegressionError(f"Benchmark executable does not exist: {executable}")
        actual_hashes = validate_input_hashes(input_paths, baseline["input_hashes"])
        actual["input_hashes"] = actual_hashes
        manifest = load_simulation_manifest(input_paths["manifest"], actual_hashes)
        validate_fixture(manifest, baseline)
        actual = make_empty_actual(actual_hashes)
        actual["generation_record"] = {name: value for name, value in manifest.items() if name != "atoms"}

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
                check=False,
                cwd=temp_dir)
            wall_time_seconds = time.perf_counter() - start_time
            log_text = completed.stdout.decode("utf-8", errors="replace")
            log_path.write_text(log_text, encoding="utf-8")
            if completed.returncode != 0:
                raise RegressionError(
                    f"Benchmark command exited with status {completed.returncode}.")
            if not temporary_database.is_file():
                raise RegressionError(
                    "Benchmark command did not create the temporary output database.")

            # Do not attach truth to an execution whose inputs changed during fitting.
            validate_input_hashes(input_paths, actual_hashes)
            paired = pair_atom_truth(read_atom_results(temporary_database), manifest)
            metrics = calculate_quality_metrics(paired)
            actual["atoms"] = paired
            actual["quality_metrics"] = metrics
            actual["diagnostics"] = {
                "maximum_absolute_offset": max(abs(atom["intercept_mdpde"]) for atom in paired),
                "fallback_charge_count": manifest["fallback_charge_count"],
            }
            scoring_status = "complete"
            actual["second_stage_summary"] = parse_second_stage_summary(log_text)
            actual["atom_cutoff_summary"] = parse_atom_cutoff_summary(log_text)
            actual["production_fitting"] = parse_final_state_certificate(log_text)

        gates = evaluate_gates(baseline, actual)
        differences = [message for gate in gates.values() for message in gate["differences"]]
    except Exception as error:  # Preserve artifacts for all benchmark failures.
        errors.append(str(error))

    log_path.write_text(log_text, encoding="utf-8")
    write_json(actual_path, actual)
    certificate = ((actual.get("production_fitting") or {}).get("certificate") or {})
    residuals = certificate.get("operator_nominal_p99", [])
    convergence_passed = (
        (actual["second_stage_summary"] or {}).get("stop_reason") == "converged" and
        certificate.get("status") == "evaluated" and certificate.get("complete") is True and
        certificate.get("qualified") is True and len(residuals) == 3 and
        all(value is not None and 0 <= value < 1e-4 for value in residuals) and
        gates["iteration_gate"]["passed"] is True and
        (actual.get("production_fitting") or {}).get("attempts", MAXIMUM_ACCEPTED_ITERATIONS + 1) <= MAXIMUM_ACCEPTED_ITERATIONS)
    passed = not errors and convergence_passed and all(gate["passed"] is True for gate in gates.values())
    report = {
        "schema_version": SCHEMA_VERSION,
        "passed": passed,
        "truth_scoring": {"status": scoring_status},
        "convergence_acceptance": {"passed": convergence_passed,
            "reason": "certified-within-25" if convergence_passed else "convergence-not-established-within-25"},
        **gates,
        "stop_reason": (actual["second_stage_summary"] or {}).get("stop_reason"),
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
