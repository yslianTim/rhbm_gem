#!/usr/bin/env python3
"""Run the FirstStage-calibrated A/B2/B4 research comparison."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path


SYNTHETIC_CASES = (
    "well-synthetic",
    "moderate-synthetic",
    "weak-synthetic",
    "noise-high-synthetic",
    "mixed-width-synthetic",
    "ac-correlated-synthetic",
)
SEEDS = (1001, 1002, 1003, 1004, 1005)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def run_group(runner: Path, dataset: str, output: Path, group: str, seed: int,
              model: Path | None = None, map_path: Path | None = None) -> tuple[dict | None, str | None]:
    command = [str(runner), dataset, str(output), group, str(seed)]
    if model is not None and map_path is not None:
        command.extend((str(model), str(map_path)))
    log_path = output.parent / f"{output.name}.log"
    output.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    if completed.returncode != 0:
        return None, f"exit={completed.returncode}; log={log_path}"
    if log_path.stat().st_size == 0:
        log_path.unlink()
    result_path = output / "group.json"
    if not result_path.is_file():
        return None, f"missing {result_path}"
    return read_json(result_path), None


def group_metrics(group: dict | None) -> dict | None:
    if group is None:
        return None
    keys = (
        "group", "runtime_convergence", "component_stop_reason", "candidate_available",
        "candidate_objective", "A_rmse", "B_rmse_angstrom", "C_rmse", "prediction_rmse",
        "truth_inside_bound_A_rmse", "truth_inside_bound_B_rmse_angstrom",
        "truth_inside_bound_C_rmse", "truth_outside_bound_count", "final_active_bound_atoms",
        "maximum_bound_violation_log_width", "projected_width_gradient_inf_norm",
        "boundary_gradient_sign_violation", "boundary_multiplier_max",
        "boundary_complementarity_residual", "independent_replay_check", "width_stationarity",
        "active_set_transition_count", "active_set_bound_hits", "active_set_bound_releases",
        "active_set_iterations", "active_set_subproblem_failures",
        "maximum_active_set_subproblem_kkt_inf_norm", "serial_1032_width_angstrom",
        "maximum_local_ac_condition", "maximum_ac_cancellation_ratio", "search_sweeps",
        "profile_evaluations", "wall_seconds", "total_wall_seconds", "process_peak_rss_bytes",
        "initial_width_projection_count", "maximum_initial_projection_angstrom",
        "truth_coverage_after_calibration",
    )
    return {key: group.get(key) for key in keys}


def audit_pair(groups: dict[str, dict | None]) -> dict:
    checks: dict[str, bool] = {}
    reasons: dict[str, str] = {}
    if any(groups.get(name) is None for name in ("A", "B2", "B4")):
        return {"status": "failed", "all_audits_pass": False,
                "checks": {"all_groups_produced_results": False},
                "failure_reasons": {name: "run-failed-or-result-missing" for name in ("A", "B2", "B4")
                                    if groups.get(name) is None}}

    a, b2, b4 = groups["A"], groups["B2"], groups["B4"]
    for field in ("input_fingerprint", "model_sha256", "map_sha256", "target_count",
                  "observation_count", "normalization_divisor", "chain_selection",
                  "pure_halo_support_exclusion", "first_stage_initialization_fingerprint",
                  "calibration", "problem_objective"):
        checks[f"paired_{field}"] = a.get(field) == b2.get(field) == b4.get(field)
    checks["A_uses_original_initialization"] = (
        a.get("effective_initialization_fingerprint") == a.get("first_stage_initialization_fingerprint"))
    checks["B2_uses_active_set_lm"] = b2.get("strategy") == "active-set-bounded-lm"
    checks["B4_uses_active_set_lm"] = b4.get("strategy") == "active-set-bounded-lm"
    calibration = a.get("calibration", {})
    p95 = calibration.get("p95")
    bmax2, bmax4 = calibration.get("bmax2"), calibration.get("bmax4")
    checks["calibration_available"] = calibration.get("status") == "Available" and p95 is not None
    checks["B2_bound_is_2P95"] = p95 is not None and bmax2 is not None and math.isclose(
        bmax2, 2 * p95, rel_tol=1e-14, abs_tol=0.0)
    checks["B4_bound_is_4P95"] = p95 is not None and bmax4 is not None and math.isclose(
        bmax4, 4 * p95, rel_tol=1e-14, abs_tol=0.0)
    for label, group, expected in (("B2", b2, bmax2), ("B4", b4, bmax4)):
        checks[f"{label}_uses_global_bound"] = expected is not None and math.isclose(
            group.get("global_width_upper_bound", math.nan), expected, rel_tol=1e-14, abs_tol=0.0)
        checks[f"{label}_candidate_available"] = group.get("candidate_available") is True
        violation = group.get("maximum_bound_violation_log_width")
        eta_max = group.get("global_eta_upper_bound")
        feasibility = 64 * sys.float_info.epsilon * max(1.0, abs(eta_max or 0.0))
        checks[f"{label}_bound_feasible"] = violation is not None and violation <= feasibility
        replay = group.get("independent_replay_check")
        checks[f"{label}_independent_replay_executed"] = isinstance(replay, dict) and replay.get("status") in (
            "Passed", "Failed")
        checks[f"{label}_independent_replay_passed"] = isinstance(replay, dict) and replay.get("status") == "Passed"
        stationarity = group.get("width_stationarity")
        checks[f"{label}_bound_aware_stationarity_executed"] = isinstance(stationarity, dict) and stationarity.get("status") in (
            "Passed", "Failed")
        checks[f"{label}_bound_aware_stationarity_passed"] = isinstance(stationarity, dict) and stationarity.get("status") == "Passed"
        checks[f"{label}_active_set_subproblems_solved"] = group.get("active_set_subproblem_failures") == 0
        checks[f"{label}_truth_not_used_for_calibration"] = group.get("truth_leakage") is False
        if not checks[f"{label}_independent_replay_passed"]:
            reasons[f"{label}_independent_replay"] = str(replay)
        if not checks[f"{label}_bound_feasible"]:
            reasons[f"{label}_bound_feasibility"] = f"violation={violation}; tolerance={feasibility}"
    checks["truth_coverage_scored_after_calibration"] = all(
        group.get("truth_leakage") is False for group in (a, b2, b4))
    return {"status": "passed" if all(checks.values()) else "failed",
            "all_audits_pass": all(checks.values()), "checks": checks, "failure_reasons": reasons}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--simulation-map", type=Path, required=True)
    parser.add_argument("--unsharpened-map", type=Path, required=True)
    parser.add_argument("--phase", choices=("synthetic", "real", "all"), default="all")
    args = parser.parse_args()
    runner = args.runner.resolve()
    root = args.output.resolve()
    if root.exists():
        raise SystemExit(f"Output path already exists: {root}")
    for path in (runner, args.model, args.simulation_map, args.unsharpened_map):
        if not path.is_file():
            raise SystemExit(f"Required input does not exist: {path}")
    root.mkdir(parents=True)
    failures: list[dict] = []
    synthetic_rows: list[dict] = []

    if args.phase in ("synthetic", "all"):
        for case in SYNTHETIC_CASES:
            for seed in SEEDS:
                print(f"Synthetic A: {case}, seed {seed}", flush=True)
                case_key = f"synthetic/{case}/seed-{seed}"
                group, error = run_group(runner, case, root / case_key / "A", "A", seed)
                if error:
                    failures.append({"case": case_key, "group": "A", "error": error})
                else:
                    synthetic_rows.append({"case": case, "seed": seed, "A": group_metrics(group),
                        "calibration": group.get("calibration"),
                        "B2_status": "CalibrationUnavailable; no FirstStage sampling/MDPDE provenance",
                        "B4_status": "CalibrationUnavailable; no FirstStage sampling/MDPDE provenance",
                        "truth_coverage2": None, "truth_coverage4": None,
                        "truth_leakage": group.get("truth_leakage")})

    real_rows: list[dict] = []
    real_cases = (
        ("6z6u-simulation", args.simulation_map),
        ("6z6u-unsharpened", args.unsharpened_map),
    ) if args.phase in ("real", "all") else ()
    for dataset, map_path in real_cases:
        print(f"Paired {dataset}: A / B2 / B4", flush=True)
        case_dir = root / dataset
        groups: dict[str, dict | None] = {}
        group_errors: dict[str, str] = {}
        for label in ("A", "B2", "B4"):
            output = case_dir / label
            group, error = run_group(runner, dataset, output, label, 0, args.model, map_path)
            groups[label] = group
            if error:
                group_errors[label] = error
                failures.append({"case": dataset, "group": label, "error": error})
        for label in ("B2", "B4"):
            group = groups[label]
            if group and group.get("initial_width_projection_count", 0) > 0:
                auxiliary = f"Aproj{2 if label == 'B2' else 4}"
                print(f"Auxiliary projected-initialization control: {dataset}, {auxiliary}", flush=True)
                _, error = run_group(runner, dataset, case_dir / auxiliary, auxiliary, 0, args.model, map_path)
                if error:
                    failures.append({"case": dataset, "group": auxiliary, "error": error})
        audit = audit_pair(groups)
        audit["run_errors"] = group_errors
        (case_dir / "pair-audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
        real_rows.append({"dataset": dataset,
            "calibration": groups["A"].get("calibration") if groups.get("A") else None,
            "truth_coverage": groups["A"].get("truth_coverage_after_calibration") if groups.get("A") else None,
            "groups": {label: group_metrics(groups[label]) for label in ("A", "B2", "B4")},
            "paired_audit": audit})
        if not audit.get("all_audits_pass"):
            failures.append({"case": dataset, "group": "A/B2/B4", "error": "paired-audit-failed"})

    real_pairs_passed = bool(real_rows) and all(
        row["paired_audit"].get("all_audits_pass") is True for row in real_rows)
    summary = {"contract": "global-firststage-p95-width-bounds-active-set-lm-v1", "phase": args.phase,
        "formal_groups": ["A", "B2", "B4"],
        "synthetic_seeds": list(SEEDS), "synthetic": synthetic_rows,
        "real_map_cases": real_rows,
        "auxiliary_projection_control": "Aproj2/Aproj4 only when their bounded initialization projects at least one width",
        "failures": failures,
        "all_executed_runs_passed": not failures,
        "real_map_paired_comparisons_passed": real_pairs_passed,
        "synthetic_A_B2_B4_comparison_complete": False}
    (root / "paired-summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"Saved compact experiment summary: {root / 'paired-summary.json'}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
