"""Shared implementation for the FixedNeighbor outer-core qualification phases.

The canonical CLI in joint_fixed_neighbor_outer_core_qualification.py is the
only public entry point; phase functions below are internal support.
"""
from __future__ import annotations

import csv
import json
import math
import statistics
import subprocess
import time
from pathlib import Path

from experiment_io import ROOT, read, sha, write
from experiment_process import RSS_LIMIT_BYTES, monitored
from experiment_provenance import source_hash

screen_DEFAULT_CASES = ("chain-512", "cube-512", "cube-1024")
screen_HISTORICAL_CORE_SIZES = (64, 128, 256)
# Keep the historical artifact reproducible while allowing every campaign to
# choose its own explicit outer-core list.
screen_CORE_SIZES = screen_HISTORICAL_CORE_SIZES
screen_POLICY = "FixedNeighbor"
screen_MODE = "--scaling-only"


def screen_finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def screen_parse_core_sizes(value):
    if isinstance(value, (tuple, list)):
        values = tuple(int(item) for item in value)
    else:
        values = tuple(int(item.strip()) for item in str(value).split(",") if item.strip())
    if not values or any(item <= 0 for item in values) or len(set(values)) != len(values):
        raise ValueError("core sizes must be a non-empty comma-separated list of unique positive integers")
    return values


def screen_case(value):
    topology, atoms = value.split("-", 1)
    atoms = int(atoms)
    if topology not in ("chain", "cube") or atoms not in (512, 1024):
        raise ValueError(f"Unsupported core-size case: {value}")
    return topology, atoms


def screen_status(process):
    if process.get("status") == "completed":
        return "completed"
    if process.get("status") == "time-limit":
        return "timeout"
    if process.get("status") in ("rss-limit", "rss-limit-observed-after-exit"):
        return "rss_limit"
    return "process_error"


def screen_partial(progress, core_atoms):
    if not progress:
        return None
    fixed = dict(progress)
    fixed.setdefault("search_converged", False)
    fixed.setdefault("search_reason", "interrupted-before-search-completion")
    fixed.setdefault("sweeps", len(progress.get("sweep_telemetry", [])))
    fixed.setdefault("confirmed_stationarity_sweep", None)
    fixed.setdefault("block_solves", sum(row.get("block_solves", 0)
                                          for row in progress.get("sweep_telemetry", [])))
    fixed.setdefault("profile_evaluations", sum(row.get("profile_evaluations", 0)
                                                 for row in progress.get("sweep_telemetry", [])))
    fixed.setdefault("accepted_local_updates", 0)
    fixed.setdefault("local_factor_seconds", None)
    fixed.setdefault("search_seconds", sum(row.get("wall_seconds", 0.0)
                                            for row in progress.get("sweep_telemetry", [])))
    fixed.setdefault("total_elapsed_seconds", fixed["search_seconds"])
    fixed.setdefault("final_global_ac_kkt", None)
    fixed.setdefault("final_raw_width_gradient_inf_norm", None)
    fixed.setdefault("maximum_local_columns", None)
    return {"topology": progress.get("topology"), "atoms": progress.get("atoms"),
            "core_atoms": core_atoms, "measurement_scope": "fixed-neighbor-search-only",
            "fixed_neighbor": fixed, "peak_rss_mb": None}


def screen_run_case(args, output_dir, topology, atoms, core_atoms):
    case = f"{topology}-{atoms}-core{core_atoms}"
    individual = output_dir / "individual-results"
    raw_path = individual / f"{case}.json"
    wrapper_path = individual / f"{case}-run.json"
    process_dir = individual / f"{case}-process"
    if wrapper_path.is_file():
        return read(wrapper_path)
    command = [str(args.build_dir / "bin" / "joint_fixed_neighbor_experiment"), screen_MODE,
               str(raw_path), topology, str(atoms), str(core_atoms)]
    print(f"Running {case} with {args.timeout:g}s cap", flush=True)
    process = monitored(command, process_dir, time.monotonic() + args.timeout,
                        rss_limit=args.rss_limit, seconds=args.timeout)
    result = read(raw_path) if raw_path.is_file() else None
    if result is None:
        progress_path = Path(str(raw_path) + ".progress.json")
        result = screen_partial(read(progress_path) if progress_path.is_file() else None, core_atoms)
    report = {"topology": topology, "atoms": atoms, "core_atoms": core_atoms,
              "policy": screen_POLICY, "measurement_scope": "fixed-neighbor-search-only",
              "command": command, "process": process, "status": screen_status(process),
              "completed": process.get("status") == "completed" and result is not None,
              "result": result}
    write(wrapper_path, report)
    return report


def screen_profile_work(fixed):
    profile = fixed.get("local_profile_work")
    if isinstance(profile, dict):
        return profile
    return (fixed.get("fixed_neighbor_work") or {}).get("local_profile_work", {})


def screen_block_geometry(result, fixed, fallback_core):
    records = fixed.get("block_telemetry") or []
    sizes = {}
    for record in records:
        block = record.get("block")
        size = record.get("core_atoms", record.get("atoms"))
        if block is not None and screen_finite(size):
            sizes.setdefault(block, int(size))
    realized_sizes = list(sizes.values())
    requested = fixed.get("outer_core_atoms", result.get("outer_core_atoms", fallback_core))
    return {
        "requested_outer_core_atoms": requested,
        "realized_outer_block_count": len(realized_sizes) or fixed.get("blocks_per_sweep"),
        "realized_outer_core_atoms_minimum": min(realized_sizes) if realized_sizes else None,
        "realized_outer_core_atoms_mean": (sum(realized_sizes) / len(realized_sizes)
                                           if realized_sizes else None),
        "realized_outer_core_atoms_maximum": max(realized_sizes) if realized_sizes else None,
    }


def screen_maximum_sweep_value(fixed, field):
    values = [sweep.get(field) for sweep in fixed.get("sweep_telemetry", [])]
    values = [value for value in values if screen_finite(value)]
    return max(values) if values else None


def screen_work_value(work, fixed, field, fallback=None):
    value = work.get(field)
    if screen_finite(value):
        return value
    value = fixed.get(field, fallback)
    return value if screen_finite(value) else fallback


def screen_ratio(numerator, denominator):
    return numerator / denominator if screen_finite(numerator) and screen_finite(denominator) and denominator > 0 else None


def screen_summarize(report):
    result = report.get("result") or {}
    fixed = result.get("fixed_neighbor", result)
    process = report.get("process") or {}
    work = fixed.get("fixed_neighbor_work") or {}
    profile = screen_profile_work(fixed)
    profile_total = profile.get("total") or {}
    search_seconds = fixed.get("search_seconds")
    sweeps = fixed.get("sweeps", len(fixed.get("sweep_telemetry", [])))
    block_solves = fixed.get("block_solves", sum(s.get("block_solves", 0)
                                                  for s in fixed.get("sweep_telemetry", [])))
    profile_evaluations = fixed.get("profile_evaluations", sum(
        s.get("profile_evaluations", 0) for s in fixed.get("sweep_telemetry", [])))
    local_search_seconds = screen_work_value(work, fixed, "local_search_seconds",
                                       fixed.get("local_factor_seconds"))
    candidate_replay_seconds = screen_work_value(work, fixed, "candidate_replay_seconds")
    sweep_replay_seconds = screen_work_value(work, fixed, "sweep_replay_seconds")
    sweep_global_state_seconds = screen_work_value(work, fixed, "sweep_global_state_seconds")
    objective = fixed.get("objective")
    maximum_cache_replay_error = screen_maximum_sweep_value(fixed, "cache_replay_error")
    maximum_objective_replay_error = screen_maximum_sweep_value(fixed, "objective_replay_error")
    cache_replay_limit = 2e-12
    objective_replay_limit = (1e-12 + 2e-12 * abs(objective)
                               if screen_finite(objective) else None)
    row = {
        "topology": report.get("topology"), "atoms": report.get("atoms"),
        "core_atoms": report.get("core_atoms", result.get("core_atoms")),
        "policy": report.get("policy", screen_POLICY), "status": report.get("status", "completed"),
        "measurement_scope": report.get("measurement_scope", "fixed-neighbor-search-only"),
        "search_converged": fixed.get("search_converged"),
        "search_reason": fixed.get("search_reason"),
        "sweeps": sweeps,
        "confirmed_stationarity_sweep": fixed.get("confirmed_stationarity_sweep"),
        "first_order_stationarity_sweep": fixed.get("first_order_stationarity_sweep"),
        "confirmation_extra_sweeps": fixed.get("confirmation_extra_sweeps"),
        "block_solves": block_solves,
        "accepted_blocks": fixed.get("accepted_blocks", sum(s.get("accepted_blocks", 0)
                                                              for s in fixed.get("sweep_telemetry", []))),
        "profile_evaluations": profile_evaluations,
        "accepted_local_updates": fixed.get("accepted_local_updates"),
        "local_work_seconds": local_search_seconds,
        "local_search_seconds": local_search_seconds,
        "profile_factor_seconds": fixed.get("profile_factor_seconds"),
        "search_seconds": search_seconds,
        "total_seconds": fixed.get("total_elapsed_seconds"),
        "peak_rss_mb": result.get("peak_rss_mb", process.get("sampled_tree_peak_rss_bytes", 0) / 1024**2
                       if process.get("sampled_tree_peak_rss_bytes") is not None else None),
        "objective": objective,
        "final_global_ac_kkt": fixed.get("final_global_ac_kkt"),
        "final_raw_width_gradient_inf_norm": fixed.get("final_raw_width_gradient_inf_norm"),
        "maximum_cache_replay_error": maximum_cache_replay_error,
        "maximum_objective_replay_error": maximum_objective_replay_error,
        "cache_replay_limit": cache_replay_limit,
        "objective_replay_limit": objective_replay_limit,
        "cache_replay_within_limit": (maximum_cache_replay_error is not None and
                                       maximum_cache_replay_error <= cache_replay_limit),
        "objective_replay_within_limit": (maximum_objective_replay_error is not None and
                                           objective_replay_limit is not None and
                                           maximum_objective_replay_error <= objective_replay_limit),
        "maximum_local_columns": fixed.get("maximum_local_columns"),
        "maximum_local_rows": fixed.get("maximum_local_rows"),
        "old_core_seconds": screen_work_value(work, fixed, "old_core_seconds"),
        "effective_response_seconds": screen_work_value(work, fixed, "effective_response_seconds"),
        "local_state_seconds": screen_work_value(work, fixed, "local_state_seconds"),
        "candidate_copy_seconds": screen_work_value(work, fixed, "candidate_copy_seconds"),
        "candidate_replay_seconds": candidate_replay_seconds,
        "cache_update_seconds": screen_work_value(work, fixed, "cache_update_seconds"),
        "sweep_replay_seconds": sweep_replay_seconds,
        "sweep_global_state_seconds": sweep_global_state_seconds,
        "derivative_prepare_seconds": profile_total.get("derivative_prepare_seconds"),
        "derivative_reduce_seconds": profile_total.get("derivative_reduce_seconds"),
        "derivative_jacobian_qr_seconds": profile_total.get("derivative_jacobian_qr_seconds"),
        "tiled_qr_assembly_copy_seconds": profile_total.get("tiled_qr_assembly_copy_seconds"),
        "tiled_qr_householder_seconds": profile_total.get("tiled_qr_householder_seconds"),
        "tiled_qr_rhs_transform_seconds": profile_total.get("tiled_qr_rhs_transform_seconds"),
    }
    row.update(screen_block_geometry(result, fixed, row["core_atoms"]))
    row.update({
        "search_seconds_per_sweep": screen_ratio(search_seconds, sweeps),
        "local_search_seconds_per_block_solve": screen_ratio(local_search_seconds, block_solves),
        "candidate_replay_seconds_per_block_solve": screen_ratio(candidate_replay_seconds, block_solves),
        "profile_evaluations_per_block_solve": screen_ratio(profile_evaluations, block_solves),
        "accepted_updates_per_sweep": screen_ratio(fixed.get("accepted_local_updates"), sweeps),
        "local_search_fraction": screen_ratio(local_search_seconds, search_seconds),
        "candidate_replay_fraction": screen_ratio(candidate_replay_seconds, search_seconds),
        "global_overhead_fraction": screen_ratio(
            (sweep_replay_seconds or 0.0) + (sweep_global_state_seconds or 0.0), search_seconds),
    })
    return row


def screen_correct(row):
    return row["status"] == "completed" and row["search_converged"] is True and \
        row.get("search_reason") == "block-stationary" and \
        screen_finite(row.get("final_global_ac_kkt")) and row["final_global_ac_kkt"] <= 1e-10 and \
        screen_finite(row.get("final_raw_width_gradient_inf_norm")) and \
        row["final_raw_width_gradient_inf_norm"] <= 1e-12 and \
        row.get("cache_replay_within_limit") is True and \
        row.get("objective_replay_within_limit") is True


def screen_analyze(rows, core_sizes=screen_CORE_SIZES):
    core_sizes = screen_parse_core_sizes(core_sizes)
    core_order = {core: index for index, core in enumerate(core_sizes)}
    groups = {}
    for row in rows:
        groups.setdefault((row["topology"], row["atoms"]), []).append(row)
    cases = []
    for key, entries in sorted(groups.items()):
        entries = sorted(entries, key=lambda row: core_order.get(row.get("core_atoms"), len(core_order)))
        correct = [row["core_atoms"] for row in entries if row["core_atoms"] in core_order and screen_correct(row)]
        observed = {row["core_atoms"] for row in entries}
        fastest = min((row for row in entries if row.get("core_atoms") in core_order and screen_correct(row) and
                       screen_finite(row.get("search_seconds"))),
                      key=lambda row: row["search_seconds"], default=None)
        cases.append({"topology": key[0], "atoms": key[1], "cores": entries,
                      "observed_core_sizes": sorted(observed & set(core_sizes), key=core_order.get),
                      "correct_core_sizes": correct,
                      "fastest_correct_core": fastest["core_atoms"] if fastest else None})
    complete = bool(cases) and all(set(case["observed_core_sizes"]) == set(core_sizes) and
                                   set(case["correct_core_sizes"]) == set(core_sizes)
                                   for case in cases)
    aggregate = {}
    for core in core_sizes:
        values = [row["search_seconds"] for case in cases for row in case["cores"]
                  if row.get("core_atoms") == core and screen_correct(row) and
                  screen_finite(row.get("search_seconds"))]
        aggregate[str(core)] = sum(values) if len(values) == len(cases) else None
    selected = min((int(core) for core, value in aggregate.items() if screen_finite(value) and value > 0),
                   key=lambda core: aggregate[str(core)], default=None) if complete else None
    return {"phase": "FixedNeighbor outer-core performance study",
            "policy": screen_POLICY, "core_sizes": list(core_sizes), "cases": cases,
            "correctness_gate": "passed" if complete else "failed",
            "aggregate_search_seconds_by_core": aggregate,
            "selected_core_size": selected,
            "selection_basis": "minimum aggregate search seconds among cores correct on every study case; resource, geometry, and attribution metrics remain reported"}


def screen_write_outputs(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = ["topology", "atoms", "core_atoms", "policy", "status", "sweeps",
               "first_order_stationarity_sweep", "confirmed_stationarity_sweep",
               "confirmation_extra_sweeps", "block_solves", "accepted_blocks",
               "profile_evaluations", "accepted_local_updates", "requested_outer_core_atoms",
               "realized_outer_block_count", "realized_outer_core_atoms_minimum",
               "realized_outer_core_atoms_mean", "realized_outer_core_atoms_maximum",
               "local_work_seconds", "local_search_seconds", "profile_factor_seconds",
               "old_core_seconds", "effective_response_seconds", "local_state_seconds",
               "candidate_copy_seconds", "candidate_replay_seconds", "cache_update_seconds",
               "sweep_replay_seconds", "sweep_global_state_seconds", "derivative_prepare_seconds",
               "derivative_reduce_seconds", "derivative_jacobian_qr_seconds",
               "tiled_qr_assembly_copy_seconds", "tiled_qr_householder_seconds",
               "tiled_qr_rhs_transform_seconds", "search_seconds", "search_seconds_per_sweep",
               "local_search_seconds_per_block_solve", "candidate_replay_seconds_per_block_solve",
               "profile_evaluations_per_block_solve", "accepted_updates_per_sweep",
               "local_search_fraction", "candidate_replay_fraction", "global_overhead_fraction",
               "total_seconds", "peak_rss_mb", "search_converged", "search_reason", "objective",
               "final_global_ac_kkt", "final_raw_width_gradient_inf_norm",
               "maximum_cache_replay_error", "maximum_objective_replay_error",
               "cache_replay_within_limit", "objective_replay_within_limit",
               "maximum_local_rows", "maximum_local_columns"]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            for row in case["cores"]:
                writer.writerow({column: row.get(column) for column in columns})


def screen_run_campaign(args):
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = [screen_case(value) for value in args.cases]
    core_sizes = screen_parse_core_sizes(args.core_sizes)
    reports = []
    for topology, atoms in cases:
        for core_atoms in core_sizes:
            reports.append(screen_run_case(args, output_dir, topology, atoms, core_atoms))
            write(output_dir / "runs.json", reports)
    rows = [screen_summarize(report) for report in reports if report.get("result") is not None]
    analysis = screen_analyze(rows, core_sizes)
    manifest = {
        "schema_version": 1, "phase": analysis["phase"],
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {"driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
                                    "support": sha(Path(__file__))},
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "backend": "EIGEN", "eigen_threads": 1, "cases": [f"{t}-{a}" for t, a in cases],
        "policy": screen_POLICY, "core_sizes": list(core_sizes),
        "measurement_scope": "fixed-neighbor-search-only",
        "resource_envelope": {"wall_seconds": args.timeout, "rss_bytes": args.rss_limit},
        "correctness_contract": ["search_converged", "A/C KKT <= 1e-10",
                                 "width gradient <= 1e-12"],
        "analysis": {"correctness_gate": analysis["correctness_gate"],
                     "selected_core_size": analysis["selected_core_size"]},
    }
    write(output_dir / "campaign-manifest.json", manifest)
    screen_write_outputs(analysis, output_dir)
    (output_dir / "README.md").write_text(
        "# FixedNeighbor outer-core performance study\n\n"
        "This search-only study keeps FixedNeighbor's intrinsic one-accepted local-update "
        "contract, forward serial Gauss-Seidel order, EIGEN backend, one Eigen thread, and "
        "frozen stationarity checks. It varies only outer core_atoms over "
        f"{', '.join(str(core) for core in core_sizes)} for "
        f"{', '.join(f'{topology}-{atoms}' for topology, atoms in cases)}.\n\n"
        "The selected core is the minimum aggregate search time among core sizes that are "
        "search-correct on every study case. Local work, profile factor telemetry, RSS, "
        "and maximum local columns remain in the machine-readable artifact.\n\n"
        f"Correctness gate: **{analysis['correctness_gate']}**. Selected core: "
        f"**{analysis['selected_core_size'] or 'none'}**.\n"
    )
    return analysis

endpoint_DEFAULT_CASES = ("chain-256", "cube-256", "chain-512", "cube-512")
endpoint_FINALIST_CORES = (12, 16, 64)
endpoint_MODE = "--neighbor-only"
endpoint_KKT_LIMIT = 1e-10
endpoint_WIDTH_GRADIENT_LIMIT = 1e-12
endpoint_REPLAY_CACHE_LIMIT = 2e-12
endpoint_REPLAY_OBJECTIVE_BASE = 1e-12
endpoint_REPLAY_OBJECTIVE_SCALE = 2e-12
endpoint_PARITY_LIMIT = 1e-10


def endpoint_case(value):
    topology, atoms = value.split("-", 1)
    atoms = int(atoms)
    if topology not in ("chain", "cube") or atoms not in (256, 512):
        raise ValueError(f"Unsupported endpoint case: {value}")
    return topology, atoms


def endpoint_check_value(fixed, name):
    value = fixed.get(name)
    return value.get("value") if isinstance(value, dict) else value


def endpoint_check_status(fixed, name):
    value = fixed.get(name)
    return value.get("status") if isinstance(value, dict) else value


def endpoint_inf_difference(lhs, rhs):
    if not isinstance(lhs, list) or not isinstance(rhs, list) or len(lhs) != len(rhs):
        return None
    return max((abs(left - right) for left, right in zip(lhs, rhs)), default=0.0)


def endpoint_scaled_ac_difference(candidate, control):
    beta = candidate.get("final_beta")
    reference_beta = control.get("final_beta")
    weights = candidate.get("final_ac_scaling_weights")
    reference_weights = control.get("final_ac_scaling_weights")
    if not all(isinstance(value, list) for value in (beta, reference_beta, weights, reference_weights)):
        return None
    if not (len(beta) == len(reference_beta) == len(weights) == len(reference_weights)):
        return None
    return max((abs(left - right) * max(weight, reference_weight)
                for left, right, weight, reference_weight in
                zip(beta, reference_beta, weights, reference_weights)), default=0.0)


def endpoint_rank_evidence(fixed):
    assessment = fixed.get("endpoint_assessment") or {}
    return {
        "projected_width_rank": (assessment.get("projected_width") or {}).get("rank"),
        "corrected_jacobian_rank": (assessment.get("corrected_jacobian") or {}).get("rank"),
        "normalized_width_rank": (assessment.get("normalized_width") or {}).get("rank"),
    }


def endpoint_run_case(args, output_dir, topology, atoms, core_atoms):
    case = f"{topology}-{atoms}-core{core_atoms}"
    individual = output_dir / "individual-results"
    raw_path = individual / f"{case}.json"
    wrapper_path = individual / f"{case}-run.json"
    process_dir = individual / f"{case}-process"
    if wrapper_path.is_file():
        return read(wrapper_path)
    command = [str(args.build_dir / "bin" / "joint_fixed_neighbor_experiment"), endpoint_MODE,
               str(raw_path), topology, str(atoms), str(core_atoms)]
    print(f"Running {case} with {args.timeout:g}s cap", flush=True)
    process = monitored(command, process_dir, time.monotonic() + args.timeout,
                        rss_limit=args.rss_limit, seconds=args.timeout)
    result = read(raw_path) if raw_path.is_file() else None
    if result is None:
        progress_path = Path(str(raw_path) + ".progress.json")
        result = screen_partial(read(progress_path) if progress_path.is_file() else None, core_atoms)
    report = {"topology": topology, "atoms": atoms, "core_atoms": core_atoms,
              "command": command, "process": process, "status": screen_status(process),
              "completed": process.get("status") == "completed" and result is not None,
              "result": result}
    write(wrapper_path, report)
    return report


def endpoint_summarize_endpoint(report):
    result = report.get("result") or {}
    fixed = result.get("fixed_neighbor", result)
    row = screen_summarize({**report, "policy": screen_POLICY,
                     "measurement_scope": "fixed-neighbor-endpoint-assessed"})
    global_kkt = endpoint_check_value(fixed, "global_kkt")
    objective = fixed.get("objective")
    objective_limit = (endpoint_REPLAY_OBJECTIVE_BASE + endpoint_REPLAY_OBJECTIVE_SCALE * abs(objective)
                       if screen_finite(objective) else None)
    maximum_cache_replay_error = screen_maximum_sweep_value(fixed, "cache_replay_error")
    maximum_objective_replay_error = screen_maximum_sweep_value(fixed, "objective_replay_error")
    row.update({
        "measurement_scope": "fixed-neighbor-endpoint-assessed",
        "global_ac_kkt": global_kkt,
        "width_gradient_inf_norm": fixed.get("width_gradient_inf_norm"),
        "endpoint_certified": fixed.get("endpoint_certified"),
        "runtime_convergence": fixed.get("runtime_convergence"),
        "assessment_inner": endpoint_check_status(fixed, "assessment_inner"),
        "assessment_gradient": endpoint_check_status(fixed, "assessment_gradient"),
        "assessment_local": endpoint_check_status(fixed, "assessment_local"),
        "assessment_identified": endpoint_check_status(fixed, "assessment_identified"),
        "endpoint_trust": ((fixed.get("endpoint_assessment") or {}).get("endpoint_trust") or {}).get("passed"),
        "maximum_cache_replay_error": maximum_cache_replay_error,
        "maximum_objective_replay_error": maximum_objective_replay_error,
        "cache_replay_within_limit": (maximum_cache_replay_error is not None and
                                       maximum_cache_replay_error <= endpoint_REPLAY_CACHE_LIMIT),
        "objective_replay_within_limit": (maximum_objective_replay_error is not None and
                                           objective_limit is not None and
                                           maximum_objective_replay_error <= objective_limit),
        "final_eta": fixed.get("final_eta"),
        "final_beta": fixed.get("final_beta"),
        "final_ac_scaling_weights": fixed.get("final_ac_scaling_weights"),
        "rank_evidence": endpoint_rank_evidence(fixed),
    })
    return row


def endpoint_endpoint_gate(row):
    return row.get("status") == "completed" and row.get("search_converged") is True and \
        row.get("search_reason") == "block-stationary" and \
        isinstance(row.get("confirmed_stationarity_sweep"), int) and \
        0 < row["confirmed_stationarity_sweep"] <= row.get("sweeps", 0) and \
        screen_finite(row.get("global_ac_kkt")) and row["global_ac_kkt"] <= endpoint_KKT_LIMIT and \
        screen_finite(row.get("width_gradient_inf_norm")) and \
        row["width_gradient_inf_norm"] <= endpoint_WIDTH_GRADIENT_LIMIT and \
        row.get("cache_replay_within_limit") is True and \
        row.get("objective_replay_within_limit") is True and \
        row.get("endpoint_certified") is True and \
        row.get("runtime_convergence") == "Passed" and \
        all(row.get(name) == "Passed" for name in
            ("assessment_inner", "assessment_gradient", "assessment_local", "assessment_identified")) and \
        row.get("endpoint_trust") is True


def endpoint_parity(candidate, control):
    objective_difference = (abs(candidate["objective"] - control["objective"])
                            if screen_finite(candidate.get("objective")) and screen_finite(control.get("objective"))
                            else None)
    eta_difference = endpoint_inf_difference(candidate.get("final_eta"), control.get("final_eta"))
    ac_difference = endpoint_scaled_ac_difference(candidate, control)
    return {
        "objective_abs_difference": objective_difference,
        "eta_inf_difference": eta_difference,
        "scaled_ac_inf_difference": ac_difference,
        "parameter_parity": eta_difference is not None and ac_difference is not None and
            eta_difference <= endpoint_PARITY_LIMIT and ac_difference <= endpoint_PARITY_LIMIT,
        "rank_evidence": {"candidate": candidate.get("rank_evidence"),
                          "control": control.get("rank_evidence")},
    }


def endpoint_analyze(rows, core_sizes=endpoint_FINALIST_CORES):
    core_sizes = screen_parse_core_sizes(core_sizes)
    groups = {}
    for row in rows:
        groups.setdefault((row["topology"], row["atoms"]), []).append(row)
    cases = []
    for key, entries in sorted(groups.items()):
        entries = sorted(entries, key=lambda row: core_sizes.index(row["core_atoms"]))
        control = next((row for row in entries if row["core_atoms"] == 64), None)
        correct = [row["core_atoms"] for row in entries if endpoint_endpoint_gate(row)]
        parity = {str(row["core_atoms"]): endpoint_parity(row, control)
                  for row in entries if control and row["core_atoms"] != 64}
        cases.append({"topology": key[0], "atoms": key[1], "cores": entries,
                      "correct_core_sizes": correct,
                      "parity_vs_core64": parity})
    complete = bool(cases) and all(set(case["correct_core_sizes"]) == set(core_sizes)
                                   for case in cases)
    qualified = [core for core in core_sizes if core != 64 and complete and
                 all(core in case["correct_core_sizes"] for case in cases)]
    return {
        "phase": "C3 FixedNeighbor outer-core full endpoint qualification",
        "core_sizes": list(core_sizes), "cases": cases,
        "qualification_gate": "passed" if complete else "failed",
        "qualified_core_sizes": qualified,
        "selection_basis": "all finalist cases pass search, replay, endpoint assessment, trust, and RuntimeConvergence; endpoint parity versus core64 is reported with existing 1e-10 parameter tolerance",
    }


def endpoint_write_outputs(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = ["topology", "atoms", "core_atoms", "status", "search_converged", "search_reason",
               "sweeps", "first_order_stationarity_sweep", "confirmed_stationarity_sweep",
               "confirmation_extra_sweeps", "search_seconds", "total_seconds", "peak_rss_mb",
               "global_ac_kkt", "width_gradient_inf_norm", "endpoint_certified",
               "runtime_convergence", "assessment_inner", "assessment_gradient", "assessment_local",
               "assessment_identified", "endpoint_trust", "maximum_cache_replay_error",
               "maximum_objective_replay_error", "cache_replay_within_limit",
               "objective_replay_within_limit", "requested_outer_core_atoms",
               "realized_outer_block_count", "realized_outer_core_atoms_minimum",
               "realized_outer_core_atoms_mean", "realized_outer_core_atoms_maximum"]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            for row in case["cores"]:
                writer.writerow({column: row.get(column) for column in columns})


def endpoint_run_campaign(args):
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = [endpoint_case(value) for value in args.cases]
    core_sizes = screen_parse_core_sizes(args.core_sizes)
    reports = []
    for topology, atoms in cases:
        for core_atoms in core_sizes:
            reports.append(endpoint_run_case(args, output_dir, topology, atoms, core_atoms))
            write(output_dir / "runs.json", reports)
    rows = [endpoint_summarize_endpoint(report) for report in reports if report.get("result") is not None]
    analysis = endpoint_analyze(rows, core_sizes)
    manifest = {
        "schema_version": 1, "phase": analysis["phase"],
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {"driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
                                    "support": sha(Path(__file__))},
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "backend": "EIGEN", "eigen_threads": 1,
        "cases": [f"{topology}-{atoms}" for topology, atoms in cases],
        "core_sizes": list(core_sizes), "local_search": "LegacyCompact",
        "block_order": "Forward",
        "maximum_sweeps": 30, "measurement_scope": "fixed-neighbor-endpoint-assessed",
        "resource_envelope": {"wall_seconds": args.timeout, "rss_bytes": args.rss_limit},
        "correctness_contract": ["block-stationary", "AC KKT <= 1e-10",
                                 "width gradient <= 1e-12", "eta confirmation <= 1e-10",
                                 "cache/objective replay within existing limits",
                                 "endpoint assessment/trust", "RuntimeConvergence"],
        "analysis": {"qualification_gate": analysis["qualification_gate"],
                     "qualified_core_sizes": analysis["qualified_core_sizes"]},
    }
    write(output_dir / "campaign-manifest.json", manifest)
    endpoint_write_outputs(analysis, output_dir)
    (output_dir / "README.md").write_text(
        "# FixedNeighbor outer-core full endpoint qualification\n\n"
        "This C3 campaign qualifies finalist outer cores "
        f"{', '.join(str(core) for core in core_sizes)} on "
        f"{', '.join(f'{topology}-{atoms}' for topology, atoms in cases)}. "
        "It keeps LegacyCompact, the intrinsic one-accepted local-update contract, "
        "Forward serial Gauss-Seidel, "
        "EIGEN, one Eigen thread, maximum_sweeps=30, and all existing thresholds.\n\n"
        "Core 64 is the historical/control comparison. Endpoint parity reports objective, eta, "
        "scaled A/C, and rank evidence; it does not require bitwise-identical trajectories.\n\n"
        f"Qualification gate: **{analysis['qualification_gate']}**. Qualified cores: "
        f"**{', '.join(str(core) for core in analysis['qualified_core_sizes']) or 'none'}**.\n"
    )
    return analysis

frontier_DEFAULT_CASES = ("chain-512", "cube-512", "chain-1024", "cube-1024")
frontier_DEFAULT_FINALIST_CORE_SIZES = (12, 16)
frontier_CONTROL_CORE_SIZE = 64
frontier_WARMUP_REPETITIONS = 1
frontier_MEASURED_REPETITIONS = 3
frontier_MIN_SEARCH_IMPROVEMENT = 0.10
frontier_RSS_MATERIAL_REGRESSION_FACTOR = 1.25
frontier_TOTAL_MATERIAL_REGRESSION_FACTOR = 1.25
frontier_MAX_SWEEPS_WARNING = 24


def frontier_finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def frontier_median(values):
    finite = [value for value in values if frontier_finite(value)]
    return statistics.median(finite) if finite else None


def frontier_maximum(values):
    finite = [value for value in values if frontier_finite(value)]
    return max(finite) if finite else None


def frontier_case_list(values):
    cases = [screen_case(value) for value in values]
    if len(set(cases)) != len(cases):
        raise ValueError("cases must be unique")
    return cases


def frontier_validate_finalists(finalists, control):
    finalists = screen_parse_core_sizes(finalists)
    if control <= 0 or control in finalists:
        raise ValueError("control core must be positive and distinct from finalists")
    return finalists


def frontier_run_case(args, output_dir, topology, atoms, core_atoms, phase,
              warmup, repetition, execution_index):
    tag = "warmup" if warmup else f"rep{repetition}"
    case = f"{topology}-{atoms}-core{core_atoms}-{tag}"
    individual = output_dir / "individual-results"
    raw_path = individual / f"{case}.json"
    wrapper_path = individual / f"{case}-run.json"
    process_dir = individual / f"{case}-process"
    if wrapper_path.is_file():
        return read(wrapper_path)
    command = [str(args.build_dir / "bin" / "joint_fixed_neighbor_experiment"), screen_MODE,
               str(raw_path), topology, str(atoms), str(core_atoms)]
    print(f"Running {case} with {args.timeout:g}s cap", flush=True)
    process = monitored(command, process_dir, time.monotonic() + args.timeout,
                        rss_limit=args.rss_limit, seconds=args.timeout)
    result = read(raw_path) if raw_path.is_file() else None
    if result is None:
        progress_path = Path(str(raw_path) + ".progress.json")
        result = screen_partial(read(progress_path) if progress_path.is_file() else None, core_atoms)
    report = {
        "topology": topology,
        "atoms": atoms,
        "core_atoms": core_atoms,
        "phase": phase,
        "warmup": warmup,
        "repetition": repetition,
        "execution_index": execution_index,
        "policy": screen_POLICY,
        "measurement_scope": "fixed-neighbor-search-only",
        "command": command,
        "process": process,
        "status": screen_status(process),
        "completed": process.get("status") == "completed" and result is not None,
        "result": result,
    }
    write(wrapper_path, report)
    return report


def frontier_interleaved_order(core_sizes, repetition):
    offset = (repetition - 1) % len(core_sizes)
    return list(core_sizes[offset:]) + list(core_sizes[:offset])


def frontier_run_repeated_case(args, output_dir, topology, atoms, core_sizes, phase, reports):
    for core_atoms in core_sizes:
        reports.append(frontier_run_case(args, output_dir, topology, atoms, core_atoms, phase,
                                  True, 0, len(reports) + 1))
        write(output_dir / "runs.json", reports)
    for repetition in range(1, frontier_MEASURED_REPETITIONS + 1):
        for core_atoms in frontier_interleaved_order(tuple(core_sizes), repetition):
            reports.append(frontier_run_case(args, output_dir, topology, atoms, core_atoms, phase,
                                     False, repetition, len(reports) + 1))
            write(output_dir / "runs.json", reports)


def frontier_aggregate(rows, topology, atoms, core_atoms):
    selected = [row for row in rows
                if row.get("topology") == topology and row.get("atoms") == atoms
                and row.get("core_atoms") == core_atoms]
    warmups = [row for row in selected if row.get("warmup") is True]
    measured = [row for row in selected if row.get("warmup") is False]
    repetitions = sorted(row.get("repetition") for row in measured)
    expected = list(range(1, frontier_MEASURED_REPETITIONS + 1))
    complete = repetitions == expected and len(measured) == frontier_MEASURED_REPETITIONS
    correct = complete and all(screen_correct(row) for row in measured)
    return {
        "topology": topology,
        "atoms": atoms,
        "core_atoms": core_atoms,
        "warmup_runs": len(warmups),
        "warmup_completed": len(warmups) == frontier_WARMUP_REPETITIONS and
        all(row.get("status") == "completed" for row in warmups),
        "measured_runs": len(measured),
        "measured_repetitions": repetitions,
        "complete": complete,
        "all_correct": correct,
        "correct_repetitions": sum(1 for row in measured if screen_correct(row)),
        "statuses": sorted({row.get("status") for row in selected}),
        "median_search_seconds": frontier_median([row.get("search_seconds") for row in measured]),
        "median_total_seconds": frontier_median([row.get("total_seconds") for row in measured]),
        "median_peak_rss_mb": frontier_median([row.get("peak_rss_mb") for row in measured]),
        "maximum_peak_rss_mb": frontier_maximum([row.get("peak_rss_mb") for row in measured]),
        "median_sweeps": frontier_median([row.get("sweeps") for row in measured]),
        "maximum_sweeps": frontier_maximum([row.get("sweeps") for row in measured]),
        "median_block_solves": frontier_median([row.get("block_solves") for row in measured]),
        "median_local_work_seconds": frontier_median([row.get("local_work_seconds") for row in measured]),
        "maximum_local_rows": frontier_maximum([row.get("maximum_local_rows") for row in measured]),
        "maximum_local_columns": frontier_maximum([row.get("maximum_local_columns") for row in measured]),
        "search_seconds": [row.get("search_seconds") for row in measured],
        "peak_rss_mb": [row.get("peak_rss_mb") for row in measured],
        "sweeps": [row.get("sweeps") for row in measured],
    }


def frontier_compare(candidate, control):
    candidate_search = candidate.get("median_search_seconds")
    control_search = control.get("median_search_seconds")
    improvement = ((control_search - candidate_search) / control_search
                   if frontier_finite(candidate_search) and frontier_finite(control_search) and control_search > 0
                   else None)
    candidate_rss = candidate.get("median_peak_rss_mb")
    control_rss = control.get("median_peak_rss_mb")
    rss_ok = (frontier_finite(candidate_rss) and frontier_finite(control_rss) and
              candidate_rss <= frontier_RSS_MATERIAL_REGRESSION_FACTOR * control_rss)
    candidate_total = candidate.get("median_total_seconds")
    control_total = control.get("median_total_seconds")
    total_ok = (frontier_finite(candidate_total) and frontier_finite(control_total) and
                candidate_total <= frontier_TOTAL_MATERIAL_REGRESSION_FACTOR * control_total)
    candidate_rows = candidate.get("maximum_local_rows")
    control_rows = control.get("maximum_local_rows")
    candidate_columns = candidate.get("maximum_local_columns")
    control_columns = control.get("maximum_local_columns")
    local_geometry_ok = True
    if frontier_finite(candidate_rows) and frontier_finite(control_rows):
        local_geometry_ok = local_geometry_ok and candidate_rows <= control_rows
    if frontier_finite(candidate_columns) and frontier_finite(control_columns):
        local_geometry_ok = local_geometry_ok and candidate_columns <= control_columns
    candidate_sweeps = candidate.get("median_sweeps")
    control_sweeps = control.get("median_sweeps")
    return {
        "candidate_core_atoms": candidate.get("core_atoms"),
        "control_core_atoms": control.get("core_atoms"),
        "candidate_all_correct": candidate.get("all_correct") is True,
        "control_all_correct": control.get("all_correct") is True,
        "search_improvement_fraction": improvement,
        "candidate_median_search_seconds": candidate_search,
        "control_median_search_seconds": control_search,
        "candidate_median_total_seconds": candidate_total,
        "control_median_total_seconds": control_total,
        "candidate_median_peak_rss_mb": candidate_rss,
        "control_median_peak_rss_mb": control_rss,
        "rss_no_material_regression": rss_ok,
        "total_no_material_regression": total_ok,
        "local_geometry_no_regression": local_geometry_ok,
        "no_material_regression": rss_ok and total_ok and local_geometry_ok,
        "candidate_median_sweeps": candidate_sweeps,
        "control_median_sweeps": control_sweeps,
        "sweeps_warning": (frontier_finite(candidate_sweeps) and candidate_sweeps > frontier_MAX_SWEEPS_WARNING),
        "maximum_candidate_sweeps_warning": (frontier_finite(candidate.get("maximum_sweeps")) and
                                               candidate["maximum_sweeps"] > frontier_MAX_SWEEPS_WARNING),
    }


def frontier_phase_report(rows, cases, core_sizes):
    report_cases = []
    aggregate_by_core = {}
    for topology, atoms in cases:
        cores = [frontier_aggregate(rows, topology, atoms, core) for core in core_sizes]
        control = next((row for row in cores if row["core_atoms"] == frontier_CONTROL_CORE_SIZE), None)
        comparisons = {}
        if control is not None:
            for row in cores:
                if row["core_atoms"] != frontier_CONTROL_CORE_SIZE:
                    comparisons[str(row["core_atoms"])] = frontier_compare(row, control)
        report_cases.append({"topology": topology, "atoms": atoms,
                             "cores": cores, "comparisons_vs_control": comparisons})
        for row in cores:
            aggregate_by_core.setdefault(str(row["core_atoms"]), []).append(row)
    aggregate = {}
    for core, entries in aggregate_by_core.items():
        values = [row.get("median_search_seconds") for row in entries]
        aggregate[core] = sum(values) if len(values) == len(cases) and all(frontier_finite(value) for value in values) else None
    return {"cases": report_cases,
            "aggregate_median_search_seconds_by_core": aggregate,
            "core_sizes": list(core_sizes)}


def frontier_phase_core(report, topology, atoms, core_atoms):
    case = next((item for item in report.get("cases", [])
                 if item.get("topology") == topology and item.get("atoms") == atoms), None)
    if case is None:
        return None
    return next((row for row in case.get("cores", [])
                 if row.get("core_atoms") == core_atoms), None)


def frontier_phase_comparison(report, topology, atoms, core_atoms):
    case = next((item for item in report.get("cases", [])
                 if item.get("topology") == topology and item.get("atoms") == atoms), None)
    return case.get("comparisons_vs_control", {}).get(str(core_atoms)) if case else None


def frontier_select_512(report, finalists):
    aggregate = report.get("aggregate_median_search_seconds_by_core", {})
    eligible = []
    for core in finalists:
        rows = [frontier_phase_core(report, topology, 512, core) for topology in ("chain", "cube")]
        if all(row is not None and row.get("all_correct") is True for row in rows):
            value = aggregate.get(str(core))
            if frontier_finite(value):
                eligible.append((value, core))
    return min(eligible)[1] if eligible else None


def frontier_analyze(rows, finalist_core_sizes=frontier_DEFAULT_FINALIST_CORE_SIZES,
            control_core_size=frontier_CONTROL_CORE_SIZE, selected_core_size=None):
    finalists = frontier_validate_finalists(finalist_core_sizes, control_core_size)
    coarse_cases = (("chain", 512), ("cube", 512))
    coarse_report = frontier_phase_report(rows, coarse_cases, (*finalists, control_core_size))
    selected = selected_core_size or frontier_select_512(coarse_report, finalists)
    large_cores = ((selected, control_core_size) if selected is not None else (control_core_size,))
    large_cases = (("chain", 1024), ("cube", 1024))
    large_report = frontier_phase_report(rows, large_cases, large_cores)
    comparisons = []
    for topology, atoms in (*coarse_cases, *large_cases):
        if selected is None:
            continue
        comparison = frontier_phase_comparison(
            coarse_report if atoms == 512 else large_report, topology, atoms, selected)
        if comparison is not None:
            comparison = dict(comparison)
            comparison.update({"topology": topology, "atoms": atoms})
            comparisons.append(comparison)
    gate_reasons = []
    if selected is None:
        gate_reasons.append("no finalist is correct on both 512-atom topologies")
    for comparison in comparisons:
        label = f"{comparison['topology']}-{comparison['atoms']}"
        if comparison.get("candidate_all_correct") is not True:
            gate_reasons.append(f"{label}: candidate measured correctness gate failed")
        if comparison.get("control_all_correct") is not True:
            gate_reasons.append(f"{label}: control measured correctness gate failed")
        improvement = comparison.get("search_improvement_fraction")
        if not frontier_finite(improvement) or improvement < frontier_MIN_SEARCH_IMPROVEMENT:
            gate_reasons.append(f"{label}: median search improvement below 10%")
        if comparison.get("no_material_regression") is not True:
            gate_reasons.append(f"{label}: material resource or total-time regression")
    expected_comparisons = 4 if selected is not None else 0
    if selected is not None and len(comparisons) != expected_comparisons:
        gate_reasons.append("frontier is missing one or more required topology/size comparisons")
    warnings = []
    for comparison in comparisons:
        if comparison.get("sweeps_warning") or comparison.get("maximum_candidate_sweeps_warning"):
            warnings.append(f"{comparison['topology']}-{comparison['atoms']}: candidate sweeps exceed 24")
    return {
        "phase": "FixedNeighbor outer-core repeated frontier",
        "policy": screen_POLICY,
        "finalist_core_sizes": list(finalists),
        "control_core_size": control_core_size,
        "warmup_repetitions": frontier_WARMUP_REPETITIONS,
        "measured_repetitions": frontier_MEASURED_REPETITIONS,
        "interleaving": "cyclic rotation of each core list for each measured repetition",
        "minimum_search_improvement_fraction": frontier_MIN_SEARCH_IMPROVEMENT,
        "rss_material_regression_factor": frontier_RSS_MATERIAL_REGRESSION_FACTOR,
        "total_material_regression_factor": frontier_TOTAL_MATERIAL_REGRESSION_FACTOR,
        "maximum_sweeps_warning": frontier_MAX_SWEEPS_WARNING,
        "coarse_512": coarse_report,
        "large_1024": large_report,
        "selected_core_size": selected,
        "comparisons": comparisons,
        "warnings": warnings,
        "qualification_gate": "passed" if not gate_reasons else "failed",
        "gate_reasons": gate_reasons,
        "selection_basis": "fastest aggregate median search time among finalists correct on both 512-atom topologies",
    }


def frontier_all_rows(reports):
    return [screen_summarize(report) | {
        "phase": report.get("phase"),
        "warmup": report.get("warmup"),
        "repetition": report.get("repetition"),
        "execution_index": report.get("execution_index"),
    } for report in reports if report.get("result") is not None]


def frontier_write_outputs(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = [
        "topology", "atoms", "core_atoms", "warmup_runs", "warmup_completed",
        "measured_runs", "measured_repetitions", "complete", "all_correct",
        "correct_repetitions", "statuses", "median_search_seconds", "median_total_seconds",
        "median_peak_rss_mb", "maximum_peak_rss_mb", "median_sweeps", "maximum_sweeps",
        "median_block_solves", "median_local_work_seconds", "maximum_local_rows",
        "maximum_local_columns", "search_seconds", "peak_rss_mb", "sweeps",
    ]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for phase_name in ("coarse_512", "large_1024"):
            for case in report[phase_name]["cases"]:
                for row in case["cores"]:
                    writer.writerow({column: row.get(column) for column in columns})


def frontier_run_campaign(args):
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = frontier_case_list(args.cases)
    finalists = frontier_validate_finalists(args.finalist_core_sizes, args.control_core_size)
    required = {("chain", 512), ("cube", 512), ("chain", 1024), ("cube", 1024)}
    if set(cases) != required:
        raise ValueError("frontier requires exactly chain/cube 512 and chain/cube 1024 cases")
    reports = []
    coarse_cores = (*finalists, args.control_core_size)
    for topology, atoms in (("chain", 512), ("cube", 512)):
        frontier_run_repeated_case(args, output_dir, topology, atoms, coarse_cores,
                           "coarse-512", reports)
    coarse_rows = frontier_all_rows(reports)
    coarse_analysis = frontier_analyze(coarse_rows, finalists, args.control_core_size)
    selected = coarse_analysis["selected_core_size"]
    if selected is not None:
        large_cores = (selected, args.control_core_size)
        for topology, atoms in (("chain", 1024), ("cube", 1024)):
            frontier_run_repeated_case(args, output_dir, topology, atoms, large_cores,
                               "large-1024", reports)
    write(output_dir / "runs.json", reports)
    analysis = frontier_analyze(frontier_all_rows(reports), finalists, args.control_core_size, selected)
    manifest = {
        "schema_version": 1,
        "phase": analysis["phase"],
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {
            "driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
            "support": sha(Path(__file__)),
        },
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "backend": "EIGEN",
        "eigen_threads": 1,
        "cases": [f"{topology}-{atoms}" for topology, atoms in cases],
        "policy": screen_POLICY,
        "finalist_core_sizes": list(finalists),
        "control_core_size": args.control_core_size,
        "measurement_scope": "fixed-neighbor-search-only",
        "resource_envelope": {"wall_seconds": args.timeout, "rss_bytes": args.rss_limit},
        "repetitions": {"warmup": frontier_WARMUP_REPETITIONS, "measured": frontier_MEASURED_REPETITIONS},
        "correctness_contract": ["search_converged", "A/C KKT <= 1e-10",
                                 "width gradient <= 1e-12", "cache/objective replay limits"],
        "analysis": {"qualification_gate": analysis["qualification_gate"],
                     "selected_core_size": analysis["selected_core_size"],
                     "warnings": analysis["warnings"]},
    }
    write(output_dir / "campaign-manifest.json", manifest)
    frontier_write_outputs(analysis, output_dir)
    (output_dir / "README.md").write_text(
        "# FixedNeighbor outer-core repeated frontier\n\n"
        "This P4 frontier keeps LegacyCompact, FixedNeighbor's intrinsic one-accepted "
        "local-update contract, Forward serial "
        "Gauss-Seidel, EIGEN, one Eigen thread, and the existing search gates. It uses "
        f"one warmup and {frontier_MEASURED_REPETITIONS} interleaved measured repetitions. The "
        f"512-atom cases compare finalists {', '.join(str(core) for core in finalists)} "
        f"against control {args.control_core_size}; the fastest correct finalist by "
        "aggregate median search time is then compared with control on both 1024-atom "
        "cases.\n\n"
        "Gate C4 requires at least 10% median search improvement for the selected core "
        "on chain and cube at both sizes, no material RSS/total-time/local-geometry "
        f"regression, and correctness on every measured run. Sweeps above {frontier_MAX_SWEEPS_WARNING} "
        "are reported as warnings.\n\n"
        f"Qualification gate: **{analysis['qualification_gate']}**. Selected core: "
        f"**{analysis['selected_core_size'] or 'none'}**.\n"
    )
    return analysis
