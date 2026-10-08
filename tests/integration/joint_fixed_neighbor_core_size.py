"""Compare selected bounded local-work performance across core sizes."""
# Internal screen-phase implementation. Use joint_fixed_neighbor_outer_core_qualification.py.
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import time
from pathlib import Path

from experiment_io import ROOT, read, sha, write
from experiment_process import RSS_LIMIT_BYTES, monitored
from experiment_provenance import source_hash


DEFAULT_CASES = ("chain-512", "cube-512", "cube-1024")
HISTORICAL_CORE_SIZES = (64, 128, 256)
# Keep the historical artifact reproducible while allowing every campaign to
# choose its own explicit outer-core list.
CORE_SIZES = HISTORICAL_CORE_SIZES
POLICY = "OneAccepted"
MODE = "--inexact-one-search"


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _parse_core_sizes(value):
    if isinstance(value, (tuple, list)):
        values = tuple(int(item) for item in value)
    else:
        values = tuple(int(item.strip()) for item in str(value).split(",") if item.strip())
    if not values or any(item <= 0 for item in values) or len(set(values)) != len(values):
        raise ValueError("core sizes must be a non-empty comma-separated list of unique positive integers")
    return values


def _case(value):
    topology, atoms = value.split("-", 1)
    atoms = int(atoms)
    if topology not in ("chain", "cube") or atoms not in (512, 1024):
        raise ValueError(f"Unsupported core-size case: {value}")
    return topology, atoms


def _status(process):
    if process.get("status") == "completed":
        return "completed"
    if process.get("status") == "time-limit":
        return "timeout"
    if process.get("status") in ("rss-limit", "rss-limit-observed-after-exit"):
        return "rss_limit"
    return "process_error"


def _partial(progress, core_atoms):
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


def _run_case(args, output_dir, topology, atoms, core_atoms):
    case = f"{topology}-{atoms}-core{core_atoms}"
    individual = output_dir / "individual-results"
    raw_path = individual / f"{case}.json"
    wrapper_path = individual / f"{case}-run.json"
    process_dir = individual / f"{case}-process"
    if wrapper_path.is_file():
        return read(wrapper_path)
    command = [str(args.build_dir / "bin" / "joint_fixed_neighbor_experiment"), MODE,
               str(raw_path), topology, str(atoms), str(core_atoms)]
    print(f"Running {case} with {args.timeout:g}s cap", flush=True)
    process = monitored(command, process_dir, time.monotonic() + args.timeout,
                        rss_limit=args.rss_limit, seconds=args.timeout)
    result = read(raw_path) if raw_path.is_file() else None
    if result is None:
        progress_path = Path(str(raw_path) + ".progress.json")
        result = _partial(read(progress_path) if progress_path.is_file() else None, core_atoms)
    report = {"topology": topology, "atoms": atoms, "core_atoms": core_atoms,
              "policy": POLICY, "measurement_scope": "fixed-neighbor-search-only",
              "command": command, "process": process, "status": _status(process),
              "completed": process.get("status") == "completed" and result is not None,
              "result": result}
    write(wrapper_path, report)
    return report


def _profile_work(fixed):
    profile = fixed.get("local_profile_work")
    if isinstance(profile, dict):
        return profile
    return (fixed.get("fixed_neighbor_work") or {}).get("local_profile_work", {})


def _block_geometry(result, fixed, fallback_core):
    records = fixed.get("block_telemetry") or []
    sizes = {}
    for record in records:
        block = record.get("block")
        size = record.get("core_atoms", record.get("atoms"))
        if block is not None and _finite(size):
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


def _maximum_sweep_value(fixed, field):
    values = [sweep.get(field) for sweep in fixed.get("sweep_telemetry", [])]
    values = [value for value in values if _finite(value)]
    return max(values) if values else None


def _work_value(work, fixed, field, fallback=None):
    value = work.get(field)
    if _finite(value):
        return value
    value = fixed.get(field, fallback)
    return value if _finite(value) else fallback


def _ratio(numerator, denominator):
    return numerator / denominator if _finite(numerator) and _finite(denominator) and denominator > 0 else None


def summarize(report):
    result = report.get("result") or {}
    fixed = result.get("fixed_neighbor", result)
    process = report.get("process") or {}
    work = fixed.get("fixed_neighbor_work") or {}
    profile = _profile_work(fixed)
    profile_total = profile.get("total") or {}
    search_seconds = fixed.get("search_seconds")
    sweeps = fixed.get("sweeps", len(fixed.get("sweep_telemetry", [])))
    block_solves = fixed.get("block_solves", sum(s.get("block_solves", 0)
                                                  for s in fixed.get("sweep_telemetry", [])))
    profile_evaluations = fixed.get("profile_evaluations", sum(
        s.get("profile_evaluations", 0) for s in fixed.get("sweep_telemetry", [])))
    local_search_seconds = _work_value(work, fixed, "local_search_seconds",
                                       fixed.get("local_factor_seconds"))
    candidate_replay_seconds = _work_value(work, fixed, "candidate_replay_seconds")
    sweep_replay_seconds = _work_value(work, fixed, "sweep_replay_seconds")
    sweep_global_state_seconds = _work_value(work, fixed, "sweep_global_state_seconds")
    objective = fixed.get("objective")
    maximum_cache_replay_error = _maximum_sweep_value(fixed, "cache_replay_error")
    maximum_objective_replay_error = _maximum_sweep_value(fixed, "objective_replay_error")
    cache_replay_limit = 2e-12
    objective_replay_limit = (1e-12 + 2e-12 * abs(objective)
                               if _finite(objective) else None)
    row = {
        "topology": report.get("topology"), "atoms": report.get("atoms"),
        "core_atoms": report.get("core_atoms", result.get("core_atoms")),
        "policy": report.get("policy", POLICY), "status": report.get("status", "completed"),
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
        "old_core_seconds": _work_value(work, fixed, "old_core_seconds"),
        "effective_response_seconds": _work_value(work, fixed, "effective_response_seconds"),
        "local_state_seconds": _work_value(work, fixed, "local_state_seconds"),
        "candidate_copy_seconds": _work_value(work, fixed, "candidate_copy_seconds"),
        "candidate_replay_seconds": candidate_replay_seconds,
        "cache_update_seconds": _work_value(work, fixed, "cache_update_seconds"),
        "sweep_replay_seconds": sweep_replay_seconds,
        "sweep_global_state_seconds": sweep_global_state_seconds,
        "derivative_prepare_seconds": profile_total.get("derivative_prepare_seconds"),
        "derivative_reduce_seconds": profile_total.get("derivative_reduce_seconds"),
        "derivative_jacobian_qr_seconds": profile_total.get("derivative_jacobian_qr_seconds"),
        "tiled_qr_assembly_copy_seconds": profile_total.get("tiled_qr_assembly_copy_seconds"),
        "tiled_qr_householder_seconds": profile_total.get("tiled_qr_householder_seconds"),
        "tiled_qr_rhs_transform_seconds": profile_total.get("tiled_qr_rhs_transform_seconds"),
    }
    row.update(_block_geometry(result, fixed, row["core_atoms"]))
    row.update({
        "search_seconds_per_sweep": _ratio(search_seconds, sweeps),
        "local_search_seconds_per_block_solve": _ratio(local_search_seconds, block_solves),
        "candidate_replay_seconds_per_block_solve": _ratio(candidate_replay_seconds, block_solves),
        "profile_evaluations_per_block_solve": _ratio(profile_evaluations, block_solves),
        "accepted_updates_per_sweep": _ratio(fixed.get("accepted_local_updates"), sweeps),
        "local_search_fraction": _ratio(local_search_seconds, search_seconds),
        "candidate_replay_fraction": _ratio(candidate_replay_seconds, search_seconds),
        "global_overhead_fraction": _ratio(
            (sweep_replay_seconds or 0.0) + (sweep_global_state_seconds or 0.0), search_seconds),
    })
    return row


def _correct(row):
    return row["status"] == "completed" and row["search_converged"] is True and \
        row.get("search_reason") == "block-stationary" and \
        _finite(row.get("final_global_ac_kkt")) and row["final_global_ac_kkt"] <= 1e-10 and \
        _finite(row.get("final_raw_width_gradient_inf_norm")) and \
        row["final_raw_width_gradient_inf_norm"] <= 1e-12 and \
        row.get("cache_replay_within_limit") is True and \
        row.get("objective_replay_within_limit") is True


def analyze(rows, core_sizes=CORE_SIZES):
    core_sizes = _parse_core_sizes(core_sizes)
    core_order = {core: index for index, core in enumerate(core_sizes)}
    groups = {}
    for row in rows:
        groups.setdefault((row["topology"], row["atoms"]), []).append(row)
    cases = []
    for key, entries in sorted(groups.items()):
        entries = sorted(entries, key=lambda row: core_order.get(row.get("core_atoms"), len(core_order)))
        correct = [row["core_atoms"] for row in entries if row["core_atoms"] in core_order and _correct(row)]
        observed = {row["core_atoms"] for row in entries}
        fastest = min((row for row in entries if row.get("core_atoms") in core_order and _correct(row) and
                       _finite(row.get("search_seconds"))),
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
                  if row.get("core_atoms") == core and _correct(row) and
                  _finite(row.get("search_seconds"))]
        aggregate[str(core)] = sum(values) if len(values) == len(cases) else None
    selected = min((int(core) for core, value in aggregate.items() if _finite(value) and value > 0),
                   key=lambda core: aggregate[str(core)], default=None) if complete else None
    return {"phase": "FixedNeighbor outer-core performance study",
            "policy": POLICY, "core_sizes": list(core_sizes), "cases": cases,
            "correctness_gate": "passed" if complete else "failed",
            "aggregate_search_seconds_by_core": aggregate,
            "selected_core_size": selected,
            "selection_basis": "minimum aggregate search seconds among cores correct on every study case; resource, geometry, and attribution metrics remain reported"}


def write_outputs(report, output_dir):
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


def run_campaign(args):
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = [_case(value) for value in args.cases]
    core_sizes = _parse_core_sizes(args.core_sizes)
    reports = []
    for topology, atoms in cases:
        for core_atoms in core_sizes:
            reports.append(_run_case(args, output_dir, topology, atoms, core_atoms))
            write(output_dir / "runs.json", reports)
    rows = [summarize(report) for report in reports if report.get("result") is not None]
    analysis = analyze(rows, core_sizes)
    manifest = {
        "schema_version": 1, "phase": analysis["phase"],
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {"driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
                                    "core_size": sha(Path(__file__))},
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "backend": "SPQR", "eigen_threads": 1, "cases": [f"{t}-{a}" for t, a in cases],
        "policy": POLICY, "core_sizes": list(core_sizes),
        "measurement_scope": "fixed-neighbor-search-only",
        "resource_envelope": {"wall_seconds": args.timeout, "rss_bytes": args.rss_limit},
        "correctness_contract": ["search_converged", "A/C KKT <= 1e-10",
                                 "width gradient <= 1e-12"],
        "analysis": {"correctness_gate": analysis["correctness_gate"],
                     "selected_core_size": analysis["selected_core_size"]},
    }
    write(output_dir / "campaign-manifest.json", manifest)
    write_outputs(analysis, output_dir)
    (output_dir / "README.md").write_text(
        "# FixedNeighbor outer-core performance study\n\n"
        "This search-only study keeps the production-compatible OneAcceptedLocalUpdate "
        "policy, forward serial Gauss-Seidel order, SPQR backend, one Eigen thread, and "
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


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=7200.0)
    parser.add_argument("--rss-limit", type=int, default=RSS_LIMIT_BYTES)
    parser.add_argument("--cases", nargs="+", default=DEFAULT_CASES)
    parser.add_argument("--core-sizes", type=_parse_core_sizes,
                        default=HISTORICAL_CORE_SIZES,
                        help="comma-separated outer core sizes (default: 64,128,256)")
    return parser


def main(argv=None):
    parser = build_parser(); args = parser.parse_args(argv)
    args.build_dir = args.build_dir.resolve(); args.output_dir = args.output_dir.resolve()
    args.core_sizes = _parse_core_sizes(args.core_sizes)
    if args.timeout <= 0 or args.rss_limit <= 0:
        parser.error("resource limits must be positive")
    if not (args.build_dir / "bin" / "joint_fixed_neighbor_experiment").is_file():
        parser.error("build joint_fixed_neighbor_experiment before running core-size study")
    for value in args.cases:
        try:
            _case(value)
        except (ValueError, TypeError) as error:
            parser.error(str(error))
    try:
        analysis = run_campaign(args)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.error(str(error))
    print(json.dumps({"correctness_gate": analysis["correctness_gate"],
                      "selected_core_size": analysis["selected_core_size"],
                      "output_dir": str(args.output_dir)}))
    return 0 if analysis["correctness_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
