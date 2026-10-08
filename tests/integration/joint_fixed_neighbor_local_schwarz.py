"""Screen local FixedNeighbor OperatorPcg Schwarz geometries on 256-atom cases."""
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


DEFAULT_CASES = ("chain-256", "cube-256")
OUTER_CORE_ATOMS = 64
LOCAL_MAX_BLOCK_ATOMS = 64
CANDIDATES = (
    {"name": "Identity", "mode": "--local-operator-identity-search", "preconditioner": "Identity",
     "core": None, "overlap": None, "max_block": None},
    {"name": "Diagonal", "mode": "--local-operator-diagonal-search", "preconditioner": "Diagonal",
     "core": None, "overlap": None, "max_block": None},
    {"name": "Schwarz-128/1", "mode": "--local-operator-schwarz-search", "preconditioner": "Schwarz",
     "core": 128, "overlap": 1, "max_block": 512},
    {"name": "Schwarz-64/0", "mode": "--local-operator-schwarz-search", "preconditioner": "Schwarz",
     "core": 64, "overlap": 0, "max_block": 64},
    {"name": "Schwarz-64/1", "mode": "--local-operator-schwarz-search", "preconditioner": "Schwarz",
     "core": 64, "overlap": 1, "max_block": 64},
    {"name": "Schwarz-32/0", "mode": "--local-operator-schwarz-search", "preconditioner": "Schwarz",
     "core": 32, "overlap": 0, "max_block": 64},
    {"name": "Schwarz-32/1", "mode": "--local-operator-schwarz-search", "preconditioner": "Schwarz",
     "core": 32, "overlap": 1, "max_block": 64},
    {"name": "Schwarz-16/0", "mode": "--local-operator-schwarz-search", "preconditioner": "Schwarz",
     "core": 16, "overlap": 0, "max_block": 64},
    {"name": "Schwarz-16/1", "mode": "--local-operator-schwarz-search", "preconditioner": "Schwarz",
     "core": 16, "overlap": 1, "max_block": 64},
    {"name": "Schwarz-8/0", "mode": "--local-operator-schwarz-search", "preconditioner": "Schwarz",
     "core": 8, "overlap": 0, "max_block": 64},
    {"name": "Schwarz-8/1", "mode": "--local-operator-schwarz-search", "preconditioner": "Schwarz",
     "core": 8, "overlap": 1, "max_block": 64},
)


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _case(value):
    topology, atoms = value.split("-", 1)
    atoms = int(atoms)
    if topology not in ("chain", "cube") or atoms != 256:
        raise ValueError(f"Unsupported local Schwarz case: {value}")
    return topology, atoms


def _status(process):
    if process.get("status") == "completed":
        return "completed"
    if process.get("status") == "time-limit":
        return "timeout"
    if process.get("status") in ("rss-limit", "rss-limit-observed-after-exit"):
        return "rss_limit"
    if process.get("status") == "not-run-budget":
        return "not_run"
    return "process_error"


def _partial(progress, topology, atoms, candidate):
    if not progress:
        return None
    fixed = dict(progress)
    fixed.setdefault("search_converged", False)
    fixed.setdefault("search_reason", "interrupted-before-search-completion")
    fixed.setdefault("sweeps", len(progress.get("sweep_telemetry", [])))
    fixed.setdefault("confirmed_stationarity_sweep", None)
    fixed.setdefault("search_seconds", sum(row.get("wall_seconds", 0.0)
                                            for row in progress.get("sweep_telemetry", [])))
    fixed.setdefault("total_elapsed_seconds", fixed["search_seconds"])
    fixed.setdefault("final_global_ac_kkt", None)
    fixed.setdefault("final_raw_width_gradient_inf_norm", None)
    fixed.setdefault("local_operator_work", {})
    return {"topology": topology, "atoms": atoms, "outer_core_atoms": OUTER_CORE_ATOMS,
            "candidate": candidate["name"], "measurement_scope": "fixed-neighbor-search-only",
            "fixed_neighbor": fixed, "peak_rss_mb": None}


def _run_case(args, output_dir, topology, atoms, candidate):
    slug = candidate["name"].lower().replace("/", "-")
    case = f"{topology}-{atoms}-{slug}"
    individual = output_dir / "individual-results"
    raw_path = individual / f"{case}.json"
    wrapper_path = individual / f"{case}-run.json"
    process_dir = individual / f"{case}-process"
    if wrapper_path.is_file():
        return read(wrapper_path)
    command = [str(args.build_dir / "bin" / "joint_fixed_neighbor_experiment"), candidate["mode"],
               str(raw_path), topology, str(atoms), str(OUTER_CORE_ATOMS),
               str(candidate["core"] if candidate["core"] is not None else OUTER_CORE_ATOMS),
               str(candidate["overlap"] if candidate["overlap"] is not None else 0),
               str(candidate["max_block"] if candidate["max_block"] is not None else LOCAL_MAX_BLOCK_ATOMS)]
    print(f"Running {case} with {args.timeout:g}s cap", flush=True)
    process = monitored(command, process_dir, time.monotonic() + args.timeout,
                        rss_limit=args.rss_limit, seconds=args.timeout)
    result = read(raw_path) if raw_path.is_file() else None
    if result is None:
        progress_path = Path(str(raw_path) + ".progress.json")
        result = _partial(read(progress_path) if progress_path.is_file() else None,
                          topology, atoms, candidate)
    report = {"topology": topology, "atoms": atoms, "outer_core_atoms": OUTER_CORE_ATOMS,
              "candidate": candidate, "measurement_scope": "fixed-neighbor-search-only",
              "command": command, "process": process, "status": _status(process),
              "completed": process.get("status") == "completed" and result is not None,
              "result": result}
    write(wrapper_path, report)
    return report


def _max_replay(fixed, key):
    values = [row.get(key) for row in fixed.get("sweep_telemetry", []) if _finite(row.get(key))]
    return max(values, default=None)


def summarize(report):
    result = report.get("result") or {}
    fixed = result.get("fixed_neighbor", result)
    work = fixed.get("local_operator_work", {})
    candidate = report.get("candidate", {})
    local_atoms = fixed.get("local_atoms", {})
    return {
        "topology": report.get("topology"), "atoms": report.get("atoms"),
        "candidate": candidate.get("name"), "preconditioner": candidate.get("preconditioner"),
        "outer_core_atoms": fixed.get("outer_core_atoms", report.get("outer_core_atoms")),
        "local_atoms_min": local_atoms.get("minimum"), "local_atoms_mean": local_atoms.get("mean"),
        "local_atoms_max": local_atoms.get("maximum"),
        "requested_schwarz_core_atoms": fixed.get("local_schwarz_core_atoms"),
        "requested_schwarz_overlap_hops": fixed.get("local_schwarz_overlap_hops"),
        "requested_schwarz_max_block_atoms": fixed.get("local_schwarz_max_block_atoms"),
        "preconditioner_partition_count": work.get("preconditioner_partition_count"),
        "minimum_core_atoms": work.get("minimum_core_atoms"),
        "mean_core_atoms": work.get("mean_core_atoms"),
        "maximum_core_atoms": work.get("maximum_core_atoms"),
        "minimum_overlap_atoms": work.get("minimum_overlap_atoms"),
        "mean_overlap_atoms": work.get("mean_overlap_atoms"),
        "maximum_overlap_atoms": work.get("maximum_overlap_atoms"),
        "minimum_realized_block_atoms": work.get("minimum_realized_block_atoms"),
        "mean_realized_block_atoms": work.get("mean_realized_block_atoms"),
        "maximum_realized_block_atoms": work.get("maximum_realized_block_atoms"),
        "mean_preconditioner_coverage_ratio": work.get("mean_preconditioner_coverage_ratio"),
        "maximum_preconditioner_coverage_ratio": work.get("maximum_preconditioner_coverage_ratio"),
        "status": report.get("status", "completed"),
        "search_converged": fixed.get("search_converged"), "search_reason": fixed.get("search_reason"),
        "sweeps": fixed.get("sweeps"), "block_solves": fixed.get("block_solves"),
        "profile_evaluations": fixed.get("profile_evaluations"),
        "search_seconds": fixed.get("search_seconds"),
        "operator_setup_seconds": work.get("operator_setup_seconds"),
        "preconditioner_setup_seconds": work.get("preconditioner_setup_seconds"),
        "partition_seconds": work.get("partition_seconds"),
        "preconditioner_model_seconds": work.get("preconditioner_model_seconds"),
        "preconditioner_factor_seconds": work.get("preconditioner_factor_seconds"),
        "inverse_seconds": work.get("inverse_seconds"), "pcg_seconds": work.get("pcg_seconds"),
        "pcg_solves": work.get("pcg_solves"), "pcg_iterations": work.get("pcg_iterations"),
        "pcg_mean_iterations": work.get("pcg_mean_iterations"),
        "pcg_max_iterations": work.get("pcg_max_iterations"),
        "damping_trials": work.get("damping_trials"),
        "peak_rss_mb": result.get("peak_rss_mb", (report.get("process") or {}).get(
            "sampled_tree_peak_rss_bytes", 0) / 1024**2),
        "objective": fixed.get("objective"),
        "final_global_ac_kkt": fixed.get("final_global_ac_kkt"),
        "final_raw_width_gradient_inf_norm": fixed.get("final_raw_width_gradient_inf_norm"),
        "maximum_cache_replay_error": _max_replay(fixed, "cache_replay_error"),
        "maximum_objective_replay_error": _max_replay(fixed, "objective_replay_error"),
        "assessment_local": None, "runtime_convergence": "NotRun",
        "measurement_scope": report.get("measurement_scope", "fixed-neighbor-search-only"),
    }


def _correct(row):
    objective = row.get("objective")
    cache_limit = 2e-12
    objective_limit = 1e-12 + 2e-12 * abs(objective) if _finite(objective) else math.inf
    return row.get("status") == "completed" and row.get("search_converged") is True and \
        _finite(row.get("final_global_ac_kkt")) and row["final_global_ac_kkt"] <= 1e-10 and \
        _finite(row.get("final_raw_width_gradient_inf_norm")) and \
        row["final_raw_width_gradient_inf_norm"] <= 1e-12 and \
        (row.get("maximum_cache_replay_error") is None or
         row["maximum_cache_replay_error"] <= cache_limit) and \
        (row.get("maximum_objective_replay_error") is None or
         row["maximum_objective_replay_error"] <= objective_limit)


def analyze(rows):
    groups = {}
    for row in rows:
        groups.setdefault((row["topology"], row["atoms"]), []).append(row)
    cases = []
    for key, entries in sorted(groups.items()):
        entries = sorted(entries, key=lambda row: next(
            index for index, item in enumerate(CANDIDATES) if item["name"] == row["candidate"]))
        correct = [row["candidate"] for row in entries if _correct(row)]
        schwarz = [row for row in entries if row.get("preconditioner") == "Schwarz" and _correct(row)]
        cases.append({"topology": key[0], "atoms": key[1], "candidates": entries,
                      "correct_candidates": correct,
                      "correct_schwarz_candidates": [row["candidate"] for row in schwarz]})
    all_candidates = {row["candidate"] for case in cases for row in case["candidates"]}
    common_schwarz = [item["name"] for item in CANDIDATES if item["preconditioner"] == "Schwarz" and
                      cases and all(item["name"] in case["correct_schwarz_candidates"] for case in cases)]
    aggregate = {name: sum(row.get("search_seconds", math.inf) for case in cases for row in case["candidates"]
                           if row["candidate"] == name and _correct(row)) for name in common_schwarz}
    retained = sorted(common_schwarz, key=lambda name: aggregate[name])[:3]
    return {
        "phase": "P0-P2 FixedNeighbor local Schwarz geometry search screen",
        "cases": cases, "candidate_count": len(CANDIDATES), "observed_candidate_count": len(all_candidates),
        "correctness_gate": "passed" if cases and all(
            all(_correct(row) for row in case["candidates"]) for case in cases) else "failed",
        "common_correct_schwarz_candidates": common_schwarz,
        "aggregate_search_seconds_by_schwarz": aggregate,
        "retained_schwarz_candidates": retained,
        "gate_s_status": "pending-full-endpoint-qualification",
        "selection_basis": "search correctness first; search seconds, PCG work, setup, RSS, and coverage are reported; endpoint stationarity is not run in this screen",
    }


def write_outputs(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = ["topology", "atoms", "candidate", "preconditioner", "outer_core_atoms",
               "local_atoms_min", "local_atoms_mean", "local_atoms_max",
               "requested_schwarz_core_atoms", "requested_schwarz_overlap_hops",
               "requested_schwarz_max_block_atoms", "preconditioner_partition_count",
               "minimum_core_atoms", "mean_core_atoms", "maximum_core_atoms",
               "minimum_overlap_atoms", "mean_overlap_atoms", "maximum_overlap_atoms",
               "minimum_realized_block_atoms", "mean_realized_block_atoms", "maximum_realized_block_atoms",
               "mean_preconditioner_coverage_ratio", "maximum_preconditioner_coverage_ratio",
               "operator_setup_seconds", "preconditioner_setup_seconds", "partition_seconds",
               "preconditioner_model_seconds", "preconditioner_factor_seconds", "inverse_seconds",
               "pcg_seconds", "pcg_solves", "pcg_iterations", "pcg_mean_iterations", "pcg_max_iterations",
               "search_seconds", "peak_rss_mb", "assessment_local", "runtime_convergence",
               "search_converged", "search_reason", "maximum_cache_replay_error",
               "maximum_objective_replay_error"]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            for row in case["candidates"]:
                writer.writerow({column: row.get(column) for column in columns})


def run_campaign(args):
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = [_case(value) for value in args.cases]
    reports = []
    for topology, atoms in cases:
        for candidate in CANDIDATES:
            reports.append(_run_case(args, output_dir, topology, atoms, candidate))
            write(output_dir / "runs.json", reports)
    rows = [summarize(report) for report in reports if report.get("result") is not None]
    analysis = analyze(rows)
    manifest = {
        "schema_version": 1, "phase": analysis["phase"],
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {"driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
                                    "campaign_runner": sha(Path(__file__))},
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "backend": "SPQR", "eigen_threads": 1, "outer_core_atoms": OUTER_CORE_ATOMS,
        "cases": [f"{topology}-{atoms}" for topology, atoms in cases],
        "candidates": list(CANDIDATES), "measurement_scope": "fixed-neighbor-search-only",
        "endpoint_claim": "not-run", "resource_envelope": {"wall_seconds": args.timeout,
        "rss_bytes": args.rss_limit}, "correctness_contract": [
            "search_converged", "A/C KKT <= 1e-10", "width gradient <= 1e-12",
            "cache replay and objective replay pass"],
        "analysis": {"correctness_gate": analysis["correctness_gate"],
                     "retained_schwarz_candidates": analysis["retained_schwarz_candidates"]},
    }
    manifest["analysis"] = analysis
    write(output_dir / "campaign-manifest.json", manifest)
    write_outputs(analysis, output_dir)
    (output_dir / "README.md").write_text(
        "# FixedNeighbor local Schwarz geometry search screen\n\n"
        "This P2 screen fixes the outer FixedNeighbor core at 64 atoms and varies only the local "
        "OperatorPcg preconditioner. It covers Identity, Diagonal, historical Schwarz 128/1/512, "
        "and Schwarz 8/16/32/64 with overlap 0/1 on chain-256 and cube-256.\n\n"
        "The process is search-only: endpoint assessment, endpoint trust, endpoint certification, "
        "and RuntimeConvergence are explicitly not run. Realized block geometry and coverage are "
        "therefore the primary confounder check; Gate S remains pending full endpoint qualification.\n\n"
        f"Screen correctness gate: **{analysis['correctness_gate']}**. Retained search candidates: "
        f"**{', '.join(analysis['retained_schwarz_candidates']) or 'none'}**.\n"
    )
    return analysis


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=7200.0)
    parser.add_argument("--rss-limit", type=int, default=RSS_LIMIT_BYTES)
    parser.add_argument("--cases", nargs="+", default=DEFAULT_CASES)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.build_dir = args.build_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    if args.timeout <= 0 or args.rss_limit <= 0:
        parser.error("resource limits must be positive")
    if not (args.build_dir / "bin" / "joint_fixed_neighbor_experiment").is_file():
        parser.error("build joint_fixed_neighbor_experiment before running local Schwarz screen")
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
                      "retained_schwarz_candidates": analysis["retained_schwarz_candidates"],
                      "output_dir": str(args.output_dir)}))
    return 0 if analysis["correctness_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
