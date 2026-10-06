"""Run the bounded P2 FixedNeighbor local-work attribution campaign."""
from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

from experiment_io import ROOT, read, sha, write
from experiment_process import RSS_LIMIT_BYTES, monitored
from experiment_provenance import source_hash
from joint_fixed_neighbor_attribution import analyze, write_outputs


DEFAULT_CASES = ("chain-256", "chain-512", "cube-256", "cube-512")


def _case(value):
    topology, atoms = value.split("-", 1)
    atoms = int(atoms)
    if topology not in ("chain", "cube") or atoms not in (256, 512):
        raise ValueError(f"Unsupported attribution case: {value}")
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


def _partial_result(progress):
    if not progress:
        return None
    if "fixed_neighbor" in progress:
        return progress
    return dict(progress, fixed_neighbor={
        "search_converged": False,
        "search_reason": "interrupted-before-search-completion",
        "sweeps": len(progress.get("sweep_telemetry", [])),
        "sweep_telemetry": progress.get("sweep_telemetry", []),
        "block_telemetry": progress.get("local_block_telemetry", []),
        "profile_factor_seconds": sum(
            block.get("profile_factor_seconds", 0.0)
            for block in progress.get("local_block_telemetry", [])),
        "local_trajectory_telemetry": True,
    })


def _run_case(args, output_dir, topology, atoms):
    case = f"{topology}-{atoms}"
    individual = output_dir / "individual-results"
    raw_path = individual / f"{case}.json"
    wrapper_path = individual / f"{case}-run.json"
    process_dir = individual / f"{case}-process"
    if wrapper_path.is_file():
        return read(wrapper_path)
    command = [str(args.build_dir / "bin" / "joint_fixed_neighbor_experiment"),
               "--attribution", str(raw_path), topology, str(atoms)]
    print(f"Running {case} local-work attribution with {args.timeout:g}s cap", flush=True)
    process = monitored(command, process_dir, time.monotonic() + args.timeout,
                        rss_limit=args.rss_limit, seconds=args.timeout)
    result = read(raw_path) if raw_path.is_file() else None
    if result is None:
        progress_path = Path(str(raw_path) + ".progress.json")
        result = _partial_result(read(progress_path) if progress_path.is_file() else None)
    report = {
        "topology": topology, "atoms": atoms,
        "measurement_scope": "fixed-neighbor-search-only",
        "command": command, "process": process,
        "status": _status(process),
        "completed": process.get("status") == "completed" and result is not None,
        "result": result,
    }
    write(wrapper_path, report)
    return report


def run_campaign(args):
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = [_case(value) for value in args.cases]
    manifest = {
        "schema_version": 1,
        "phase": "P2 FixedNeighbor local-work attribution",
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {
            "attribution_driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
            "attribution_analyzer": sha(Path(__file__).with_name("joint_fixed_neighbor_attribution.py")),
            "campaign_runner": sha(Path(__file__)),
        },
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "backend": "SPQR", "eigen_threads": 1,
        "cases": [f"{topology}-{atoms}" for topology, atoms in cases],
        "measurement_scope": "fixed-neighbor-search-only",
        "resource_envelope": {"wall_seconds": args.timeout, "rss_bytes": args.rss_limit},
        "local_trajectory_fields": [
            "profile_evaluation", "accepted", "accepted_update", "local_objective_before",
            "local_objective_after", "objective_reduction", "eta_change_inf",
            "gradient_inf_norm", "factor_seconds", "cumulative_factor_seconds",
        ],
        "route_policy": "serial forward FixedNeighbor, 128 atom core, full local LegacyCompact search",
    }
    reports = [_run_case(args, output_dir, topology, atoms) for topology, atoms in cases]
    results = [report["result"] for report in reports if report.get("result") is not None]
    decision = args.decision
    analysis = analyze(results, decision)
    manifest["analysis"] = {"over_solving_decision": analysis["over_solving_decision"]}
    write(output_dir / "campaign-manifest.json", manifest)
    write_outputs(analysis, output_dir)
    write(output_dir / "runs.json", reports)
    (output_dir / "README.md").write_text(
        "# FixedNeighbor local-work attribution\n\n"
        "This P2 campaign records every local LegacyCompact profile trial for chain/cube "
        "256 and 512 atom workloads. It uses separate processes, one Eigen thread, a "
        f"{args.rss_limit / 1024**3:g} GiB RSS limit, and a {args.timeout:g} second diagnostic cap. "
        "The campaign is search-only; it does not alter block acceptance or endpoint thresholds.\n\n"
        f"Over-solving decision recorded by the analyzer: **{analysis['over_solving_decision']}**.\n"
    )
    return analysis


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=7200.0)
    parser.add_argument("--rss-limit", type=int, default=RSS_LIMIT_BYTES)
    parser.add_argument("--decision", choices=("yes", "no", "insufficient evidence"))
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
        parser.error("build joint_fixed_neighbor_experiment before running attribution")
    for value in args.cases:
        try:
            _case(value)
        except (ValueError, TypeError) as error:
            parser.error(str(error))
    try:
        analysis = run_campaign(args)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.error(str(error))
    print(json.dumps({"over_solving_decision": analysis["over_solving_decision"],
                      "output_dir": str(args.output_dir)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
