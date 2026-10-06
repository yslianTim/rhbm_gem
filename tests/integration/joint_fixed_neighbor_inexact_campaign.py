"""Run Full, OneAccepted, and TwoAccepted FixedNeighbor endpoint comparisons."""
from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

from experiment_io import ROOT, read, sha, write
from experiment_process import RSS_LIMIT_BYTES, monitored
from experiment_provenance import source_hash
from joint_fixed_neighbor_inexact import POLICIES, analyze, write_outputs, summarize


DEFAULT_CASES = ("chain-256", "chain-512", "cube-256", "cube-512")
MODES = {"Full": "--full-attribution", "OneAccepted": "--inexact-one", "TwoAccepted": "--inexact-two"}


def _case(value):
    topology, atoms = value.split("-", 1)
    atoms = int(atoms)
    if topology not in ("chain", "cube") or atoms not in (256, 512):
        raise ValueError(f"Unsupported inexact case: {value}")
    return topology, atoms


def _status(process):
    if process.get("status") == "completed":
        return "completed"
    if process.get("status") == "time-limit":
        return "timeout"
    if process.get("status") in ("rss-limit", "rss-limit-observed-after-exit"):
        return "rss_limit"
    return "process_error"


def _partial(progress, policy):
    if not progress:
        return None
    if "fixed_neighbor" in progress:
        return progress
    return dict(progress, fixed_neighbor={
        "search_converged": False,
        "search_reason": "interrupted-before-endpoint-completion",
        "sweeps": len(progress.get("sweep_telemetry", [])),
        "sweep_telemetry": progress.get("sweep_telemetry", []),
        "block_telemetry": progress.get("local_block_telemetry", []),
        "runtime_convergence": "NotRun",
        "endpoint_certified": False,
        "search_seconds": sum(row.get("wall_seconds", 0.0)
                               for row in progress.get("sweep_telemetry", [])),
    })


def _run_case(args, output_dir, topology, atoms, policy):
    case = f"{topology}-{atoms}"
    individual = output_dir / "individual-results"
    raw_path = individual / f"{case}-{policy}.json"
    wrapper_path = individual / f"{case}-{policy}-run.json"
    process_dir = individual / f"{case}-{policy}-process"
    if wrapper_path.is_file():
        return read(wrapper_path)
    command = [str(args.build_dir / "bin" / "joint_fixed_neighbor_experiment"), MODES[policy],
               str(raw_path), topology, str(atoms)]
    print(f"Running {case} {policy} with {args.timeout:g}s cap", flush=True)
    process = monitored(command, process_dir, time.monotonic() + args.timeout,
                        rss_limit=args.rss_limit, seconds=args.timeout)
    result = read(raw_path) if raw_path.is_file() else None
    if result is None:
        progress_path = Path(str(raw_path) + ".progress.json")
        result = _partial(read(progress_path) if progress_path.is_file() else None, policy)
    report = {"topology": topology, "atoms": atoms, "policy": policy,
              "measurement_scope": "fixed-neighbor-full-endpoint", "command": command,
              "process": process, "status": _status(process),
              "completed": process.get("status") == "completed" and result is not None,
              "result": result}
    write(wrapper_path, report)
    return report


def run_campaign(args):
    output_dir = args.output_dir.resolve(); output_dir.mkdir(parents=True, exist_ok=True)
    cases = [_case(value) for value in args.cases]
    manifest = {
        "schema_version": 1, "phase": "P3 bounded FixedNeighbor local-work comparison",
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {
            "inexact_driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
            "inexact_analyzer": sha(Path(__file__).with_name("joint_fixed_neighbor_inexact.py")),
            "campaign_runner": sha(Path(__file__)),
        },
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "backend": "SPQR", "eigen_threads": 1,
        "cases": [f"{topology}-{atoms}" for topology, atoms in cases],
        "policies": list(POLICIES), "measurement_scope": "fixed-neighbor-full-endpoint",
        "resource_envelope": {"wall_seconds": args.timeout, "rss_bytes": args.rss_limit},
        "policy_definitions": {
            "Full": "FullLocalSearch",
            "OneAccepted": "at most one trusted accepted local update per block visit",
            "TwoAccepted": "at most two trusted accepted local updates per block visit",
        },
        "correctness_contract": ["inner", "gradient", "local", "identified",
                                 "endpoint_trust", "RuntimeConvergence"],
    }
    reports = []
    for topology, atoms in cases:
        for policy in POLICIES:
            reports.append(_run_case(args, output_dir, topology, atoms, policy))
            write(output_dir / "runs.json", reports)
    rows = [summarize(report["result"], report["policy"], report.get("process"))
            for report in reports if report.get("result") is not None]
    analysis = analyze(rows)
    manifest["analysis"] = {"correctness_gate": analysis["correctness_gate"],
                             "selected_policy": analysis["selected_policy"]}
    write(output_dir / "campaign-manifest.json", manifest)
    write_outputs(analysis, output_dir)
    (output_dir / "README.md").write_text(
        "# FixedNeighbor bounded local-work comparison\n\n"
        "This P3 campaign compares FullLocalSearch with one and two trusted accepted "
        "local updates per block visit. All cases use full endpoint assessment and the "
        f"existing {args.rss_limit / 1024**3:g} GiB RSS / {args.timeout:g} second diagnostic envelope.\n\n"
        f"Correctness gate: **{analysis['correctness_gate']}**. Selected policy: "
        f"**{analysis['selected_policy'] or 'none'}**.\n"
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
    parser = build_parser(); args = parser.parse_args(argv)
    args.build_dir = args.build_dir.resolve(); args.output_dir = args.output_dir.resolve()
    if args.timeout <= 0 or args.rss_limit <= 0:
        parser.error("resource limits must be positive")
    if not (args.build_dir / "bin" / "joint_fixed_neighbor_experiment").is_file():
        parser.error("build joint_fixed_neighbor_experiment before running inexact comparison")
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
                      "selected_policy": analysis["selected_policy"],
                      "output_dir": str(args.output_dir)}))
    return 0 if analysis["correctness_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
