"""Compare selected bounded local-work performance across core sizes."""
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
CORE_SIZES = (64, 128, 256)
POLICY = "OneAccepted"
MODE = "--inexact-one-search"


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


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


def summarize(report):
    result = report.get("result") or {}
    fixed = result.get("fixed_neighbor", result)
    process = report.get("process") or {}
    return {
        "topology": report.get("topology"), "atoms": report.get("atoms"),
        "core_atoms": report.get("core_atoms", result.get("core_atoms")),
        "policy": report.get("policy", POLICY), "status": report.get("status", "completed"),
        "measurement_scope": report.get("measurement_scope", "fixed-neighbor-search-only"),
        "search_converged": fixed.get("search_converged"),
        "search_reason": fixed.get("search_reason"),
        "sweeps": fixed.get("sweeps"),
        "confirmed_stationarity_sweep": fixed.get("confirmed_stationarity_sweep"),
        "block_solves": fixed.get("block_solves", sum(s.get("block_solves", 0)
                                                        for s in fixed.get("sweep_telemetry", []))),
        "profile_evaluations": fixed.get("profile_evaluations", sum(
            s.get("profile_evaluations", 0) for s in fixed.get("sweep_telemetry", []))),
        "accepted_local_updates": fixed.get("accepted_local_updates"),
        "local_work_seconds": fixed.get("local_factor_seconds"),
        "profile_factor_seconds": fixed.get("profile_factor_seconds"),
        "search_seconds": fixed.get("search_seconds"),
        "total_seconds": fixed.get("total_elapsed_seconds"),
        "peak_rss_mb": result.get("peak_rss_mb", process.get("sampled_tree_peak_rss_bytes", 0) / 1024**2
                               if process.get("sampled_tree_peak_rss_bytes") is not None else None),
        "objective": fixed.get("objective"),
        "final_global_ac_kkt": fixed.get("final_global_ac_kkt"),
        "final_raw_width_gradient_inf_norm": fixed.get("final_raw_width_gradient_inf_norm"),
        "maximum_local_columns": fixed.get("maximum_local_columns"),
    }


def _correct(row):
    return row["status"] == "completed" and row["search_converged"] is True and \
        _finite(row.get("final_global_ac_kkt")) and row["final_global_ac_kkt"] <= 1e-10 and \
        _finite(row.get("final_raw_width_gradient_inf_norm")) and \
        row["final_raw_width_gradient_inf_norm"] <= 1e-12


def analyze(rows):
    groups = {}
    for row in rows:
        groups.setdefault((row["topology"], row["atoms"]), []).append(row)
    cases = []
    for key, entries in sorted(groups.items()):
        entries = sorted(entries, key=lambda row: CORE_SIZES.index(row["core_atoms"]))
        correct = [row["core_atoms"] for row in entries if _correct(row)]
        fastest = min((row for row in entries if _correct(row) and _finite(row.get("search_seconds"))),
                      key=lambda row: row["search_seconds"], default=None)
        cases.append({"topology": key[0], "atoms": key[1], "cores": entries,
                      "correct_core_sizes": correct,
                      "fastest_correct_core": fastest["core_atoms"] if fastest else None})
    complete = bool(cases) and all(set(case["correct_core_sizes"]) == set(CORE_SIZES)
                                   for case in cases)
    aggregate = {str(core): sum(row["search_seconds"] for case in cases for row in case["cores"]
                                if row["core_atoms"] == core and _correct(row) and
                                _finite(row.get("search_seconds"))) for core in CORE_SIZES}
    selected = min((int(core) for core, value in aggregate.items() if value > 0),
                   key=lambda core: aggregate[str(core)], default=None) if complete else None
    return {"phase": "P5 FixedNeighbor core-size performance study",
            "policy": POLICY, "core_sizes": list(CORE_SIZES), "cases": cases,
            "correctness_gate": "passed" if complete else "failed",
            "aggregate_search_seconds_by_core": aggregate,
            "selected_core_size": selected,
            "selection_basis": "minimum aggregate search seconds among cores correct on every study case; resource and peak RSS remain reported"}


def write_outputs(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = ["topology", "atoms", "core_atoms", "policy", "status", "sweeps",
               "confirmed_stationarity_sweep", "block_solves", "profile_evaluations",
               "accepted_local_updates", "local_work_seconds", "profile_factor_seconds",
               "search_seconds", "total_seconds", "peak_rss_mb", "search_converged",
               "search_reason", "objective", "final_global_ac_kkt",
               "final_raw_width_gradient_inf_norm", "maximum_local_columns"]
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
    reports = []
    for topology, atoms in cases:
        for core_atoms in CORE_SIZES:
            reports.append(_run_case(args, output_dir, topology, atoms, core_atoms))
            write(output_dir / "runs.json", reports)
    rows = [summarize(report) for report in reports if report.get("result") is not None]
    analysis = analyze(rows)
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
        "policy": POLICY, "core_sizes": list(CORE_SIZES),
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
        "# FixedNeighbor core-size performance study\n\n"
        "This P5 study keeps the P4-selected OneAcceptedLocalUpdate policy, forward "
        "serial Gauss-Seidel order, SPQR backend, one Eigen thread, and frozen stationarity "
        "checks. It varies only core_atoms over 64, 128, and 256 for chain-512, cube-512, "
        "and cube-1024. All runs are search-only; no full endpoint claim is made for 1024.\n\n"
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
    return parser


def main(argv=None):
    parser = build_parser(); args = parser.parse_args(argv)
    args.build_dir = args.build_dir.resolve(); args.output_dir = args.output_dir.resolve()
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
