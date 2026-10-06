"""Qualify the selected bounded local-work policy at 1024 atoms."""
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


DEFAULT_CASES = ("chain-1024", "cube-1024")
POLICY = "OneAccepted"
MODE = "--inexact-one-search"
BASELINE_DIR = ROOT / "docs/developer/figures/joint-fixed-neighbor-scaling-r1/individual-results"


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _case(value):
    topology, atoms = value.split("-", 1)
    atoms = int(atoms)
    if topology not in ("chain", "cube") or atoms != 1024:
        raise ValueError(f"Unsupported qualification case: {value}")
    return topology, atoms


def _status(process):
    if process.get("status") == "completed":
        return "completed"
    if process.get("status") == "time-limit":
        return "timeout"
    if process.get("status") in ("rss-limit", "rss-limit-observed-after-exit"):
        return "rss_limit"
    return "process_error"


def _partial(progress):
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
    return {"topology": progress.get("topology"), "atoms": progress.get("atoms"),
            "core_atoms": 128, "measurement_scope": "fixed-neighbor-search-only",
            "fixed_neighbor": fixed, "peak_rss_mb": None}


def _run_candidate(args, output_dir, topology, atoms):
    case = f"{topology}-{atoms}"
    individual = output_dir / "individual-results"
    raw_path = individual / f"{case}-{POLICY}.json"
    wrapper_path = individual / f"{case}-{POLICY}-run.json"
    process_dir = individual / f"{case}-{POLICY}-process"
    if wrapper_path.is_file():
        return read(wrapper_path)
    command = [str(args.build_dir / "bin" / "joint_fixed_neighbor_experiment"), MODE,
               str(raw_path), topology, str(atoms)]
    print(f"Running {case} {POLICY} with {args.timeout:g}s cap", flush=True)
    process = monitored(command, process_dir, time.monotonic() + args.timeout,
                        rss_limit=args.rss_limit, seconds=args.timeout)
    result = read(raw_path) if raw_path.is_file() else None
    if result is None:
        progress_path = Path(str(raw_path) + ".progress.json")
        result = _partial(read(progress_path) if progress_path.is_file() else None)
    report = {"topology": topology, "atoms": atoms, "policy": POLICY,
              "measurement_scope": "fixed-neighbor-search-only", "command": command,
              "process": process, "status": _status(process),
              "completed": process.get("status") == "completed" and result is not None,
              "result": result}
    write(wrapper_path, report)
    return report


def _baseline(topology, atoms):
    path = BASELINE_DIR / f"{topology}-{atoms}.json"
    if not path.is_file():
        raise FileNotFoundError(f"missing FullLocalSearch baseline: {path}")
    result = read(path)
    return {"topology": topology, "atoms": atoms, "policy": "Full",
            "measurement_scope": result.get("measurement_scope", "fixed-neighbor-search-only"),
            "status": "completed", "completed": True, "result": result,
            "source_file": str(path.relative_to(ROOT))}


def _factor_seconds(fixed):
    """Return the comparable local profile-work field from both campaigns."""
    return fixed.get("local_factor_seconds")


def summarize(report):
    result = report.get("result") or {}
    fixed = result.get("fixed_neighbor", result)
    process = report.get("process") or {}
    return {
        "topology": report.get("topology"), "atoms": report.get("atoms"),
        "policy": report.get("policy"), "status": report.get("status", "completed"),
        "measurement_scope": report.get("measurement_scope",
                                         result.get("measurement_scope", "fixed-neighbor-search-only")),
        "source_file": report.get("source_file"),
        "search_converged": fixed.get("search_converged"),
        "search_reason": fixed.get("search_reason"),
        "sweeps": fixed.get("sweeps"),
        "confirmed_stationarity_sweep": fixed.get("confirmed_stationarity_sweep"),
        "block_solves": fixed.get("block_solves", sum(s.get("block_solves", 0)
                                                        for s in fixed.get("sweep_telemetry", []))),
        "profile_evaluations": fixed.get("profile_evaluations", sum(
            s.get("profile_evaluations", 0) for s in fixed.get("sweep_telemetry", []))),
        "accepted_local_updates": fixed.get("accepted_local_updates"),
        "local_work_seconds": _factor_seconds(fixed),
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


def _search_gate(row):
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
        entries = sorted(entries, key=lambda row: (0 if row["policy"] == "Full" else 1))
        full = next((row for row in entries if row["policy"] == "Full"), None)
        candidate = next((row for row in entries if row["policy"] == POLICY), None)
        comparable = bool(full and candidate and _finite(full.get("local_work_seconds")) and
                          _finite(candidate.get("local_work_seconds")) and
                          _finite(full.get("search_seconds")) and
                          _finite(candidate.get("search_seconds")))
        if comparable:
            candidate["local_work_reduction"] = 1.0 - candidate["local_work_seconds"] / full["local_work_seconds"]
            candidate["search_reduction"] = 1.0 - candidate["search_seconds"] / full["search_seconds"]
            candidate["search_speedup"] = full["search_seconds"] / candidate["search_seconds"]
        cases.append({"topology": key[0], "atoms": key[1], "full": full,
                      "candidate": candidate, "candidate_search_gate": bool(candidate and _search_gate(candidate)),
                      "comparable": comparable,
                      "material_local_work_reduction": bool(comparable and candidate["local_work_reduction"] > 0),
                      "material_search_reduction": bool(comparable and candidate["search_reduction"] > 0)})
    gate = bool(cases) and all(case["candidate_search_gate"] and case["comparable"] and
                               case["material_local_work_reduction"] and
                               case["material_search_reduction"] for case in cases)
    return {"phase": "P4 bounded FixedNeighbor local-work qualification",
            "candidate_policy": POLICY, "cases": cases,
            "qualification_gate": "passed" if gate else "failed",
            "selected_policy": POLICY if gate else None,
            "selection_basis": "smallest P3-correct bounded budget with positive measured local-work and search reductions at 1024; no preset percentage threshold",
            "measurement_note": "local_work_seconds is the comparable legacy local profile-work field; profile_factor_seconds is reported when available but the pre-P3 Full baseline does not contain it."}


def write_outputs(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = ["topology", "atoms", "policy", "status", "sweeps", "confirmed_stationarity_sweep",
               "block_solves", "profile_evaluations", "accepted_local_updates", "local_work_seconds",
               "profile_factor_seconds", "search_seconds", "total_seconds", "peak_rss_mb",
               "search_converged", "search_reason", "final_global_ac_kkt",
               "final_raw_width_gradient_inf_norm", "local_work_reduction", "search_reduction",
               "search_speedup", "source_file"]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            for row in (case["full"], case["candidate"]):
                if row is not None:
                    writer.writerow({column: row.get(column) for column in columns})


def run_campaign(args):
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = [_case(value) for value in args.cases]
    baseline_manifest = ROOT / "docs/developer/figures/joint-fixed-neighbor-scaling-r1/campaign-manifest.json"
    reports = []
    for topology, atoms in cases:
        reports.append(_baseline(topology, atoms))
        reports.append(_run_candidate(args, output_dir, topology, atoms))
        write(output_dir / "runs.json", reports)
    rows = [summarize(report) for report in reports if report.get("result") is not None]
    analysis = analyze(rows)
    manifest = {
        "schema_version": 1, "phase": analysis["phase"],
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {"driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
                                    "qualification": sha(Path(__file__))},
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "backend": "SPQR", "eigen_threads": 1, "cases": [f"{t}-{a}" for t, a in cases],
        "candidate_policy": POLICY, "measurement_scope": "fixed-neighbor-search-only",
        "resource_envelope": {"wall_seconds": args.timeout, "rss_bytes": args.rss_limit},
        "full_baseline_campaign": str(baseline_manifest.relative_to(ROOT)),
        "full_baseline_manifest_sha256": sha(baseline_manifest),
        "analysis": {"qualification_gate": analysis["qualification_gate"],
                     "selected_policy": analysis["selected_policy"]},
    }
    write(output_dir / "campaign-manifest.json", manifest)
    write_outputs(analysis, output_dir)
    (output_dir / "README.md").write_text(
        "# FixedNeighbor bounded local-work qualification\n\n"
        "This P4 campaign runs the P3-selected OneAcceptedLocalUpdate policy at "
        "chain-1024 and cube-1024 in search-only mode. FullLocalSearch values are "
        "the existing matched search-only frontier from "
        "`joint-fixed-neighbor-scaling-r1`; the candidate uses the same SPQR backend, "
        "128-atom core, forward order, one Eigen thread, and frozen stationarity checks.\n\n"
        "`local_work_seconds` is the comparable legacy local profile-work field. "
        "The candidate also reports the newer `profile_factor_seconds` telemetry; the "
        "pre-P3 Full baseline does not contain that field. No full endpoint claim is "
        "made for these 1024 points; endpoint correctness was gated separately in P3.\n\n"
        f"Qualification gate: **{analysis['qualification_gate']}**. Selected policy: "
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
    parser = build_parser()
    args = parser.parse_args(argv)
    args.build_dir = args.build_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    if args.timeout <= 0 or args.rss_limit <= 0:
        parser.error("resource limits must be positive")
    if not (args.build_dir / "bin" / "joint_fixed_neighbor_experiment").is_file():
        parser.error("build joint_fixed_neighbor_experiment before running qualification")
    for value in args.cases:
        try:
            _case(value)
        except (ValueError, TypeError) as error:
            parser.error(str(error))
    try:
        analysis = run_campaign(args)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.error(str(error))
    print(json.dumps({"qualification_gate": analysis["qualification_gate"],
                      "selected_policy": analysis["selected_policy"],
                      "output_dir": str(args.output_dir)}))
    return 0 if analysis["qualification_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
