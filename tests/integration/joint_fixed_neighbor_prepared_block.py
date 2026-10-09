"""Check prepared FixedNeighbor structural invariants on small deterministic cases."""
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


CORE_ATOMS = 12
SMALL_CASE_ATOMS = 24
DEFAULT_CASES = tuple(f"{topology}-{SMALL_CASE_ATOMS}" for topology in ("chain", "cube"))
POLICY = "FixedNeighbor"
MODE = "--scaling-only"


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _case(value):
    topology, atoms = value.split("-", 1)
    atoms = int(atoms)
    if topology not in ("chain", "cube") or atoms != SMALL_CASE_ATOMS:
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
    fixed.setdefault("accepted_local_updates", None)
    for key in ("prepared_block_count", "block_preparations", "domain_preparations",
                "mapping_preparations", "symbolic_factorizations", "symbolic_reuses",
                "numeric_factorizations", "fresh_workspace_symbolic_factorizations"):
        fixed.setdefault(key, None)
    fixed.setdefault("final_global_ac_kkt", None)
    fixed.setdefault("final_raw_width_gradient_inf_norm", None)
    fixed.setdefault("search_seconds", sum(row.get("wall_seconds", 0.0)
                                             for row in progress.get("sweep_telemetry", [])))
    fixed.setdefault("total_elapsed_seconds", fixed["search_seconds"])
    return {"topology": progress.get("topology"), "atoms": progress.get("atoms"),
            "core_atoms": CORE_ATOMS, "measurement_scope": "fixed-neighbor-search-only",
            "fixed_neighbor": fixed, "peak_rss_mb": None}


def _run_case(args, output_dir, topology, atoms):
    case = f"{topology}-{atoms}"
    individual = output_dir / "individual-results"
    raw_path = individual / f"{case}.json"
    wrapper_path = individual / f"{case}-run.json"
    process_dir = individual / f"{case}-process"
    if wrapper_path.is_file():
        return read(wrapper_path)
    command = [str(args.build_dir / "bin" / "joint_fixed_neighbor_experiment"), MODE,
               str(raw_path), topology, str(atoms), str(CORE_ATOMS)]
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


def summarize(report):
    result = report.get("result") or report
    fixed = result.get("fixed_neighbor", result)
    process = report.get("process") or {}
    keys = ("prepared_block_count", "block_preparations", "domain_preparations",
            "mapping_preparations", "symbolic_factorizations", "symbolic_reuses",
            "numeric_factorizations", "fresh_workspace_symbolic_factorizations")
    return {
        "topology": report.get("topology"), "atoms": report.get("atoms"),
        "policy": report.get("policy", POLICY), "status": report.get("status", "completed"),
        "measurement_scope": report.get("measurement_scope", "fixed-neighbor-search-only"),
        "core_atoms": result.get("core_atoms", CORE_ATOMS),
        "search_converged": fixed.get("search_converged"),
        "search_reason": fixed.get("search_reason"),
        "sweeps": fixed.get("sweeps"),
        "confirmed_stationarity_sweep": fixed.get("confirmed_stationarity_sweep"),
        "blocks_per_sweep": fixed.get("blocks_per_sweep"),
        "block_solves": fixed.get("block_solves"),
        "profile_evaluations": fixed.get("profile_evaluations"),
        "accepted_local_updates": fixed.get("accepted_local_updates"),
        "final_objective": fixed.get("objective"),
        "final_global_ac_kkt": fixed.get("final_global_ac_kkt"),
        "final_width_gradient_inf_norm": fixed.get("final_raw_width_gradient_inf_norm"),
        "search_seconds": fixed.get("search_seconds"),
        "assessment_seconds": fixed.get("assessment_seconds"),
        "total_seconds": fixed.get("total_elapsed_seconds"),
        "peak_rss_mb": result.get("peak_rss_mb", process.get("sampled_tree_peak_rss_bytes", 0) / 1024**2
                               if process.get("sampled_tree_peak_rss_bytes") is not None else None),
        **{key: fixed.get(key) for key in keys},
    }


def _numerical_gate(row):
    confirmed = row.get("confirmed_stationarity_sweep")
    return row["status"] == "completed" and row["measurement_scope"] == "fixed-neighbor-search-only" and \
        row.get("search_converged") is True and isinstance(confirmed, int) and confirmed > 0 and \
        isinstance(row.get("sweeps"), int) and confirmed <= row["sweeps"] and \
        _finite(row.get("final_global_ac_kkt")) and row["final_global_ac_kkt"] <= 1e-10 and \
        _finite(row.get("final_width_gradient_inf_norm")) and \
        row["final_width_gradient_inf_norm"] <= 1e-12


def _preparation_gate(row):
    blocks = row.get("blocks_per_sweep")
    return row.get("core_atoms") == CORE_ATOMS and isinstance(blocks, int) and blocks > 0 and \
        row.get("prepared_block_count") == blocks and \
        row.get("block_preparations") == blocks and \
        row.get("domain_preparations") == blocks and \
        row.get("mapping_preparations") == blocks


def _symbolic_gate(row):
    symbolic = row.get("symbolic_factorizations")
    reuses = row.get("symbolic_reuses")
    numeric = row.get("numeric_factorizations")
    fresh = row.get("fresh_workspace_symbolic_factorizations")
    return all(isinstance(value, int) for value in (symbolic, reuses, numeric, fresh)) and \
        symbolic > 0 and reuses > 0 and numeric >= symbolic and fresh == numeric and \
        symbolic + reuses == numeric


def analyze(rows, expected=None):
    expected = set(expected or ((_case(value) for value in DEFAULT_CASES)))
    by_case = {(row["topology"], row["atoms"]): row for row in rows}
    cases = []
    failures = []
    for key in sorted(expected):
        row = by_case.get(key)
        reasons = []
        if row is None:
            reasons.append("missing-case")
        else:
            if not _numerical_gate(row):
                reasons.append("numerical-gate")
            if not _preparation_gate(row):
                reasons.append("prepared-block-gate")
            if not _symbolic_gate(row):
                reasons.append("symbolic-reuse-gate")
            if reasons:
                failures.append({"topology": key[0], "atoms": key[1], "reasons": reasons})
        cases.append({"topology": key[0], "atoms": key[1], "passed": not reasons,
                      "reasons": reasons, "metrics": row})
    return {
        "phase": "prepared FixedNeighbor structural contract",
        "policy": {"core_atoms": CORE_ATOMS, "maximum_sweeps": 30,
                   "block_order": "forward",
                   "local_search": "LegacyCompact", "backend": "SPQR", "eigen_threads": 1},
        "measurement_scope": "fixed-neighbor-search-only",
        "expected_cases": [{"topology": topology, "atoms": atoms} for topology, atoms in sorted(expected)],
        "cases": cases, "failures": failures,
        "symbolic_attribution": "fresh workspace symbolic count is the observed numeric factorization count; persistent count and reuse are measured by SparseWork",
        "wall_time_gate": "not-run",
        "qualification_gate": "passed" if not failures and len(cases) == len(expected) else "failed",
    }


def write_outputs(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    write(output_dir / "analysis.json", report)
    columns = ["topology", "atoms", "status", "sweeps", "confirmed_stationarity_sweep",
               "block_solves", "profile_evaluations", "accepted_local_updates", "search_seconds",
               "assessment_seconds", "total_seconds", "peak_rss_mb", "final_global_ac_kkt",
               "final_width_gradient_inf_norm", "prepared_block_count", "block_preparations",
               "domain_preparations", "mapping_preparations", "symbolic_factorizations",
               "symbolic_reuses", "numeric_factorizations", "fresh_workspace_symbolic_factorizations"]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            if case["metrics"] is not None:
                writer.writerow({column: case["metrics"].get(column) for column in columns})


def run_campaign(args):
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = [_case(value) for value in args.cases]
    reports = []
    for topology, atoms in cases:
        reports.append(_run_case(args, output_dir, topology, atoms))
        write(output_dir / "runs.json", reports)
    rows = [summarize(report) for report in reports if report.get("result") is not None]
    analysis = analyze(rows, cases)
    manifest = {
        "schema_version": 1, "phase": analysis["phase"],
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {"driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
                                    "qualification": sha(Path(__file__))},
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "cases": [f"{topology}-{atoms}" for topology, atoms in cases],
        "policy": analysis["policy"], "measurement_scope": analysis["measurement_scope"],
        "resource_envelope": {"wall_seconds": args.timeout, "rss_bytes": args.rss_limit},
        "analysis": {"qualification_gate": analysis["qualification_gate"],
                     "wall_time_gate": analysis["wall_time_gate"]},
    }
    write(output_dir / "campaign-manifest.json", manifest)
    write_outputs(analysis, output_dir)
    (output_dir / "README.md").write_text(
        "# Prepared FixedNeighbor structural contract\n\n"
        "This small deterministic check exercises the production FixedNeighbor policy on chain/cube 24 "
        "with a 12-atom core, producing two prepared blocks. It checks preparation counts and "
        "symbolic reuse invariants; it is not a performance campaign and makes no large-case claim.\n\n"
        f"Qualification gate: **{analysis['qualification_gate']}**. Wall-time gate: **{analysis['wall_time_gate']}**.\n"
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
                      "wall_time_gate": analysis["wall_time_gate"],
                      "output_dir": str(args.output_dir)}))
    return 0 if analysis["qualification_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
