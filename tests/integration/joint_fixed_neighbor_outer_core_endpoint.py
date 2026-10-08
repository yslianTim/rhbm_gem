"""Qualify finalist FixedNeighbor outer cores with full endpoint assessment."""
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
from joint_fixed_neighbor_core_size import (
    _finite, _maximum_sweep_value, _parse_core_sizes, _partial, _status, summarize,
)


DEFAULT_CASES = ("chain-256", "cube-256", "chain-512", "cube-512")
FINALIST_CORES = (12, 16, 64)
MODE = "--inexact-one-endpoint"
KKT_LIMIT = 1e-10
WIDTH_GRADIENT_LIMIT = 1e-12
REPLAY_CACHE_LIMIT = 2e-12
REPLAY_OBJECTIVE_BASE = 1e-12
REPLAY_OBJECTIVE_SCALE = 2e-12
PARITY_LIMIT = 1e-10


def _case(value):
    topology, atoms = value.split("-", 1)
    atoms = int(atoms)
    if topology not in ("chain", "cube") or atoms not in (256, 512):
        raise ValueError(f"Unsupported endpoint case: {value}")
    return topology, atoms


def _check_value(fixed, name):
    value = fixed.get(name)
    return value.get("value") if isinstance(value, dict) else value


def _check_status(fixed, name):
    value = fixed.get(name)
    return value.get("status") if isinstance(value, dict) else value


def _inf_difference(lhs, rhs):
    if not isinstance(lhs, list) or not isinstance(rhs, list) or len(lhs) != len(rhs):
        return None
    return max((abs(left - right) for left, right in zip(lhs, rhs)), default=0.0)


def _scaled_ac_difference(candidate, control):
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


def _rank_evidence(fixed):
    assessment = fixed.get("endpoint_assessment") or {}
    return {
        "projected_width_rank": (assessment.get("projected_width") or {}).get("rank"),
        "corrected_jacobian_rank": (assessment.get("corrected_jacobian") or {}).get("rank"),
        "normalized_width_rank": (assessment.get("normalized_width") or {}).get("rank"),
    }


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
              "command": command, "process": process, "status": _status(process),
              "completed": process.get("status") == "completed" and result is not None,
              "result": result}
    write(wrapper_path, report)
    return report


def summarize_endpoint(report):
    result = report.get("result") or {}
    fixed = result.get("fixed_neighbor", result)
    row = summarize({**report, "policy": "OneAccepted",
                     "measurement_scope": "fixed-neighbor-endpoint-assessed"})
    global_kkt = _check_value(fixed, "global_kkt")
    objective = fixed.get("objective")
    objective_limit = (REPLAY_OBJECTIVE_BASE + REPLAY_OBJECTIVE_SCALE * abs(objective)
                       if _finite(objective) else None)
    maximum_cache_replay_error = _maximum_sweep_value(fixed, "cache_replay_error")
    maximum_objective_replay_error = _maximum_sweep_value(fixed, "objective_replay_error")
    row.update({
        "measurement_scope": "fixed-neighbor-endpoint-assessed",
        "global_ac_kkt": global_kkt,
        "width_gradient_inf_norm": fixed.get("width_gradient_inf_norm"),
        "endpoint_certified": fixed.get("endpoint_certified"),
        "runtime_convergence": fixed.get("runtime_convergence"),
        "assessment_inner": _check_status(fixed, "assessment_inner"),
        "assessment_gradient": _check_status(fixed, "assessment_gradient"),
        "assessment_local": _check_status(fixed, "assessment_local"),
        "assessment_identified": _check_status(fixed, "assessment_identified"),
        "endpoint_trust": ((fixed.get("endpoint_assessment") or {}).get("endpoint_trust") or {}).get("passed"),
        "maximum_cache_replay_error": maximum_cache_replay_error,
        "maximum_objective_replay_error": maximum_objective_replay_error,
        "cache_replay_within_limit": (maximum_cache_replay_error is not None and
                                       maximum_cache_replay_error <= REPLAY_CACHE_LIMIT),
        "objective_replay_within_limit": (maximum_objective_replay_error is not None and
                                           objective_limit is not None and
                                           maximum_objective_replay_error <= objective_limit),
        "final_eta": fixed.get("final_eta"),
        "final_beta": fixed.get("final_beta"),
        "final_ac_scaling_weights": fixed.get("final_ac_scaling_weights"),
        "rank_evidence": _rank_evidence(fixed),
    })
    return row


def _endpoint_gate(row):
    return row.get("status") == "completed" and row.get("search_converged") is True and \
        row.get("search_reason") == "block-stationary" and \
        isinstance(row.get("confirmed_stationarity_sweep"), int) and \
        0 < row["confirmed_stationarity_sweep"] <= row.get("sweeps", 0) and \
        _finite(row.get("global_ac_kkt")) and row["global_ac_kkt"] <= KKT_LIMIT and \
        _finite(row.get("width_gradient_inf_norm")) and \
        row["width_gradient_inf_norm"] <= WIDTH_GRADIENT_LIMIT and \
        row.get("cache_replay_within_limit") is True and \
        row.get("objective_replay_within_limit") is True and \
        row.get("endpoint_certified") is True and \
        row.get("runtime_convergence") == "Passed" and \
        all(row.get(name) == "Passed" for name in
            ("assessment_inner", "assessment_gradient", "assessment_local", "assessment_identified")) and \
        row.get("endpoint_trust") is True


def _parity(candidate, control):
    objective_difference = (abs(candidate["objective"] - control["objective"])
                            if _finite(candidate.get("objective")) and _finite(control.get("objective"))
                            else None)
    eta_difference = _inf_difference(candidate.get("final_eta"), control.get("final_eta"))
    ac_difference = _scaled_ac_difference(candidate, control)
    return {
        "objective_abs_difference": objective_difference,
        "eta_inf_difference": eta_difference,
        "scaled_ac_inf_difference": ac_difference,
        "parameter_parity": eta_difference is not None and ac_difference is not None and
            eta_difference <= PARITY_LIMIT and ac_difference <= PARITY_LIMIT,
        "rank_evidence": {"candidate": candidate.get("rank_evidence"),
                          "control": control.get("rank_evidence")},
    }


def analyze(rows, core_sizes=FINALIST_CORES):
    core_sizes = _parse_core_sizes(core_sizes)
    groups = {}
    for row in rows:
        groups.setdefault((row["topology"], row["atoms"]), []).append(row)
    cases = []
    for key, entries in sorted(groups.items()):
        entries = sorted(entries, key=lambda row: core_sizes.index(row["core_atoms"]))
        control = next((row for row in entries if row["core_atoms"] == 64), None)
        correct = [row["core_atoms"] for row in entries if _endpoint_gate(row)]
        parity = {str(row["core_atoms"]): _parity(row, control)
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


def write_outputs(report, output_dir):
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
    rows = [summarize_endpoint(report) for report in reports if report.get("result") is not None]
    analysis = analyze(rows, core_sizes)
    manifest = {
        "schema_version": 1, "phase": analysis["phase"],
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {"driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
                                    "endpoint": sha(Path(__file__))},
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "backend": "SPQR", "eigen_threads": 1,
        "cases": [f"{topology}-{atoms}" for topology, atoms in cases],
        "core_sizes": list(core_sizes), "local_search": "LegacyCompact",
        "local_work": "OneAcceptedUpdate", "block_order": "Forward",
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
    write_outputs(analysis, output_dir)
    (output_dir / "README.md").write_text(
        "# FixedNeighbor outer-core full endpoint qualification\n\n"
        "This C3 campaign qualifies finalist outer cores "
        f"{', '.join(str(core) for core in core_sizes)} on "
        f"{', '.join(f'{topology}-{atoms}' for topology, atoms in cases)}. "
        "It keeps LegacyCompact, OneAcceptedUpdate, Forward serial Gauss-Seidel, "
        "SPQR, one Eigen thread, maximum_sweeps=30, and all existing thresholds.\n\n"
        "Core 64 is the production control. Endpoint parity reports objective, eta, "
        "scaled A/C, and rank evidence; it does not require bitwise-identical trajectories.\n\n"
        f"Qualification gate: **{analysis['qualification_gate']}**. Qualified cores: "
        f"**{', '.join(str(core) for core in analysis['qualified_core_sizes']) or 'none'}**.\n"
    )
    return analysis


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=7200.0)
    parser.add_argument("--rss-limit", type=int, default=RSS_LIMIT_BYTES)
    parser.add_argument("--cases", nargs="+", default=DEFAULT_CASES)
    parser.add_argument("--core-sizes", type=_parse_core_sizes, default=FINALIST_CORES)
    return parser


def main(argv=None):
    parser = build_parser(); args = parser.parse_args(argv)
    args.build_dir = args.build_dir.resolve(); args.output_dir = args.output_dir.resolve()
    args.core_sizes = _parse_core_sizes(args.core_sizes)
    if args.timeout <= 0 or args.rss_limit <= 0:
        parser.error("resource limits must be positive")
    if not (args.build_dir / "bin" / "joint_fixed_neighbor_experiment").is_file():
        parser.error("build joint_fixed_neighbor_experiment before running endpoint qualification")
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
                      "qualified_core_sizes": analysis["qualified_core_sizes"],
                      "output_dir": str(args.output_dir)}))
    return 0 if analysis["qualification_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
