"""Run a matched fresh-workspace/persistent-workspace FixedNeighbor timing campaign."""
from __future__ import annotations

import argparse
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


DEFAULT_CASES = tuple(
    f"{topology}-{atoms}"
    for topology in ("chain", "cube")
    for atoms in (256, 512, 1024)
)
CORE_ATOMS = 64
CONTROL = "fresh-workspace"
TREATMENT = "persistent-workspace"
MODES = {CONTROL: "--matched-control", TREATMENT: "--inexact-one-search"}


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _case(value):
    topology, atoms = value.split("-", 1)
    atoms = int(atoms)
    if topology not in ("chain", "cube") or atoms not in (256, 512, 1024):
        raise ValueError(f"Unsupported timing case: {value}")
    return topology, atoms


def _status(process):
    if process.get("status") == "completed":
        return "completed"
    if process.get("status") == "time-limit":
        return "timeout"
    if process.get("status") in ("rss-limit", "rss-limit-observed-after-exit"):
        return "rss_limit"
    return "process_error"


def _fixed(report):
    result = report.get("result") or {}
    return result.get("fixed_neighbor", result)


def _work(report):
    return _fixed(report).get("fixed_neighbor_work", {})


def _run_case(args, output_dir, topology, atoms, variant, repetition, warmup):
    case = f"{topology}-{atoms}"
    label = f"{case}-{variant}-{'warmup' if warmup else 'measurement'}-{repetition}"
    individual = output_dir / "individual-results"
    raw_path = individual / f"{label}.json"
    wrapper_path = individual / f"{label}-run.json"
    if wrapper_path.is_file():
        return read(wrapper_path)
    command = [str(args.build_dir / "bin" / "joint_fixed_neighbor_experiment"), MODES[variant],
               str(raw_path), topology, str(atoms), str(CORE_ATOMS)]
    print(f"Running {label} with {args.timeout:g}s cap", flush=True)
    process_dir = individual / f"{label}-process"
    process = monitored(command, process_dir, time.monotonic() + args.timeout,
                        rss_limit=args.rss_limit, seconds=args.timeout)
    result = read(raw_path) if raw_path.is_file() else None
    report = {"topology": topology, "atoms": atoms, "variant": variant,
              "repetition": repetition, "warmup": warmup, "command": command,
              "process": process, "status": _status(process),
              "completed": process.get("status") == "completed" and result is not None,
              "result": result}
    write(wrapper_path, report)
    return report


def _summary(report):
    fixed = _fixed(report)
    result = report.get("result") or {}
    process = report.get("process") or {}
    work = _work(report)
    return {
        "topology": report.get("topology"), "atoms": report.get("atoms"),
        "variant": report.get("variant"), "repetition": report.get("repetition"),
        "warmup": report.get("warmup", False), "status": report.get("status"),
        "workspace_mode": result.get("workspace_mode", fixed.get("workspace_mode")),
        "measurement_scope": result.get("measurement_scope"),
        "search_converged": fixed.get("search_converged"),
        "search_reason": fixed.get("search_reason"),
        "sweeps": fixed.get("sweeps"),
        "first_order_stationarity_sweep": fixed.get("first_order_stationarity_sweep"),
        "confirmed_stationarity_sweep": fixed.get("confirmed_stationarity_sweep"),
        "block_solves": fixed.get("block_solves"),
        "profile_evaluations": fixed.get("profile_evaluations"),
        "accepted_local_updates": fixed.get("accepted_local_updates"),
        "objective": fixed.get("objective"),
        "final_global_ac_kkt": fixed.get("final_global_ac_kkt"),
        "final_width_gradient_inf_norm": fixed.get("final_raw_width_gradient_inf_norm"),
        "search_seconds": fixed.get("search_seconds"),
        "total_seconds": fixed.get("total_elapsed_seconds"),
        "peak_rss_mb": result.get("peak_rss_mb", process.get("sampled_tree_peak_rss_bytes", 0) / 1024**2
                               if process.get("sampled_tree_peak_rss_bytes") is not None else None),
        "symbolic_factorizations": fixed.get("symbolic_factorizations"),
        "symbolic_reuses": fixed.get("symbolic_reuses"),
        "numeric_factorizations": fixed.get("numeric_factorizations"),
        "matrix_preparation_seconds": fixed.get("matrix_preparation_seconds"),
        "symbolic_seconds": fixed.get("symbolic_seconds"),
        "numeric_seconds": fixed.get("numeric_seconds"),
        "full_candidate_replays": work.get("full_candidate_replays"),
        "candidate_state_full_copies": work.get("candidate_state_full_copies"),
        "candidate_replay_seconds": work.get("candidate_replay_seconds"),
        "candidate_copy_seconds": work.get("candidate_copy_seconds"),
        "fixed_neighbor_work": work,
        "final_eta": fixed.get("final_eta"), "final_beta": fixed.get("final_beta"),
    }


def _median(rows, key):
    values = [row.get(key) for row in rows if _finite(row.get(key))]
    return statistics.median(values) if values else None


def _scaled_difference(lhs, rhs):
    if not isinstance(lhs, list) or not isinstance(rhs, list) or len(lhs) != len(rhs) or not lhs:
        return math.inf
    return max((abs(left - right) / (1.0 + max(abs(left), abs(right)))
                for left, right in zip(lhs, rhs)), default=math.inf)


def _trajectory(row):
    result = row.get("result") or {}
    fixed = result.get("fixed_neighbor", result)
    blocks = fixed.get("block_telemetry") or []
    return [(block.get("sweep"), block.get("block"), block.get("status"),
             block.get("accepted"), block.get("accepted_updates")) for block in blocks]


def compare(control, treatment):
    reasons = []
    if control.get("status") != "completed" or treatment.get("status") != "completed":
        reasons.append("process-not-completed")
    for key in ("search_converged", "search_reason", "sweeps", "first_order_stationarity_sweep",
                "confirmed_stationarity_sweep", "block_solves", "profile_evaluations",
                "accepted_local_updates"):
        if control.get(key) != treatment.get(key):
            reasons.append(f"{key}-mismatch")
    for key, threshold in (("final_global_ac_kkt", 1e-10), ("final_width_gradient_inf_norm", 1e-12)):
        for label, row in (("control", control), ("treatment", treatment)):
            if not _finite(row.get(key)) or row[key] > threshold:
                reasons.append(f"{label}-{key}-gate")
    if not _finite(control.get("objective")) or not _finite(treatment.get("objective")):
        reasons.append("objective-unavailable")
        objective_difference = math.inf
    else:
        objective_difference = abs(control["objective"] - treatment["objective"])
        if objective_difference > 1e-12:
            reasons.append("objective-mismatch")
    eta_difference = _scaled_difference(control.get("final_eta"), treatment.get("final_eta"))
    beta_difference = _scaled_difference(control.get("final_beta"), treatment.get("final_beta"))
    if eta_difference > 1e-10:
        reasons.append("eta-endpoint-mismatch")
    if beta_difference > 1e-10:
        reasons.append("beta-endpoint-mismatch")
    if _trajectory(control) != _trajectory(treatment):
        reasons.append("block-trajectory-mismatch")
    return {
        "passed": not reasons, "reasons": reasons,
        "objective_difference": objective_difference,
        "eta_endpoint_difference": eta_difference,
        "beta_endpoint_difference": beta_difference,
        "endpoint_certification": "not-run-search-only",
    }


def _performance_class(improvement):
    if not _finite(improvement):
        return "not-run"
    if improvement >= 0.15:
        return "material"
    if improvement >= 0.05:
        return "measurable-secondary"
    if improvement >= -0.05:
        return "not-major"
    return "regression"


def _case_report(rows, warmup, measurements):
    measured = [row for row in rows if not row.get("warmup")]
    controls = [row for row in measured if row.get("variant") == CONTROL]
    treatments = [row for row in measured if row.get("variant") == TREATMENT]
    by_repetition = {}
    for row in measured:
        by_repetition.setdefault(row.get("repetition"), {})[row.get("variant")] = row
    comparisons = [compare(pair[CONTROL], pair[TREATMENT]) for pair in by_repetition.values()
                   if CONTROL in pair and TREATMENT in pair]
    control_search = _median(controls, "search_seconds")
    treatment_search = _median(treatments, "search_seconds")
    control_total = _median(controls, "total_seconds")
    treatment_total = _median(treatments, "total_seconds")
    improvement = ((control_search - treatment_search) / control_search
                   if _finite(control_search) and _finite(treatment_search) and control_search > 0 else math.nan)
    candidate_fractions = [
        (row.get("candidate_replay_seconds", 0.0) + row.get("candidate_copy_seconds", 0.0)) /
        row["search_seconds"]
        for row in treatments
        if _finite(row.get("search_seconds")) and row["search_seconds"] > 0
    ]
    return {
        "topology": rows[0].get("topology") if rows else None,
        "atoms": rows[0].get("atoms") if rows else None,
        "warmup_runs": warmup, "measurement_runs": measurements,
        "control": {"search_seconds_median": control_search, "total_seconds_median": control_total,
                     "peak_rss_mb_median": _median(controls, "peak_rss_mb"),
                     "symbolic_factorizations_median": _median(controls, "symbolic_factorizations"),
                     "symbolic_reuses_median": _median(controls, "symbolic_reuses"),
                     "numeric_factorizations_median": _median(controls, "numeric_factorizations"),
                     "matrix_preparation_seconds_median": _median(controls, "matrix_preparation_seconds"),
                     "symbolic_seconds_median": _median(controls, "symbolic_seconds"),
                     "numeric_seconds_median": _median(controls, "numeric_seconds")},
        "treatment": {"search_seconds_median": treatment_search, "total_seconds_median": treatment_total,
                       "peak_rss_mb_median": _median(treatments, "peak_rss_mb"),
                       "symbolic_factorizations_median": _median(treatments, "symbolic_factorizations"),
                       "symbolic_reuses_median": _median(treatments, "symbolic_reuses"),
                       "numeric_factorizations_median": _median(treatments, "numeric_factorizations"),
                       "matrix_preparation_seconds_median": _median(treatments, "matrix_preparation_seconds"),
                       "symbolic_seconds_median": _median(treatments, "symbolic_seconds"),
                       "numeric_seconds_median": _median(treatments, "numeric_seconds"),
                       "candidate_replay_copy_fraction_median": (
                           statistics.median(candidate_fractions) if candidate_fractions else None)},
        "relative_search_improvement": improvement,
        "relative_total_improvement": ((control_total - treatment_total) / control_total
                                        if _finite(control_total) and _finite(treatment_total) and control_total > 0
                                        else math.nan),
        "performance_class": _performance_class(improvement),
        "numerical_gate": "passed" if comparisons and all(item["passed"] for item in comparisons) else "failed",
        "comparisons": comparisons,
    }


def analyze(reports, expected=None, warmup=1, measurements=3):
    expected = set(expected or (_case(value) for value in DEFAULT_CASES))
    rows = [_summary(report) for report in reports]
    grouped = {(row["topology"], row["atoms"]): [] for row in rows}
    for row in rows:
        grouped.setdefault((row["topology"], row["atoms"]), []).append(row)
    cases = [_case_report(grouped[key], warmup, measurements) if grouped.get(key) else
             {"topology": key[0], "atoms": key[1], "numerical_gate": "failed",
              "performance_class": "not-run", "comparisons": []}
             for key in sorted(expected)]
    valid_large = [case for case in cases if case.get("atoms") == 1024]
    material_large = [case for case in valid_large if case.get("performance_class") == "material"]
    regressions = [case for case in valid_large if case.get("performance_class") == "regression"]
    all_numerical = all(case.get("numerical_gate") == "passed" for case in cases)
    if all_numerical and len(material_large) >= 2 and not regressions:
        wall_time_gate = "passed"
    elif all_numerical:
        wall_time_gate = "not-material"
    else:
        wall_time_gate = "failed"
    return {
        "phase": "matched FixedNeighbor prepared-workspace timing",
        "measurement_scope": "fixed-neighbor-search-only",
        "policy": {"core_atoms": CORE_ATOMS, "maximum_sweeps": 30, "block_order": "forward",
                    "local_work": "OneAcceptedUpdate", "local_search": "LegacyCompact",
                    "backend": "SPQR", "eigen_threads": 1},
        "control": CONTROL, "treatment": TREATMENT,
        "expected_cases": [{"topology": topology, "atoms": atoms} for topology, atoms in sorted(expected)],
        "cases": cases,
        "wall_time_gate": wall_time_gate,
        "qualification_gate": "passed" if all_numerical else "failed",
        "interpretation": {"material": ">=15%", "measurable-secondary": "5-15%",
                            "not-major": "<5%", "regression": "more than 5% slower"},
    }


def write_outputs(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    write(output_dir / "analysis.json", report)
    columns = ["topology", "atoms", "numerical_gate", "performance_class",
               "control_search_seconds_median", "treatment_search_seconds_median",
               "relative_search_improvement", "relative_total_improvement",
               "treatment_candidate_replay_copy_fraction_median",
               "control_symbolic_factorizations_median", "control_symbolic_reuses_median",
               "treatment_symbolic_factorizations_median", "treatment_symbolic_reuses_median"]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            writer.writerow({
                "topology": case.get("topology"), "atoms": case.get("atoms"),
                "numerical_gate": case.get("numerical_gate"),
                "performance_class": case.get("performance_class"),
                "control_search_seconds_median": (case.get("control") or {}).get("search_seconds_median"),
                "treatment_search_seconds_median": (case.get("treatment") or {}).get("search_seconds_median"),
                "relative_search_improvement": case.get("relative_search_improvement"),
                "relative_total_improvement": case.get("relative_total_improvement"),
                "treatment_candidate_replay_copy_fraction_median":
                    (case.get("treatment") or {}).get("candidate_replay_copy_fraction_median"),
                "control_symbolic_factorizations_median":
                    (case.get("control") or {}).get("symbolic_factorizations_median"),
                "control_symbolic_reuses_median": (case.get("control") or {}).get("symbolic_reuses_median"),
                "treatment_symbolic_factorizations_median":
                    (case.get("treatment") or {}).get("symbolic_factorizations_median"),
                "treatment_symbolic_reuses_median":
                    (case.get("treatment") or {}).get("symbolic_reuses_median"),
            })


def run_campaign(args):
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = [_case(value) for value in args.cases]
    reports = []
    for repetition in range(args.warmup):
        for topology, atoms in cases:
            for variant in (CONTROL, TREATMENT):
                reports.append(_run_case(args, output_dir, topology, atoms, variant, repetition, True))
    for repetition in range(args.measurements):
        for topology, atoms in cases:
            for variant in (CONTROL, TREATMENT):
                reports.append(_run_case(args, output_dir, topology, atoms, variant, repetition, False))
        write(output_dir / "runs.json", reports)
    report = analyze(reports, cases, args.warmup, args.measurements)
    manifest = {
        "schema_version": 1, "phase": report["phase"],
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {"driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
                                    "timing": sha(Path(__file__))},
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "cases": [f"{topology}-{atoms}" for topology, atoms in cases],
        "warmup_runs": args.warmup, "measurements": args.measurements,
        "policy": report["policy"], "measurement_scope": report["measurement_scope"],
        "resource_envelope": {"wall_seconds": args.timeout, "rss_bytes": args.rss_limit},
        "analysis": {"qualification_gate": report["qualification_gate"],
                     "wall_time_gate": report["wall_time_gate"]},
    }
    write(output_dir / "campaign-manifest.json", manifest)
    write_outputs(report, output_dir)
    (output_dir / "README.md").write_text(
        "# Matched FixedNeighbor workspace timing\n\n"
        "Control uses a fresh block-local LinearWorkspace on every block visit. Treatment uses the "
        "persistent prepared-block workspace and SPQR symbolic reuse. Both are search-only and use "
        f"{args.warmup} warmup run(s) plus {args.measurements} measured run(s).\n\n"
        f"Numerical qualification: **{report['qualification_gate']}**. "
        f"Wall-time gate: **{report['wall_time_gate']}**.\n"
    )
    return report


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=7200.0)
    parser.add_argument("--rss-limit", type=int, default=RSS_LIMIT_BYTES)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--measurements", type=int, default=3)
    parser.add_argument("--cases", nargs="+", default=DEFAULT_CASES)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.build_dir = args.build_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    if args.warmup < 1 or args.measurements < 3 or args.timeout <= 0 or args.rss_limit <= 0:
        parser.error("warmup must be >=1, measurements must be >=3, and resource limits must be positive")
    if not (args.build_dir / "bin" / "joint_fixed_neighbor_experiment").is_file():
        parser.error("build joint_fixed_neighbor_experiment before running timing campaign")
    for value in args.cases:
        try:
            _case(value)
        except (ValueError, TypeError) as error:
            parser.error(str(error))
    try:
        report = run_campaign(args)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.error(str(error))
    print(json.dumps({"qualification_gate": report["qualification_gate"],
                      "wall_time_gate": report["wall_time_gate"],
                      "output_dir": str(args.output_dir)}))
    return 0 if report["qualification_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
