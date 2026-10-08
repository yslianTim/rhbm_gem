"""Measure repeated fixed-neighbor search performance for qualified outer cores."""
# Internal frontier-phase implementation. Use joint_fixed_neighbor_outer_core_qualification.py.
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
from joint_fixed_neighbor_core_size import (
    MODE,
    POLICY,
    _case,
    _correct,
    _parse_core_sizes,
    _partial,
    _status,
    summarize,
)


DEFAULT_CASES = ("chain-512", "cube-512", "chain-1024", "cube-1024")
DEFAULT_FINALIST_CORE_SIZES = (12, 16)
CONTROL_CORE_SIZE = 64
WARMUP_REPETITIONS = 1
MEASURED_REPETITIONS = 3
MIN_SEARCH_IMPROVEMENT = 0.10
RSS_MATERIAL_REGRESSION_FACTOR = 1.25
TOTAL_MATERIAL_REGRESSION_FACTOR = 1.25
MAX_SWEEPS_WARNING = 24


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _median(values):
    finite = [value for value in values if _finite(value)]
    return statistics.median(finite) if finite else None


def _maximum(values):
    finite = [value for value in values if _finite(value)]
    return max(finite) if finite else None


def _case_list(values):
    cases = [_case(value) for value in values]
    if len(set(cases)) != len(cases):
        raise ValueError("cases must be unique")
    return cases


def _validate_finalists(finalists, control):
    finalists = _parse_core_sizes(finalists)
    if control <= 0 or control in finalists:
        raise ValueError("control core must be positive and distinct from finalists")
    return finalists


def _run_case(args, output_dir, topology, atoms, core_atoms, phase,
              warmup, repetition, execution_index):
    tag = "warmup" if warmup else f"rep{repetition}"
    case = f"{topology}-{atoms}-core{core_atoms}-{tag}"
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
    report = {
        "topology": topology,
        "atoms": atoms,
        "core_atoms": core_atoms,
        "phase": phase,
        "warmup": warmup,
        "repetition": repetition,
        "execution_index": execution_index,
        "policy": POLICY,
        "measurement_scope": "fixed-neighbor-search-only",
        "command": command,
        "process": process,
        "status": _status(process),
        "completed": process.get("status") == "completed" and result is not None,
        "result": result,
    }
    write(wrapper_path, report)
    return report


def _interleaved_order(core_sizes, repetition):
    offset = (repetition - 1) % len(core_sizes)
    return list(core_sizes[offset:]) + list(core_sizes[:offset])


def _run_repeated_case(args, output_dir, topology, atoms, core_sizes, phase, reports):
    for core_atoms in core_sizes:
        reports.append(_run_case(args, output_dir, topology, atoms, core_atoms, phase,
                                  True, 0, len(reports) + 1))
        write(output_dir / "runs.json", reports)
    for repetition in range(1, MEASURED_REPETITIONS + 1):
        for core_atoms in _interleaved_order(tuple(core_sizes), repetition):
            reports.append(_run_case(args, output_dir, topology, atoms, core_atoms, phase,
                                     False, repetition, len(reports) + 1))
            write(output_dir / "runs.json", reports)


def _aggregate(rows, topology, atoms, core_atoms):
    selected = [row for row in rows
                if row.get("topology") == topology and row.get("atoms") == atoms
                and row.get("core_atoms") == core_atoms]
    warmups = [row for row in selected if row.get("warmup") is True]
    measured = [row for row in selected if row.get("warmup") is False]
    repetitions = sorted(row.get("repetition") for row in measured)
    expected = list(range(1, MEASURED_REPETITIONS + 1))
    complete = repetitions == expected and len(measured) == MEASURED_REPETITIONS
    correct = complete and all(_correct(row) for row in measured)
    return {
        "topology": topology,
        "atoms": atoms,
        "core_atoms": core_atoms,
        "warmup_runs": len(warmups),
        "warmup_completed": len(warmups) == WARMUP_REPETITIONS and
        all(row.get("status") == "completed" for row in warmups),
        "measured_runs": len(measured),
        "measured_repetitions": repetitions,
        "complete": complete,
        "all_correct": correct,
        "correct_repetitions": sum(1 for row in measured if _correct(row)),
        "statuses": sorted({row.get("status") for row in selected}),
        "median_search_seconds": _median([row.get("search_seconds") for row in measured]),
        "median_total_seconds": _median([row.get("total_seconds") for row in measured]),
        "median_peak_rss_mb": _median([row.get("peak_rss_mb") for row in measured]),
        "maximum_peak_rss_mb": _maximum([row.get("peak_rss_mb") for row in measured]),
        "median_sweeps": _median([row.get("sweeps") for row in measured]),
        "maximum_sweeps": _maximum([row.get("sweeps") for row in measured]),
        "median_block_solves": _median([row.get("block_solves") for row in measured]),
        "median_local_work_seconds": _median([row.get("local_work_seconds") for row in measured]),
        "maximum_local_rows": _maximum([row.get("maximum_local_rows") for row in measured]),
        "maximum_local_columns": _maximum([row.get("maximum_local_columns") for row in measured]),
        "search_seconds": [row.get("search_seconds") for row in measured],
        "peak_rss_mb": [row.get("peak_rss_mb") for row in measured],
        "sweeps": [row.get("sweeps") for row in measured],
    }


def _compare(candidate, control):
    candidate_search = candidate.get("median_search_seconds")
    control_search = control.get("median_search_seconds")
    improvement = ((control_search - candidate_search) / control_search
                   if _finite(candidate_search) and _finite(control_search) and control_search > 0
                   else None)
    candidate_rss = candidate.get("median_peak_rss_mb")
    control_rss = control.get("median_peak_rss_mb")
    rss_ok = (_finite(candidate_rss) and _finite(control_rss) and
              candidate_rss <= RSS_MATERIAL_REGRESSION_FACTOR * control_rss)
    candidate_total = candidate.get("median_total_seconds")
    control_total = control.get("median_total_seconds")
    total_ok = (_finite(candidate_total) and _finite(control_total) and
                candidate_total <= TOTAL_MATERIAL_REGRESSION_FACTOR * control_total)
    candidate_rows = candidate.get("maximum_local_rows")
    control_rows = control.get("maximum_local_rows")
    candidate_columns = candidate.get("maximum_local_columns")
    control_columns = control.get("maximum_local_columns")
    local_geometry_ok = True
    if _finite(candidate_rows) and _finite(control_rows):
        local_geometry_ok = local_geometry_ok and candidate_rows <= control_rows
    if _finite(candidate_columns) and _finite(control_columns):
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
        "sweeps_warning": (_finite(candidate_sweeps) and candidate_sweeps > MAX_SWEEPS_WARNING),
        "maximum_candidate_sweeps_warning": (_finite(candidate.get("maximum_sweeps")) and
                                               candidate["maximum_sweeps"] > MAX_SWEEPS_WARNING),
    }


def _phase_report(rows, cases, core_sizes):
    report_cases = []
    aggregate_by_core = {}
    for topology, atoms in cases:
        cores = [_aggregate(rows, topology, atoms, core) for core in core_sizes]
        control = next((row for row in cores if row["core_atoms"] == CONTROL_CORE_SIZE), None)
        comparisons = {}
        if control is not None:
            for row in cores:
                if row["core_atoms"] != CONTROL_CORE_SIZE:
                    comparisons[str(row["core_atoms"])] = _compare(row, control)
        report_cases.append({"topology": topology, "atoms": atoms,
                             "cores": cores, "comparisons_vs_control": comparisons})
        for row in cores:
            aggregate_by_core.setdefault(str(row["core_atoms"]), []).append(row)
    aggregate = {}
    for core, entries in aggregate_by_core.items():
        values = [row.get("median_search_seconds") for row in entries]
        aggregate[core] = sum(values) if len(values) == len(cases) and all(_finite(value) for value in values) else None
    return {"cases": report_cases,
            "aggregate_median_search_seconds_by_core": aggregate,
            "core_sizes": list(core_sizes)}


def _phase_core(report, topology, atoms, core_atoms):
    case = next((item for item in report.get("cases", [])
                 if item.get("topology") == topology and item.get("atoms") == atoms), None)
    if case is None:
        return None
    return next((row for row in case.get("cores", [])
                 if row.get("core_atoms") == core_atoms), None)


def _phase_comparison(report, topology, atoms, core_atoms):
    case = next((item for item in report.get("cases", [])
                 if item.get("topology") == topology and item.get("atoms") == atoms), None)
    return case.get("comparisons_vs_control", {}).get(str(core_atoms)) if case else None


def _select_512(report, finalists):
    aggregate = report.get("aggregate_median_search_seconds_by_core", {})
    eligible = []
    for core in finalists:
        rows = [_phase_core(report, topology, 512, core) for topology in ("chain", "cube")]
        if all(row is not None and row.get("all_correct") is True for row in rows):
            value = aggregate.get(str(core))
            if _finite(value):
                eligible.append((value, core))
    return min(eligible)[1] if eligible else None


def analyze(rows, finalist_core_sizes=DEFAULT_FINALIST_CORE_SIZES,
            control_core_size=CONTROL_CORE_SIZE, selected_core_size=None):
    finalists = _validate_finalists(finalist_core_sizes, control_core_size)
    coarse_cases = (("chain", 512), ("cube", 512))
    coarse_report = _phase_report(rows, coarse_cases, (*finalists, control_core_size))
    selected = selected_core_size or _select_512(coarse_report, finalists)
    large_cores = ((selected, control_core_size) if selected is not None else (control_core_size,))
    large_cases = (("chain", 1024), ("cube", 1024))
    large_report = _phase_report(rows, large_cases, large_cores)
    comparisons = []
    for topology, atoms in (*coarse_cases, *large_cases):
        if selected is None:
            continue
        comparison = _phase_comparison(
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
        if not _finite(improvement) or improvement < MIN_SEARCH_IMPROVEMENT:
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
        "policy": POLICY,
        "finalist_core_sizes": list(finalists),
        "control_core_size": control_core_size,
        "warmup_repetitions": WARMUP_REPETITIONS,
        "measured_repetitions": MEASURED_REPETITIONS,
        "interleaving": "cyclic rotation of each core list for each measured repetition",
        "minimum_search_improvement_fraction": MIN_SEARCH_IMPROVEMENT,
        "rss_material_regression_factor": RSS_MATERIAL_REGRESSION_FACTOR,
        "total_material_regression_factor": TOTAL_MATERIAL_REGRESSION_FACTOR,
        "maximum_sweeps_warning": MAX_SWEEPS_WARNING,
        "coarse_512": coarse_report,
        "large_1024": large_report,
        "selected_core_size": selected,
        "comparisons": comparisons,
        "warnings": warnings,
        "qualification_gate": "passed" if not gate_reasons else "failed",
        "gate_reasons": gate_reasons,
        "selection_basis": "fastest aggregate median search time among finalists correct on both 512-atom topologies",
    }


def _all_rows(reports):
    return [summarize(report) | {
        "phase": report.get("phase"),
        "warmup": report.get("warmup"),
        "repetition": report.get("repetition"),
        "execution_index": report.get("execution_index"),
    } for report in reports if report.get("result") is not None]


def write_outputs(report, output_dir):
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


def run_campaign(args):
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = _case_list(args.cases)
    finalists = _validate_finalists(args.finalist_core_sizes, args.control_core_size)
    required = {("chain", 512), ("cube", 512), ("chain", 1024), ("cube", 1024)}
    if set(cases) != required:
        raise ValueError("frontier requires exactly chain/cube 512 and chain/cube 1024 cases")
    reports = []
    coarse_cores = (*finalists, args.control_core_size)
    for topology, atoms in (("chain", 512), ("cube", 512)):
        _run_repeated_case(args, output_dir, topology, atoms, coarse_cores,
                           "coarse-512", reports)
    coarse_rows = _all_rows(reports)
    coarse_analysis = analyze(coarse_rows, finalists, args.control_core_size)
    selected = coarse_analysis["selected_core_size"]
    if selected is not None:
        large_cores = (selected, args.control_core_size)
        for topology, atoms in (("chain", 1024), ("cube", 1024)):
            _run_repeated_case(args, output_dir, topology, atoms, large_cores,
                               "large-1024", reports)
    write(output_dir / "runs.json", reports)
    analysis = analyze(_all_rows(reports), finalists, args.control_core_size, selected)
    manifest = {
        "schema_version": 1,
        "phase": analysis["phase"],
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {
            "driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
            "core_size": sha(ROOT / "tests/integration/joint_fixed_neighbor_core_size.py"),
            "frontier": sha(Path(__file__)),
        },
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "backend": "SPQR",
        "eigen_threads": 1,
        "cases": [f"{topology}-{atoms}" for topology, atoms in cases],
        "policy": POLICY,
        "finalist_core_sizes": list(finalists),
        "control_core_size": args.control_core_size,
        "measurement_scope": "fixed-neighbor-search-only",
        "resource_envelope": {"wall_seconds": args.timeout, "rss_bytes": args.rss_limit},
        "repetitions": {"warmup": WARMUP_REPETITIONS, "measured": MEASURED_REPETITIONS},
        "correctness_contract": ["search_converged", "A/C KKT <= 1e-10",
                                 "width gradient <= 1e-12", "cache/objective replay limits"],
        "analysis": {"qualification_gate": analysis["qualification_gate"],
                     "selected_core_size": analysis["selected_core_size"],
                     "warnings": analysis["warnings"]},
    }
    write(output_dir / "campaign-manifest.json", manifest)
    write_outputs(analysis, output_dir)
    (output_dir / "README.md").write_text(
        "# FixedNeighbor outer-core repeated frontier\n\n"
        "This P4 frontier keeps LegacyCompact, OneAcceptedUpdate, Forward serial "
        "Gauss-Seidel, SPQR, one Eigen thread, and the existing search gates. It uses "
        f"one warmup and {MEASURED_REPETITIONS} interleaved measured repetitions. The "
        f"512-atom cases compare finalists {', '.join(str(core) for core in finalists)} "
        f"against control {args.control_core_size}; the fastest correct finalist by "
        "aggregate median search time is then compared with control on both 1024-atom "
        "cases.\n\n"
        "Gate C4 requires at least 10% median search improvement for the selected core "
        "on chain and cube at both sizes, no material RSS/total-time/local-geometry "
        f"regression, and correctness on every measured run. Sweeps above {MAX_SWEEPS_WARNING} "
        "are reported as warnings.\n\n"
        f"Qualification gate: **{analysis['qualification_gate']}**. Selected core: "
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
    parser.add_argument("--finalist-core-sizes", type=_parse_core_sizes,
                        default=DEFAULT_FINALIST_CORE_SIZES)
    parser.add_argument("--control-core-size", type=int, default=CONTROL_CORE_SIZE)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.build_dir = args.build_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.finalist_core_sizes = _validate_finalists(args.finalist_core_sizes, args.control_core_size)
    if args.timeout <= 0 or args.rss_limit <= 0:
        parser.error("resource limits must be positive")
    if not (args.build_dir / "bin" / "joint_fixed_neighbor_experiment").is_file():
        parser.error("build joint_fixed_neighbor_experiment before running frontier study")
    try:
        analysis = run_campaign(args)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.error(str(error))
    print(json.dumps({"qualification_gate": analysis["qualification_gate"],
                      "selected_core_size": analysis["selected_core_size"],
                      "output_dir": str(args.output_dir)}))
    return 0 if analysis["qualification_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
