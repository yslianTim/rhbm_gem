"""Compare the P1 route frontier with the qualified FixedNeighbor benchmark policy."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from experiment_io import ROOT, read, sha, write
from experiment_provenance import source_hash
from joint_route_frontier import compare_group, record_from_report


DEFAULT_CASES = tuple(
    f"{topology}-{atoms}"
    for atoms in (512, 768, 1024)
    for topology in ("chain", "cube")
)
BASELINE_ROUTES = ("LegacyCompact", "OperatorPcg", "FixedNeighbor")
OPTIMIZED_ROUTE = "FixedNeighborOptimized"
ROUTES = (*BASELINE_ROUTES, OPTIMIZED_ROUTE)


def _case(value):
    topology, size = value.split("-", 1)
    atoms = int(size)
    if topology not in ("chain", "cube") or atoms not in (512, 768, 1024):
        raise ValueError(f"Unsupported optimized frontier case: {value}")
    return topology, atoms


def _baseline_record(baseline_dir, case, route, scope="search-only"):
    topology, atoms = _case(case)
    individual = baseline_dir / "individual-results" / scope
    formal_path = individual / f"{case}-{route}-formal.json"
    if not formal_path.is_file():
        raise ValueError(f"Missing P1 formal report: {formal_path}")
    formal = record_from_report(read(formal_path), topology, atoms, route, scope, 600., formal_path)
    record = {
        "topology": topology,
        "atoms": atoms,
        "rows": formal.get("rows"),
        "parameter_count": formal.get("parameter_count"),
        "route": route,
        "measurement_scope": scope,
        "input_sha256": formal.get("input_sha256"),
        "formal_result": formal["result"],
        "diagnostic_result": None,
        "source": "P1 baseline",
        "report_path": formal["report_path"],
    }
    if formal["result"].get("completed"):
        record["diagnostic_result"] = dict(formal["result"], source="formal-envelope")
        return record
    diagnostic_path = individual / f"{case}-{route}-diagnostic.json"
    if diagnostic_path.is_file():
        diagnostic = record_from_report(read(diagnostic_path), topology, atoms, route, scope,
                                        7200., diagnostic_path)
        record["diagnostic_result"] = dict(diagnostic["result"], source="diagnostic-extended-envelope")
        record["diagnostic_report_path"] = diagnostic["report_path"]
    else:
        record["diagnostic_result"] = {
            "status": "not-run",
            "completed": False,
            "reason": "P1 extended diagnostic was deferred",
            "source": "P1 baseline",
        }
    return record


def _candidate_command(args, case, timeout, output):
    return [
        sys.executable, str(ROOT / "tests/integration/joint_benchmark.py"),
        "--profile", "search", "--case", case,
        "--build-dir", str(args.build_dir), "--output", str(output),
        "--timeout", str(timeout), "--rss-limit", str(args.rss_limit),
        "--preconditioner", "fixed-neighbor", "--repeat", "1", "--warmup", "0",
        "--operator-rank", "auto", "--operator-rank-seconds", "120",
        "--operator-rank-work-entries", "100000000", "--operator-rank-workspace-mib", "256",
        "--schwarz-core-atoms", "128", "--schwarz-overlap-hops", "1",
        "--schwarz-max-block-atoms", "512", "--schwarz-storage-mib", "512",
        "--schwarz-scratch-mib", "256", "--fixed-local-work", args.fixed_local_work,
        "--fixed-core-atoms", str(args.fixed_core_atoms),
    ]


def _run_candidate(args, case, timeout, report_path):
    if report_path.is_file():
        report = read(report_path)
    else:
        command = _candidate_command(args, case, timeout, report_path)
        print(f"Running {case} {OPTIMIZED_ROUTE} with {timeout:g}s cap", flush=True)
        completed = subprocess.run(command, cwd=ROOT, check=False)
        if not report_path.is_file():
            raise RuntimeError(f"Benchmark runner did not write {report_path} (exit {completed.returncode})")
        report = read(report_path)
    topology, atoms = _case(case)
    return record_from_report(report, topology, atoms, OPTIMIZED_ROUTE, "search-only",
                              timeout, report_path)


def _candidate_record(args, output_dir, case):
    individual = output_dir / "individual-results" / "search-only"
    formal_path = individual / f"{case}-{OPTIMIZED_ROUTE}-formal.json"
    formal = _run_candidate(args, case, args.formal_timeout, formal_path)
    record = {
        "topology": formal["topology"],
        "atoms": formal["atoms"],
        "rows": formal.get("rows"),
        "parameter_count": formal.get("parameter_count"),
        "route": OPTIMIZED_ROUTE,
        "measurement_scope": "search-only",
        "input_sha256": formal.get("input_sha256"),
        "formal_result": formal["result"],
        "diagnostic_result": None,
        "source": "P6 optimized candidate",
        "report_path": formal["report_path"],
    }
    if formal["result"].get("completed"):
        record["diagnostic_result"] = dict(formal["result"], source="formal-envelope")
        return record
    if args.skip_diagnostics:
        record["diagnostic_result"] = {
            "status": "not-run",
            "completed": False,
            "reason": "P6 extended diagnostic was deferred",
            "source": "deferred-by-policy",
        }
        return record
    diagnostic_path = individual / f"{case}-{OPTIMIZED_ROUTE}-diagnostic.json"
    diagnostic = _run_candidate(args, case, args.diagnostic_timeout, diagnostic_path)
    if record.get("rows") is None:
        record["rows"] = diagnostic.get("rows")
    if record.get("parameter_count") is None:
        record["parameter_count"] = diagnostic.get("parameter_count")
    record["diagnostic_result"] = dict(diagnostic["result"], source="diagnostic-extended-envelope")
    record["diagnostic_report_path"] = diagnostic["report_path"]
    return record


def _classify(groups):
    wins = {}
    resource_advantage = False
    for group in groups:
        winner = group.get("fastest_route")
        if winner:
            wins[winner] = wins.get(winner, 0) + 1
        formal = {row["route"]: row.get("formal_envelope") or {}
                  for row in group.get("routes", [])}
        optimized_pass = formal.get(OPTIMIZED_ROUTE, {}).get("completed") is True
        competing_fail = all(formal.get(route, {}).get("completed") is not True
                             for route in ("LegacyCompact", "OperatorPcg"))
        resource_advantage |= optimized_pass and competing_fail
    completed_groups = sum(1 for group in groups if group.get("fastest_route"))
    if resource_advantage:
        return "FixedNeighborOptimized-memory-advantaged"
    if wins.get(OPTIMIZED_ROUTE, 0) > completed_groups / 2:
        return "FixedNeighborOptimized-time-dominant"
    if wins.get("OperatorPcg", 0) >= wins.get("LegacyCompact", 0) and wins.get("OperatorPcg", 0):
        return "OperatorPcg-dominant"
    if wins:
        return "mixed-by-topology"
    return "no-measured-benefit"


def _analyze(records):
    grouped = {}
    for record in records:
        grouped.setdefault((record["measurement_scope"], record["topology"], record["atoms"]), []).append(record)
    groups = [compare_group(sorted(rows, key=lambda row: ROUTES.index(row["route"])))
              for _, rows in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][2], item[0][1]))]
    first_failures = {}
    for route in ROUTES:
        for topology in ("chain", "cube"):
            rows = sorted((row for row in records if row["route"] == route and
                           row["topology"] == topology), key=lambda item: item["atoms"])
            first = next((row for row in rows if not (row.get("formal_result") or {}).get("completed")), None)
            first_failures[f"{topology}:{route}"] = first["atoms"] if first else None
    return {
        "phase": "P6 optimized FixedNeighbor matched route frontier",
        "route_position": _classify(groups),
        "promotion_decision": "deferred",
        "promotion_reason": ("OperatorPcg remained the fastest route and completed inside the formal "
                              "envelope at all six matched points. FixedNeighborOptimized used less RSS, "
                              "but no formal resource-survival case required a promotion to fallback."),
        "groups": groups,
        "first_formal_search_failure_atoms": first_failures,
    }


def _write_outputs(output_dir, manifest, records):
    output_dir.mkdir(parents=True, exist_ok=True)
    analysis = _analyze(records)
    write(output_dir / "campaign-manifest.json", manifest)
    write(output_dir / "analysis.json", analysis)
    columns = [
        "topology", "atoms", "rows", "parameter_count", "route", "source", "measurement_scope",
        "formal_status", "formal_pass", "formal_search_seconds", "formal_total_seconds",
        "formal_peak_rss_mb", "formal_failure_stage", "formal_failure_reason",
        "diagnostic_status", "diagnostic_completed", "diagnostic_search_seconds",
        "diagnostic_total_seconds", "diagnostic_peak_rss_mb", "diagnostic_reason",
        "objective", "factorization_count", "factor_seconds", "profile_evaluations",
        "block_solves", "sweeps",
    ]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for record in records:
            formal = record.get("formal_result") or {}
            diagnostic = record.get("diagnostic_result") or {}
            writer.writerow({
                "topology": record["topology"], "atoms": record["atoms"],
                "rows": record.get("rows"), "parameter_count": record.get("parameter_count"),
                "route": record["route"], "source": record.get("source"),
                "measurement_scope": record["measurement_scope"],
                "formal_status": formal.get("status"), "formal_pass": formal.get("completed"),
                "formal_search_seconds": formal.get("search_seconds"),
                "formal_total_seconds": formal.get("total_seconds"),
                "formal_peak_rss_mb": formal.get("peak_rss_mb"),
                "formal_failure_stage": formal.get("failure_stage"),
                "formal_failure_reason": formal.get("failure_reason"),
                "diagnostic_status": diagnostic.get("status"),
                "diagnostic_completed": diagnostic.get("completed"),
                "diagnostic_search_seconds": diagnostic.get("search_seconds"),
                "diagnostic_total_seconds": diagnostic.get("total_seconds"),
                "diagnostic_peak_rss_mb": diagnostic.get("peak_rss_mb"),
                "diagnostic_reason": diagnostic.get("reason"),
                "objective": diagnostic.get("objective"),
                "factorization_count": diagnostic.get("factorization_count"),
                "factor_seconds": diagnostic.get("factor_seconds"),
                "profile_evaluations": diagnostic.get("profile_evaluations"),
                "block_solves": diagnostic.get("block_solves"),
                "sweeps": diagnostic.get("sweeps"),
            })
    (output_dir / "README.md").write_text(
        "# Optimized FixedNeighbor matched route frontier\n\n"
        "This P6 campaign reuses the P1 LegacyCompact, OperatorPcg, and original Full/128-core "
        "FixedNeighbor reports and runs the qualified OneAccepted/64-core FixedNeighbor policy "
        "on the same chain and cube cases. All measurements are search-only with one Eigen thread, "
        "the SPQR backend, and a 4 GiB RSS limit. The P1 records are linked by the manifest; they "
        "are not rerun in this campaign.\n\n"
        f"Formal envelope: {manifest['formal_envelope']['wall_seconds']} seconds. "
        f"Diagnostic envelope: {manifest['diagnostic_envelope']['wall_seconds']} seconds.\n\n"
        f"Route position: **{analysis['route_position']}**.\n"
        "\nPromotion gate: **deferred**. OperatorPcg remained fastest and passed the formal "
        "envelope at every matched point; the optimized FixedNeighbor route had lower RSS, "
        "but this campaign did not show a formal resource-survival case that warrants an "
        "internal-route promotion.\n"
    )
    return analysis


def run_campaign(args):
    output_dir = args.output_dir.resolve()
    baseline_dir = args.baseline_dir.resolve()
    baseline_manifest_path = baseline_dir / "campaign-manifest.json"
    if not baseline_manifest_path.is_file():
        raise ValueError(f"Missing P1 campaign manifest: {baseline_manifest_path}")
    baseline_manifest = read(baseline_manifest_path)
    records = []
    manifest = {
        "schema_version": 1,
        "phase": "P6 optimized FixedNeighbor matched route frontier",
        "source_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        "source_sha256": source_hash(ROOT),
        "benchmark_source_hashes": {
            "optimized_frontier": sha(Path(__file__)),
            "shared_runner": sha(ROOT / "tests/integration/joint_benchmark.py"),
            "fixed_neighbor_driver": sha(ROOT / "tests/experiments/joint_fixed_neighbor.cpp"),
        },
        "baseline_artifact": str(baseline_dir.relative_to(ROOT)) if baseline_dir.is_relative_to(ROOT)
        else str(baseline_dir),
        "baseline_source_revision": baseline_manifest.get("source_revision"),
        "branch": subprocess.run(["git", "branch", "--show-current"], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        "backend": "SPQR", "eigen_threads": 1,
        "cases": list(args.cases), "routes": list(ROUTES),
        "measurement_scope": "search-only",
        "formal_envelope": {"wall_seconds": args.formal_timeout, "rss_bytes": args.rss_limit},
        "diagnostic_envelope": {"wall_seconds": args.diagnostic_timeout, "rss_bytes": args.rss_limit},
        "diagnostic_retry_policy": "deferred" if args.skip_diagnostics else "retry every formal non-pass",
        "optimized_policy": {"local_work": "OneAcceptedUpdate", "core_atoms": args.fixed_core_atoms,
                             "block_order": "forward"},
        "promotion_gate": {"decision": "deferred",
                            "reason": "OperatorPcg is fastest and formally completes at all six matched points; optimized FixedNeighbor's RSS advantage does not yet demonstrate resource survival."},
        "route_policy": {
            "LegacyCompact": "P1 baseline global LegacyCompact search",
            "OperatorPcg": "P1 baseline global OperatorPcg with Schwarz",
            "FixedNeighbor": "P1 baseline Full FixedNeighbor with 128 atom core",
            "FixedNeighborOptimized": "P6 OneAcceptedUpdate FixedNeighbor with 64 atom core",
        },
    }
    cache = (args.build_dir / "CMakeCache.txt").read_text()
    if "RHBM_GEM_JOINT_SPARSE_BACKEND:STRING=SPQR" not in cache:
        raise ValueError("Optimized frontier requires the SPQR benchmark build")
    for case in args.cases:
        for route in BASELINE_ROUTES:
            records.append(_baseline_record(baseline_dir, case, route))
        records.append(_candidate_record(args, output_dir, case))
        _write_outputs(output_dir, manifest, records)
    manifest["analysis"] = _analyze(records)
    _write_outputs(output_dir, manifest, records)
    return manifest["analysis"]


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--baseline-dir", type=Path,
                        default=ROOT / "docs/developer/figures/joint-fixed-neighbor-route-frontier-r1")
    parser.add_argument("--formal-timeout", type=float, default=600.)
    parser.add_argument("--diagnostic-timeout", type=float, default=7200.)
    parser.add_argument("--rss-limit", type=int, default=4 * 1024**3)
    parser.add_argument("--skip-diagnostics", action="store_true")
    parser.add_argument("--fixed-local-work", choices=("one",), default="one")
    parser.add_argument("--fixed-core-atoms", type=int, default=64)
    parser.add_argument("--cases", nargs="+", default=DEFAULT_CASES)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.build_dir = args.build_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.baseline_dir = args.baseline_dir.resolve()
    if args.formal_timeout <= 0 or args.diagnostic_timeout < args.formal_timeout or args.rss_limit <= 0:
        parser.error("resource limits must be positive and diagnostic timeout must cover formal timeout")
    if args.fixed_core_atoms <= 0:
        parser.error("fixed core atoms must be positive")
    for case in args.cases:
        try:
            _case(case)
        except (ValueError, IndexError) as error:
            parser.error(str(error))
    try:
        analysis = run_campaign(args)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.error(str(error))
    print(json.dumps({"route_position": analysis["route_position"],
                      "output_dir": str(args.output_dir)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
