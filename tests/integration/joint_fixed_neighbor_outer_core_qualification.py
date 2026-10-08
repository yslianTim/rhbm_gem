"""Run one phase of the FixedNeighbor outer-core qualification pipeline."""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from experiment_io import ROOT, read, sha, write
from experiment_process import RSS_LIMIT_BYTES
import joint_fixed_neighbor_outer_core_support as support


PHASES = ("screen", "endpoint", "frontier")


def analyze(phase, rows, *, core_sizes=None, finalist_core_sizes=None,
            control_core_size=support.frontier_CONTROL_CORE_SIZE, selected_core_size=None):
    """Apply the phase-specific correctness gate through one public API."""
    if phase == "screen":
        return support.screen_analyze(rows, core_sizes or support.screen_CORE_SIZES)
    if phase == "endpoint":
        return support.endpoint_analyze(rows, core_sizes or support.endpoint_FINALIST_CORES)
    if phase == "frontier":
        return support.frontier_analyze(
            rows, finalist_core_sizes or support.frontier_DEFAULT_FINALIST_CORE_SIZES,
            control_core_size, selected_core_size)
    raise ValueError(f"unknown qualification phase: {phase}")


def _cases(phase, values):
    defaults = {
        "screen": support.screen_DEFAULT_CASES,
        "endpoint": support.endpoint_DEFAULT_CASES,
        "frontier": support.frontier_DEFAULT_CASES,
    }
    return list(values or defaults[phase])


def _phase_args(args):
    phase = args.phase
    values = _cases(phase, args.cases)
    common = {
        "build_dir": args.build_dir,
        "output_dir": args.output_dir,
        "timeout": args.timeout,
        "rss_limit": args.rss_limit,
        "cases": values,
    }
    if phase == "screen":
        return argparse.Namespace(
            **common,
            core_sizes=args.core_sizes or support.screen_HISTORICAL_CORE_SIZES,
        )
    if phase == "endpoint":
        return argparse.Namespace(
            **common,
            core_sizes=args.core_sizes or support.endpoint_FINALIST_CORES,
        )
    finalists = args.core_sizes or args.finalist_core_sizes
    return argparse.Namespace(
        **common,
        finalist_core_sizes=finalists,
        control_core_size=args.control_core_size,
    )


def _patch_manifest(args, phase):
    path = args.output_dir / "campaign-manifest.json"
    manifest = read(path)
    manifest["qualification_driver"] = {
        "path": str(Path(__file__).relative_to(ROOT)),
        "phase": phase,
        "command_phase": ["--phase", phase],
    }
    manifest.setdefault("benchmark_source_hashes", {})["qualification_driver"] = sha(Path(__file__))
    write(path, manifest)


def run_phase(args):
    phase_args = _phase_args(args)
    runners = {
        "screen": support.screen_run_campaign,
        "endpoint": support.endpoint_run_campaign,
        "frontier": support.frontier_run_campaign,
    }
    analysis = runners[args.phase](phase_args)
    _patch_manifest(args, args.phase)
    return analysis


def _gate(analysis):
    return analysis.get("correctness_gate", analysis.get("qualification_gate"))


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=PHASES, required=True,
                        help="screen, endpoint, or repeated frontier")
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=7200.0)
    parser.add_argument("--rss-limit", type=int, default=RSS_LIMIT_BYTES)
    parser.add_argument("--cases", nargs="+",
                        help="explicit phase-compatible topology-atom cases")
    parser.add_argument("--core-sizes", type=support.screen_parse_core_sizes,
                        help="explicit comma-separated candidate cores")
    parser.add_argument("--finalist-core-sizes", type=support.screen_parse_core_sizes,
                        default=support.frontier_DEFAULT_FINALIST_CORE_SIZES,
                        help="frontier finalists; --core-sizes is an alias")
    parser.add_argument("--control-core-size", type=int,
                        default=support.frontier_CONTROL_CORE_SIZE)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.build_dir = args.build_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    if args.timeout <= 0 or args.rss_limit <= 0:
        parser.error("resource limits must be positive")
    if args.control_core_size <= 0:
        parser.error("control core size must be positive")
    if not (args.build_dir / "bin" / "joint_fixed_neighbor_experiment").is_file():
        parser.error("build joint_fixed_neighbor_experiment before running qualification")
    try:
        phase_args = _phase_args(args)
        validators = {
            "screen": support.screen_case,
            "endpoint": support.endpoint_case,
            "frontier": support.screen_case,
        }
        for value in phase_args.cases:
            validators[args.phase](value)
        analysis = run_phase(args)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.error(str(error))
    print(json.dumps({"phase": args.phase, "gate": _gate(analysis),
                      "output_dir": str(args.output_dir)}))
    return 0 if _gate(analysis) == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
