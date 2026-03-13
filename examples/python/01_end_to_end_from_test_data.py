#!/usr/bin/env python3
"""Run a minimal Python pipeline: simulate -> analyze -> dump."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import rhbm_gem_module as rgm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model",
        type=Path,
        default=None,
        help="Input model path (.cif/.pdb). Default: auto-detect bundled sample model.",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path("example_output/python_pipeline"),
        help="Working directory for map/db/dump outputs.",
    )
    parser.add_argument(
        "--database",
        type=Path,
        default=None,
        help="SQLite database path. Default: <workdir>/demo.sqlite",
    )
    parser.add_argument(
        "--saved-key-tag",
        default="demo_model",
        help="Key tag used by PotentialAnalysisCommand for saving model object.",
    )
    parser.add_argument(
        "--blurring-width",
        default="1.50",
        help="Comma-separated blurring widths for map simulation.",
    )
    parser.add_argument(
        "--sampling-size",
        type=int,
        default=200,
        help="Sampling size for potential analysis (smaller value is faster for onboarding).",
    )
    return parser.parse_args()


def find_default_model(script_path: Path) -> Path:
    source_tree_candidate = script_path.parents[2] / "tests" / "data" / "test_model.cif"
    installed_candidate = script_path.parent / "data" / "test_model.cif"
    for candidate in (source_tree_candidate, installed_candidate):
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "Could not find sample model. Checked:\n"
        f"  - {source_tree_candidate}\n"
        f"  - {installed_candidate}\n"
        "Pass --model to specify your own input file."
    )


def ensure_execute(report, step_name: str) -> None:
    if report.prepared and report.executed:
        return

    details = "\n".join(
        f"[{issue.phase}/{issue.level}] {issue.option_name}: {issue.message}"
        for issue in report.validation_issues
    )
    if not details:
        details = "(no validation issue details)"
    raise RuntimeError(
        f"{step_name} failed (prepared={report.prepared}, executed={report.executed}).\n{details}"
    )


def main() -> int:
    args = parse_args()
    script_path = Path(__file__).resolve()
    model_path = args.model.resolve() if args.model else find_default_model(script_path)

    workdir = args.workdir.resolve()
    maps_dir = workdir / "maps"
    dump_dir = workdir / "dump"
    maps_dir.mkdir(parents=True, exist_ok=True)
    dump_dir.mkdir(parents=True, exist_ok=True)

    database_path = args.database.resolve() if args.database else (workdir / "demo.sqlite")

    print("[1/3] Map simulation")
    simulation_request = rgm.MapSimulationRequest()
    simulation_request.model_file_path = str(model_path)
    simulation_request.common.folder_path = str(maps_dir)
    simulation_request.blurring_width_list = args.blurring_width
    ensure_execute(rgm.RunMapSimulation(simulation_request), "RunMapSimulation")

    map_files = sorted(maps_dir.glob("sim_map_*.map"))
    if not map_files:
        raise RuntimeError(f"No map output generated under: {maps_dir}")
    map_path = map_files[0]
    print(f"  model: {model_path}")
    print(f"  map:   {map_path}")

    print("[2/3] Potential analysis")
    analysis_request = rgm.PotentialAnalysisRequest()
    analysis_request.common.database_path = str(database_path)
    analysis_request.model_file_path = str(model_path)
    analysis_request.map_file_path = str(map_path)
    analysis_request.saved_key_tag = args.saved_key_tag
    analysis_request.sampling_size = args.sampling_size
    ensure_execute(rgm.RunPotentialAnalysis(analysis_request), "RunPotentialAnalysis")
    print(f"  database: {database_path}")

    print("[3/3] Result dump")
    dump_request = rgm.ResultDumpRequest()
    dump_request.common.database_path = str(database_path)
    dump_request.common.folder_path = str(dump_dir)
    dump_request.model_key_tag_list = args.saved_key_tag
    dump_request.printer_choice = rgm.PrinterType.GAUS_ESTIMATES
    ensure_execute(rgm.RunResultDump(dump_request), "RunResultDump")

    dump_outputs = sorted(path.name for path in dump_dir.glob("*"))
    print("  dump outputs:")
    for name in dump_outputs:
        print(f"    - {name}")

    print("")
    print("Pipeline completed successfully.")
    print(f"Workdir:  {workdir}")
    print(f"Maps dir: {maps_dir}")
    print(f"Dump dir: {dump_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - entrypoint UX
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
