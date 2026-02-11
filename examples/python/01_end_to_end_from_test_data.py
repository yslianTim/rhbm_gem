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


def ensure_execute(success: bool, step_name: str) -> None:
    if not success:
        raise RuntimeError(f"{step_name} failed (Execute() returned false)")


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
    simulator = rgm.MapSimulationCommand()
    simulator.SetModelFilePath(str(model_path))
    simulator.SetFolderPath(str(maps_dir))
    simulator.SetBlurringWidthList(args.blurring_width)
    ensure_execute(simulator.Execute(), "MapSimulationCommand")

    map_files = sorted(maps_dir.glob("sim_map_*.map"))
    if not map_files:
        raise RuntimeError(f"No map output generated under: {maps_dir}")
    map_path = map_files[0]
    print(f"  model: {model_path}")
    print(f"  map:   {map_path}")

    print("[2/3] Potential analysis")
    analyzer = rgm.PotentialAnalysisCommand()
    analyzer.SetDatabasePath(str(database_path))
    analyzer.SetModelFilePath(str(model_path))
    analyzer.SetMapFilePath(str(map_path))
    analyzer.SetSavedKeyTag(args.saved_key_tag)
    analyzer.SetSamplingSize(args.sampling_size)
    ensure_execute(analyzer.Execute(), "PotentialAnalysisCommand")
    print(f"  database: {database_path}")

    print("[3/3] Result dump")
    dumper = rgm.ResultDumpCommand()
    dumper.SetDatabasePath(str(database_path))
    dumper.SetFolderPath(str(dump_dir))
    dumper.SetModelKeyTagList(args.saved_key_tag)
    dumper.SetPrinterChoice(rgm.PrinterType.GAUS_ESTIMATES)
    ensure_execute(dumper.Execute(), "ResultDumpCommand")

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
