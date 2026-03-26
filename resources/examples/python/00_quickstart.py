#!/usr/bin/env python3
"""Run the Python quickstart on one example model/map pair."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib import error, request

import rhbm_gem_module as rgm


DEFAULT_WORKDIR = Path("example_output/python_quickstart")
MODEL_URL = "https://files.rcsb.org/download/6Z6U.cif"
MAP_URL = (
    "https://drive.usercontent.google.com/download"
    "?id=1MA7jQRZAKc0hGMLom_I2RAPcEqDc6ZSh&confirm=xxx"
)
MODEL_NAME = "6Z6U.cif"
MAP_NAME = "emd_11103_additional.map"
DATABASE_NAME = "database.sqlite"
KEY_TAG = "6z6u"
JOBS = 4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "work_dir",
        nargs="?",
        type=Path,
        help="Working directory for downloads and generated outputs.",
    )
    parser.add_argument(
        "--workdir",
        dest="workdir_option",
        type=Path,
        default=None,
        help="Working directory for downloads and generated outputs.",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Reuse existing files under <workdir>/data instead of downloading them again.",
    )
    args = parser.parse_args()
    if args.work_dir is not None and args.workdir_option is not None:
        parser.error("provide either positional work_dir or --workdir, not both")
    return args


def resolve_workdir(args: argparse.Namespace) -> Path:
    selected = args.workdir_option or args.work_dir or DEFAULT_WORKDIR
    return selected.resolve()


def ensure_execute(report: rgm.ExecutionReport, step_name: str) -> None:
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


def require_file(path: Path) -> None:
    if path.is_file():
        return
    raise FileNotFoundError(
        f"Required input file not found: {path}\n"
        "Prepare the file manually or rerun without --skip-download."
    )


def download_file(url: str, destination: Path) -> None:
    temp_path = destination.with_suffix(destination.suffix + ".tmp")
    req = request.Request(url, headers={"User-Agent": "RHBM-GEM Python quickstart"})
    try:
        with request.urlopen(req) as response, temp_path.open("wb") as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
        temp_path.replace(destination)
    except (error.URLError, error.HTTPError, OSError) as exc:
        if temp_path.exists():
            temp_path.unlink()
        raise RuntimeError(
            f"Failed to download {url} to {destination}: {exc}\n"
            "Place the file under <workdir>/data manually, then rerun with --skip-download."
        ) from exc


def prepare_inputs(data_dir: Path, skip_download: bool) -> tuple[Path, Path]:
    model_path = data_dir / MODEL_NAME
    map_path = data_dir / MAP_NAME

    if skip_download:
        print(f"Skipping downloads; expecting existing files under {data_dir}")
    else:
        print(f"Downloading data files under {data_dir}")
        download_file(MODEL_URL, model_path)
        download_file(MAP_URL, map_path)

    require_file(model_path)
    require_file(map_path)
    return model_path, map_path


def main() -> int:
    args = parse_args()
    workdir = resolve_workdir(args)
    data_dir = workdir / "data"
    database_path = data_dir / DATABASE_NAME

    workdir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    print(f"Working directory: {workdir}")
    print(f"Data directory:    {data_dir}")
    print("")

    print("[1/3] Prepare data")
    model_path, map_path = prepare_inputs(data_dir, args.skip_download)
    print(f"  model: {model_path}")
    print(f"  map:   {map_path}")

    print("[2/3] Potential analysis")
    analysis_request = rgm.PotentialAnalysisRequest()
    analysis_request.common.database_path = str(database_path)
    analysis_request.common.thread_size = JOBS
    analysis_request.model_file_path = str(model_path)
    analysis_request.map_file_path = str(map_path)
    analysis_request.saved_key_tag = KEY_TAG
    analysis_request.training_alpha_flag = False
    ensure_execute(rgm.RunPotentialAnalysis(analysis_request), "RunPotentialAnalysis")
    print(f"  database: {database_path}")
    print(f"  key tag:  {KEY_TAG}")
    print(f"  jobs:     {JOBS}")

    print("[3/3] Result dump")
    dump_request = rgm.ResultDumpRequest()
    dump_request.common.database_path = str(database_path)
    dump_request.common.folder_path = str(workdir)
    dump_request.model_key_tag_list = KEY_TAG
    dump_request.printer_choice = rgm.PrinterType.GAUS_ESTIMATES
    ensure_execute(rgm.RunResultDump(dump_request), "RunResultDump")

    dump_outputs = sorted(path.name for path in workdir.iterdir() if path.is_file())
    print("  dump outputs:")
    for name in dump_outputs:
        print(f"    - {name}")

    print("")
    print("Quickstart completed successfully.")
    print(f"Workdir:      {workdir}")
    print(f"Data dir:     {data_dir}")
    print(f"Database:     {database_path}")
    print(f"Result dumps: {workdir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - entrypoint UX
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
