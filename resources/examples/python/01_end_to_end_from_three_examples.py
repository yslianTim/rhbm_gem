#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
#
# Script:
#   01_end_to_end_from_three_examples.py
#
# Purpose:
#   Download three example model/map pairs (PDB-6Z6U/EMD-11103, PDB-6DRV/EMD-8908, PDB-9GXM/EMD-51667),
#   run `potential_analysis`, and then export three estimation results with `result_dump`
#   through the Python bindings.
#
# Audience / Platform:
#   New RHBM-GEM Python users on macOS, Linux, or Windows after the bindings are installed.
#
# Prerequisites:
#   Import `rhbm_gem_module` successfully before running this script.
#   Use `--skip-download` for offline reruns after staging files under `<workdir>/data`.
#
# Inputs / Outputs:
#   Download and read inputs from `<workdir>/data`.
#   Write the SQLite database under `<workdir>/data` and generated outputs under `<workdir>`.
#
# References:
#   docs/user/getting-started.md (Python Examples)
#   resources/README.md
#   resources/examples/cli/01_end_to_end_from_three_examples.sh
#   https://www.rcsb.org/structure/6Z6U
#   https://www.rcsb.org/structure/6DRV
#   https://www.rcsb.org/structure/9GXM
"""Run the Python end-to-end workflow on three example model/map pairs."""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib import error, request

import rhbm_gem_module as rgm


DEFAULT_WORKDIR = Path("example_output/python_three_examples")
DATABASE_NAME = "database.sqlite"
JOBS = 4


@dataclass(frozen=True)
class ExampleSpec:
    key_tag: str
    model_name: str
    model_url: str
    map_name: str
    map_url: str
    compressed_map_name: str | None = None


EXAMPLE_SPECS = (
    ExampleSpec(
        key_tag="6z6u",
        model_name="6Z6U.cif",
        model_url="https://files.rcsb.org/download/6Z6U.cif",
        map_name="emd_11103_additional.map",
        map_url=(
            "https://drive.usercontent.google.com/download"
            "?id=1MA7jQRZAKc0hGMLom_I2RAPcEqDc6ZSh&confirm=xxx"
        ),
    ),
    ExampleSpec(
        key_tag="6drv",
        model_name="6DRV.cif",
        model_url="https://files.rcsb.org/download/6DRV.cif",
        map_name="emd_8908_additional.map",
        map_url="https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-8908/other/emd_8908_additional.map.gz",
        compressed_map_name="emd_8908_additional.map.gz",
    ),
    ExampleSpec(
        key_tag="9gxm",
        model_name="9GXM.cif",
        model_url="https://files.rcsb.org/download/9GXM.cif",
        map_name="emd_51667_additional_1.map",
        map_url="https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-51667/other/emd_51667_additional_1.map.gz",
        compressed_map_name="emd_51667_additional_1.map.gz",
    ),
)


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


def ensure_execute(result: rgm.CommandResult, step_name: str) -> None:
    if result.outcome == rgm.CommandOutcome.Succeeded:
        return

    details = "\n".join(
        f"[{issue.phase}/{issue.level}] {issue.option_name}: {issue.message}"
        for issue in result.issues
    )
    if not details:
        details = "(no validation issue details)"
    raise RuntimeError(
        f"{step_name} failed (outcome={result.outcome}).\n{details}"
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
    req = request.Request(url, headers={"User-Agent": "RHBM-GEM Python example"})
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


def extract_gzip(source: Path, destination: Path) -> None:
    temp_path = destination.with_suffix(destination.suffix + ".tmp")
    try:
        with gzip.open(source, "rb") as compressed, temp_path.open("wb") as handle:
            shutil.copyfileobj(compressed, handle)
        temp_path.replace(destination)
    except OSError as exc:
        if temp_path.exists():
            temp_path.unlink()
        raise RuntimeError(f"Failed to extract {source} to {destination}: {exc}") from exc
    finally:
        if source.exists():
            source.unlink()


def prepare_inputs(data_dir: Path, skip_download: bool) -> list[tuple[ExampleSpec, Path, Path]]:
    staged_inputs: list[tuple[ExampleSpec, Path, Path]] = []

    if skip_download:
        print(f"Skipping downloads; expecting existing files under {data_dir}")
    else:
        print(f"Downloading data files under {data_dir}")

    for spec in EXAMPLE_SPECS:
        model_path = data_dir / spec.model_name
        map_path = data_dir / spec.map_name

        if not skip_download:
            download_file(spec.model_url, model_path)
            if spec.compressed_map_name is None:
                download_file(spec.map_url, map_path)
            else:
                compressed_map_path = data_dir / spec.compressed_map_name
                download_file(spec.map_url, compressed_map_path)
                extract_gzip(compressed_map_path, map_path)

        require_file(model_path)
        require_file(map_path)
        staged_inputs.append((spec, model_path, map_path))

    return staged_inputs


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
    staged_inputs = prepare_inputs(data_dir, args.skip_download)
    for spec, model_path, map_path in staged_inputs:
        print(f"  {spec.key_tag}:")
        print(f"    model: {model_path}")
        print(f"    map:   {map_path}")

    print("[2/3] Potential analysis")
    key_tags: list[str] = []
    for spec, model_path, map_path in staged_inputs:
        analysis_request = rgm.PotentialAnalysisRequest()
        analysis_request.database_path = str(database_path)
        analysis_request.job_count = JOBS
        analysis_request.model_file_path = str(model_path)
        analysis_request.map_file_path = str(map_path)
        analysis_request.saved_key_tag = spec.key_tag
        analysis_request.training_alpha_flag = False
        ensure_execute(rgm.RunPotentialAnalysis(analysis_request), "RunPotentialAnalysis")
        key_tags.append(spec.key_tag)
        print(f"  key tag:  {spec.key_tag}")
        print(f"  model:    {model_path.name}")
        print(f"  map:      {map_path.name}")
        print(f"  jobs:     {JOBS}")

    print("[3/3] Result dump")
    dump_request = rgm.ResultDumpRequest()
    dump_request.database_path = str(database_path)
    dump_request.output_dir = str(workdir)
    dump_request.model_key_tag_list = key_tags
    dump_request.printer_choice = rgm.PrinterType.GAUS_ESTIMATES
    ensure_execute(rgm.RunResultDump(dump_request), "RunResultDump")

    dump_outputs = sorted(path.name for path in workdir.iterdir() if path.is_file())
    print(f"  key tags: {dump_request.model_key_tag_list}")
    print("  dump outputs:")
    for name in dump_outputs:
        print(f"    - {name}")

    print("")
    print("Three-example workflow completed successfully.")
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
