#!/usr/bin/env python3
"""Determinism checks for the build-gated exposure case runner."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import subprocess


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def load_analyzer():
    path = (PROJECT_ROOT / "resources" / "tools" / "developer" /
            "analyze_counterfactual_convergence.py")
    spec = importlib.util.spec_from_file_location("exposure_runner_analyzer", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_case(executable: Path, seed: int, threads: int) -> tuple[str, dict]:
    completed = subprocess.run([
        str(executable),
        "--case-id", f"determinism-{seed}-{threads}",
        "--family", "natural", "--topology", "c-c",
        "--level", "0", "--replica", "0", "--seed", str(seed),
        "--threads", str(threads), "--noise-sigma", "0.01",
    ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
    text = completed.stdout.decode("utf-8", errors="strict")
    truth = "\n".join(
        line for line in text.splitlines()
        if "Convergence exposure truth:" in line)
    return truth, load_analyzer().parse_log(text)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    args = parser.parse_args()

    truth_a, parsed_a = run_case(args.executable, 410000, 1)
    truth_b, parsed_b = run_case(args.executable, 410000, 1)
    truth_c, _ = run_case(args.executable, 410001, 1)
    _, parsed_parallel = run_case(args.executable, 410000, 4)

    if truth_a != truth_b:
        raise AssertionError("identical seeds did not reproduce scenario truth")
    if truth_a == truth_c:
        raise AssertionError("different seeds did not alter the scenario")
    if parsed_a["trajectory_records"] != parsed_b["trajectory_records"]:
        raise AssertionError("identical seeds did not reproduce the predicate sequence")
    if parsed_a["trajectory_records"] != parsed_parallel["trajectory_records"]:
        raise AssertionError("one-thread and four-thread predicate sequences differ")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
