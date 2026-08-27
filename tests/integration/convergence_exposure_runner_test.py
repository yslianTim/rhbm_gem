#!/usr/bin/env python3
"""Determinism checks for the build-gated exposure case runner."""

from __future__ import annotations

import argparse
import importlib.util
import math
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


def run_case(
    executable: Path,
    seed: int,
    threads: int,
    topology: str = "c-c",
    noise_sigma: float = 0.01,
) -> tuple[str, dict]:
    completed = subprocess.run([
        str(executable),
        "--case-id", f"determinism-{seed}-{threads}",
        "--family", "natural", "--topology", topology,
        "--level", "0", "--replica", "0", "--seed", str(seed),
        "--threads", str(threads), "--noise-sigma", str(noise_sigma),
    ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
    text = completed.stdout.decode("utf-8", errors="strict")
    truth = "\n".join(
        line.split(", amplitude=", 1)[1]
        for line in text.splitlines()
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
    if truth_a != truth_c:
        raise AssertionError("natural reference truth changed with the noise seed")
    if parsed_a["trajectory_records"] != parsed_b["trajectory_records"]:
        raise AssertionError("identical seeds did not reproduce the predicate sequence")
    if parsed_a["trajectory_records"] != parsed_parallel["trajectory_records"]:
        raise AssertionError("one-thread and four-thread predicate sequences differ")
    if parsed_a["shadow_policy_checkpoints"] != parsed_b["shadow_policy_checkpoints"]:
        raise AssertionError("identical seeds did not reproduce shadow policies")
    if (parsed_a["shadow_policy_checkpoints"] !=
            parsed_parallel["shadow_policy_checkpoints"]):
        raise AssertionError("one-thread and four-thread shadow policies differ")
    if parsed_a["audit_terminal"] != parsed_b["audit_terminal"]:
        raise AssertionError("identical seeds did not reproduce terminal state")

    for topology in ("unk-c", "unk-n", "unk-o"):
        _, parsed = run_case(
            args.executable, 410000, 1, topology=topology, noise_sigma=0.0)
        production = next(
            (record for record in parsed["checkpoints"]
             if record["policy"] == "production"),
            None)
        if production is None:
            raise AssertionError(
                f"{topology} positive control did not reach production convergence")
        objective = [float(value) for value in production["objective"].split("/")]
        if not all(math.isfinite(value) for value in objective):
            raise AssertionError(f"{topology} produced a non-finite objective")
        production_try = production["try"]
        trajectory = next(
            record for record in parsed["trajectory_records"]
            if record["try"] == production_try)
        blockers = [
            int(float(value)) for value in trajectory["blockers"].split("/")]
        if any(blockers):
            raise AssertionError(
                f"{topology} production checkpoint retained safety blockers")
        for atom in parsed["atoms"][(production["experiment"], "production")]:
            parameters = [
                float(atom[name]) for name in ("amplitude", "width", "offset")]
            if not all(math.isfinite(value) for value in parameters):
                raise AssertionError(f"{topology} produced non-finite parameters")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
