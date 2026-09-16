#!/usr/bin/env python3
"""Verify the production switch against a preserved fold-168 run, without policy overrides."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import time
from types import SimpleNamespace

import fold_168_regression as fold
import mdpde_experiment as experiment


def validate_solves(records, enabled):
    refinements = []
    for row in records:
        refinement = row.get("refinement")
        if row["native_status"] == 0 or not enabled:
            if refinement is not None:
                raise RuntimeError("A native successful or disabled solve was refined.")
            if row["qualification"] != (1 if row["native_status"] == 0 else 0):
                raise RuntimeError("Native qualification changed without refinement.")
        elif refinement is None:
            raise RuntimeError("A failed solve has no refinement evidence.")
        if refinement is None:
            continue
        if not refinement["accepted"] or row["qualification"] != 2:
            raise RuntimeError("The preserved fold has an unexpected refinement rejection.")
        if not (refinement["candidate_equation_evaluations"] <= 128
                and refinement["candidate_residual"] <= 1e-8
                and refinement["reference_residual"] <= 1e-10
                and refinement["reference_stop"] == "fresh-residual"
                and max(refinement["relative_coordinate_difference"]) <= 1e-6
                and refinement["weight_max_difference"] <= 1e-6
                and refinement["floor_masks_equal"]):
            raise RuntimeError("Accepted refinement lacks equation/branch qualification.")
        refinements.append(refinement)
    if enabled and (len(records) != 3192 or len(refinements) != 20):
        raise RuntimeError("Incomplete or changed production solve population.")
    if not enabled and len(records) != 5544:
        raise RuntimeError("Incomplete or changed legacy solve population.")
    return {"shape_calls": len(records), "accepted": len(refinements),
            "candidate_equation_evaluations": sum(r["candidate_equation_evaluations"] for r in refinements),
            "reference_updates": sum(r["reference_updates"] for r in refinements),
            "maximum_candidate_equations": max((r["candidate_equation_evaluations"] for r in refinements), default=0)}


def run(args):
    executable, output = args.executable.resolve(), args.output.resolve()
    if output.exists():
        raise RuntimeError("Use a fresh output directory.")
    baseline = fold.load_baseline(experiment.ROOT / "tests/benchmarks/fold_168_simulation_baseline.json")
    inputs = {"model": args.model.resolve(), "map": args.map.resolve(),
              "manifest": Path(str(args.map.resolve()) + ".simulation.json")}
    hashes = fold.validate_input_hashes(inputs, baseline["input_hashes"])
    before = experiment.provenance(executable, experiment.ROOT)
    output.mkdir(parents=True)
    command = fold.build_command(executable, output / "database.sqlite", inputs["model"], inputs["map"])
    command[command.index("-j") + 1] = str(args.jobs)
    command[command.index("-v") + 1] = str(args.verbosity)
    if args.legacy:
        command += ["--second-stage-failed-only-refinement", "false"]
    # Enabled runs deliberately omit the option to test its production default.
    env = {k: v for k, v in os.environ.items() if not k.startswith("RHBM_TEST_")}
    if args.capture:
        env["RHBM_TEST_SOLVER_CAPTURE_DIR"] = str(output / "solver-failures")
    start = time.perf_counter()
    with (output / "run.log").open("w") as log:
        code = subprocess.run(command, env=env, cwd=output, stdout=log, stderr=subprocess.STDOUT).returncode
    experiment.write(output / "execution.json", {"command": command, "returncode": code,
                     "seconds": time.perf_counter()-start, "capture": args.capture})
    experiment.write(output / "provenance.json", before)
    if code or before != experiment.provenance(executable, experiment.ROOT):
        raise RuntimeError("Execution failed or sources/binaries changed during the run.")
    fold.validate_input_hashes(inputs, hashes)
    actual = experiment.score(output, inputs["model"], inputs["map"])
    experiment.verify(SimpleNamespace(reference=args.reference, runs=[output], output=output / "comparison.json"))
    summary = {"enabled": not args.legacy, "exact_reference_match": True,
               "summary": actual["second_stage_summary"], "certificate": actual["production_fitting"]}
    def recoveries(directory):
        records = [json.loads(line.split("payload=", 1)[1]) for line in (directory / "run.log").read_text().splitlines()
                   if "Second-stage audit: schema=2, payload=" in line]
        return [(r["attempt"], r["recovery"]) for r in records if r.get("recovery", {}).get("attempted")]
    recovery = recoveries(output)
    if recovery:
        if recovery != recoveries(args.reference):
            raise RuntimeError("Recovery trials differ from the preserved reference.")
        summary["recovery_accepted"] = sum(r["accepted"] for _, r in recovery)
    if args.capture:
        records = [json.loads(line) for line in (output / "solver-failures/shape-solves.jsonl").read_text().splitlines()]
        summary["work"] = validate_solves(records, not args.legacy)
    experiment.write(output / "validation.json", summary)
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "model", "map", "output", "reference"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--legacy", action="store_true")
    parser.add_argument("--capture", action="store_true")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--verbosity", type=int, default=3)
    run(parser.parse_args())
