"""Shared MDPDE capture, scoring, provenance, and comparison support."""
from __future__ import annotations
import json
import os
import sqlite3
import subprocess
import time
from pathlib import Path
from types import SimpleNamespace

import fold_168_support as fold
from experiment_io import ROOT, write_json as write
from experiment_provenance import source_and_build_provenance as provenance


def score(directory, model, map_path):
    baseline = fold.load_baseline(ROOT / "tests/benchmarks/fold_168_simulation_baseline.json")
    inputs = {"model": model, "map": map_path, "manifest": Path(str(map_path) + ".simulation.json")}
    hashes = fold.validate_input_hashes(inputs, baseline["input_hashes"])
    manifest = fold.load_simulation_manifest(inputs["manifest"], hashes)
    fold.validate_fixture(manifest, baseline)
    actual = fold.make_empty_actual(hashes)
    actual["generation_record"] = {k: v for k, v in manifest.items() if k != "atoms"}
    actual["atoms"] = fold.pair_atom_truth(fold.read_atom_results(directory / "database.sqlite"), manifest)
    actual["quality_metrics"] = fold.calculate_quality_metrics(actual["atoms"])
    log = (directory / "run.log").read_text()
    actual["second_stage_summary"] = fold.parse_second_stage_summary(log)
    actual["atom_cutoff_summary"] = fold.parse_atom_cutoff_summary(log)
    actual["production_fitting"] = fold.parse_final_state_certificate(log)
    gates = fold.evaluate_gates(baseline, actual)
    write(directory / "actual.json", actual)
    write(directory / "report.json", {"schema_version": 7, "truth_scoring": "complete", **gates,
        "convergence_acceptance": actual["second_stage_summary"]["stop_reason"] == "converged"
        and actual["production_fitting"]["certificate"]["qualified"]
        and actual["production_fitting"]["certificate"]["complete"]
        and actual["production_fitting"]["attempts"] <= 25
        and all(x is not None and 0 <= x < 1e-4 for x in actual["production_fitting"]["certificate"]["operator_nominal_p99"]),
        "stop_reason": actual["second_stage_summary"]["stop_reason"]})
    return actual


def capture(args):
    directory = args.output.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    executable = args.executable.resolve()
    before = provenance(executable, args.source.resolve())
    if args.score_existing:
        score(directory, args.model.resolve(), args.map.resolve())
        write(directory / "provenance.json", before)
        return
    if (directory / "database.sqlite").exists():
        raise RuntimeError("Refusing to overwrite an existing experiment database.")
    baseline = fold.load_baseline(ROOT / "tests/benchmarks/fold_168_simulation_baseline.json")
    fold.validate_input_hashes({"model": args.model, "map": args.map,
        "manifest": Path(str(args.map) + ".simulation.json")}, baseline["input_hashes"])
    command = fold.build_command(executable, directory / "database.sqlite", args.model.resolve(), args.map.resolve())
    # Historical binaries predate the production switch; they are already native-only.
    help_text = subprocess.check_output([str(executable), "potential_analysis", "--help"], text=True)
    if "--second-stage-failed-only-refinement" in help_text:
        command += ["--second-stage-failed-only-refinement", "false"]
    command[command.index("-v") + 1] = str(args.verbosity)
    command[command.index("-j") + 1] = str(args.jobs)
    env = os.environ.copy()
    if not args.no_capture:
        env["RHBM_TEST_SOLVER_CAPTURE_DIR"] = str(directory / "solver-failures")
    else:
        env.pop("RHBM_TEST_SOLVER_CAPTURE_DIR", None)
    start = time.perf_counter()
    with (directory / "run.log").open("w") as log:
        status = subprocess.run(command, env=env, cwd=directory, stdout=log, stderr=subprocess.STDOUT).returncode
    write(directory / "execution.json", {"command": command, "returncode": status,
        "seconds": time.perf_counter() - start, "capture_enabled": not args.no_capture})
    after = provenance(executable, args.source.resolve())
    if before != after:
        raise RuntimeError("Sources or executable changed during capture.")
    write(directory / "provenance.json", before)
    if status:
        raise RuntimeError(f"Fitting exited {status}.")
    actual = score(directory, args.model.resolve(), args.map.resolve())
    if not args.no_capture:
        populations = {}
        for path in (directory / "solver-failures").glob("shape-*.txt.json"):
            metadata = json.loads(path.read_text())
            populations.setdefault(metadata["phase"],set()).update(metadata["indices"])
        required = ["final"]
        if actual["second_stage_summary"]["stop_reason"] == "recovery-failed":
            required.append("recovery-current")
        for phase in required:
            if populations.get(phase) != set(range(len(actual["atoms"]))):
                raise RuntimeError(f"Incomplete {phase} shape capture population.")
        write(directory / "capture-populations.json", {k: sorted(v) for k,v in populations.items()})
    print(json.dumps(actual["second_stage_summary"]))


def verify(args):
    def read(directory):
        with sqlite3.connect(f"file:{directory.resolve() / 'database.sqlite'}?mode=ro", uri=True) as db:
            return db.execute("SELECT * FROM model_atom_local_potential ORDER BY 1,2,3").fetchall()
    reference = read(args.reference)
    result = {str(run): read(run) == reference for run in args.runs}
    expected = json.loads((args.reference / "actual.json").read_text())
    summaries = {}
    for run in args.runs:
        actual = json.loads((run / "actual.json").read_text())
        summaries[str(run)] = all(actual[k] == expected[k] for k in
            ["quality_metrics","second_stage_summary","production_fitting","atom_cutoff_summary"])
    write(args.output, {"reference": str(args.reference), "rows": len(reference),
        "exact_database_records": result, "exact_summary_and_certificate": summaries})
    if not all(result.values()) or not all(summaries.values()):
        raise RuntimeError("Numerical neutrality comparison failed.")
