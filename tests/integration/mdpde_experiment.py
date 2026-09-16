#!/usr/bin/env python3
"""Offline experiment orchestration. Never changes fitting or regression gates."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import subprocess
import time

import fold_168_regression as fold

ROOT = Path(__file__).resolve().parents[2]


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def provenance(executable, source):
    build = executable.parent.parent
    files = [executable, build / "src/librhbm_gem.dylib", build / "CMakeCache.txt"]
    links = {str(p): subprocess.check_output(["otool", "-L", str(p)], text=True)
             for p in files[:2] if p.is_file()}
    dependencies = {Path(line.strip().split(" (")[0]) for output in links.values()
                    for line in output.splitlines()[1:] if line.strip().startswith("/")}
    eigen_metadata = []
    if files[2].is_file():
        for line in files[2].read_text().splitlines():
            if line.startswith("Eigen3_DIR:PATH="):
                eigen_metadata = sorted(Path(line.split("=",1)[1]).glob("*.cmake"))
    source_files = sorted(p for folder in ("src", "include", "tests", "cmake")
                          for p in (source / folder).rglob("*")
                          if p.is_file() and p.suffix in {".cpp", ".hpp", ".py", ".cmake", ".txt", ".json"})
    source_files += [source / "CMakeLists.txt"]
    return {"artifacts": {str(p): fold.sha256_file(p) for p in files if p.is_file()},
            "source_files": {str(p.relative_to(source)): fold.sha256_file(p) for p in source_files},
            "linked_libraries": links,
            "dependency_hashes": {str(p): fold.sha256_file(p) for p in sorted(dependencies) + eigen_metadata if p.is_file()}}


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


def stats(values):
    values = sorted(v for v in values if math.isfinite(v))
    if not values:
        return {"n": 0}
    absolute = sorted(map(abs, values))
    return {"n": len(values), "bias": sum(values)/len(values),
        "rmse": math.sqrt(sum(x*x for x in values)/len(values)),
        "p99": absolute[math.ceil(.99*len(values))-1], "max": absolute[-1]}


def experiment(args):
    output = args.output.resolve()
    if output.exists():
        raise RuntimeError("Use a fresh directory for experiment evidence.")
    output.mkdir(parents=True)
    executable = args.executable.resolve()
    before = provenance(executable,args.source.resolve())
    captures = args.captures.resolve()
    inputs = {str(p):fold.sha256_file(p) for p in sorted(captures.rglob("*")) if p.is_file()}
    if args.operation == "solve":
        command = [str(executable),"solve",str(captures),str(output)]
    else:
        if not args.manifest or not args.map:
            raise RuntimeError("Forward experiment requires --manifest and --map.")
        for p in [args.manifest.resolve(),args.map.resolve()]: inputs[str(p)] = fold.sha256_file(p)
        command = [str(executable),"forward",str(args.manifest.resolve()),str(args.map.resolve()),str(captures),str(output)]
    start = time.perf_counter()
    with (output / "run.log").open("w") as log:
        status = subprocess.run(command,stdout=log,stderr=subprocess.STDOUT).returncode
    write(output / "execution.json",{"command":command,"returncode":status,"seconds":time.perf_counter()-start})
    write(output / "provenance.json",before)
    write(output / "input-hashes.json",inputs)
    if status or before != provenance(executable,args.source.resolve()):
        raise RuntimeError("Experiment failed or source/executable changed during execution.")
    if any(fold.sha256_file(Path(p)) != digest for p,digest in inputs.items()):
        raise RuntimeError("Experiment inputs changed during execution.")


def summarize(args):
    output = args.output; output.mkdir(parents=True, exist_ok=True)
    rows = []
    for file in sorted(args.solver.glob("shape-*.json")):
        data = json.loads(file.read_text())
        for m in data["methods"]:
            rows.append({"fixture": data["fixture"], **{k:v for k,v in m.items() if k != "trace"}})
    if not rows:
        raise RuntimeError("No solver results; refusing an empty experiment summary.")
    write(output / "solver-summary.json", rows)
    with (output / "solver-summary.csv").open("w") as out:
        fields = ["fixture","method","equation_pass","reference_pass","residual_inf","variance",
                  "amplitude","width","condition","floor_count","linear_solves","equation_evaluations","seconds"]
        writer = csv.DictWriter(out, fields, extrasaction="ignore"); writer.writeheader(); writer.writerows(rows)
    groups = {}
    differences = {"kernel": ("generator","estimator"), "interpolation": ("grid","generator"),
        "roundtrip": ("roundtrip","grid"), "fixed_map": ("original","roundtrip"), "total": ("original","estimator")}
    with (args.forward / "samples.csv").open() as inp:
        samples = list(csv.DictReader(inp))
        for row in samples:
            radius = float(row["distance"])
            labels = ["all", "signal" if radius <= 1 else "tail" if 1.2 <= radius <= 2 else "other",
                      "crosses-cutoff" if row["crosses_cutoff"] == "1" else "smooth-stencil",
                      "map-boundary" if row["map_boundary"] == "1" else "map-interior"]
            for metric, (left,right) in differences.items():
                delta = float(row[left])-float(row[right])
                if not math.isfinite(delta):
                    continue
                for label in labels:
                    key = (row["case"],row["h"],metric,label)
                    group = groups.setdefault(key, [[],[]]); group[0].append(delta); group[1].append(delta/float(row["peak"]))
    # Hold membership fixed across grid spacings before interpreting refinement rates.
    memberships = {}
    for row in samples:
        key = (row["case"],row["serial_id"],row["sample"])
        memberships.setdefault(key,set()).add(row["crosses_cutoff"])
    for row in samples:
        membership = memberships[(row["case"],row["serial_id"],row["sample"])]
        if len(membership) != 1:
            continue
        label = "common-crosses-cutoff" if membership == {"1"} else "common-smooth-stencil"
        for metric,(left,right) in differences.items():
            delta = float(row[left])-float(row[right])
            if math.isfinite(delta):
                key = (row["case"],row["h"],metric,label)
                group = groups.setdefault(key,[[],[]]);group[0].append(delta);group[1].append(delta/float(row["peak"]))
    forward = [{"case": k[0], "h": float(k[1]), "metric": k[2], "group": k[3],
        **stats(v[0]), "peak_normalized": stats(v[1])} for k,v in sorted(groups.items())]
    write(output / "forward-summary.json", forward)
    fits = [json.loads(line) for line in (args.forward / "fits.jsonl").read_text().splitlines()]
    write(output / "forward-fit-summary.json", fits)
    files = [p for folder in [args.solver,args.forward] for p in folder.rglob("*") if p.is_file()]
    write(output / "artifact-index.json", {str(p): fold.sha256_file(p) for p in sorted(files)})
    if args.plots:
        plot_results(args.solver, forward, output)


def plot_results(solver, forward, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    files = sorted(solver.glob("shape-maximum-iterations*.json"))
    fig, axes = plt.subplots(max(1,len(files)),1,figsize=(9,3.6*max(1,len(files))),squeeze=False)
    for ax, file in zip(axes[:,0],files):
        data = json.loads(file.read_text())
        for method in data["methods"]:
            name = method["method"]
            if name not in {"production","normal-production-init","qr-production-init","svd-production-init","root-ols","root-endpoint"}:
                continue
            values = [(p.get("equation") or p).get("residual_inf") for p in method.get("trace",[])]
            ax.semilogy(range(1,len(values)+1),[max(v,1e-17) if v is not None else float("nan") for v in values],label=name)
        ax.axhline(1e-8,color="grey",linestyle="--",linewidth=.8,label="equation criterion")
        ax.set(title=file.stem.replace("shape-maximum-iterations-",""),xlabel="Fixed-point iteration / root equation evaluation",ylabel="Fresh scaled residual (max norm)")
        ax.grid(alpha=.2); ax.legend(fontsize=7,ncol=2)
    fig.tight_layout(); fig.savefig(output/"solver-residuals.png",dpi=170); fig.savefig(output/"solver-residuals.pdf"); plt.close(fig)
    cases = sorted({r["case"] for r in forward if r["case"] != "fold"})
    fig, axes = plt.subplots(2,1,figsize=(9,8))
    for case in cases:
        for ax, group in zip(axes,["common-smooth-stencil","common-crosses-cutoff"]):
            points = sorted((r["h"],r["rmse"]) for r in forward if r["case"] == case and r["group"] == group and r["metric"] == "interpolation")
            if points:
                ax.loglog(*zip(*points),marker="o",label=case)
            ax.set(title=group,xlabel="Grid spacing (angstrom)",ylabel="Interpolation RMSE")
            ax.grid(alpha=.2); ax.legend(fontsize=8,ncol=3)
    fig.tight_layout(); fig.savefig(output/"forward-discrepancy.png",dpi=170); fig.savefig(output/"forward-discrepancy.pdf"); plt.close(fig)


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode",required=True)
    p = sub.add_parser("capture")
    for name in ["executable","model","map","output"]: p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--source",type=Path,default=ROOT)
    p.add_argument("--verbosity",type=int,default=3); p.add_argument("--jobs",type=int,default=4)
    p.add_argument("--score-existing",action="store_true"); p.add_argument("--no-capture",action="store_true")
    p = sub.add_parser("summarize")
    for name in ["solver","forward","output"]: p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--plots",action="store_true")
    p = sub.add_parser("verify")
    p.add_argument("--reference",type=Path,required=True); p.add_argument("--runs",type=Path,nargs="+",required=True)
    p.add_argument("--output",type=Path,required=True)
    p = sub.add_parser("experiment")
    for name in ["executable","captures","output"]: p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--operation",choices=["solve","forward"],required=True)
    p.add_argument("--source",type=Path,default=ROOT)
    p.add_argument("--manifest",type=Path); p.add_argument("--map",type=Path)
    args = parser.parse_args(); {"capture":capture,"summarize":summarize,"verify":verify,"experiment":experiment}[args.mode](args)


if __name__ == "__main__":
    main()
