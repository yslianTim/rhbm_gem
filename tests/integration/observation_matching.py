#!/usr/bin/env python3
"""Reproducible observation-matched experiments; no production policy changes."""
from __future__ import annotations

import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import fold_168_regression as fold
import mdpde_experiment as experiment

ROOT = Path(__file__).resolve().parents[2]
PARAMETERS = ("A", "B", "C")


def read(path):
    return json.loads(Path(path).read_text())


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def run(args):
    output = args.output.resolve()
    require(not output.exists(), "Use a fresh output directory; historical evidence is immutable.")
    baseline = fold.load_baseline(ROOT / "tests/benchmarks/fold_168_simulation_baseline.json")
    inputs = {"model": args.model.resolve(), "map": args.map.resolve(), "manifest": args.manifest.resolve()}
    hashes = fold.validate_input_hashes(inputs, baseline["input_hashes"])
    manifest = fold.load_simulation_manifest(inputs["manifest"], hashes)
    fold.validate_fixture(manifest, baseline)
    captures = args.captures.resolve()
    files = sorted(p for p in captures.rglob("*") if p.is_file())
    require(files, "Missing capture evidence.")
    input_hashes = {str(p): fold.sha256_file(p) for p in list(inputs.values()) + files}
    executable = args.executable.resolve()
    before = experiment.provenance(executable, args.source.resolve())
    output.mkdir(parents=True)
    experiment.write(output / "input-hashes.json", input_hashes)
    experiment.write(output / "provenance.json", before)
    experiment.write(output / "environment.json", {"python": sys.version, "platform": platform.platform(),
        "jobs": args.jobs, "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=args.source, text=True).strip()})
    command = [str(executable), "matched", str(inputs["manifest"]), str(inputs["map"]), str(captures), str(output)]
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.jobs)
    start = time.perf_counter()
    with (output / "run.log").open("w") as log:
        result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
    experiment.write(output / "execution.json", {"command": command, "returncode": result.returncode,
        "seconds": time.perf_counter() - start, "jobs": args.jobs})
    require(before == experiment.provenance(executable, args.source.resolve()), "Source or binary changed during experiment.")
    require(all(fold.sha256_file(Path(p)) == h for p, h in input_hashes.items()), "Input changed during experiment.")
    require(result.returncode == 0, f"Experiment failed; inspect {output / 'run.log'}.")
    summarize(output)


def load_fits(directory):
    index = read(directory / "fit-index.json")
    require(index["targets"] == 183 and index["fits"] == 2928, "Incomplete experiment matrix.")
    rows = [row for i in range(index["targets"]) for row in read(directory / "fits" / f"{i}.json")]
    require(len(rows) == index["fits"], "Missing fits.")
    keys = {(r["case"], r["serial_id"], r["h"], r["layer"], r["prediction"], r["parameters"]) for r in rows}
    require(len(keys) == len(rows), "Duplicate fit identity.")
    for layer in ("analytic", "grid", "header", "file"):
        for prediction in ("analytic", "matched"):
            for parameters in ("AB", "ABC"):
                group = [r for r in rows if r["case"] == "fold" and r["layer"] == layer
                         and r["prediction"] == prediction and r["parameters"] == parameters]
                require({r["serial_id"] for r in group} == set(range(1, 169)), "Incomplete fold population.")
                require(all(r["n"] == (100 if parameters == "AB" else 200) for r in group), "Changed membership.")
    return rows


def fit_statistics(rows):
    groups = collections.defaultdict(list)
    for row in rows:
        groups[(row["case"], row["h"], row["layer"], row["prediction"], row["parameters"])].append(row)
    results = []
    for key, members in sorted(groups.items()):
        summary = dict(zip(("case", "h", "layer", "prediction", "parameters"), key))
        summary.update(total=len(members), qualified=sum(r["qualified"] for r in members),
                       reasons=dict(collections.Counter(r["reason"] for r in members)))
        for i, parameter in enumerate(PARAMETERS):
            values = [r[parameter] - r["truth"][i] for r in members if r[parameter] is not None]
            differences = [r[parameter] - r["reference"][parameter] for r in members
                           if r[parameter] is not None and r["reference"][parameter] is not None]
            summary[parameter] = {**experiment.stats(values), "numerical_difference": experiment.stats(differences),
                                  "fixed": parameter == "C" and key[-1] == "AB"}
        summary["condition"] = experiment.stats(r["condition"] for r in members if r["condition"] is not None)
        summary["evaluations"] = sum(r["evaluations"] for r in members)
        summary["reference_evaluations"] = sum(r["reference_evaluations"] for r in members)
        summary["linear_solves"] = sum(r["linear_solves"] for r in members)
        results.append(summary)
    return results


def paired_results(rows):
    pairs = collections.defaultdict(dict)
    for r in rows:
        pairs[(r["case"], r["h"], r["serial_id"], r["layer"], r["parameters"])][r["prediction"]] = r
    result = []
    for key, pair in sorted(pairs.items()):
        require(set(pair) == {"analytic", "matched"}, "Incomplete prediction pair.")
        a, m = pair["analytic"], pair["matched"]
        require(a["membership"] == m["membership"], "Prediction modes use different samples.")
        entry = dict(zip(("case", "h", "serial_id", "layer", "parameters"), key))
        entry["both_qualified"] = a["qualified"] and m["qualified"]
        for i, p in enumerate(PARAMETERS):
            if p == "C" and key[-1] == "AB":
                continue
            if a[p] is None or m[p] is None:
                entry[p] = {"status": "unavailable"}
                continue
            ea, em = a[p] - a["truth"][i], m[p] - m["truth"][i]
            uncertainty = abs(a[p] - a["reference"][p]) + abs(m[p] - m["reference"][p])
            gain = abs(ea) - abs(em)
            entry[p] = {"analytic_error": ea, "matched_error": em, "absolute_error_gain": gain,
                        "numerical_difference": uncertainty,
                        "status": "improved" if gain > uncertainty else "regressed" if gain < -uncertainty else "unresolved"}
        result.append(entry)
    return result


def integration_decision(statistics):
    evidence = []
    for parameters in ("AB", "ABC"):
        pair = {r["prediction"]: r for r in statistics
                if r["case"] == "fold" and r["layer"] == "file" and r["parameters"] == parameters}
        require(set(pair) == {"analytic", "matched"}, "Missing file-level accuracy comparison.")
        a, m = pair["analytic"], pair["matched"]
        for p in parameters:
            ar, mr = a[p].get("rmse"), m[p].get("rmse")
            uncertainty = a[p]["numerical_difference"].get("rmse", math.inf) + m[p]["numerical_difference"].get("rmse", math.inf)
            qualified = a["qualified"] == m["qualified"] == a["total"] == m["total"] == 168
            available = ar is not None and mr is not None and a[p]["n"] == m[p]["n"] == 168
            evidence.append({"parameters": parameters, "parameter": p, "analytic_rmse": ar, "matched_rmse": mr,
                "ratio": mr/ar if available and ar else None, "numerical_difference": uncertainty if math.isfinite(uncertainty) else None,
                "all_qualified": qualified, "improvement_resolved": available and ar-mr > uncertainty})
    return {"recommend_peeling_experiment": all(e["all_qualified"] and e["improvement_resolved"] for e in evidence),
            "evidence": evidence, "production_quality_gate": "unchanged-uncalibrated"}


def forward_statistics(directory):
    with (directory / "samples.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == 35220 and all(r["passed"] == "1" for r in rows), "Missing or failed forward samples.")
    memberships = collections.defaultdict(set)
    for r in rows:
        for geometry in ("grid", "file"):
            memberships[(r["case"], r["serial_id"], r["sample"], geometry)].add(r["crosses_" + geometry])
    groups = collections.defaultdict(list)
    for r in rows:
        radius = float(r["distance"])
        for metric, left, right in [("interpolation", "grid", "analytic"), ("geometry", "header", "grid"),
                                    ("quantization", "roundtrip", "header"), ("matched-file", "matched_header", "original"),
                                    ("compact-double", "matched_grid", "grid"), ("compact-float", "matched_float", "roundtrip")]:
            geometry = "grid" if metric in ("interpolation", "compact-double", "geometry") else "file"
            labels = ["all", "signal" if radius <= 1 else "tail" if 1.2 <= radius <= 2 else "other",
                      "crosses-cutoff" if r["crosses_" + geometry] == "1" else "smooth-stencil",
                      "map-boundary" if r["boundary_" + geometry] == "1" else "map-interior"]
            membership = memberships[(r["case"], r["serial_id"], r["sample"], geometry)]
            if len(membership) == 1:
                labels.append("common-crosses-cutoff" if membership == {"1"} else "common-smooth-stencil")
            for label in labels:
                groups[(r["case"], float(r["h"]), metric, label)].append(float(r[left])-float(r[right]))
    return [{**dict(zip(("case", "h", "metric", "group"), key)), **experiment.stats(values)}
            for key, values in sorted(groups.items())]


def control_checks(rows):
    controls = [r for r in rows if (r["layer"] == "analytic" and r["prediction"] == "analytic")
                or (r["layer"] in ("grid", "header") and r["prediction"] == "matched")]
    failures = []
    maximum = 0.0
    for r in controls:
        error = max(abs(r[p] - r["truth"][i]) / max(1.0, abs(r["truth"][i]))
                    if r[p] is not None else math.inf for i, p in enumerate(r["parameters"]))
        maximum = max(maximum, error)
        if not r["qualified"] or error > 1e-6:
            failures.append({k: r[k] for k in ("case", "h", "serial_id", "layer", "prediction", "parameters", "reason")})
    return {"count": len(controls), "passed": bool(controls) and not failures, "failures": failures,
            "maximum_scaled_parameter_error": maximum if math.isfinite(maximum) else None, "tolerance": 1e-6}


def summarize(directory):
    require(read(directory / "forward-status.json")["passed"], "Forward stage did not pass.")
    checks = read(directory / "forward-checks.json")
    require(len(checks) == 16 and all(c["passed"] for c in checks), "Incomplete forward cases.")
    rows = load_fits(directory)
    controls = control_checks(rows)
    statistics = fit_statistics(rows)
    pairs = paired_results(rows)
    summary = {"schema_version": 1, "technical_complete": controls["passed"], "controls": controls, "forward": forward_statistics(directory),
               "fit_statistics": statistics, "decision": integration_decision(statistics),
               "fit_count": len(rows), "failure_count": sum(not r["qualified"] for r in rows),
               "forward_checks": checks}
    experiment.write(directory / "results.json", summary)
    require(controls["passed"], "Noiseless self-consistency controls did not recover identifiable parameters.")
    experiment.write(directory / "paired-results.json", pairs)
    with (directory / "fits.csv").open("w") as stream:
        fields = ["case", "serial_id", "h", "layer", "prediction", "parameters", "n", "A", "B", "C", "loss",
                  "qualified", "reason", "stationarity", "condition", "reference_difference", "evaluations", "reference_evaluations", "seconds"]
        writer = csv.DictWriter(stream, fields, extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)
    files = sorted(p for p in directory.rglob("*") if p.is_file() and p.name != "artifact-index.json")
    experiment.write(directory / "artifact-index.json", {str(p.relative_to(directory)): fold.sha256_file(p) for p in files})
    print(json.dumps(summary["decision"], indent=2))


def compare(left, right, output):
    a, b = load_fits(left), load_fits(right)
    require(len(a) == len(b), "Different fit populations.")
    differences = []
    for i, (x, y) in enumerate(zip(a, b)):
        x = {k: v for k, v in x.items() if k != "seconds"}
        y = {k: v for k, v in y.items() if k != "seconds"}
        if x != y:
            differences.append(i)
    sample_equal = (left / "samples.csv").read_bytes() == (right / "samples.csv").read_bytes()
    checks_equal = read(left / "forward-checks.json") == read(right / "forward-checks.json")
    result = {"passed": not differences and sample_equal and checks_equal,
              "fit_differences": differences, "samples_equal": sample_equal, "checks_equal": checks_equal}
    experiment.write(output, result)
    require(result["passed"], "j1/j4 evidence differs.")


def plots(directory, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    result = read(directory / "results.json")
    output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.5), constrained_layout=True)
    for axis, p in zip(axes, PARAMETERS):
        rows = [r for r in result["fit_statistics"] if r["case"] == "fold" and r["layer"] == "file" and r["parameters"] == "ABC"]
        axis.bar([r["prediction"] for r in rows], [r[p]["rmse"] for r in rows], color=["#b56b45", "#246b86"])
        axis.set_yscale("log"); axis.set_title(f"{p} RMSE"); axis.grid(axis="y", alpha=.2)
    fig.suptitle("Fold-168 · fixed truth neighbors · raw least squares")
    for extension in ("png", "pdf"):
        fig.savefig(output / f"parameter-errors.{extension}", dpi=180)
    plt.close(fig)
    fig, axis = plt.subplots(figsize=(7, 4), constrained_layout=True)
    data = [r for r in result["forward"] if r["case"] == "fold" and r["group"] == "all"
            and r["metric"] in ("interpolation", "geometry", "quantization", "matched-file", "compact-double")]
    axis.bar([r["metric"] for r in data], [r["rmse"] for r in data], color="#246b86")
    axis.set_yscale("log"); axis.set_ylabel("Response RMSE"); axis.tick_params(axis="x", rotation=20)
    axis.set_title("Forward discrepancy decomposition"); axis.grid(axis="y", alpha=.2)
    for extension in ("png", "pdf"):
        fig.savefig(output / f"forward-errors.{extension}", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    run_parser = sub.add_parser("run")
    for name in ("executable", "model", "map", "manifest", "captures", "output"):
        run_parser.add_argument("--" + name, type=Path, required=True)
    run_parser.add_argument("--source", type=Path, default=ROOT)
    run_parser.add_argument("--jobs", type=int, choices=(1, 4), default=4)
    summarize_parser = sub.add_parser("summarize")
    summarize_parser.add_argument("directory", type=Path)
    compare_parser = sub.add_parser("compare")
    for name in ("left", "right", "output"):
        compare_parser.add_argument("--" + name, type=Path, required=True)
    plot_parser = sub.add_parser("plots")
    plot_parser.add_argument("directory", type=Path)
    plot_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.operation == "run": run(args)
    elif args.operation == "summarize": summarize(args.directory)
    elif args.operation == "compare": compare(args.left, args.right, args.output)
    else: plots(args.directory, args.output)


if __name__ == "__main__":
    main()
