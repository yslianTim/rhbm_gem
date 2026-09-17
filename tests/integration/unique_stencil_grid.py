#!/usr/bin/env python3
"""Unique stencil voxels, fixed checkpoint B, common-alpha joint MDPDE (testing only)."""
from __future__ import annotations

import argparse
from array import array
import csv
import math
import os
from pathlib import Path
import platform
import resource
import struct
import subprocess
import sys
import time

import estimated_neighbor_sweep as sweep
import fold_168_regression as fold
import matched_joint_ac as ac
import mdpde_experiment as experiment
from observation_matching import read, require

ROOT = Path(__file__).resolve().parents[2]
ALPHAS = [0, .1, .5, 1]
BUDGET = 100
ROWS = 602995
HISTORY = ROOT / "docs/developer/figures/matched-joint-ac/fits.csv"


def run(args):
    output = args.output.resolve()
    require(not output.exists(), "Use a fresh output directory; historical evidence is immutable.")
    inputs = {k: getattr(args, k).resolve() for k in ("model", "map", "manifest")}
    baseline = fold.load_baseline(ROOT / "tests/benchmarks/fold_168_simulation_baseline.json")
    hashes = fold.validate_input_hashes(inputs, baseline["input_hashes"])
    manifest = fold.load_simulation_manifest(inputs["manifest"], hashes)
    fold.validate_fixture(manifest, baseline)
    index, evidence = sweep.prepare_index(args, manifest, hashes)
    ids = [s["id"] for s in index["states"]]
    require(args.state is None or args.state in ids, "Unknown selected state.")
    index["selected_state_ids"] = [args.state] if args.state else ids
    provenance = experiment.provenance(args.executable.resolve(), args.source.resolve())
    input_hashes = {str(p): fold.sha256_file(p) for p in list(inputs.values()) + evidence + [HISTORY]}
    output.mkdir(parents=True)
    for name, value in (("state-index", index), ("input-hashes", input_hashes), ("provenance", provenance)):
        experiment.write(output / f"{name}.json", value)
    experiment.write(output / "scoring-truth.json", {str(a["serial_id"]):
        [float(a["element"]), manifest["settings"]["blurring_width"], a["charge_used"]] for a in manifest["atoms"]})
    (output / "historical-fits.csv").write_bytes(HISTORY.read_bytes())
    experiment.write(output / "environment.json", {"python": sys.version, "platform": platform.platform(),
        "jobs": 1, "iteration_budget": BUDGET, "refinement_budget": BUDGET,
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=args.source, text=True).strip()})
    command = [str(args.executable.resolve()), "unique-stencil-grid", str(inputs["manifest"]), str(inputs["map"]),
               str(output / "state-index.json"), str(output)]
    env = dict(os.environ, OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    start = time.perf_counter()
    with (output / "run.log").open("w") as log:
        status = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    # ru_maxrss covers child processes of this fresh runner. Record the platform unit explicitly.
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    experiment.write(output / "execution.json", {"command": command, "returncode": status,
        "seconds": time.perf_counter()-start, "jobs": 1, "peak_child_rss_bytes": peak if sys.platform == "darwin" else peak*1024})
    require(provenance == experiment.provenance(args.executable.resolve(), args.source.resolve()), "Source or binary changed during run.")
    require(all(fold.sha256_file(Path(p)) == h for p, h in input_hashes.items()), "Experiment input changed.")
    require(status == 0, f"Unique-grid run failed; inspect {output / 'run.log'}.")
    summarize(output)


def validate_dataset(directory, dataset):
    require(dataset.get("schema_version") == 1 and dataset.get("experiment") == "unique-stencil-grid", "Unsupported dataset.")
    require(len(dataset["atoms"]) == 168 and len(dataset["samples"]) == 33600 and dataset["row_count"] == ROWS
            and dataset["slot_count"] == 33600*64 and dataset["alphas"] == ALPHAS, "Incomplete grid population.")
    require(fold.sha256_file(directory / "voxels.csv") == dataset["voxel_table_sha256"]
            and fold.sha256_file(directory / "slots.csv") == dataset["slot_table_sha256"], "Dataset table hash mismatch.")
    indices, values, multiplicities = array("Q"), array("d"), array("I")
    nx, ny, nz = dataset["grid_size"]
    with (directory / "voxels.csv").open() as stream:
        for p, row in enumerate(csv.DictReader(stream)):
            i = int(row["index"])
            require(int(row["row"]) == p and 0 <= i < nx*ny*nz and (not indices or i > indices[-1]), "Duplicate or invalid voxel index.")
            xyz = [i % nx, i // nx % ny, i // (nx*ny)]
            expected = [xyz[k]*dataset["generation_spacing"][k]+dataset["generation_origin"][k] for k in range(3)]
            # C++ may fuse index*spacing+origin; Python rounds the product first.
            require(all(abs(float(row[key])-expected[k]) <= 16*sys.float_info.epsilon*
                        max(1, abs(dataset["generation_origin"][k]), abs(xyz[k]*dataset["generation_spacing"][k]))
                        for k, key in enumerate(("x", "y", "z"))), "Generation voxel geometry mismatch.")
            observed, reference, bound = (float(row[k]) for k in ("observed", "reference_double", "quantization_bound"))
            require(all(sweep.finite(v) for v in (observed, reference, bound)) and int(row["multiplicity"]) > 0, "Invalid voxel value.")
            require(struct.unpack("f", struct.pack("f", reference))[0] == observed and abs(observed-reference) == bound,
                    "Voxel quantization check failed.")
            indices.append(i); values.append(observed); multiplicities.append(int(row["multiplicity"]))
    require(len(indices) == ROWS, "Missing voxels.")
    counts = array("I", [0])*ROWS
    origin, spacing, size = dataset["sampling_origin"], dataset["sampling_spacing"], dataset["grid_size"]
    sampled = array("d", [0])*33600
    total = 0
    with (directory / "slots.csv").open() as stream:
        for total, row in enumerate(csv.DictReader(stream), 1):
            p, k, v = (int(row[key]) for key in ("sample", "slot", "voxel_row"))
            require((p, k) == divmod(total-1, 64) and p < 33600 and 0 <= v < ROWS, "Incomplete slot sequence.")
            if k == 0:
                position = dataset["samples"][p]["position"]
                outside = any(position[d] < origin[d]-.5*spacing[d] or position[d] > origin[d]+(size[d]-.5)*spacing[d] for d in range(3))
                cell = [-1]*3 if outside else [math.floor((position[d]-origin[d])/spacing[d]) for d in range(3)]
                coefficients = []
                for d in range(3):
                    t = (position[d]-origin[d]-cell[d]*spacing[d])/spacing[d]
                    t2, t3 = t*t, t*t*t
                    coefficients.append([-.5*t3+t2-.5*t, 1.5*t3-2.5*t2+1, -1.5*t3+2*t2+.5*t, .5*t3-.5*t2])
            offsets = [k//16, k//4 % 4, k % 4]
            xyz = [min(size[d]-1, max(0, cell[d]+offsets[d]-1)) for d in range(3)]
            require(indices[v] == xyz[0]+nx*(xyz[1]+ny*xyz[2]), "Stencil union membership mismatch.")
            coefficient = float(row["coefficient"])
            expected = coefficients[0][offsets[0]]*coefficients[1][offsets[1]]*coefficients[2][offsets[2]]
            require(abs(coefficient-expected) <= 64*sys.float_info.epsilon*max(1, abs(expected), *map(abs, cell)),
                    "Interpolation coefficient mismatch.")
            sampled[p] += coefficient*values[v]; counts[v] += 1
    require(total == 33600*64 and counts == multiplicities, "Slot multiplicity mismatch.")
    require(all(abs(y-s["response"]) <= 512*sys.float_info.epsilon*max(1, abs(s["response"]))
                for y, s in zip(sampled, dataset["samples"])), "Original sample replay mismatch.")
    return values


def validate_fit(fit, context, alpha, rows=ROWS):
    require(fit.get("schema_version") == 2 and fit.get("experiment") == "unique-stencil-grid", "Unsupported fit schema.")
    require(fit["rows"] == rows and fit["columns"] == 2*len(context["state"]) and fit["alpha"] == alpha, "Fit population mismatch.")
    require(fit["iteration_budget"] == fit["refinement_budget"] == BUDGET, "Changed iteration budget.")
    require(fit.get("linear_solver") == "sparse-qr" and fit.get("sparse_row_reduction") == "householder-qr-1024"
            and 0 < fit.get("design_nonzeros", 0) <= rows*fit["columns"],
            "Missing exact sparse design evidence.")
    require(fit["blocks"] == [{"scope": "global", "alpha": alpha, "rows": rows, "lambda": 1}], "Invalid global block.")
    initial = [value for a in context["state"] for value in (a[0], a[2])]
    require(fit["initial_beta"] == initial and len(fit["abc"]) == len(context["state"]), "Changed initialization or assembly.")
    for k, (abc, saved) in enumerate(zip(fit["abc"], context["state"])):
        require(struct.pack("d", abc[1]) == struct.pack("d", saved[1]), "Width changed.")
        require([abc[0], abc[2]] == fit["beta"][2*k:2*k+2] and all(sweep.finite(v) for v in abc) and abc[0] >= 0,
                "Infeasible or inconsistent assembly.")
    branches = fit["branches"]
    require(not branches or [b["seed"] for b in branches] == ["checkpoint", "constrained-ls"], "Missing initialization branch.")
    for b in branches:
        for phase, trace in (("primary", "trace"), ("reference", "reference_trace")):
            require(0 <= b[phase]["iterations"] <= BUDGET and len(b[trace]) <= BUDGET+1, "Iteration budget exceeded.")
        if b["qualified"]:
            require(b["primary"]["stop"] == b["reference"]["stop"] == "stationary"
                    and b["primary"]["stationarity"] <= 1e-8 and b["reference"]["stationarity"] <= 1e-10
                    and b["svd_reason"] == "solved" and b["coefficient_difference"] <= 1e-6
                    and b["log_variance_difference"] <= 1e-6 and b["weighted_spectrum"]["rank"] == fit["columns"],
                    "Missing qualified numerical evidence.")
    if fit["qualified"]:
        require(branches and branches[fit["selected_seed"]]["qualified"] and fit["reason"] == "qualified"
                and len(fit["uncertainty"]) == fit["columns"] and all(sweep.finite(v) and v >= 0 for v in fit["uncertainty"]),
                "Unverified selected fit.")
    if branches:
        require(isinstance(fit["selected_seed"], int) and 0 <= fit["selected_seed"] < len(branches), "Invalid selected branch.")
        endpoint = branches[fit["selected_seed"]]["primary"]
        require(fit["beta"] == endpoint["beta"] and fit["variances"] == endpoint["variances"], "Selected endpoint mismatch.")
        if not fit["qualified"]:
            require(fit["reason"] == ("reference-unverified" if endpoint["stop"] == "stationary" else endpoint["stop"]),
                    "Unqualified endpoint classification mismatch.")


def residual_statistics(directory, name, values, samples, expected, fit=None):
    output = {}
    for suffix, observations, label, key in (("", values, "row", "grid_rmse"),
                                            ("-samples", [s["response"] for s in samples], "sample", "matched_sample_rmse")):
        ss, count = 0.0, 0
        with (directory / "residuals" / f"{name}{suffix}.csv").open() as stream:
            for count, row in enumerate(csv.DictReader(stream), 1):
                require(int(row[label]) == count-1 and count <= len(observations), "Residual row mismatch.")
                prediction, residual = float(row["prediction"]), float(row["residual"])
                require(sweep.finite(prediction) and sweep.finite(residual)
                        and observations[count-1]-prediction == residual, "Residual replay mismatch.")
                if not suffix and fit and row["weight"]:
                    variance = fit["variances"][0]
                    # Eigen's y-X*beta and separately materialized X*beta can
                    # round differently. Propagate a scaled residual bound
                    # through exp(-alpha*r^2/(2*v)); tiny v amplifies this gap.
                    delta = 64*sys.float_info.epsilon*max(1, abs(observations[count-1]), abs(prediction))
                    low = math.exp(-.5*fit["alpha"]*((abs(residual)+delta)**2/variance))
                    high = math.exp(-.5*fit["alpha"]*(max(0, abs(residual)-delta)**2/variance))
                    require(low*(1-8*sys.float_info.epsilon) <= float(row["weight"]) <= high*(1+8*sys.float_info.epsilon),
                            "Endpoint weight replay mismatch.")
                ss += residual*residual
        require(count == len(observations), "Missing residual rows.")
        output[key] = math.sqrt(ss/count)
        require(math.isclose(output[key], expected[key], rel_tol=1e-11, abs_tol=1e-14), "Residual aggregate mismatch.")
    return output


def summarize(directory):
    directory = directory.resolve()
    index, states, dataset = (read(directory / f"{n}.json") for n in ("fit-index", "state-index", "dataset"))
    selected = states["selected_state_ids"]
    require(index.get("schema_version") == 1 and index["states"] == states["states"] and len(index["states"]) == 8,
            "State index mismatch.")
    require(selected and len(set(selected)) == len(selected) and set(selected) <= {s["id"] for s in index["states"]}
            and index["selected_state_ids"] == selected, "Selected state mismatch.")
    expected = [f"{s['id']}-alpha-{k}" for s in states["states"] if s["id"] in selected for k in range(4)]
    require(index["cases"] == expected and index["joint_fits"] == len(expected) and index["alphas"] == ALPHAS
            and index["iteration_budget"] == index["refinement_budget"] == BUDGET, "Incomplete experiment matrix.")
    forward = read(directory / "forward-status.json")
    require(forward["passed"] and forward["voxel_count"] == ROWS and forward["sample_count"] == 33600, "Missing forward verification.")
    values = validate_dataset(directory, dataset)
    truth = read(directory / "scoring-truth.json")
    fit_rows, estimates, statistics, inputs = [], [], [], []
    for state in states["states"]:
        if state["id"] not in selected: continue
        sid = state["id"]; context = read(Path(state["context"]))
        require(fold.sha256_file(Path(state["context"])) == state["sha256"], "Checkpoint hash changed.")
        require([a["identity"] for a in context["atoms"]] == dataset["atoms"], "Atom identity mismatch.")
        input_row = read(directory / f"{sid}-input.json")
        require(input_row["abc"] == context["state"], "Input state changed.")
        residual_statistics(directory, sid+"-input", values, dataset["samples"], input_row)
        inputs.append({"state_id": sid, **input_row})
        for k, alpha in enumerate(ALPHAS):
            name = f"{sid}-alpha-{k}"; fit = read(directory / "fits" / f"{name}.json")
            require(fit["state_id"] == sid, "Fit state mismatch.")
            validate_fit(fit, context, alpha)
            residual_statistics(directory, name, values, dataset["samples"], fit["residuals"], fit)
            chosen = fit["branches"][fit["selected_seed"]] if fit["branches"] else {}
            diagnostic = chosen.get("endpoint_diagnostics", {})
            fit_rows.append({"state_id": sid, "alpha": alpha, "qualified": fit["qualified"], "reason": fit["reason"],
                "branch_sensitive": fit.get("branch_sensitive", False), "stationarity": fit.get("stationarity"),
                "variance": fit.get("variances", [None])[0], "objective": fit.get("objective"),
                "primary_iterations": chosen.get("primary", {}).get("iterations"),
                "reference_iterations": chosen.get("reference", {}).get("iterations"),
                "condition": diagnostic.get("condition"), "rank": diagnostic.get("rank"),
                "minimum_singular": diagnostic.get("minimum_singular"), "effective_n": diagnostic.get("effective_n"),
                "maximum_row_share": diagnostic.get("maximum_row_share"), "top_one_percent_share": diagnostic.get("top_one_percent_share"),
                "seconds": fit.get("seconds"), **fit["residuals"]})
            errors = {p: [] for p in ("A", "B", "C")}; input_errors = {p: [] for p in errors}
            for a, identity in enumerate(dataset["atoms"]):
                t = truth[str(identity["serial_id"])]; abc = fit["abc"][a]; saved = context["state"][a]
                row = {"state_id": sid, "alpha": alpha, "serial_id": identity["serial_id"], "qualified": fit["qualified"]}
                for q, p in enumerate(errors):
                    error, before = abc[q]-t[q], saved[q]-t[q]
                    errors[p].append(error); input_errors[p].append(before)
                    row.update({p: abc[q], f"truth_{p}": t[q], f"input_{p}": saved[q], f"error_{p}": error})
                    if p != "B":
                        uncertainty = fit["uncertainty"][2*a+(p == "C")] if fit["qualified"] else None
                        row[f"pair_{p}"] = ac.pair_status(saved[q], abc[q], t[q], 0, uncertainty, fit["qualified"]) if fit["qualified"] else "unqualified"
                estimates.append(row)
            for p in errors:
                for mode, err in (("input", input_errors[p]), ("joint", errors[p])):
                    statistics.append({"state_id": sid, "alpha": alpha, "parameter": p, "mode": mode,
                        "qualified": fit["qualified"] if mode == "joint" else None,
                        "below_0_01": sum(abs(e) < .01 for e in err), **sweep.population_stats(err, 168)})
    result = {"schema_version": 1, "experiment": "unique-stencil-grid", "states": states["states"], "selected_state_ids": selected,
        "alphas": ALPHAS, "row_count": ROWS, "iteration_budget": BUDGET, "refinement_budget": BUDGET,
        "fits": fit_rows, "statistics": statistics, "inputs": inputs, "forward": forward,
        "historical_comparison": "descriptive only: observation space, alpha/variance structure, and iteration budgets differ"}
    experiment.write(directory / "results.json", result)
    for name, rows in (("fits", fit_rows), ("estimates", estimates), ("statistics", statistics)): ac.csv_write(directory / f"{name}.csv", rows)
    artifact_index(directory)
    return result


def artifact_index(directory):
    experiment.write(directory / "artifact-index.json", {"files": {str(p.relative_to(directory)): fold.sha256_file(p)
        for p in sorted(directory.rglob("*")) if p.is_file() and p.name != "artifact-index.json"}})


def scientific(value):
    if isinstance(value, dict): return {k: scientific(v) for k, v in value.items() if k not in ("seconds", "peak_child_rss_bytes")}
    if isinstance(value, list): return [scientific(v) for v in value]
    return value


def compare(left, right, output, state):
    li, ri = read(left / "fit-index.json"), read(right / "fit-index.json")
    require(state in li["selected_state_ids"] and state in ri["selected_state_ids"], "State was not executed in both runs.")
    require(li["states"] == ri["states"], "Source states differ.")
    names = ["dataset.json", "forward-status.json", "input-hashes.json", "provenance.json", f"{state}-input.json"]
    cases = [f"{state}-alpha-{k}" for k in range(4)]
    require(all(c in li["cases"] and c in ri["cases"] for c in cases), "Missing alpha case.")
    names += [f"fits/{c}.json" for c in cases]
    differences = [n for n in names if scientific(read(left/n)) != scientific(read(right/n))]
    tables = ["voxels.csv", "slots.csv"]+[f"residuals/{c}{suffix}.csv" for c in [state+"-input"]+cases for suffix in ("", "-samples")]
    differences += [n for n in tables if fold.sha256_file(left/n) != fold.sha256_file(right/n)]
    for key in ("fits", "statistics", "inputs"):
        a, b = ([r for r in read(p/"results.json")[key] if r["state_id"] == state] for p in (left, right))
        if scientific(a) != scientific(b): differences.append("results/"+key)
    experiment.write(output, {"passed": not differences, "state": state, "compared_fits": 4, "differences": differences,
        "not_compared_states": [s["id"] for s in li["states"] if s["id"] != state], "left": str(left.resolve()), "right": str(right.resolve())})
    require(not differences, f"Reproducibility comparison failed: {differences}")


def plots(directory, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    result = read(directory / "results.json"); output.mkdir(parents=True, exist_ok=True)
    ids = result["selected_state_ids"]
    labels = [s.replace("baseline-", "base ").replace("failed-only-", "refined ") for s in ids]
    colors = ["#606c80", "#187c8c", "#b46b24", "#964568"]
    def save(fig, name):
        for ext in ("png", "pdf"): fig.savefig(output / f"{name}.{ext}", dpi=180)
        plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, p in zip(axes, ("A", "C")):
        for k, alpha in enumerate(ALPHAS):
            rows = [next(r for r in result["statistics"] if r["state_id"] == s and r["alpha"] == alpha and r["parameter"] == p and r["mode"] == "joint") for s in ids]
            ax.plot(range(len(ids)), [r["rmse"] for r in rows], "o-", color=colors[k], label=f"alpha={alpha}")
            for i, r in enumerate(rows):
                if not r["qualified"]: ax.plot(i, r["rmse"], "x", color="black", markersize=9)
        before = [next(r["rmse"] for r in result["statistics"] if r["state_id"] == s and r["parameter"] == p and r["mode"] == "input") for s in ids]
        ax.plot(range(len(ids)), before, "--", color="black", label="Checkpoint")
        ax.set_title(f"{p} absolute RMSE"); ax.set_xticks(range(len(ids)), labels, rotation=55, ha="right"); ax.grid(alpha=.2); ax.set_ylim(bottom=0)
    axes[0].legend(); fig.suptitle("Unique stencil grid | fixed B | x = unqualified endpoint"); save(fig, "parameter-errors")
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), constrained_layout=True)
    for ax, key, title in zip(axes, ("stationarity", "condition", "effective_n"), ("Stationarity", "Weighted condition", "Effective sample size")):
        for k, alpha in enumerate(ALPHAS):
            rows = [next(r for r in result["fits"] if r["state_id"] == s and r["alpha"] == alpha) for s in ids]
            ax.plot(range(len(ids)), [r[key] if r[key] is not None else float("nan") for r in rows], "o-", color=colors[k], label=f"alpha={alpha}")
        if key == "stationarity": ax.axhline(1e-8, ls="--", color="black", label="Primary threshold"); ax.set_yscale("log")
        ax.set_title(title); ax.set_xticks(range(len(ids)), labels, rotation=55, ha="right"); ax.grid(alpha=.2)
    axes[0].legend(); fig.suptitle("Unique stencil grid | numerical diagnostics"); save(fig, "solver-diagnostics")


def main():
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest="operation", required=True)
    runner = sub.add_parser("run")
    for name in ("executable", "model", "map", "manifest", "baseline-run", "refined-run", "output"):
        runner.add_argument("--"+name, type=Path, required=True)
    runner.add_argument("--source", type=Path, default=ROOT); runner.add_argument("--state")
    summary = sub.add_parser("summarize"); summary.add_argument("directory", type=Path)
    comparison = sub.add_parser("compare")
    for name in ("left", "right", "output"): comparison.add_argument("--"+name, type=Path, required=True)
    comparison.add_argument("--state", required=True)
    plot = sub.add_parser("plots"); plot.add_argument("directory", type=Path); plot.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.operation == "run": run(args)
    elif args.operation == "summarize": summarize(args.directory)
    elif args.operation == "compare": compare(args.left, args.right, args.output, args.state)
    else: plots(args.directory, args.output)


if __name__ == "__main__": main()
