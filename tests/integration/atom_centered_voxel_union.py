#!/usr/bin/env python3
"""Experiment B: exact atom-centered voxel union, fixed B, global MDPDE."""
from __future__ import annotations

import argparse
from array import array
import csv
import math
import os
from pathlib import Path
import platform
import resource
import signal
import struct
import subprocess
import sys
import time

import estimated_neighbor_sweep as sweep
import fold_168_regression as fold
import matched_joint_ac as ac
import mdpde_experiment as experiment
import unique_stencil_grid as grid
from observation_matching import read, require

ROOT = Path(__file__).resolve().parents[2]
NAME = "atom-centered-voxel-union"
ROWS = 3768656
ALPHAS = [0, .1, .5, 1]
RADIUS = 2.5


def expected_cases(index):
    require(len(index["selected_state_ids"]) == 1 and index["alphas"] and
            index["alphas"] == [a for a in ALPHAS if a in index["alphas"]], "Invalid requested cases.")
    return [f"{index['selected_state_ids'][0]}-alpha-{ALPHAS.index(a)}" for a in index["alphas"]]


def run(args):
    output = args.output.resolve()
    require(not output.exists(), "Use a fresh output directory.")
    inputs = {k: getattr(args, k).resolve() for k in ("model", "map", "manifest")}
    baseline = fold.load_baseline(ROOT / "tests/benchmarks/fold_168_simulation_baseline.json")
    hashes = fold.validate_input_hashes(inputs, baseline["input_hashes"])
    manifest = fold.load_simulation_manifest(inputs["manifest"], hashes)
    fold.validate_fixture(manifest, baseline)
    index, evidence = sweep.prepare_index(args, manifest, hashes)
    require(args.state in [s["id"] for s in index["states"]], "Unknown checkpoint.")
    index.update(selected_state_ids=[args.state], alphas=ALPHAS if args.alpha is None else [args.alpha],
                 reference_a=str(args.reference_a.resolve()))
    reference_files = [args.reference_a.resolve()/"fits"/f"{args.state}-alpha-{k}.json" for k in range(4)]
    require(all(p.is_file() for p in reference_files), "Missing four paired Experiment A fits.")
    source = args.source.resolve(); executable = args.executable.resolve()
    provenance = experiment.provenance(executable, source)
    input_hashes = {str(p): fold.sha256_file(p) for p in list(inputs.values())+evidence+reference_files}
    output.mkdir(parents=True)
    for name, value in (("state-index", index), ("input-hashes", input_hashes), ("provenance", provenance)):
        experiment.write(output/f"{name}.json", value)
    experiment.write(output/"scoring-truth.json", {str(a["serial_id"]):
        [float(a["element"]), manifest["settings"]["blurring_width"], a["charge_used"]] for a in manifest["atoms"]})
    experiment.write(output/"environment.json", {"python": sys.version, "platform": platform.platform(), "jobs": 1,
        "iteration_budget": 100, "refinement_budget": 100, "radius": RADIUS,
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()})
    command = [str(executable), NAME, str(inputs["manifest"]), str(inputs["map"]), str(output/"state-index.json"), str(output)]
    env = dict(os.environ, OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    start = time.perf_counter(); interrupted = False
    with (output/"run.log").open("w") as log:
        child = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        try:
            status = child.wait()
        except KeyboardInterrupt:
            interrupted = True; child.send_signal(signal.SIGINT); status = child.wait()
        finally:
            peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
            experiment.write(output/"execution.json", {"command": command, "returncode": child.returncode,
                "seconds": time.perf_counter()-start, "jobs": 1, "user_interrupted": interrupted,
                "peak_child_rss_bytes": peak if sys.platform == "darwin" else peak*1024})
    require(provenance == experiment.provenance(executable, source), "Source or binary changed during run.")
    require(all(fold.sha256_file(Path(p)) == h for p, h in input_hashes.items()), "Experiment input changed.")
    if status != 0:
        completed = read(output/"completion.json").get("cases", []) if (output/"completion.json").exists() else []
        experiment.write(output/"partial-status.json", {"complete": False, "completed": completed,
            "missing": [c for c in expected_cases(index) if c not in completed],
            "reason": "user-interrupted" if interrupted else "technical-failure"})
    require(status == 0, f"Incomplete execution; retained completed cases in {output}.")
    summarize(output)


def sphere_membership(dataset):
    import numpy as np
    size = np.array(dataset["grid_size"]); origin = np.array(dataset["generation_origin"])
    spacing = np.array(dataset["generation_spacing"]); total = int(np.prod(size)); radius = dataset["radius"]
    coverage = np.zeros(total, dtype=np.uint16); nearest = np.full(total, np.inf)
    for atom in dataset["atoms"]:
        center = np.array(atom["position"])
        lo = np.maximum(np.floor((center-radius-origin)/spacing).astype(int), 0)
        hi = np.minimum(np.floor((center+radius-origin)/spacing).astype(int), size-1)
        xx, yy, zz = [np.arange(lo[k], hi[k]+1) for k in range(3)]
        dx, dy, dz = [origin[k]+v*spacing[k]-center[k] for k, v in enumerate((xx, yy, zz))]
        squared = dx[None, None, :]**2+dy[None, :, None]**2+dz[:, None, None]**2
        use = squared <= radius*radius
        ids = (xx[None, None, :]+size[0]*(yy[None, :, None]+size[1]*zz[:, None, None]))[use]
        coverage[ids] += 1; nearest[ids] = np.minimum(nearest[ids], squared[use])
    indices = np.flatnonzero(coverage)
    return indices, coverage[indices], np.sqrt(nearest[indices])


def validate_dataset(directory, dataset, expected_rows=ROWS):
    import numpy as np
    require(dataset.get("experiment") == NAME and dataset.get("radius") == RADIUS
            and dataset.get("membership_geometry") == "generation", "Invalid union geometry.")
    require(dataset["row_count"] == expected_rows and dataset["alphas"] == ALPHAS, "Incomplete voxel population.")
    indices, coverage, nearest = sphere_membership(dataset)
    require(len(indices) == expected_rows, "Independent sphere union population differs.")
    require(fold.sha256_file(directory/"voxels.csv") == dataset["voxel_table_sha256"] and
            fold.sha256_file(directory/"slots.csv") == dataset["slot_table_sha256"], "Dataset hash mismatch.")
    values = array("d"); flags = np.zeros(expected_rows, dtype=bool)
    nx, ny, nz = dataset["grid_size"]; origin = dataset["generation_origin"]; spacing = dataset["generation_spacing"]
    with (directory/"voxels.csv").open() as stream:
        for p, row in enumerate(csv.DictReader(stream)):
            require(p < expected_rows and int(row["row"]) == p and int(row["index"]) == indices[p], "Duplicate or incorrect voxel index.")
            require(int(row["multiplicity"]) == coverage[p] and
                    abs(float(row["nearest_distance"])-nearest[p]) <= 2e-13, "Sphere membership diagnostic differs.")
            i = int(indices[p]); xyz = [i % nx, i//nx % ny, i//(nx*ny)]
            for k, key in enumerate(("x", "y", "z")):
                expected = origin[k]+xyz[k]*spacing[k]
                require(abs(float(row[key])-expected) <= 16*sys.float_info.epsilon*max(1, abs(origin[k]), abs(xyz[k]*spacing[k])), "Voxel geometry differs.")
            y, reference, bound = (float(row[k]) for k in ("observed", "reference_double", "quantization_bound"))
            require(all(math.isfinite(v) for v in (y, reference, bound)) and
                    struct.unpack("f", struct.pack("f", reference))[0] == y and abs(y-reference) == bound, "Quantization differs.")
            values.append(y); require(row["in_stencil"] in ("0", "1"), "Invalid stencil flag."); flags[p] = row["in_stencil"] == "1"
    require(len(values) == expected_rows, "Missing voxels.")
    seen = np.zeros(expected_rows, dtype=bool); sampled = array("d", [0])*len(dataset["samples"])
    origin, spacing = dataset["sampling_origin"], dataset["sampling_spacing"]; size = dataset["grid_size"]; total = 0
    with (directory/"slots.csv").open() as stream:
        for total, row in enumerate(csv.DictReader(stream), 1):
            p, k, v = (int(row[key]) for key in ("sample", "slot", "voxel_row"))
            require((p, k) == divmod(total-1, 64) and p < len(sampled) and 0 <= v < expected_rows, "Incomplete slot mapping.")
            if k == 0:
                point = dataset["samples"][p]["position"]
                outside = any(point[d] < origin[d]-.5*spacing[d] or point[d] > origin[d]+(size[d]-.5)*spacing[d] for d in range(3))
                cell = [-1]*3 if outside else [math.floor((point[d]-origin[d])/spacing[d]) for d in range(3)]
                coefficients = []
                for d in range(3):
                    t = (point[d]-origin[d]-cell[d]*spacing[d])/spacing[d]; t2, t3 = t*t, t*t*t
                    coefficients.append([-.5*t3+t2-.5*t, 1.5*t3-2.5*t2+1, -1.5*t3+2*t2+.5*t, .5*t3-.5*t2])
            offsets = [k//16, k//4 % 4, k % 4]
            xyz = [min(size[d]-1, max(0, cell[d]+offsets[d]-1)) for d in range(3)]
            require(indices[v] == xyz[0]+nx*(xyz[1]+ny*xyz[2]), "Stencil voxel mismatch.")
            coefficient = float(row["coefficient"])
            expected = coefficients[0][offsets[0]]*coefficients[1][offsets[1]]*coefficients[2][offsets[2]]
            require(abs(coefficient-expected) <= 64*sys.float_info.epsilon*max(1, abs(expected), *map(abs, cell)), "Stencil coefficient differs.")
            sampled[p] += coefficient*values[v]; seen[v] = True
    require(total == dataset["slot_count"] == 64*len(sampled) and np.array_equal(seen, flags)
            and int(flags.sum()) == dataset["stencil_rows"], "A subset mapping differs.")
    require(all(abs(y-s["response"]) <= 512*sys.float_info.epsilon*max(1, abs(s["response"]))
                for y, s in zip(sampled, dataset["samples"])), "Sample replay differs.")
    return values, flags, nearest


def validate_fit(fit, context, alpha, rows=ROWS):
    grid.validate_fit(fit, context, alpha, rows, NAME)
    require(fit.get("svd_preconditioner") == "independent-tsqr-8192", "Missing independent SVD path.")


def validate_regions(directory, name, flags, nearest, expected):
    sums = [[0, 0., 0., 0.] for _ in range(8)]
    with (directory/"residuals"/f"{name}.csv").open() as stream:
        for p, row in enumerate(csv.DictReader(stream)):
            r = float(row["residual"]); w = float(row["weight"]) if row["weight"] else 0
            for k in (0, 1 if flags[p] else 2, 3+min(4, int(nearest[p]/.5))):
                sums[k][0] += 1; sums[k][1] += r*r; sums[k][2] += w; sums[k][3] += w*w
    require(len(expected) == 8, "Missing region evidence.")
    for actual, reference in zip(sums, expected):
        n, ss, mass, square = actual
        require(n == reference["rows"], "Region row count differs.")
        for key, value in (("rmse", math.sqrt(ss/n) if n else None),
                           ("weight_share", mass/sums[0][2] if sums[0][2] else None),
                           ("effective_n", mass*mass/square if square else None)):
            require(value == reference[key] if value is None else math.isclose(value, reference[key], rel_tol=1e-10, abs_tol=1e-12), "Region diagnostic differs.")


def summarize(directory):
    directory = directory.resolve(); index = read(directory/"state-index.json"); completion = read(directory/"completion.json")
    names = expected_cases(index); fit_index = read(directory/"fit-index.json")
    require(completion["complete"] and completion["cases"] == names and fit_index["cases"] == names
            and fit_index["alphas"] == index["alphas"], "Missing requested cases.")
    dataset = read(directory/"dataset.json"); forward = read(directory/"forward-status.json")
    require(len(dataset["atoms"]) == 168 and len(dataset["samples"]) == 33600 and forward["passed"]
            and forward["voxel_count"] == ROWS, "Incomplete forward evidence.")
    values, flags, nearest = validate_dataset(directory, dataset)
    sid = index["selected_state_ids"][0]; state = next(s for s in index["states"] if s["id"] == sid)
    require(fold.sha256_file(Path(state["context"])) == state["sha256"], "Checkpoint changed.")
    context = read(Path(state["context"])); truth = read(directory/"scoring-truth.json")
    require([a["identity"] for a in context["atoms"]] == dataset["atoms"], "Atom identities differ.")
    replay = read(directory/f"{sid}-a-replay.json")
    require(replay["passed"] and [c["alpha"] for c in replay["cases"]] == index["alphas"]
            and all(c["passed"] for c in replay["cases"]), "Missing A replay evidence.")
    input_row = read(directory/f"{sid}-input.json"); require(input_row["abc"] == context["state"], "Changed input B.")
    grid.residual_statistics(directory, sid+"-input", values, dataset["samples"], input_row)
    validate_regions(directory, sid+"-input", flags, nearest, input_row["regions"])
    fits, statistics, estimates, comparisons, regions = [], [], [], [], []
    for name, alpha in zip(names, index["alphas"]):
        fit = read(directory/"fits"/f"{name}.json"); require(fit["state_id"] == sid, "Wrong fit checkpoint.")
        validate_fit(fit, context, alpha)
        grid.residual_statistics(directory, name, values, dataset["samples"], fit["residuals"], fit)
        validate_regions(directory, name, flags, nearest, fit["residuals"]["regions"])
        chosen = fit["branches"][fit["selected_seed"]] if fit["branches"] else {}
        diagnostics = chosen.get("endpoint_diagnostics", {})
        row = {"state_id": sid, "alpha": alpha, "qualified": fit["qualified"], "reason": fit["reason"],
            "stationarity": fit.get("stationarity"), "variance": fit.get("variances", [None])[0],
            "objective": fit.get("objective"), "branch_sensitive": fit.get("branch_sensitive", False),
            "primary_iterations": chosen.get("primary", {}).get("iterations"),
            "reference_iterations": chosen.get("reference", {}).get("iterations"),
            **{k: diagnostics.get(k) for k in ("rank", "condition", "minimum_singular", "effective_n", "effective_fraction", "maximum_row_share", "top_one_percent_share")},
            "grid_rmse": fit["residuals"]["grid_rmse"], "matched_sample_rmse": fit["residuals"]["matched_sample_rmse"], "seconds": fit.get("seconds")}
        fits.append(row)
        reference = read(directory/f"{name}-reference-a.json")
        source = Path(index["reference_a"])/"fits"/f"{name}.json"
        require(fold.sha256_file(source) == reference["source_sha256"], "A reference changed.")
        old = read(source)
        require(reference["abc"] == old["abc"] and reference["qualified"] == old["qualified"], "A projected coefficients changed.")
        grid.residual_statistics(directory, name+"-reference-a", values, dataset["samples"], reference)
        validate_regions(directory, name+"-reference-a", flags, nearest, reference["regions"])
        for method, result in (("A", reference), ("B", dict(fit["residuals"], qualified=fit["qualified"]))):
            regions += [dict(r, method=method, alpha=alpha, qualified=result["qualified"]) for r in result["regions"]]
            comparisons.append({"method": method, "alpha": alpha, "qualified": result["qualified"],
                "grid_rmse": result["grid_rmse"], "matched_sample_rmse": result["matched_sample_rmse"]})
        for method, abc, qualified in (("input", context["state"], None), ("A", reference["abc"], reference["qualified"]), ("B", fit["abc"], fit["qualified"])):
            errors = {p: [] for p in ("A", "B", "C")}
            for atom, identity in enumerate(dataset["atoms"]):
                target = truth[str(identity["serial_id"])]
                row = {"method": method, "alpha": alpha, "serial_id": identity["serial_id"], "qualified": qualified}
                for k, parameter in enumerate(errors):
                    error = abc[atom][k]-target[k]; errors[parameter].append(error)
                    row.update({parameter: abc[atom][k], f"truth_{parameter}": target[k], f"error_{parameter}": error})
                estimates.append(row)
            for parameter, errors in errors.items():
                statistics.append({"method": method, "alpha": alpha, "parameter": parameter, "qualified": qualified,
                    "below_0_01": sum(abs(e) < .01 for e in errors), **sweep.population_stats(errors, 168)})
    result = {"experiment": NAME, "state_id": sid, "alphas": index["alphas"], "radius": RADIUS,
        "row_count": ROWS, "stencil_rows": int(flags.sum()), "added_rows": ROWS-int(flags.sum()),
        "fits": fits, "statistics": statistics, "comparisons": comparisons, "regions": regions,
        "input": input_row, "forward": forward, "a_replay": replay}
    experiment.write(directory/"results.json", result)
    for name, rows in (("fits", fits), ("statistics", statistics), ("estimates", estimates), ("comparisons", comparisons), ("regions", regions)):
        ac.csv_write(directory/f"{name}.csv", rows)
    grid.artifact_index(directory)
    return result


def scientific(value):
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items() if not k.endswith("seconds") and k != "peak_child_rss_bytes"}
    if isinstance(value, list): return [scientific(v) for v in value]
    return value


def compare(left, right, output):
    li, ri = read(left/"state-index.json"), read(right/"state-index.json")
    require(li["selected_state_ids"] == ri["selected_state_ids"] and li["states"] == ri["states"], "Checkpoint mismatch.")
    require(ri["alphas"] == [0] and 0 in li["alphas"], "Expected independent alpha=0 replay.")
    sid = ri["selected_state_ids"][0]; case = f"{sid}-alpha-0"
    for directory, index in ((left, li), (right, ri)):
        complete = read(directory/"completion.json")
        require(complete["complete"] and complete["cases"] == expected_cases(index), "Missing completed cases.")
    names = ["dataset.json", "forward-status.json", "input-hashes.json", "provenance.json", "scoring-truth.json",
             f"{sid}-input.json", f"fits/{case}.json", f"{case}-reference-a.json"]
    differences = [name for name in names if scientific(read(left/name)) != scientific(read(right/name))]
    tables = ["voxels.csv", "slots.csv"]+[f"residuals/{name}{suffix}.csv" for name in (sid+"-input", case, case+"-reference-a") for suffix in ("", "-samples")]
    differences += [name for name in tables if fold.sha256_file(left/name) != fold.sha256_file(right/name)]
    for key in ("fits", "statistics", "comparisons", "regions"):
        a, b = ([r for r in read(p/"results.json")[key] if r["alpha"] == 0] for p in (left, right))
        if scientific(a) != scientific(b): differences.append(f"results/{key}")
    a, b = ([r for r in read(p/f"{sid}-a-replay.json")["cases"] if r["alpha"] == 0] for p in (left, right))
    if a != b: differences.append("a-replay")
    experiment.write(output, {"passed": not differences, "compared_alphas": [0], "not_repeated_alphas": [.1, .5, 1], "differences": differences})
    require(not differences, f"Scientific replay differs: {differences}")


def plots(directory, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    result = read(directory/"results.json"); output.mkdir(parents=True, exist_ok=True)
    def save(fig, name):
        for extension in ("png", "pdf"): fig.savefig(output/f"{name}.{extension}", dpi=180)
        plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    for axis, parameter in zip(axes, ("A", "C")):
        for method, color in (("A", "#596b8a"), ("B", "#16857b")):
            rows = [r for r in result["statistics"] if r["method"] == method and r["parameter"] == parameter]
            axis.plot([r["alpha"] for r in rows], [r["rmse"] for r in rows], "o-", color=color, label=f"Experiment {method}")
            for r in rows:
                if not r["qualified"]: axis.plot(r["alpha"], r["rmse"], "x", color="black", markersize=9)
        axis.set(title=f"{parameter} RMSE", xlabel="Common alpha", ylim=(0, None)); axis.grid(alpha=.2)
    axes[0].legend(); fig.suptitle("Fixed B | baseline-best-28 | x = unqualified"); save(fig, "parameter-errors")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    for axis, region, title in zip(axes, ("stencil", "added", "union"), ("Original A voxels", "Added voxels", "Complete B union")):
        for method, color in (("A", "#596b8a"), ("B", "#16857b")):
            rows = [r for r in result["regions"] if r["method"] == method and r["region"] == region]
            axis.plot([r["alpha"] for r in rows], [r["rmse"] for r in rows], "o-", color=color, label=f"Experiment {method}")
            for r in rows:
                if not r["qualified"]: axis.plot(r["alpha"], r["rmse"], "x", color="black", markersize=9)
        axis.set(title=title, xlabel="Common alpha", ylabel="Residual RMSE", ylim=(0, None)); axis.grid(alpha=.2)
    axes[0].legend(); fig.suptitle("A and B evaluated on identical regions | x = unqualified"); save(fig, "common-domain-residuals")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    fits = result["fits"]
    for axis, key, title in zip(axes[:2], ("stationarity", "effective_fraction"), ("Stationarity", "ESS / voxel count")):
        axis.plot([r["alpha"] for r in fits], [r[key] for r in fits], "o-", color="#16857b"); axis.set(title=title, xlabel="Common alpha"); axis.grid(alpha=.2)
    axes[0].axhline(1e-8, color="black", ls="--"); axes[0].set_yscale("log")
    for alpha in result["alphas"]:
        rows = [r for r in result["regions"] if r["method"] == "B" and r["alpha"] == alpha and r["region"].startswith("distance-")]
        axes[2].plot([.25, .75, 1.25, 1.75, 2.25], [r["weight_share"] for r in rows], "o-", label=f"alpha={alpha}")
    axes[2].set(title="Weight mass by nearest-atom distance", xlabel="Distance (angstrom)", ylabel="Share of all weights"); axes[2].legend(); axes[2].grid(alpha=.2)
    save(fig, "solver-and-spatial-diagnostics")


def main():
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest="operation", required=True)
    run_parser = sub.add_parser("run")
    for name in ("executable", "model", "map", "manifest", "baseline-run", "refined-run", "output"):
        run_parser.add_argument("--"+name, type=Path, required=True)
    run_parser.add_argument("--source", type=Path, default=ROOT)
    run_parser.add_argument("--state", default="baseline-best-28")
    run_parser.add_argument("--alpha", type=float, choices=ALPHAS)
    run_parser.add_argument("--reference-a", type=Path, default=ROOT/"build/unique-stencil-grid/final")
    summary = sub.add_parser("summarize"); summary.add_argument("directory", type=Path)
    comparison = sub.add_parser("compare")
    for name in ("left", "right", "output"): comparison.add_argument("--"+name, type=Path, required=True)
    plot = sub.add_parser("plots"); plot.add_argument("directory", type=Path); plot.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.operation == "run": run(args)
    elif args.operation == "summarize": summarize(args.directory)
    elif args.operation == "compare": compare(args.left, args.right, args.output)
    else: plots(args.directory, args.output)


if __name__ == "__main__": main()
