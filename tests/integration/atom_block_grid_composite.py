#!/usr/bin/env python3
"""Experiment C: overlapping atom blocks, fixed B, checkpoint alpha composite MDPDE."""
from __future__ import annotations

import argparse
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
import atom_centered_voxel_union as union
from observation_matching import read, require

ROOT = Path(__file__).resolve().parents[2]
NAME = "atom-block-grid-composite"
ROWS = 3768656
ALPHAS = [0, .1, .5, 1]
RADIUS = 2.5


def expected_cases(index):
    require(len(index["selected_state_ids"]) == 1 and index.get("alpha_source") == "checkpoint", "Invalid requested cases.")
    return [index["selected_state_ids"][0]+"-checkpoint-alpha"]


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
    index.update(selected_state_ids=[args.state], alpha_source="checkpoint",
                 reference_a=str(args.reference_a.resolve()), reference_b=str(args.reference_b.resolve()))
    reference_files = [directory.resolve()/"fits"/f"{args.state}-alpha-{k}.json"
                       for directory in (args.reference_a, args.reference_b) for k in range(4)]
    reference_files += [args.reference_b.resolve()/name for name in ("dataset.json", "voxels.csv", "slots.csv")]
    require(all(p.is_file() for p in reference_files), "Missing paired Experiment A/B evidence.")
    source = args.source.resolve(); executable = args.executable.resolve()
    provenance = experiment.provenance(executable, source)
    input_hashes = {str(p): fold.sha256_file(p) for p in list(inputs.values())+evidence+reference_files}
    output.mkdir(parents=True)
    for name, value in (("state-index", index), ("input-hashes", input_hashes), ("provenance", provenance)):
        experiment.write(output/f"{name}.json", value)
    experiment.write(output/"scoring-truth.json", {str(a["serial_id"]):
        [float(a["element"]), manifest["settings"]["blurring_width"], a["charge_used"]] for a in manifest["atoms"]})
    experiment.write(output/"environment.json", {"python": sys.version, "platform": platform.platform(), "jobs": 1,
        "iteration_budget": 100, "refinement_budget": 100, "radius": RADIUS, "independent_repeat": False,
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


def geometry_blocks(dataset):
    import numpy as np
    size = np.array(dataset["grid_size"]); origin = np.array(dataset["generation_origin"])
    spacing = np.array(dataset["generation_spacing"]); radius = dataset["radius"]
    for atom in dataset["atoms"]:
        center = np.array(atom["position"])
        lo = np.maximum(np.floor((center-radius-origin)/spacing).astype(int), 0)
        hi = np.minimum(np.floor((center+radius-origin)/spacing).astype(int), size-1)
        xx, yy, zz = [np.arange(lo[k], hi[k]+1) for k in range(3)]
        dx, dy, dz = [origin[k]+v*spacing[k]-center[k] for k, v in enumerate((xx, yy, zz))]
        square = dx[None, None, :]**2+dy[None, :, None]**2+dz[:, None, None]**2
        use = square <= radius*radius
        ids = (xx[None, None, :]+size[0]*(yy[None, :, None]+size[1]*zz[:, None, None]))[use]
        yield ids, square[use]


def validate_memberships(directory, dataset):
    import numpy as np
    indices, coverage, _ = union.sphere_membership(dataset)
    require(fold.sha256_file(directory/"memberships.csv") == dataset["membership_table_sha256"], "Membership hash differs.")
    blocks = []; count = 0; actual_coverage = np.zeros(len(indices), dtype=np.uint16)
    with (directory/"memberships.csv").open() as stream:
        reader = csv.DictReader(stream)
        for block, (ids, square) in enumerate(geometry_blocks(dataset)):
            rows = np.searchsorted(indices, ids)
            require(np.array_equal(indices[rows], ids), "Membership outside union.")
            for row, voxel in zip(rows, ids):
                entry = next(reader, None)
                require(entry is not None and (int(entry["block"]), int(entry["row"]), int(entry["index"])) ==
                        (block, int(row), int(voxel)), "Duplicate, missing or incorrect membership.")
            actual_coverage[rows] += 1; count += len(ids); blocks.append((rows, square))
        require(next(reader, None) is None, "Extra membership.")
    require(count == dataset["membership_count"] and np.array_equal(actual_coverage, coverage), "Membership coverage differs.")
    return blocks


def validate_fit(fit, context, blocks, rows=ROWS):
    require(fit.get("schema_version") == 2 and fit.get("experiment") == NAME and fit.get("alpha_source") == "checkpoint", "Unsupported composite fit.")
    require(fit["rows"] == rows and fit["columns"] == 2*len(context["state"]), "Fit population differs.")
    require(fit["iteration_budget"] == fit["refinement_budget"] == 100, "Changed iteration budget.")
    require(fit.get("linear_solver") == "sparse-qr" and fit.get("sparse_row_reduction") == "householder-qr-1024" and
            fit.get("svd_preconditioner") == "independent-tsqr-8192" and 0 < fit.get("design_nonzeros", 0) <= rows*fit["columns"], "Missing sparse/SVD evidence.")
    require(len(fit["blocks"]) == len(blocks) == len(context["atoms"]), "Missing atom blocks.")
    for i, ((members, _), metadata, atom) in enumerate(zip(blocks, fit["blocks"], context["atoms"])):
        require(metadata == {"owner": i, "alpha": atom["alpha"], "rows": len(members), "lambda": 1/len(blocks)}, "Block alpha or membership changed.")
        require(struct.pack("d", metadata["alpha"]) == struct.pack("d", atom["alpha"]), "Alpha changed.")
    require(fit["membership_count"] == sum(len(rows) for rows, _ in blocks), "Membership count differs.")
    require(fit["initial_beta"] == [v for a in context["state"] for v in (a[0], a[2])] and len(fit["abc"]) == len(blocks), "Changed initialization.")
    for i, (abc, saved) in enumerate(zip(fit["abc"], context["state"])):
        require(struct.pack("d", abc[1]) == struct.pack("d", saved[1]), "Width changed.")
        require([abc[0], abc[2]] == fit["beta"][2*i:2*i+2] and abc[0] >= 0 and all(math.isfinite(v) for v in abc), "Infeasible assembly.")
    branches = fit["branches"]
    require([b["seed"] for b in branches] == ["checkpoint", "constrained-ls"], "Missing initialization branch.")
    for b in branches:
        require(len(b["initial_variances"]) == len(blocks), "Missing initial variances.")
        for phase, trace in (("primary", "trace"), ("reference", "reference_trace")):
            endpoint = b[phase]; require(0 <= endpoint["iterations"] <= 100 and len(b[trace]) <= 101, "Iteration budget exceeded.")
            require(len(endpoint["variances"]) == len(blocks), "Missing endpoint variances.")
            if phase == "primary" or b["primary"]["stop"] == "stationary":
                require(len(b[trace]) == endpoint["iterations"]+1 and b[trace][-1]["beta"] == endpoint["beta"], "Incomplete trajectory.")
            else: require(not b[trace] and endpoint["iterations"] == 0, "Unqualified endpoint was refined.")
        if b["qualified"]:
            require(b["primary"]["stop"] == b["reference"]["stop"] == "stationary" and
                    b["primary"]["stationarity"] <= 1e-8 and b["reference"]["stationarity"] <= 1e-10 and
                    b.get("svd_reason") == "solved" and b["coefficient_difference"] <= 1e-6 and
                    b["log_variance_difference"] <= 1e-6 and b["weighted_spectrum"]["rank"] == fit["columns"], "Missing qualified evidence.")
    require(type(fit["selected_seed"]) is int and fit["selected_seed"] in (0, 1), "Invalid selected branch.")
    chosen = branches[fit["selected_seed"]]; endpoint = chosen["primary"]
    require(fit["beta"] == endpoint["beta"] and fit["variances"] == endpoint["variances"], "Selected endpoint differs.")
    if fit["qualified"]: require(chosen["qualified"] and fit["reason"] == "qualified", "Unverified selected fit.")
    else: require(not any(b["qualified"] for b in branches) and fit["selected_seed"] == 0 and
                  fit["reason"] == ("reference-unverified" if endpoint["stop"] == "stationary" else endpoint["stop"]), "Incorrect failure classification.")


def basis_columns(blocks, context):
    import numpy as np
    for (rows, square), abc in zip(blocks, context["state"]):
        r = np.sqrt(square); width = abc[1]
        gaussian = (2*math.pi*width*width)**-1.5*np.exp(-(r*r)/(2*width*width))
        charge = np.full(len(r), math.sqrt(2/math.pi)/width)
        noncentral = r >= 1e-5
        charge[noncentral] = np.fromiter((math.erf(float(d)/width/math.sqrt(2))/float(d) if d <= 2.5 else 0
                                        for d in r[noncentral]), dtype=float, count=int(noncentral.sum()))
        yield rows, gaussian, charge


def composite_evidence(residual, variances, alphas, blocks, columns, beta):
    import numpy as np
    n = np.array([len(p) for p, _ in blocks]); v = np.asarray(variances); a = np.asarray(alphas)
    logq = np.log1p(a)-math.log(len(blocks))-np.log(n)-np.log(v)-.5*a*np.log(2*math.pi*v)
    q = np.exp(logq-logq.max()); omega = np.zeros(len(residual)); prefactors = np.zeros(len(residual))
    scales = []; objectives = []; details = []; qv = float(np.sum(n*q*v))
    for i, (rows, _) in enumerate(blocks):
        r = residual[rows]; t = r*r/v[i]; w = np.ones(len(rows)) if a[i] == 0 else np.exp(-.5*a[i]*t)
        omega[rows] += q[i]*w; prefactors[rows] += q[i]
        equation = float(np.mean(w*(t-1))+a[i]*(1+a[i])**-1.5); scales.append(equation)
        if a[i] == 0: objective = .5*(math.log(2*math.pi*v[i])+float(np.mean(t)))
        else:
            logc = -.5*a[i]*math.log(2*math.pi*v[i])
            objective = math.exp(logc)*(math.expm1(-.5*math.log1p(a[i]))-(1+a[i])/a[i]*float(np.mean(np.expm1(-.5*a[i]*t))))-math.expm1(logc)/a[i]
        objectives.append(objective)
        details.append({"log_prefactor": float(logq[i]), "scaled_prefactor": float(q[i]), "scale_equation": equation,
                        "denominator": float(w.sum()-len(rows)*a[i]*(1+a[i])**-1.5),
                        "effective_n": float(w.sum()**2/(w@w)) if w@w else None, "mass": float(q[i]*w.sum())})
    equations = []
    for i, (rows, gaussian, charge) in enumerate(columns):
        for k, column in enumerate((gaussian, charge)):
            value = float(np.sum(omega[rows]*column*residual[rows]))/math.sqrt(qv*float(np.sum(prefactors[rows]*column*column)))
            equations.append(max(0, value) if k == 0 and beta[2*i] == 0 else value)
    for d in details: d["linear_weight_share"] = d.pop("mass")/float(omega.sum())
    return {"objective": sum(objectives)/len(objectives), "stationarity": max(map(abs, equations+scales)),
            "scaled_equations": equations+scales, "blocks": details, "omega": omega}


def validate_aggregate_weights(table, variances, alphas, blocks):
    import numpy as np
    v = np.asarray(variances); a = np.asarray(alphas); n = np.array([len(rows) for rows, _ in blocks])
    logq = np.log1p(a)-math.log(len(blocks))-np.log(n)-np.log(v)-.5*a*np.log(2*math.pi*v)
    q = np.exp(logq-logq.max()); low = np.zeros(len(table)); high = np.zeros(len(table))
    # Eigen's fused y-X*beta and y-(X*beta) can round differently. Propagate
    # the existing grid residual roundoff bound through each exponential,
    # then sum the membership bounds using the unchanged block prefactors.
    delta = 64*np.finfo(float).eps*np.maximum(1, np.maximum(abs(table[:, 1]), abs(table[:, 1]+table[:, 2])))
    for i, (rows, _) in enumerate(blocks):
        r = abs(table[rows, 2]); d = delta[rows]
        low[rows] += q[i]*np.exp(-.5*a[i]*(r+d)**2/v[i])
        high[rows] += q[i]*np.exp(-.5*a[i]*np.maximum(0, r-d)**2/v[i])
    subnormal = 16*np.nextafter(0., 1.)
    require(np.all(table[:, 3] >= low*(1-1e-12)-subnormal) and
            np.all(table[:, 3] <= high*(1+1e-12)+subnormal), "Aggregate weight replay differs.")
    return float(np.max(high-low))


def validate_numerics(directory, dataset, context, fit, values, blocks):
    import numpy as np
    y = np.asarray(values); columns = list(basis_columns(blocks, context)); alphas = [a["alpha"] for a in context["atoms"]]
    checks = []
    for branch in fit["branches"]:
        for phase in ("primary", "reference"):
            if phase == "reference" and not branch["reference_trace"]: continue
            endpoint = branch[phase]; beta = endpoint["beta"]; prediction = np.zeros(len(y))
            for i, (rows, gaussian, charge) in enumerate(columns): prediction[rows] += beta[2*i]*gaussian+beta[2*i+1]*charge
            if endpoint["objective"] is None:
                require(endpoint["stop"] == "exact-fit-boundary", "Unverified nonfinite endpoint.")
                owner = endpoint["failure_owner"]; rows = blocks[owner][0]
                row_square = np.zeros(len(y))
                for rr, gaussian, charge in columns: row_square[rr] += gaussian*gaussian+charge*charge
                bound = 128*np.finfo(float).eps*(np.linalg.norm(y[rows])+math.sqrt(float(row_square[rows].sum()))*np.linalg.norm(beta))
                require(np.linalg.norm((y-prediction)[rows]) <= bound, "Exact-fit boundary differs.")
                checks.append({"seed": branch["seed"], "phase": phase, "verified_failure": endpoint["stop"]})
                continue
            e = composite_evidence(y-prediction, endpoint["variances"], alphas, blocks, columns, beta)
            denominator_failure = endpoint["stop"] == "invalid-denominator"
            if denominator_failure:
                require(any(b["denominator"] <= 0 for b in e["blocks"]) and endpoint["stationarity"] is None,
                        "Invalid denominator classification differs.")
            require(math.isclose(e["objective"], endpoint["objective"], rel_tol=1e-10, abs_tol=1e-10) and
                    (denominator_failure or abs(e["stationarity"]-endpoint["stationarity"]) <= 1e-9), "Independent composite objective/stationarity differs.")
            require(np.allclose(e["scaled_equations"][-len(blocks):] if denominator_failure else e["scaled_equations"],
                                endpoint["scaled_equations"][-len(blocks):] if denominator_failure else endpoint["scaled_equations"], rtol=1e-7, atol=1e-9), "Composite equations differ.")
            if phase == "primary":
                diagnostics = branch["composite_diagnostics"]
                require(diagnostics["available"] and len(diagnostics["blocks"]) == len(blocks), "Missing block diagnostics.")
                for actual, saved in zip(e["blocks"], diagnostics["blocks"]):
                    for key, value in actual.items():
                        require(value == saved[key] if value is None else math.isclose(value, saved[key], rel_tol=1e-7, abs_tol=1e-9), f"Block diagnostic differs: {key}")
                require(abs(sum(b["linear_weight_share"] for b in diagnostics["blocks"])-1) < 1e-10, "Block shares double counted.")
            checks.append({"seed": branch["seed"], "phase": phase, "objective_difference": abs(e["objective"]-endpoint["objective"]),
                           "stationarity_difference": None if denominator_failure else abs(e["stationarity"]-endpoint["stationarity"])})
    name = fit["state_id"]+"-checkpoint-alpha"
    table = np.loadtxt(directory/"residuals"/f"{name}.csv", delimiter=",", skiprows=1)
    require(table.shape == (len(y), 4) and np.array_equal(table[:, 0], np.arange(len(y))), "Missing or duplicate residual rows.")
    require(np.all(np.abs(y-table[:, 1]-table[:, 2]) <= 64*np.finfo(float).eps*np.maximum(1, np.abs(y))), "Residual replay differs.")
    e = composite_evidence(table[:, 2], fit["variances"], alphas, blocks, columns, fit["beta"])
    maximum_weight_bound_width = validate_aggregate_weights(table, fit["variances"], alphas, blocks)
    # Independently verify RSS_i/n_i initialization for both coefficient starts.
    for branch in fit["branches"]:
        prediction = np.zeros(len(y)); beta = branch["initial_beta"]
        for i, (rows, gaussian, charge) in enumerate(columns): prediction[rows] += beta[2*i]*gaussian+beta[2*i+1]*charge
        variances = [float(np.mean((y[rows]-prediction[rows])**2)) for rows, _ in blocks]
        require(np.allclose(variances, branch["initial_variances"], rtol=1e-10, atol=1e-20), "Initial block variance differs.")
    return {"passed": True, "endpoints": checks, "aggregate_weight_maximum_difference": float(np.max(np.abs(e["omega"]-table[:, 3]))),
            "maximum_propagated_weight_bound_width": maximum_weight_bound_width}


def summarize(directory):
    directory = directory.resolve(); index = read(directory/"state-index.json"); names = expected_cases(index)
    completion = read(directory/"completion.json"); fit_index = read(directory/"fit-index.json")
    require(completion["complete"] and completion["cases"] == names and fit_index["cases"] == names, "Missing requested cases.")
    dataset = read(directory/"dataset.json"); sid = index["selected_state_ids"][0]
    state = next(s for s in index["states"] if s["id"] == sid); context = read(Path(state["context"]))
    require(fold.sha256_file(Path(state["context"])) == state["sha256"], "Checkpoint changed.")
    require([a["identity"] for a in context["atoms"]] == dataset["atoms"] and
            dataset["alphas"] == [a["alpha"] for a in context["atoms"]], "Identity or alpha changed.")
    values, flags, nearest = union.validate_dataset(directory, dataset, ROWS, NAME, dataset["alphas"])
    blocks = validate_memberships(directory, dataset)
    require(dataset["membership_count"] == 10995198, "Wrong complete membership population.")
    forward = read(directory/"forward-status.json"); require(forward["passed"] and forward["voxel_count"] == ROWS, "Missing forward evidence.")
    replay = read(directory/f"{sid}-b-replay.json")
    require(replay["passed"] and [c["alpha"] for c in replay["cases"]] == ALPHAS and all(c["passed"] for c in replay["cases"]), "Missing B replay evidence.")
    fit = read(directory/"fits"/f"{names[0]}.json"); validate_fit(fit, context, blocks)
    numeric = validate_numerics(directory, dataset, context, fit, values, blocks)
    experiment.write(directory/"numerical-validation.json", numeric)
    grid.residual_statistics(directory, names[0], values, dataset["samples"], fit["residuals"])
    union.validate_regions(directory, names[0], flags, nearest, fit["residuals"]["regions"])
    historical = read(directory/f"{sid}-historical.json")
    require([(r["method"], r["alpha"]) for r in historical] == [(m, a) for a in ALPHAS for m in ("B", "A")], "Missing paired history.")
    for reference in historical:
        method, alpha = reference["method"], reference["alpha"]; k = ALPHAS.index(alpha)
        saved = read(Path(index[f"reference_{method.lower()}"])/"fits"/f"{sid}-alpha-{k}.json")
        require(reference["abc"] == saved["abc"] and reference["qualified"] == saved["qualified"], "Historical endpoint changed.")
        name = f"{sid}-reference-{method.lower()}-{k}"
        grid.residual_statistics(directory, name, values, dataset["samples"], reference)
        union.validate_regions(directory, name, flags, nearest, reference["regions"])
    truth = read(directory/"scoring-truth.json"); statistics = []; estimates = []; comparisons = []; regions = []
    own = dict(fit["residuals"], method="C", alpha="checkpoint", qualified=fit["qualified"], reason=fit["reason"], abc=fit["abc"])
    for result in historical+[own]:
        label = {k: result[k] for k in ("method", "alpha", "qualified", "reason")}
        comparisons.append(dict(label, grid_rmse=result["grid_rmse"], matched_sample_rmse=result["matched_sample_rmse"]))
        regions += [dict(r, **label) for r in result["regions"]]
        errors = {p: [] for p in ("A", "B", "C")}
        for atom, abc in zip(dataset["atoms"], result["abc"]):
            target = truth[str(atom["serial_id"])]; row = dict(label, serial_id=atom["serial_id"])
            for k, parameter in enumerate(errors):
                error = abc[k]-target[k]; errors[parameter].append(error)
                row.update({parameter: abc[k], f"truth_{parameter}": target[k], f"error_{parameter}": error})
            estimates.append(row)
        for parameter, error in errors.items():
            statistics.append(dict(label, parameter=parameter, below_0_01=sum(abs(v)<.01 for v in error), **sweep.population_stats(error, 168)))
    chosen = fit["branches"][fit["selected_seed"]]
    fit_row = {"state_id": sid, **{k: fit[k] for k in ("qualified", "reason", "stationarity", "objective", "branch_sensitive")},
               **chosen["endpoint_diagnostics"], "primary_iterations": chosen["primary"]["iterations"],
               "reference_iterations": chosen["reference"]["iterations"], **{k: fit["residuals"][k] for k in ("grid_rmse", "matched_sample_rmse")}}
    block_rows = [dict(b, seed=branch["seed"], serial_id=dataset["atoms"][b["owner"]]["serial_id"])
                  for branch in fit["branches"] for b in branch["composite_diagnostics"]["blocks"]]
    result = {"experiment": NAME, "state_id": sid, "alpha_source": "checkpoint", "independent_repeat": False,
              "row_count": ROWS, "membership_count": dataset["membership_count"], "fits": [fit_row],
              "statistics": statistics, "comparisons": comparisons, "regions": regions, "blocks": block_rows,
              "composite_diagnostics": chosen["composite_diagnostics"], "forward": forward, "b_replay": replay, "numeric_validation": numeric}
    experiment.write(directory/"results.json", result)
    for name, rows in (("fits", [fit_row]), ("statistics", statistics), ("estimates", estimates), ("comparisons", comparisons), ("regions", regions), ("blocks", block_rows)):
        ac.csv_write(directory/f"{name}.csv", [{k:v for k,v in r.items() if not isinstance(v, (list, dict))} for r in rows])
    grid.artifact_index(directory)
    return result


def compare(left, right, output):
    li, ri = read(left/"state-index.json"), read(right/"state-index.json")
    require(expected_cases(li) == expected_cases(ri) and li["states"] == ri["states"], "Checkpoint mismatch.")
    case = expected_cases(li)[0]
    for directory in (left, right):
        status = read(directory/"completion.json"); require(status["complete"] and status["cases"] == [case], "Incomplete comparison.")
    files = ["dataset.json", "forward-status.json", "results.json", f"fits/{case}.json"]
    differences = [p for p in files if union.scientific(read(left/p)) != union.scientific(read(right/p))]
    differences += [p for p in ("voxels.csv", "memberships.csv", "slots.csv") if fold.sha256_file(left/p) != fold.sha256_file(right/p)]
    experiment.write(output, {"passed": not differences, "differences": differences})
    require(not differences, "Scientific outputs differ.")


def plots(directory, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    result = read(directory/"results.json"); output.mkdir(parents=True, exist_ok=True)
    def save(fig, name):
        for ext in ("png", "pdf"): fig.savefig(output/f"{name}.{ext}", dpi=180)
        plt.close(fig)
    labels = [f"{m} / {a:g}" for m in ("A", "B") for a in ALPHAS]+["C / atom alpha"]
    order = [(m,a) for m in ("A", "B") for a in ALPHAS]+[("C","checkpoint")]
    fig, axes = plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
    for axis, parameter in zip(axes,("A","C")):
        rows = [next(r for r in result["statistics"] if (r["method"],r["alpha"])==key and r["parameter"]==parameter) for key in order]
        axis.bar(range(9),[r["rmse"] for r in rows],color=["#596b8a"]*4+["#16857b"]*4+["#bd6c34"])
        for i,r in enumerate(rows):
            if not r["qualified"]: axis.plot(i,r["rmse"],"kx")
        axis.set(title=f"{parameter} RMSE",xticks=range(9),xticklabels=labels); axis.tick_params(axis="x",rotation=60)
    fig.suptitle("Fixed B | baseline-best-28 | x = unqualified"); save(fig,"parameter-errors")
    fig, axes = plt.subplots(1,3,figsize=(14,4.8),constrained_layout=True)
    for axis, region in zip(axes,("stencil","added","union")):
        rows = [next(r for r in result["regions"] if (r["method"],r["alpha"])==key and r["region"]==region) for key in order]
        axis.bar(range(9),[r["rmse"] for r in rows],color=["#596b8a"]*4+["#16857b"]*4+["#bd6c34"])
        for i,r in enumerate(rows):
            if not r["qualified"]: axis.plot(i,r["rmse"],"kx")
        axis.set(title=region+" voxels",xticks=range(9),xticklabels=labels); axis.tick_params(axis="x",rotation=60)
    fig.suptitle("Identical observation domains | x = unqualified"); save(fig,"common-domain-residuals")
    fit = read(directory/"fits"/f"{result['state_id']}-checkpoint-alpha.json")
    fig, axes = plt.subplots(1,3,figsize=(14,4.8),constrained_layout=True)
    for branch in fit["branches"]:
        axes[0].semilogy([r["iterations"] for r in branch["trace"]],[r["stationarity"] for r in branch["trace"]],label=branch["seed"])
        b=branch["composite_diagnostics"]["blocks"]
        axes[1].plot(range(1,len(b)+1),[r["linear_weight_share"] for r in b],label=branch["seed"])
        axes[2].plot([.25,.75,1.25,1.75,2.25],[r["weight_share"] for r in branch["composite_diagnostics"]["membership_radial"]],"o-",label=branch["seed"])
    axes[0].axhline(1e-8,color="black",ls="--"); axes[0].set(title="Stationarity",xlabel="Primary iteration"); axes[0].legend()
    axes[1].set(title="Block mass shares",xlabel="Atom index (1-based)",ylabel="Share of aggregate mass")
    axes[2].set(title="Membership mass by block-center distance",xlabel="Distance (angstrom)",ylabel="Share of aggregate mass")
    save(fig,"solver-and-block-diagnostics")


def main():
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest="operation",required=True)
    runner=sub.add_parser("run")
    for name in ("executable","model","map","manifest","baseline-run","refined-run","output"):
        runner.add_argument("--"+name,type=Path,required=True)
    runner.add_argument("--source",type=Path,default=ROOT); runner.add_argument("--state",default="baseline-best-28")
    runner.add_argument("--reference-a",type=Path,default=ROOT/"build/unique-stencil-grid/final")
    runner.add_argument("--reference-b",type=Path,default=ROOT/"build/atom-centered-voxel-union/final")
    summary=sub.add_parser("summarize"); summary.add_argument("directory",type=Path)
    comparison=sub.add_parser("compare")
    for name in ("left","right","output"): comparison.add_argument("--"+name,type=Path,required=True)
    plot=sub.add_parser("plots"); plot.add_argument("directory",type=Path); plot.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.operation=="run": run(args)
    elif args.operation=="summarize": summarize(args.directory)
    elif args.operation=="compare": compare(args.left,args.right,args.output)
    else: plots(args.directory,args.output)


if __name__=="__main__": main()
