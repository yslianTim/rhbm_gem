#!/usr/bin/env python3
"""Four fixed-B, mean-only LS controls on the fold-168 0.30 A voxel union."""
from __future__ import annotations

import argparse
import csv
import hashlib
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

for _thread_variable in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_thread_variable] = "1"
import numpy as np
import atom_centered_voxel_union as union
import fold_168_regression as fold
import mdpde_experiment as experiment
from observation_matching import read, require

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/benchmarks/fixed_b_oracle.json"
CASES = ["true-b-double", "true-b-float32", "checkpoint-b-double", "checkpoint-b-float32"]
IDENTITY = ["serial_id", "chain_id", "sequence_id", "component_id", "atom_id", "alternate_indicator", "position"]


def checkpoint_widths(manifest, context):
    require(context["schema_version"] == 1 and len(context["atoms"]) == len(context["state"]) == 168,
            "Invalid checkpoint population.")
    captured = {a["identity"]["serial_id"]: (k, a) for k, a in enumerate(context["atoms"])}
    require(len(captured) == 168, "Duplicate checkpoint identity.")
    widths = []
    for atom in manifest["atoms"]:
        require(atom["serial_id"] in captured, "Missing checkpoint identity.")
        k, actual = captured[atom["serial_id"]]
        require(actual["index"] == k and all(atom[key] == actual["identity"][key] for key in IDENTITY),
                "Checkpoint identity mismatch.")
        b = context["state"][k][1]
        require(math.isfinite(b) and b > 0, "Invalid checkpoint B.")
        widths.append(b)
    return widths


def validate_inputs(paths, fixture):
    hashes = fold.validate_input_hashes({k: paths[k] for k in fixture["input_hashes"]}, fixture["input_hashes"])
    manifest = fold.load_simulation_manifest(paths["manifest"], hashes)
    for key, value in fixture["generation_contract"].items():
        require(manifest[key] == value, f"Changed oracle generation contract: {key}.")
    require(fold.sha256_file(paths["checkpoint"]) == fixture["checkpoint_sha256"], "Checkpoint SHA-256 mismatch.")
    require([a["preparation_index"] for a in manifest["atoms"]] == list(range(168)), "Changed preparation order.")
    return manifest, checkpoint_widths(manifest, read(paths["checkpoint"]))


def run(args, operation="fixed-b-oracle", summarize_result=None):
    output = args.output.resolve()
    require(not output.exists(), "Use a fresh output directory.")
    paths = {k: getattr(args, k).resolve() for k in ("model", "map", "manifest", "checkpoint")}
    validate_inputs(paths, read(FIXTURE))
    executable = args.executable.resolve()
    before = experiment.provenance(executable, ROOT)
    hashes = {str(p): fold.sha256_file(p) for p in [*paths.values(), FIXTURE]}
    output.mkdir(parents=True)
    for name, data in (("inputs", {k: str(v) for k, v in paths.items()}), ("input-hashes", hashes),
                       ("provenance", before), ("fixture", read(FIXTURE)),
                       ("environment", {"python": sys.version, "platform": platform.platform(), "jobs": 1,
                        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()})):
        experiment.write(output / f"{name}.json", data)
    command = [str(executable), operation, *(str(paths[k]) for k in ("manifest", "map", "checkpoint")), str(output)]
    env = dict(os.environ, OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    start = time.perf_counter(); interrupted = False
    with (output / "run.log").open("w") as log:
        child = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        try:
            status = child.wait()
        except KeyboardInterrupt:
            interrupted = True; child.send_signal(signal.SIGINT); status = child.wait()
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    experiment.write(output / "execution.json", {"command": command, "returncode": status,
        "seconds": time.perf_counter()-start, "jobs": 1, "user_interrupted": interrupted,
        "peak_child_rss_bytes": peak if sys.platform == "darwin" else peak*1024})
    stable = before == experiment.provenance(executable, ROOT) and all(fold.sha256_file(Path(p)) == h for p, h in hashes.items())
    experiment.write(output / "execution-status.json", {"inputs_and_sources_stable": stable,
        "reason": "user-interrupted" if interrupted else "technical-failure" if status or not stable else "completed"})
    require(status == 0 and stable, "Incomplete execution or changed provenance; retained completed cases.")
    (summarize_result or summarize)(output)


def load_table(path):
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")


def validate_quantization(y64, y32, q):
    require(all(np.isfinite(v).all() for v in (y64, y32, q)) and np.array_equal(y64.astype(np.float32).astype(float), y32)
            and np.array_equal(y32-y64, q), "Quantization differs.")


def validate_dataset(directory, dataset, fixture, paths, manifest, experiment_name="fixed-b-oracle"):
    ids, coverage, nearest = union.sphere_membership(dataset)
    require(dataset["experiment"] == experiment_name and dataset["radius"] == 2.5 and
            dataset["membership_geometry"] == "generation", "Invalid oracle geometry.")
    require(len(ids) == dataset["row_count"] == fixture["expected_rows"] and
            int(coverage.sum()) == dataset["memberships"] == fixture["expected_memberships"], "Wrong voxel population.")
    require(hashlib.sha256(ids.astype("<u8").tobytes()).hexdigest() == fixture["voxel_ids_sha256_le_u64"], "Wrong voxel IDs.")
    require(fold.sha256_file(directory / "voxels.csv") == dataset["voxel_table_sha256"], "Voxel table hash mismatch.")
    table = load_table(directory / "voxels.csv")
    require(len(table) == len(ids) and np.array_equal(table["row"], np.arange(len(ids))) and
            np.array_equal(table["index"], ids) and np.array_equal(table["multiplicity"], coverage), "Duplicate or incorrect voxel index.")
    nx, ny, _ = dataset["grid_size"]
    xyz = np.column_stack((ids % nx, ids//nx % ny, ids//(nx*ny)))
    points = np.asarray(dataset["generation_origin"]) + xyz*np.asarray(dataset["generation_spacing"])
    require(np.allclose(np.column_stack([table[k] for k in ("x", "y", "z")]), points, rtol=0, atol=2e-14) and
            np.allclose(table["nearest_distance"], nearest, rtol=0, atol=2e-13), "Voxel geometry differs.")
    y64, y32, q = (table[k] for k in ("reference_double", "observed", "quantization_delta"))
    validate_quantization(y64, y32, q)
    with paths["map"].open("rb") as stream:
        header = stream.read(1024)
        require(struct.unpack_from("<3i", header) == tuple(dataset["grid_size"]) and
                struct.unpack_from("<i", header, 12)[0] == 2 and struct.unpack_from("<3i", header, 64) == (1, 2, 3), "Unsupported map header.")
        stream.seek(1024+struct.unpack_from("<i", header, 92)[0]); values = np.fromfile(stream, dtype="<f4")
    require(len(values) == math.prod(dataset["grid_size"]) and np.array_equal(values[ids].astype(float), y32), "Map values differ.")
    settings = manifest["settings"]
    require(all(dataset[a] == settings[b] for a, b in (("grid_size", "grid_size"),
            ("generation_origin", "origin"), ("generation_spacing", "grid_spacing"))) and
            dataset["atoms"] == [{k: a[k] for k in dataset["atoms"][0]} for a in manifest["atoms"]], "Manifest geometry/identities differ.")
    return points, y64, y32


def columns(points, atoms, widths):
    result = []
    for atom, b in zip(atoms, widths):
        delta = points-np.asarray(atom["position"])
        square = delta[:, 0]**2+delta[:, 1]**2+delta[:, 2]**2
        indices = np.flatnonzero(square <= 2.5**2)
        r = np.sqrt(square[indices]); g = (2*math.pi*b*b)**(-1.5)*np.exp(-(r*r)/(2*b*b))
        charge = np.full(len(r), math.sqrt(2/math.pi)/b)
        use = r >= 1e-5
        charge[use] = np.fromiter((math.erf(v/b/math.sqrt(2))/v for v in r[use]), float, count=int(use.sum()))
        result.extend(((indices, g), (indices, charge)))
    return result


def predict(basis, beta, n):
    out = np.zeros(n)
    for (ids, values), coefficient in zip(basis, beta): out[ids] += values*coefficient
    return out


def scaled_difference(a, b):
    return float(np.max(np.abs(a-b)/(1+np.maximum(np.abs(a), np.abs(b)))))


def validate_fit(directory, name, fit, basis, y, widths):
    require(fit["case"] == name and fit["experiment"] == "fixed-b-oracle" and fit["alpha"] == 0 and
            fit["b"] == widths and fit["row_count"] == len(y), "Changed fixed B/case contract.")
    require(fit["linear_solver"] == "sparse-qr-householder-1024" and
            fit["reference_solver"] == "independent-tsqr-8192-svd" and
            fit["kkt_tolerance"] == fit["reference_tolerance"] == 1e-10, "Changed solver contract.")
    for label in ("primary", "reference"):
        endpoint = fit[label]; beta = np.asarray(endpoint["beta"])
        require(beta.shape == (336,) and np.isfinite(beta).all(), "Invalid coefficients.")
        residual = predict(basis, beta, len(y))-y
        norms = np.array([np.linalg.norm(v) for _, v in basis]); scale = max(1, np.linalg.norm(y))
        u = norms*beta/scale
        g = np.array([v@residual[i] for i, v in basis])/norms/scale
        projection = u-g; projection[::2] = np.maximum(0, projection[::2])
        kkt = float(np.max(np.abs(u-projection))); feasible = bool(np.all(beta[::2] >= 0))
        require(abs(kkt-endpoint["projected_kkt"]) <= 1e-13 and endpoint["feasible"] == feasible and
                endpoint["kkt_passed"] == (feasible and endpoint["projected_kkt"] <= 1e-10), "KKT certificate differs.")
        require(endpoint["active_atoms"] == np.flatnonzero(beta[::2] == 0).tolist(), "Active constraints differ.")
    primary, reference = (np.asarray(fit[k]["beta"]) for k in ("primary", "reference"))
    difference = scaled_difference(primary, reference)
    require(math.isclose(difference, fit["scaled_reference_difference"], rel_tol=1e-12, abs_tol=1e-16), "Reference difference differs.")
    qualified = all(fit[k]["valid"] and fit[k]["kkt_passed"] for k in ("primary", "reference")) and difference <= 1e-10 and fit["design_spectrum"]["rank"] == 336
    require(fit["qualified"] == qualified, "Unverified numerical qualification.")
    path = directory / "residuals" / f"{name}.csv"
    require(fold.sha256_file(path) == fit["residual_sha256"], "Residual hash mismatch.")
    table = load_table(path)
    require(len(table) == len(y) and np.array_equal(table["row"], np.arange(len(y))) and
            np.allclose(table["prediction"], predict(basis, primary, len(y)), rtol=2e-13, atol=2e-12) and
            np.array_equal(table["residual"], y-table["prediction"]), "Residual replay differs.")
    rss = float(table["residual"]@table["residual"])
    expected = {"rss": rss, "objective": rss/2, "residual_scale": rss/len(y), "residual_rmse": math.sqrt(rss/len(y)),
                "residual_max": float(np.max(np.abs(table["residual"]))), "relative_residual": math.sqrt(rss)/max(1, np.linalg.norm(y))}
    require(all(math.isclose(fit["primary"][k], v, rel_tol=1e-10, abs_tol=1e-28) for k, v in expected.items()), "Residual statistics differ.")


def statistics(values):
    out = experiment.stats(values)
    out["below_0_01"] = int(np.sum(np.abs(values) < .01))
    return out


def write_csv(path, rows):
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def summarize(directory):
    directory = directory.resolve(); start = time.perf_counter()
    fixture = read(directory / "fixture.json"); paths = {k: Path(p) for k, p in read(directory / "inputs.json").items()}
    require(fixture == read(FIXTURE), "Changed frozen fixture.")
    manifest, checkpoint = validate_inputs(paths, fixture)
    require(all(fold.sha256_file(Path(p)) == h for p, h in read(directory / "input-hashes.json").items()), "Changed input hashes.")
    completion = read(directory / "completion.json")
    require(completion["complete"] and completion["cases"] == CASES, "Missing completed case.")
    dataset = read(directory / "dataset.json"); forward = read(directory / "forward-status.json")
    require(forward["passed"] and forward["voxel_count"] == fixture["expected_rows"] and forward["generator_jobs"] == 1
            and dataset["checkpoint_b"] == checkpoint, "Invalid forward/B evidence.")
    points, y64, y32 = validate_dataset(directory, dataset, fixture, paths, manifest)
    truth = np.array([v for a in manifest["atoms"] for v in (a["element"], a["charge_used"])])
    require(read(directory / "scoring-truth.json") == {"beta": truth.tolist(), "width": .5}, "Truth mismatch.")
    fits = {}; estimates = []; summaries = []; decompositions = []
    for prefix, widths in (("true-b", [.5]*168), ("checkpoint-b", checkpoint)):
        basis = columns(points, manifest["atoms"], widths)
        if prefix == "true-b":
            expected = predict(basis, truth, len(y64))
            absolute = predict([(i, np.abs(v)) for i, v in basis], np.abs(truth), len(y64))
            require(np.all(np.abs(expected-y64) <= 512*np.finfo(float).eps*np.maximum(1, absolute)), "Independent double forward mismatch.")
        for precision, y in (("double", y64), ("float32", y32)):
            name = f"{prefix}-{precision}"; fit = read(directory / "fits" / f"{name}.json")
            validate_fit(directory, name, fit, basis, y, widths); fits[name] = fit
            error = np.asarray(fit["primary"]["beta"])-truth
            summaries.append({"case": name, "qualified": fit["qualified"], "A": statistics(error[::2]), "C": statistics(error[1::2]),
                **{k: fit["primary"][k] for k in ("residual_rmse", "residual_max", "projected_kkt", "active_atoms")},
                "scaled_reference_difference": fit["scaled_reference_difference"], **fit["design_spectrum"]})
            for a, atom in enumerate(manifest["atoms"]):
                estimates.append({"case": name, "serial_id": atom["serial_id"], "A": fit["primary"]["beta"][2*a], "B": widths[a],
                    "C": fit["primary"]["beta"][2*a+1], "A_error": error[2*a], "B_error": widths[a]-.5, "C_error": error[2*a+1]})
    beta = {k: np.asarray(v["primary"]["beta"]) for k, v in fits.items()}
    resolution = {k: np.abs(np.asarray(v["reference"]["beta"])-beta[k]) for k, v in fits.items()}
    td, tf, cd, cf = (beta[k] for k in CASES)
    effects = {"oracle_error": td-truth, "quantization_true_b": tf-td, "width_double": cd-td,
               "quantization_checkpoint_b": cf-cd, "interaction": (cf-cd)-(tf-td)}
    envelopes = {"quantization_true_b": resolution[CASES[0]]+resolution[CASES[1]],
                 "quantization_checkpoint_b": resolution[CASES[2]]+resolution[CASES[3]],
                 "interaction": sum(resolution.values())}
    effect_summary = {}
    for name, values in effects.items():
        effect_summary[name] = {"A": statistics(values[::2]), "C": statistics(values[1::2])}
        for k, value in enumerate(values):
            bound = envelopes[name][k] if name in envelopes else None
            decompositions.append({"effect": name, "serial_id": manifest["atoms"][k//2]["serial_id"], "parameter": "A" if k%2 == 0 else "C",
                "delta": value, "qr_svd_resolution": bound, "resolved": abs(value) > bound if bound is not None else None})
        if name in envelopes:
            for offset, key in enumerate(("A", "C")):
                effect_summary[name][key]["resolved_count"] = int(np.sum(np.abs(values[offset::2]) > envelopes[name][offset::2]))
    require(np.allclose(effects["oracle_error"]+effects["quantization_true_b"]+effects["width_double"]+effects["interaction"], cf-truth,
                        rtol=1e-12, atol=1e-14), "Error decomposition does not close.")
    recovery = {"scaled_coefficient_error": scaled_difference(td, truth), "relative_residual": fits[CASES[0]]["primary"]["relative_residual"]}
    recovery["passed"] = recovery["scaled_coefficient_error"] <= 1e-10 and recovery["relative_residual"] <= 1e-12
    results = {"experiment": "fixed-b-oracle", "rows": len(y64), "complete": True, "cases": summaries,
        "oracle_recovery": recovery, "effects": effect_summary, "quantization": statistics(y32-y64),
        "next_stage_ready": recovery["passed"] and all(f["qualified"] for f in fits.values())}
    experiment.write(directory / "results.json", results)
    write_csv(directory / "estimates.csv", estimates); write_csv(directory / "decomposition.csv", decompositions)
    write_csv(directory / "comparison.csv", [{**{k: v for k, v in row.items() if k not in ("A", "C", "active_atoms")},
        "active_constraints": len(row["active_atoms"]), **{f"{parameter}_{k}": v for parameter in ("A", "C") for k, v in row[parameter].items()}} for row in summaries])
    experiment.write(directory / "independent-validation.json", {"passed": True, "seconds": time.perf_counter()-start,
        "checks": ["full ROI IDs/coverage/geometry", "map bytes/float32 quantization", "analytic forward", "fixed B identity", "independent KKT/prediction", "raw residual statistics", "vector decomposition"]})
    print({"cases": len(fits), "qualified": sum(f["qualified"] for f in fits.values()), "oracle_recovery": recovery, "next_stage_ready": results["next_stage_ready"]})
    return results


def scientific(value):
    if isinstance(value, list): return [scientific(v) for v in value]
    if isinstance(value, dict): return {k: scientific(v) for k, v in value.items() if not k.endswith("seconds")}
    return value


def compare(left, right, output):
    for directory in (left, right): summarize(directory)
    files = ["dataset.json", "forward-status.json", "results.json", "completion.json", "scoring-truth.json", "fixture.json", "provenance.json", "input-hashes.json"]
    files += [f"fits/{name}.json" for name in CASES]
    differences = [name for name in files if scientific(read(left/name)) != scientific(read(right/name))]
    for name in ["voxels.csv", "estimates.csv", "decomposition.csv", "comparison.csv", *[f"residuals/{c}.csv" for c in CASES]]:
        if fold.sha256_file(left/name) != fold.sha256_file(right/name): differences.append(name)
    experiment.write(output, {"passed": not differences, "cases": CASES, "differences": differences,
        "excluded": ["timing", "resource usage", "execution/output directory"], "comparison": "exact scientific JSON and CSV bytes"})
    require(not differences, f"Independent repeat differs: {differences}")


def plots(directory, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output.mkdir(parents=True, exist_ok=True)
    table = load_table(directory/"estimates.csv")
    figure, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for row, key in enumerate(("A_error", "C_error")):
        for col, prefix in enumerate(("true-b", "checkpoint-b")):
            ax = axes[row, col]
            for precision, style in (("double", "-"), ("float32", "--")):
                selected = table[table["case"] == f"{prefix}-{precision}"]
                ax.plot(selected["serial_id"], selected[key], style, label=precision, linewidth=1)
            ax.set(xlabel="Atom serial ID", ylabel=key.replace("_", " "), title=prefix)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0)); ax.legend(); ax.grid(alpha=.2)
    for suffix in ("png", "pdf"): figure.savefig(output/f"parameter-errors.{suffix}", dpi=170)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__); commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("run")
    for key in ("executable", "model", "map", "manifest", "checkpoint", "output"): p.add_argument(f"--{key}", type=Path, required=True)
    p = commands.add_parser("summarize"); p.add_argument("directory", type=Path)
    p = commands.add_parser("compare")
    for key in ("left", "right", "output"): p.add_argument(f"--{key}", type=Path, required=True)
    p = commands.add_parser("plots"); p.add_argument("directory", type=Path); p.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "run": run(args)
    elif args.command == "summarize": summarize(args.directory)
    elif args.command == "compare": compare(args.left, args.right, args.output)
    else: plots(args.directory, args.output)


if __name__ == "__main__": main()
