#!/usr/bin/env python3
"""Eight fixed-ROI joint A/C/B variable-projection least-squares experiments."""
from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path
import time

import fixed_b_oracle as fixed  # Sets single-thread environment before NumPy import.
import numpy as np
from observation_matching import read, require

ROOT = fixed.ROOT
STARTS = ["checkpoint", "narrower", "wider", "mixed"]
CASES = [f"{start}-{precision}" for precision in ("double", "float32") for start in STARTS]
SETTINGS = {"factor": .1, "ftol": 1e-14, "xtol": 1e-12, "gtol": 1e-12,
            "profile_budget": 200, "accepted_update_budget": 100}
CONTROLS = ROOT / "docs/developer/figures/fixed-b-oracle"


def initial_widths(checkpoint, atoms, start):
    require(start in STARTS, "Unknown initialization.")
    if start == "checkpoint": return np.asarray(checkpoint).copy()
    factors = [.8 if start == "narrower" or (start == "mixed" and a["serial_id"] % 2) else 1.2 for a in atoms]
    return np.asarray(checkpoint)*factors


def raw_certificate(points, atoms, y, endpoint, column_builder=fixed.columns):
    beta = np.asarray(endpoint["beta"]); widths = np.asarray(endpoint["b"])
    require(beta.shape == (2*len(atoms),) and widths.shape == (len(atoms),) and
            np.isfinite(beta).all() and np.isfinite(widths).all() and np.all(widths > 0), "Invalid joint endpoint.")
    basis = column_builder(points, atoms, widths); prediction = fixed.predict(basis, beta, len(y)); residual = prediction-y
    norms = np.array([np.linalg.norm(v) for _, v in basis]); scale = max(1., np.linalg.norm(y))
    require(np.all(norms > 0), "Zero design column.")
    u = norms*beta/scale
    g = np.array([v@residual[ids] for ids, v in basis])/norms/scale
    projected = u-g; projected[::2] = np.maximum(0, projected[::2])
    gradient = []
    for a, (atom, b) in enumerate(zip(atoms, widths)):
        ids, gaussian = basis[2*a]
        delta = points[ids]-np.asarray(atom["position"]); r = np.sqrt(np.sum(delta*delta, axis=1))
        center = math.sqrt(2/math.pi)/b
        dg = gaussian*(r*r/(b*b)-3)
        dk = np.where(r < 1e-5, -center, -center*np.exp(-r*r/(2*b*b)))
        gradient.append(float((beta[2*a]*dg+beta[2*a+1]*dk)@residual[ids]/scale/scale))
    return {"prediction": prediction, "residual": residual, "projected_kkt": float(np.max(np.abs(u-projected))),
            "feasible": bool(np.all(beta[::2] >= 0)), "b_gradient": np.asarray(gradient),
            "active_atoms": np.flatnonzero(beta[::2] == 0).tolist()}


def qualification(fit):
    if not all(fit[k]["valid"] for k in ("primary", "reference")): return False
    if "qualification_checks" not in fit: return False
    inner = all(fit[k]["feasible"] and fit[k]["kkt_passed"] and fit[k]["projected_kkt"] <= 1e-10
                for k in ("primary", "reference"))
    inner = inner and fit["scaled_reference_difference"] <= 1e-10 and fit["design_spectrum"]["rank"] == len(fit["primary"]["beta"])
    gradient = all(fit[k]["b_gradient_inf"] <= 1e-12 for k in ("primary", "reference"))
    local = fit["local_correction_inf"] <= 1e-10
    identified = fit["width_spectrum"]["rank"] == fit["profile_jacobian_spectrum"]["rank"] == len(fit["primary"]["b"])
    checks = fit["derivative_checks"]
    verified = len(checks) == 6 and all(c["passed"] and c["same_active_face"] and c["plus_valid"] and c["minus_valid"]
                   and c["relative_l2_difference"] is not None and c["relative_l2_difference"] <= 1e-6 for c in checks)
    require(fit["derivative_verified"] == verified and fit["qualification_checks"] == {
        "inner": inner, "b_gradient": gradient, "local_correction": local, "identified": identified, "derivative": verified},
        "Joint qualification checks differ.")
    return inner and gradient and local and identified and verified


def validate_fit(directory, name, fit, points, atoms, y, expected_initial, certificate_fn=raw_certificate, numeric_check=require):
    require(fit["case"] == name and fit["experiment"] == "joint-abc-profile" and fit["alpha"] == 0 and
            fit["execution_complete"] and fit["row_count"] == len(y) and fit["initial_b"] == expected_initial.tolist(),
            "Changed joint case/initialization contract.")
    require(fit["settings"] == SETTINGS and fit["linear_solver"] == "sparse-qr-householder-1024" and
            fit["reference_solver"] == "independent-tsqr-8192-svd", "Changed joint solver contract.")
    require(len(fit["trials"]) == fit["profile_evaluations"] <= 200 and fit["accepted_updates"] <= 100 and
            fit["endpoint_evaluations"] == 2 and fit["directional_evaluations"] in (0, 12), "Incorrect evaluation budget.")
    accepted = [t for t in fit["trials"] if t["accepted"]]
    require(len(accepted) == fit["accepted_updates"]+int(fit.get("initial_accepted", fit["initial"]["valid"])), "Incomplete accepted trajectory.")
    require(all(t["valid"] and t["kkt_passed"] for t in accepted), "Accepted failed inner solve.")
    require(all(accepted[k]["rss"] >= accepted[k+1]["rss"] for k in range(len(accepted)-1)), "Nondecreasing accepted RSS.")
    require(fit["joint_qualified"] == qualification(fit), "Unverified joint qualification.")
    if not fit["primary"]["valid"]: return None
    for label in ("initial", "primary", "reference"):
        endpoint = fit[label]
        if not endpoint["valid"]: continue
        independent = certificate_fn(points, atoms, y, endpoint)
        require(independent["feasible"] == endpoint["feasible"] and independent["active_atoms"] == endpoint["active_atoms"],
                "Independent feasibility/active constraints differ.")
        numeric_check(abs(independent["projected_kkt"]-endpoint["projected_kkt"]) <= 1e-13, "Independent A/C KKT differs.")
        numeric_check(np.allclose(independent["b_gradient"], endpoint["b_gradient"], rtol=2e-9, atol=1e-13) and
                abs(np.max(np.abs(endpoint["b_gradient"]))-endpoint["b_gradient_inf"]) <= 1e-16, "Independent B gradient differs.")
        if label == "primary": primary = independent
    beta = np.asarray(fit["primary"]["beta"]); reference = np.asarray(fit["reference"]["beta"])
    if fit["reference"]["valid"]:
        require(math.isclose(fixed.scaled_difference(beta, reference), fit["scaled_reference_difference"], rel_tol=1e-12, abs_tol=1e-16),
                "Reference coefficient difference differs.")
    residual_path = directory / "residuals" / f"{name}.csv"
    require(fixed.fold.sha256_file(residual_path) == fit["residual_sha256"], "Residual hash mismatch.")
    table = fixed.load_table(residual_path)
    require(len(table) == len(y) and np.array_equal(table["row"], np.arange(len(y))) and
            np.array_equal(table["residual"], table["prediction"]-y), "Residual row contract differs.")
    numeric_check(np.allclose(primary["prediction"], table["prediction"], rtol=2e-13, atol=2e-12), "Independent prediction differs.")
    rss = float(table["residual"]@table["residual"])
    numeric_check(math.isclose(rss, fit["primary"]["rss"], rel_tol=1e-9, abs_tol=1e-24), "Raw residual RSS differs.")
    if "width_spectrum" in fit:
        require(fit["width_spectrum"]["active_face_only"] == bool(fit["primary"]["active_atoms"]), "Incorrect active-face interpretation.")
        weak = fixed.load_table(directory / "weak-directions" / f"{name}.csv")
        require(np.array_equal(weak["serial_id"], [a["serial_id"] for a in atoms]) and
                np.array_equal(np.column_stack([weak[k] for k in ("weakest", "second", "third")]).T,
                               fit["width_spectrum"]["weak_directions"]), "Width directions differ.")
    return np.column_stack((beta[::2], fit["primary"]["b"], beta[1::2]))


def compare_starts(fits, abc):
    pairs = []; representatives = {}; consistent = True
    for precision in ("double", "float32"):
        eligible = [name for name in CASES if name.endswith(precision) and fits[name]["joint_qualified"]]
        representatives[precision] = min(eligible, key=lambda n: fits[n]["primary"]["rss"]) if eligible else None
        consistent &= len(eligible) == 4
        for left, right in itertools.combinations(eligible, 2):
            difference = fixed.scaled_difference(abc[left], abc[right])
            log_difference = float(np.max(np.abs(np.log(abc[left][:, 1])-np.log(abc[right][:, 1]))))
            passed = difference <= 1e-8 and log_difference <= 1e-8; consistent &= passed
            pairs.append({"left": left, "right": right, "scaled_abc_difference": difference,
                          "maximum_log_b_difference": log_difference, "passed": passed})
    return bool(consistent), pairs, representatives


def controls(directory, dataset, fits, abc):
    require(read(CONTROLS/"fixture.json") == read(fixed.FIXTURE), "Fixed-B controls use a different fixture.")
    require(read(CONTROLS/"dataset.json")["voxel_table_sha256"] == dataset["voxel_table_sha256"], "Fixed-B observations/ROI differ.")
    rows = []
    for name, values in abc.items():
        precision = name.rsplit("-", 1)[1]
        for prefix in ("true-b", "checkpoint-b"):
            control_name = f"{prefix}-{precision}"; path = CONTROLS/"fits"/f"{control_name}.json"; control = read(path)
            beta = np.asarray(control["primary"]["beta"])
            previous = np.column_stack((beta[::2], control["b"], beta[1::2]))
            rows.append({"case": name, "control": control_name, "control_sha256": fixed.fold.sha256_file(path),
                         "joint_rss": fits[name]["primary"]["rss"], "fixed_rss": control["primary"]["rss"],
                         "scaled_abc_difference": fixed.scaled_difference(values, previous)})
    if rows: fixed.write_csv(directory/"fixed-b-comparison.csv", rows)
    return rows


def summarize(directory):
    directory = directory.resolve(); clock = time.perf_counter()
    fixture = read(directory/"fixture.json"); paths = {k: Path(p) for k, p in read(directory/"inputs.json").items()}
    require(fixture == read(fixed.FIXTURE), "Changed frozen fixture.")
    manifest, checkpoint = fixed.validate_inputs(paths, fixture)
    require(all(fixed.fold.sha256_file(Path(p)) == h for p, h in read(directory/"input-hashes.json").items()), "Changed input hashes.")
    completion = read(directory/"completion.json")
    require(completion["complete"] and completion["execution_complete"] and completion["cases"] == CASES, "Missing completed joint case.")
    dataset = read(directory/"dataset.json"); forward = read(directory/"forward-status.json")
    require(forward["passed"] and forward["generator_jobs"] == 1 and forward["voxel_count"] == fixture["expected_rows"] and
            dataset["checkpoint_b"] == checkpoint, "Invalid forward/B evidence.")
    points, y64, y32 = fixed.validate_dataset(directory, dataset, fixture, paths, manifest, "joint-abc-profile")
    atoms = manifest["atoms"]; truth = np.array([[a["element"], .5, a["charge_used"]] for a in atoms])
    truth_beta = truth[:, [0, 2]].ravel()
    require(read(directory/"scoring-truth.json") == {"beta": truth_beta.tolist(), "width": .5}, "Truth mismatch.")
    truth_basis = fixed.columns(points, atoms, [.5]*len(atoms))
    absolute = fixed.predict([(i, np.abs(v)) for i, v in truth_basis], np.abs(truth_beta), len(y64))
    require(np.all(np.abs(fixed.predict(truth_basis, truth_beta, len(y64))-y64) <=
                   512*np.finfo(float).eps*np.maximum(1, absolute)), "Independent double forward mismatch.")
    fits = {}; abc = {}; summaries = []; estimates = []
    for name in CASES:
        start, precision = name.rsplit("-", 1); y = y64 if precision == "double" else y32
        fit = read(directory/"fits"/f"{name}.json"); fits[name] = fit
        values = validate_fit(directory, name, fit, points, atoms, y, initial_widths(checkpoint, atoms, start))
        summary = {"case": name, "execution_complete": fit["execution_complete"], "joint_qualified": fit["joint_qualified"],
                   "oracle_recovered": False, "float32_accuracy_passed": False, "qualification_failure": fit["qualification_failure"]}
        if values is not None:
            abc[name] = values; error = values-truth; scaled = fixed.scaled_difference(values, truth)
            summary.update({key: fixed.statistics(error[:, k]) for k, key in enumerate(("A", "B", "C"))})
            summary.update({"scaled_truth_error": scaled, "oracle_recovered": precision == "double" and scaled <= 1e-10 and fit["primary"]["relative_residual"] <= 1e-12,
                            "float32_accuracy_passed": precision == "float32" and bool(np.all(np.abs(error) < .01)),
                            **{k: fit["primary"][k] for k in ("rss", "residual_rmse", "residual_max", "relative_residual", "projected_kkt", "b_gradient_inf")},
                            **{k: fit.get(k) for k in ("scaled_reference_difference", "local_correction_inf", "profile_evaluations", "accepted_updates")}})
            for a, atom in enumerate(atoms):
                estimates.append({"case": name, "serial_id": atom["serial_id"], **{key: values[a, k] for k, key in enumerate(("A", "B", "C"))},
                                  **{f"{key}_error": error[a, k] for k, key in enumerate(("A", "B", "C"))}})
        summaries.append(summary)
    consistent, pairs, representatives = compare_starts(fits, abc)
    qualified = all(f["joint_qualified"] for f in fits.values())
    recovered = all(s["oracle_recovered"] for s in summaries if s["case"].endswith("double"))
    accurate = all(s["float32_accuracy_passed"] for s in summaries if s["case"].endswith("float32"))
    precision_difference = None
    if all(representatives.values()):
        delta = abc[representatives["float32"]]-abc[representatives["double"]]
        precision_difference = {key: fixed.statistics(delta[:, k]) for k, key in enumerate(("A", "B", "C"))}
    results = {"experiment": "joint-abc-profile", "rows": len(y64), "execution_complete": True,
               "joint_qualified": qualified, "oracle_recovered": recovered, "float32_accuracy_passed": accurate,
               "multistart_consistent": consistent, "overall_passed": qualified and recovered and accurate and consistent,
               "representatives": representatives, "cases": summaries, "pairwise_comparisons": pairs,
               "float32_minus_double": precision_difference, "fixed_b_comparisons": controls(directory, dataset, fits, abc)}
    fixed.experiment.write(directory/"results.json", results)
    if estimates: fixed.write_csv(directory/"estimates.csv", estimates)
    # Flatten common summary fields; failed cases remain present with empty diagnostics.
    fields = sorted(set().union(*(s.keys() for s in summaries))-{"A", "B", "C"})
    flat = [{**{k: row.get(k) for k in fields}, **{f"{parameter}_{k}": row.get(parameter, {}).get(k)
              for parameter in ("A", "B", "C") for k in ("bias", "rmse", "p99", "max", "below_0_01")}} for row in summaries]
    fixed.write_csv(directory/"comparison.csv", flat)
    if pairs: fixed.write_csv(directory/"multistart.csv", pairs)
    fixed.experiment.write(directory/"independent-validation.json", {"passed": True, "seconds": time.perf_counter()-clock,
        "checks": ["complete original ROI and frozen geometry", "generator/analytic forward and float32 map bytes", "checkpoint identity and initialization",
                   "original-row prediction, A/C KKT and B gradient", "raw residuals and per-atom statistics", "joint qualification and multistart comparisons"]})
    raw = {str(p.relative_to(directory)): fixed.fold.sha256_file(p) for p in sorted(directory.rglob("*"))
           if p.is_file() and p.name not in ("raw-artifact-index.json", "reproducibility.json")}
    fixed.experiment.write(directory/"raw-artifact-index.json", raw)
    print({k: results[k] for k in ("execution_complete", "joint_qualified", "oracle_recovered", "multistart_consistent", "overall_passed")})
    return results


def scientific(value):
    if isinstance(value, list): return [scientific(v) for v in value]
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items() if not k.endswith("seconds") and k != "process_peak_rss_bytes"}
    return value


def compare(left, right, output):
    for directory in (left, right): summarize(directory)
    names = ["dataset.json", "forward-status.json", "results.json", "completion.json", "scoring-truth.json", "fixture.json",
             "provenance.json", "input-hashes.json", *[f"fits/{name}.json" for name in CASES]]
    differences = [name for name in names if scientific(read(left/name)) != scientific(read(right/name))]
    csv_files = sorted(str(p.relative_to(left)) for p in left.rglob("*.csv"))
    other_csv_files = sorted(str(p.relative_to(right)) for p in right.rglob("*.csv"))
    if csv_files != other_csv_files: differences.append("CSV inventory")
    for name in set(csv_files) & set(other_csv_files):
        if fixed.fold.sha256_file(left/name) != fixed.fold.sha256_file(right/name): differences.append(name)
    fixed.experiment.write(output, {"passed": not differences, "cases": CASES, "differences": differences,
        "excluded": ["timing", "resource usage", "execution/output directory"], "comparison": "exact scientific JSON and CSV bytes"})
    require(not differences, f"Independent repeat differs: {differences}")


def plots(directory, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output.mkdir(parents=True, exist_ok=True); table = fixed.load_table(directory/"estimates.csv")
    figure, axes = plt.subplots(3, 2, figsize=(12, 9), constrained_layout=True)
    for row, key in enumerate(("A_error", "B_error", "C_error")):
        for col, precision in enumerate(("double", "float32")):
            ax = axes[row, col]
            for start in STARTS:
                selected = table[table["case"] == f"{start}-{precision}"]
                ax.plot(selected["serial_id"], selected[key], label=start, linewidth=.8, alpha=.8)
            ax.set(xlabel="Atom serial ID", ylabel=key.replace("_", " "), title=precision)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0)); ax.grid(alpha=.2)
            if row == 0: ax.legend(fontsize=8)
    for suffix in ("png", "pdf"): figure.savefig(output/f"parameter-errors.{suffix}", dpi=170)
    plt.close(figure)
    figure, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for col, precision in enumerate(("double", "float32")):
        for start in STARTS:
            fit = read(directory/"fits"/f"{start}-{precision}.json")
            accepted = [t for t in fit["trials"] if t["accepted"]]
            axes[0, col].semilogy([t["accepted_update"] for t in accepted], [max(1e-32, t["rss"]) for t in accepted], "o-", markersize=3, label=start)
            if "width_spectrum" in fit:
                axes[1, col].semilogy(range(1, 169), fit["width_spectrum"]["singular_values"], label=start)
        axes[0, col].set(xlabel="Accepted update", ylabel="RSS", title=precision); axes[0, col].legend(fontsize=8)
        axes[1, col].set(xlabel="Singular value index", ylabel="Projected width sensitivity / observation scale", title=precision)
        for ax in axes[:, col]: ax.grid(alpha=.2)
    for suffix in ("png", "pdf"): figure.savefig(output/f"convergence-spectrum.{suffix}", dpi=170)
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
    if args.command == "run": fixed.run(args, "joint-abc-profile", summarize)
    elif args.command == "summarize": summarize(args.directory)
    elif args.command == "compare": compare(args.left, args.right, args.output)
    else: plots(args.directory, args.output)


if __name__ == "__main__": main()
