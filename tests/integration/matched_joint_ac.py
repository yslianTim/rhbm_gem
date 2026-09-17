#!/usr/bin/env python3
"""Fixed-checkpoint-width matched component A/C MDPDE experiment (testing only)."""
from __future__ import annotations

import argparse
import collections
import csv
import os
from pathlib import Path
import platform
import struct
import subprocess
import sys
import time

import estimated_neighbor_sweep as sweep
import fold_168_regression as fold
import mdpde_experiment as experiment
from observation_matching import read, require

ROOT = Path(__file__).resolve().parents[2]
HISTORY = ROOT / "docs/developer/figures/estimated-neighbor-sweep/fits.csv"


def run(args):
    output = args.output.resolve()
    require(not output.exists(), "Use a fresh output directory; historical evidence is immutable.")
    inputs = {k: getattr(args, k).resolve() for k in ("model", "map", "manifest")}
    baseline = fold.load_baseline(ROOT / "tests/benchmarks/fold_168_simulation_baseline.json")
    hashes = fold.validate_input_hashes(inputs, baseline["input_hashes"])
    manifest = fold.load_simulation_manifest(inputs["manifest"], hashes)
    fold.validate_fixture(manifest, baseline)
    index, evidence = sweep.prepare_index(args, manifest, hashes)
    all_ids = [state["id"] for state in index["states"]]
    require(args.state is None or args.state in all_ids, "Unknown selected state.")
    index["selected_state_ids"] = [args.state] if args.state else all_ids
    mappings = [[atom["alpha"] for atom in read(Path(state["context"]))["atoms"]] for state in index["states"]]
    require(all(all(sweep.finite(a) and a >= 0 for a in values) and values == mappings[0] for values in mappings),
            "Invalid or changed checkpoint alpha mapping.")
    before = experiment.provenance(args.executable.resolve(), args.source.resolve())
    input_hashes = {str(p): fold.sha256_file(p) for p in list(inputs.values()) + evidence + [HISTORY]}
    output.mkdir(parents=True)
    experiment.write(output / "state-index.json", index)
    experiment.write(output / "input-hashes.json", input_hashes)
    experiment.write(output / "provenance.json", before)
    experiment.write(output / "environment.json", {"python": sys.version, "platform": platform.platform(),
        "jobs": args.jobs, "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=args.source, text=True).strip()})
    experiment.write(output / "scoring-truth.json", {str(a["serial_id"]):
        [float(a["element"]), manifest["settings"]["blurring_width"], a["charge_used"]] for a in manifest["atoms"]})
    # Retain the exact historical comparator, separate from every fit input.
    (output / "historical-fits.csv").write_bytes(HISTORY.read_bytes())
    command = [str(args.executable.resolve()), "matched-joint-ac", str(inputs["manifest"]), str(inputs["map"]),
               str(output / "state-index.json"), str(output)]
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(args.jobs), VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    start = time.perf_counter()
    with (output / "run.log").open("w") as stream:
        result = subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT)
    experiment.write(output / "execution.json", {"command": command, "returncode": result.returncode,
        "seconds": time.perf_counter()-start, "jobs": args.jobs})
    require(before == experiment.provenance(args.executable.resolve(), args.source.resolve()), "Source or binary changed during run.")
    require(all(fold.sha256_file(Path(p)) == h for p, h in input_hashes.items()), "Experiment input changed.")
    require(result.returncode == 0, f"Joint AC run failed; inspect {output / 'run.log'}.")
    summarize(output)


def csv_write(path, rows):
    require(bool(rows), f"No rows for {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def pair_status(left, right, truth, left_uncertainty, right_uncertainty, qualified):
    if not all(sweep.finite(x) for x in (left, right, truth, left_uncertainty, right_uncertainty)):
        return "unavailable"
    if not qualified:
        return "unqualified"
    gain = abs(left-truth)-abs(right-truth)
    uncertainty = left_uncertainty+right_uncertainty
    return "improved" if gain > uncertainty else "regressed" if gain < -uncertainty else "unresolved"


def validate_case(case, dataset, fits, context):
    members = case["members"]
    component = case["component"]
    require(members == dataset["component_atoms"][component], "Component membership mismatch.")
    require(case["row_count"] == len(dataset["component_rows"][component]), "Sample population mismatch.")
    require(len(fits) == len(members)+1, "Missing frozen or joint fits.")
    require(fits[0]["mode"] == "joint" and fits[0]["columns"] == 2*len(members), "Invalid joint design.")
    expected_blocks = [{"owner": owner, "alpha": context["atoms"][owner]["alpha"],
                        "rows": sum(dataset["rows"][r]["owner"] == owner for r in dataset["component_rows"][component])}
                       for owner in members]
    expected_blocks = [b for b in expected_blocks if b["rows"]]
    for block in expected_blocks: block["lambda"] = 1 / len(expected_blocks)
    require(case["blocks"] == expected_blocks, "Block identity or alpha mismatch.")
    for k, fit in enumerate(fits):
        require(fit["blocks"] == case["blocks"] and fit["state_id"] == case["state_id"]
                and fit["component"] == component and fit["rows"] == case["row_count"], "Mismatched fit domain.")
        if k:
            require(fit["mode"] == "frozen" and fit["atom_index"] == members[k-1] and fit["columns"] == 2,
                    "Frozen target mismatch.")
        require(fit.get("schema_version") == 2, "Unsupported fit schema.")
        require(len(fit.get("variances", [])) == len(expected_blocks) or not fit["qualified"], "Missing block scales.")
        require(len(fit["beta"]) == fit["columns"], "Missing fitted coefficients.")
        require(len(fit["uncertainty"]) == fit["columns"] or not fit["qualified"], "Missing qualified numerical evidence.")
    require([r["mode"] for r in case["assembled"]] == ["input", "frozen", "joint"], "Missing assembled state.")
    for row in case["assembled"]:
        require(len(row["abc"]) == len(members), "Missing assembled member.")
        for k, a in enumerate(members):
            require(struct.pack("d", row["abc"][k][1]) == struct.pack("d", context["state"][a][1]), "Width changed.")
            if row["mode"] == "input":
                require(row["abc"][k] == context["state"][a], "Input checkpoint changed.")
            else:
                fit, offset = (fits[k+1], 0) if row["mode"] == "frozen" else (fits[0], 2*k)
                require([row["abc"][k][0], row["abc"][k][2]] == fit["beta"][offset:offset+2], "Assembly differs from fit.")


def sample_statistics(path, case, dataset):
    accum = collections.defaultdict(list)
    seen = set()
    expected = set(dataset["component_rows"][case["component"]])
    with path.open() as stream:
        for r in csv.DictReader(stream):
            index = int(r["row"])
            require(index not in seen and index in expected, "Invalid sample membership.")
            seen.add(index)
            metadata = dataset["rows"][index]
            require(int(r["owner"]) == metadata["owner"] and int(r["sample"]) == metadata["sample"]
                    and float(r["observed"]) == metadata["response"], "Sample identity mismatch.")
            groups = ("all", "signal" if int(r["signal"]) else "tail",
                      "cutoff-crossing" if int(r["cutoff_crossing"]) else "smooth",
                      "boundary" if int(r["boundary"]) else "interior")
            for mode in ("input", "frozen", "joint"):
                for metric in ("residual", "neighbor_error"):
                    value = float(r[f"{mode}_{metric}"])
                    for group in groups: accum[(mode, metric, group)].append(value)
    require(seen == expected, "Incomplete samples.")
    return [{"state_id": case["state_id"], "component": case["component"],
             "mode": mode, "metric": metric, "group": group, **sweep.population_stats(values, len(values))}
            for (mode, metric, group), values in sorted(accum.items())]


def summarize(directory):
    directory = directory.resolve()
    index, dataset = read(directory / "fit-index.json"), read(directory / "dataset.json")
    states = read(directory / "state-index.json")["states"]
    require(index["states"] == states and len(states) == 8, "State index mismatch.")
    require(index.get("schema_version") == 2, "Unsupported experiment schema.")
    selected_ids = read(directory / "state-index.json")["selected_state_ids"]
    require(index["selected_state_ids"] == selected_ids and selected_ids and len(set(selected_ids)) == len(selected_ids)
            and set(selected_ids) <= {s["id"] for s in states}, "Selected state index mismatch.")
    require(len(dataset["atoms"]) == 168 and len(dataset["rows"]) == 33600, "Incomplete population.")
    all_atoms = [a for group in dataset["component_atoms"] for a in group]
    all_rows = [a for group in dataset["component_rows"] for a in group]
    require(sorted(all_atoms) == list(range(168)) and sorted(all_rows) == list(range(33600)), "Invalid component partition.")
    for atoms, rows in zip(dataset["component_atoms"], dataset["component_rows"]):
        for p in rows:
            require(dataset["rows"][p]["owner"] in atoms and set(dataset["rows"][p]["contributors"]) <= set(atoms), "Cross-component dependency.")
    expected = {f"{s}-{c}" for s in range(len(states)) if states[s]["id"] in selected_ids
                for c in range(len(dataset["component_atoms"]))}
    require(set(index["cases"]) == expected and len(index["cases"]) == len(expected), "Incomplete experiment matrix.")
    require(index["joint_fits"] == len(expected) and index["frozen_fits"] == 168*len(selected_ids), "Incorrect fit counts.")
    truth = read(directory / "scoring-truth.json")
    with (directory / "historical-fits.csv").open() as stream:
        historical = {(r["state_id"], int(r["serial_id"])): r for r in csv.DictReader(stream)
                      if r["parameters"] == "ABC" and r["prediction"] == "matched"}
    contexts = {s["id"]: read(Path(s["context"])) for s in states}
    require(all([a["alpha"] for a in c["atoms"]] == dataset["alphas"] for c in contexts.values()), "Dataset alpha mapping changed.")
    estimates, pairs, fit_rows, sample_stats, cross_rows, common, block_rows = [], [], [], [], [], [], []
    for name in index["cases"]:
        case = read(directory / f"case-{name}.json")
        fits = [read(directory / "fits" / f"{name}-{k}.json") for k in range(len(case["members"])+1)]
        validate_case(case, dataset, fits, contexts[case["state_id"]])
        state_id = case["state_id"]
        for fit in fits:
            branches = fit["branches"]
            selected = branches[fit["selected_seed"]] if branches else {}
            weights = selected.get("weights", {})
            for branch in branches:
                for d in branch.get("block_diagnostics", []):
                    block_rows.append({"state_id": state_id, "component": case["component"], "mode": fit["mode"],
                        "atom_index": fit.get("atom_index"), "seed": branch["seed"], "qualified": branch["qualified"],
                        "selected": branch is selected, **d})
            fit_rows.append({"state_id": state_id, "component": case["component"],
                "mode": fit["mode"], "atom_index": fit.get("atom_index"), "qualified": fit["qualified"],
                "reason": fit["reason"], "failure_owner": fit.get("failure_owner"), "branch_sensitive": fit.get("branch_sensitive", False),
                "stationarity": fit.get("stationarity"), "objective": fit.get("objective"),
                "active_amplitudes": fit.get("active_amplitudes"),
                "rank": fit.get("design_spectrum", {}).get("rank"),
                "condition": fit.get("design_spectrum", {}).get("condition"),
                "minimum_singular": fit.get("design_spectrum", {}).get("minimum_singular"),
                "weighted_rank": selected.get("weighted_spectrum", {}).get("rank"),
                "weighted_condition": selected.get("weighted_spectrum", {}).get("condition"),
                "weighted_minimum_singular": selected.get("weighted_spectrum", {}).get("minimum_singular"),
                "weight_minimum": weights.get("minimum"), "weight_mean": weights.get("mean"),
                "weight_underflow_count": weights.get("underflow_count"), "effective_n": weights.get("effective_n"),
                "linear_solves": fit["initial_linear_solves"] + sum(b["primary"]["linear_solves"]
                    + b["reference"]["linear_solves"] + b.get("svd_linear_solves", 0) for b in branches),
                "iterations": sum(b["primary"]["iterations"]+b["reference"]["iterations"] for b in branches),
                "seconds": fit.get("seconds", 0)})
        assembled = {r["mode"]: r for r in case["assembled"]}
        for mode, row in assembled.items():
            common.append({"state_id": state_id, "component": case["component"], "mode": mode,
                "common_scale_status": row["common_scale_status"], "common_scale_objective": row["common_scale_objective"],
                "residual_rmse": row["residual_rmse"]})
        for k, atom_index in enumerate(case["members"]):
            identity = dataset["atoms"][atom_index]; serial = identity["serial_id"]
            t = truth[str(serial)]
            uncertainties = {"input": [0, 0], "joint": fits[0]["uncertainty"][2*k:2*k+2], "frozen": fits[k+1]["uncertainty"]}
            qualified = {"input": True, "joint": fits[0]["qualified"], "frozen": fits[k+1]["qualified"]}
            for mode in ("input", "frozen", "joint"):
                for p, parameter in enumerate(("A", "C")):
                    i = 2*p
                    value = assembled[mode]["abc"][k][i]
                    uncertainty = uncertainties[mode][p] if len(uncertainties[mode]) == 2 else None
                    estimates.append({"state_id": state_id, "component": case["component"],
                        "serial_id": serial, "chain_id": identity["chain_id"], "sequence_id": identity["sequence_id"],
                        "component_id": identity["component_id"], "atom_id": identity["atom_id"],
                        "mode": mode, "parameter": parameter, "value": value, "truth": t[i], "error": value-t[i],
                        "uncertainty": uncertainty, "qualified": qualified[mode],
                        "B": assembled[mode]["abc"][k][1], "B_error": assembled[mode]["abc"][k][1]-t[1]})
            statuses = {}
            for left, right in (("input", "frozen"), ("input", "joint"), ("frozen", "joint")):
                for p, parameter in enumerate(("A", "C")):
                    lv, rv = assembled[left]["abc"][k][2*p], assembled[right]["abc"][k][2*p]
                    lu = uncertainties[left][p] if len(uncertainties[left]) == 2 else None
                    ru = uncertainties[right][p] if len(uncertainties[right]) == 2 else None
                    status = pair_status(lv, rv, t[2*p], lu, ru, qualified[left] and qualified[right])
                    statuses[(right, parameter)] = status if left == "input" else statuses.get((right, parameter))
                    pairs.append({"state_id": state_id, "serial_id": serial, "parameter": parameter,
                        "comparison": f"{left}-to-{right}", "left_error": lv-t[2*p], "right_error": rv-t[2*p],
                        "absolute_error_gain": abs(lv-t[2*p])-abs(rv-t[2*p]),
                        "numerical_uncertainty": lu+ru if sweep.finite(lu) and sweep.finite(ru) else None, "status": status})
            old = historical[(state_id, serial)]
            old_status = pair_status(float(old["input_C"]), float(old["C"]), t[2], 0,
                                    abs(float(old["C"])-float(old["reference_C"])), old["qualified"] == "True")
            cross_rows.append({"state_id": state_id, "serial_id": serial,
                "historical_C_status": old_status, "frozen_C_status": statuses[("frozen", "C")],
                "joint_C_status": statuses[("joint", "C")]})
        sample_stats.extend(sample_statistics(directory / "samples" / f"{name}.csv", case, dataset))
    grouped = collections.defaultdict(list)
    for row in estimates: grouped[(row["state_id"], row["mode"], row["parameter"])].append(row)
    statistics = [{"state_id": s, "mode": m, "parameter": p,
        "qualified_atoms": sum(r["qualified"] for r in rows),
        **sweep.population_stats([r["error"] for r in rows], 168)} for (s, m, p), rows in grouped.items()]
    counts = collections.defaultdict(collections.Counter)
    for row in pairs: counts[(row["state_id"], row["parameter"], row["comparison"])][row["status"]] += 1
    cross_counts = collections.defaultdict(collections.Counter)
    for row in cross_rows:
        if row["historical_C_status"] == "regressed":
            cross_counts[row["state_id"]][row["joint_C_status"]] += 1
    worst = []
    for s in states:
        for parameter in ("A", "C"):
            candidates = [r for r in pairs if r["state_id"] == s["id"] and r["parameter"] == parameter
                          and r["comparison"] == "input-to-joint"]
            worst.extend(sorted(candidates, key=lambda r: r["absolute_error_gain"])[:5])
    result = {"schema_version": 2, "states": [s for s in states if s["id"] in selected_ids],
        "selected_state_ids": selected_ids, "alphas": dataset["alphas"], "alpha_source": "checkpoint-owner",
        "component_sizes": [len(c) for c in dataset["component_atoms"]], "widths_unchanged": True,
        "joint_fits": index["joint_fits"], "frozen_fits": index["frozen_fits"],
        "fit_status": {mode: dict(collections.Counter(r["reason"] for r in fit_rows if r["mode"] == mode)) for mode in ("joint", "frozen")},
        "statistics": statistics, "sample_statistics": sample_stats, "common_objectives": common,
        "pair_counts": [dict(zip(("state_id", "parameter", "comparison"), key), counts=dict(value)) for key, value in counts.items()],
        "historical_regression_cross_counts": [{"state_id": k, "joint_status": dict(v)} for k, v in cross_counts.items()],
        "worst_cases": worst, "serial_100": [r for r in pairs if r["serial_id"] == 100],
        "maximum_design_difference": max(read(directory / f"case-{name}.json")["maximum_design_difference"] for name in index["cases"]),
        "fit_work": {"linear_solves": sum(r["linear_solves"] for r in fit_rows), "iterations": sum(r["iterations"] for r in fit_rows),
                     "worker_seconds": sum(r["seconds"] for r in fit_rows)}}
    experiment.write(directory / "results.json", result)
    for name, rows in (("estimates", estimates), ("pairs", pairs), ("fits", fit_rows), ("historical-cross", cross_rows), ("blocks", block_rows)):
        if rows: csv_write(directory / f"{name}.csv", rows)
    sweep.artifact_index(directory)
    return result


def scientific(value):
    if isinstance(value, dict): return {k: scientific(v) for k, v in value.items() if k not in ("seconds", "worker_seconds")}
    if isinstance(value, list): return [scientific(v) for v in value]
    return value


def compare(left, right, output, state):
    li, ri = read(left / "fit-index.json"), read(right / "fit-index.json")
    require(state in li["selected_state_ids"] and state in ri["selected_state_ids"], "Requested state was not executed in both runs.")
    require(li["states"] == ri["states"], "Source state index mismatch.")
    names = ["dataset.json", "forward-status.json", "input-hashes.json", "provenance.json"]
    lc = [name for name in li["cases"] if read(left / f"case-{name}.json")["state_id"] == state]
    rc = [name for name in ri["cases"] if read(right / f"case-{name}.json")["state_id"] == state]
    require(lc == rc and lc, "Requested state's component coverage differs.")
    for name in lc:
        names.append(f"case-{name}.json")
        names.extend(f"fits/{name}-{k}.json" for k in range(len(read(left / f"case-{name}.json")["members"])+1))
    differences = [name for name in names if not (right / name).is_file() or scientific(read(left / name)) != scientific(read(right / name))]
    for name in lc:
        if fold.sha256_file(left / "samples" / f"{name}.csv") != fold.sha256_file(right / "samples" / f"{name}.csv"):
            differences.append(f"samples/{name}.csv")
    lr, rr = read(left / "results.json"), read(right / "results.json")
    for key in ("statistics", "sample_statistics", "common_objectives", "pair_counts", "historical_regression_cross_counts", "worst_cases", "serial_100"):
        if scientific([r for r in lr[key] if r["state_id"] == state]) != scientific([r for r in rr[key] if r["state_id"] == state]):
            differences.append(f"results/{key}")
    result = {"passed": not differences, "state": state, "compared_states": [state],
              "not_compared_states": [s["id"] for s in li["states"] if s["id"] != state],
              "json_files": len(names), "differences": differences,
              "left": str(left.resolve()), "right": str(right.resolve())}
    experiment.write(output, result)
    require(result["passed"], f"Thread comparison failed: {differences[:5]}")


def plots(directory, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    result = read(directory / "results.json")
    output.mkdir(parents=True, exist_ok=True)
    states = [s["id"] for s in result["states"]]
    labels = [s.replace("baseline-", "base ").replace("failed-only-", "refined ") for s in states]
    colors = {"input": "#777777", "frozen": "#b56b45", "joint": "#246b86"}
    def save(fig, name):
        for extension in ("png", "pdf"): fig.savefig(output / f"{name}.{extension}", dpi=180)
        plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    for ax, p in zip(axes, ("A", "C")):
        for offset, (mode, color) in enumerate(colors.items()):
            values = [next(r["rmse"] for r in result["statistics"] if r["state_id"] == s
                           and r["mode"] == mode and r["parameter"] == p) for s in states]
            ax.bar([i+(offset-1)*.25 for i in range(len(states))], values, .25, label=mode, color=color)
        ax.set_xticks(range(len(states)), labels, rotation=55, ha="right"); ax.set_title(f"{p} RMSE"); ax.set_ylim(bottom=0); ax.grid(axis="y", alpha=.2)
    axes[0].legend(); fig.suptitle("Fixed B · owner alpha · all endpoints (including unqualified)")
    save(fig, "parameter-errors")
    with (directory / "blocks.csv").open() as stream:
        blocks = [r for r in csv.DictReader(stream) if r["mode"] == "joint" and r["selected"] == "True"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    for state, label in zip(states, labels):
        rows = [r for r in blocks if r["state_id"] == state]
        for ax, key in zip(axes, ("variance", "effective_n", "linear_weight_share")):
            ax.scatter([float(r["alpha"]) for r in rows], [float(r[key]) for r in rows], s=12, alpha=.5, label=label)
            ax.set_xlabel("Checkpoint owner alpha"); ax.set_title(key.replace("_", " ")); ax.grid(alpha=.2)
    axes[0].set_yscale("log"); axes[2].set_yscale("log"); axes[1].legend(fontsize=6)
    fig.suptitle("Joint endpoint block diagnostics (qualification reported separately)")
    save(fig, "block-diagnostics")
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), constrained_layout=True, sharey=True)
    for ax, (left, right) in zip(axes, (("input", "frozen"), ("input", "joint"), ("frozen", "joint"))):
        counts = [next(r["counts"] for r in result["pair_counts"] if r["state_id"] == s
                       and r["parameter"] == "C" and r["comparison"] == f"{left}-to-{right}") for s in states]
        bottom = [0]*len(states)
        for status, color in (("improved", "#246b86"), ("regressed", "#b56b45"), ("unresolved", "#bbb"), ("unqualified", "#80689d"), ("unavailable", "#333")):
            values = [c.get(status, 0) for c in counts]; ax.bar(range(len(states)), values, bottom=bottom, label=status, color=color)
            bottom = [a+b for a, b in zip(bottom, values)]
        ax.set_xticks(range(len(states)), labels, rotation=55, ha="right", fontsize=8)
        ax.set_title(f"{left.title()} → {right.title()}"); ax.set_ylim(0, 220); ax.set_yticks((0, 56, 112, 168))
    axes[0].set_ylabel("Atoms (full denominator = 168)"); axes[0].legend(fontsize=7, ncol=2, loc="upper left")
    fig.suptitle("C paired classifications · fixed checkpoint B · owner alpha"); save(fig, "charge-pairs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    runner = sub.add_parser("run")
    for name in ("executable", "model", "map", "manifest", "baseline-run", "refined-run", "output"):
        runner.add_argument("--"+name, type=Path, required=True)
    runner.add_argument("--source", type=Path, default=ROOT); runner.add_argument("--jobs", type=int, choices=(1, 4), default=4)
    runner.add_argument("--state", help="Execute one retained state; default: all eight.")
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


if __name__ == "__main__":
    main()
