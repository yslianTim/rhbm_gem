#!/usr/bin/env python3
"""Frozen estimated-neighbor paired sweeps on eight captured fold-168 states."""
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
from observation_matching import read, require

ROOT = Path(__file__).resolve().parents[2]
PARAMETERS = ("A", "B", "C")
IDENTITY = ("serial_id", "chain_id", "sequence_id", "component_id", "atom_id", "alternate_indicator", "position")
# These are the actual, previously retained contexts, not newly selected endpoints.
CHECKPOINTS = (
    ("baseline", "final", 32, 28, "final-32-33", "d8bc49e3d9346e7ea35106fecbba0e722c08e023989efa253d24870b3c6660c0"),
    ("baseline", "recovery-current", 32, None, "recovery-current-32-32", "15bbb9c42d42613c0a7b22642bd961eaf413969c52a4489ec6eb71b8badaf384"),
    ("failed-only", "final", 14, 6, "final-14-19", "ee545f1635159351658604ad35f20003c9c5f72cc15951a9c8cae67d51ac8788"),
    ("failed-only", "recovery-current", 10, None, "recovery-current-10-10", "e8bd12767c989b59dc68714d2917c4ab33bd3cd2fabb10b2f6850dbb528bfdef"),
    ("failed-only", "recovery-current", 11, None, "recovery-current-11-12", "10edf32b32c26dcc0a789da2408df254691a2c36c966904d862e9646b3a5a55e"),
    ("failed-only", "recovery-current", 12, None, "recovery-current-12-14", "54f21b0d77985c26e616cbf34dff09145974a597d485b79b8a96f5c1db168029"),
    ("failed-only", "recovery-current", 13, None, "recovery-current-13-16", "d09fcf311b4ecd217fe1fa85db8b884edac69d60f651d770d59ab6f1594228b6"),
    ("failed-only", "recovery-current", 14, None, "recovery-current-14-18", "20f5b3eb5cc600d17f23e324970fcc6048ae1f37660e8a92818f3211cca65fe8"),
)


def state_id(policy, phase, attempt, best):
    return f"{policy}-best-{best}" if phase == "final" else f"{policy}-recovery-{attempt}"


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def validate_context(context, manifest, reference=None):
    require(context["schema_version"] == 1 and len(context["state"]) == len(context["atoms"]) == 168,
            "Incomplete checkpoint atom population.")
    truth = {a["serial_id"]: a for a in manifest["atoms"]}
    seen = set()
    for i, (atom, abc) in enumerate(zip(context["atoms"], context["state"])):
        identity = atom["identity"]
        serial = identity["serial_id"]
        require(serial in truth and serial not in seen and atom["index"] == i, "Missing or duplicate atom identity.")
        seen.add(serial)
        require(all(identity[k] == truth[serial][k] for k in IDENTITY), "Atom identity mismatch.")
        require(len(abc) == 3 and all(finite(v) for v in abc) and abc[0] > 0 and abc[1] > 0,
                "Nonfinite or invalid checkpoint parameters.")
        samples = atom["samples"]
        require(len(samples) == 200 and sum(s["distance"] <= 1 for s in samples) == 100, "Changed sample membership.")
        for sample in samples:
            require(finite(sample["response"]) and finite(sample["distance"]) and sample["distance"] >= 0
                    and len(sample["position"]) == 3 and all(finite(v) for v in sample["position"])
                    and isinstance(sample["selected"], bool), "Nonfinite or invalid sample.")
        if reference is not None:
            other = reference["atoms"][i]
            require(identity == other["identity"] and samples == other["samples"], "State samples or identities differ.")
    require(seen == set(truth), "Checkpoint and manifest populations differ.")


def prepare_index(args, manifest, input_hashes):
    runs = {"baseline": args.baseline_run.resolve(), "failed-only": args.refined_run.resolve()}
    states, evidence, reference = [], [], None
    for policy, phase, attempt, best, name, expected_hash in CHECKPOINTS:
        run = runs[policy]
        path = run / "solver-failures" / "contexts" / f"{name}.json"
        require(path.is_file(), f"Missing checkpoint: {path}")
        require(fold.sha256_file(path) == expected_hash, f"Checkpoint hash differs from retained evidence: {path}")
        context = read(path)
        require(context["phase"] == phase and context["attempt"] == attempt and context["operator_id"] == name,
                f"Checkpoint role mismatch: {path}")
        validate_context(context, manifest, reference)
        if reference is None:
            reference = context
        actual_path = run / "actual.json"
        actual = read(actual_path)
        require(actual["input_hashes"] == input_hashes, f"Source run used different simulation inputs: {run}")
        if phase == "final":
            require(str(actual["second_stage_summary"]["best_iteration"]) == str(best)
                    and actual["second_stage_summary"]["final_state_source"] == "best-audit", "Unexpected saved best state.")
            saved = {a["serial_id"]: a for a in actual["atoms"]}
            require(len(saved) == 168, "Incomplete saved final results.")
            for atom, abc in zip(context["atoms"], context["state"]):
                row = saved[atom["identity"]["serial_id"]]
                require(all(row[k] == atom["identity"][k] for k in IDENTITY), "Saved final identity mismatch.")
                require(abc == [row[k] for k in ("amplitude_mdpde", "width_mdpde", "intercept_mdpde")],
                        "Final capture differs from saved best parameters.")
        states.append({"id": state_id(policy, phase, attempt, best), "policy": policy, "phase": phase,
                       "attempt": attempt, "best_iteration": best, "context": str(path), "sha256": expected_hash,
                       "context_semantics": "best parameters with final certification context" if best else "recovery entry state"})
        evidence.extend([path, actual_path, run / "provenance.json", run / "execution.json"])
        if policy == "failed-only":
            evidence.append(run / "policy.json")
    return {"schema_version": 1, "states": states}, sorted(set(evidence))


def artifact_index(directory):
    files = sorted(p for p in directory.rglob("*") if p.is_file() and p.name != "artifact-index.json")
    experiment.write(directory / "artifact-index.json", {str(p.relative_to(directory)): fold.sha256_file(p) for p in files})


def run(args):
    output = args.output.resolve()
    require(not output.exists(), "Use a fresh output directory; historical evidence is immutable.")
    inputs = {k: getattr(args, k).resolve() for k in ("model", "map", "manifest")}
    baseline = fold.load_baseline(ROOT / "tests/benchmarks/fold_168_simulation_baseline.json")
    hashes = fold.validate_input_hashes(inputs, baseline["input_hashes"])
    manifest = fold.load_simulation_manifest(inputs["manifest"], hashes)
    fold.validate_fixture(manifest, baseline)
    index, evidence = prepare_index(args, manifest, hashes)
    before = experiment.provenance(args.executable.resolve(), args.source.resolve())
    input_hashes = {str(p): fold.sha256_file(p) for p in list(inputs.values()) + evidence}
    output.mkdir(parents=True)
    experiment.write(output / "state-index.json", index)
    experiment.write(output / "input-hashes.json", input_hashes)
    experiment.write(output / "provenance.json", before)
    experiment.write(output / "environment.json", {"python": sys.version, "platform": platform.platform(), "jobs": args.jobs,
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=args.source, text=True).strip()})
    # Scoring parameters are in a separate artifact and never passed to the fit API.
    experiment.write(output / "scoring-truth.json", {str(a["serial_id"]):
        [float(a["element"]), manifest["settings"]["blurring_width"], a["charge_used"]] for a in manifest["atoms"]})
    command = [str(args.executable.resolve()), "matched-sweep", str(inputs["manifest"]), str(inputs["map"]),
               str(output / "state-index.json"), str(output)]
    env = os.environ.copy(); env["OMP_NUM_THREADS"] = str(args.jobs)
    start = time.perf_counter()
    with (output / "run.log").open("w") as stream:
        result = subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT)
    experiment.write(output / "execution.json", {"command": command, "returncode": result.returncode,
                     "seconds": time.perf_counter()-start, "jobs": args.jobs})
    require(before == experiment.provenance(args.executable.resolve(), args.source.resolve()), "Source or binary changed during sweep.")
    require(all(fold.sha256_file(Path(p)) == h for p, h in input_hashes.items()), "Input changed during sweep.")
    require(result.returncode == 0, f"Sweep failed; inspect {output / 'run.log'}.")
    summarize(output)


def load_fits(directory):
    index = read(directory / "fit-index.json")
    expected_ids = {state_id(*c[:4]) for c in CHECKPOINTS}
    require(index["targets"] == 1344 and index["fits"] == 5376 and len(index["states"]) == 8
            and {s["id"] for s in index["states"]} == expected_ids, "Incomplete sweep matrix.")
    require(index["states"] == read(directory / "state-index.json")["states"], "Output state index differs from input.")
    rows = [r for i in range(index["targets"]) for r in read(directory / "fits" / f"{i}.json")]
    expected = {(s, serial, mode, fit) for s in expected_ids for serial in range(1, 169)
                for mode in ("analytic", "matched") for fit in ("AB", "ABC")}
    keys = [(r["state_id"], r["serial_id"], r["prediction"], r["parameters"]) for r in rows]
    require(len(rows) == len(set(keys)) == 5376 and set(keys) == expected, "Missing or duplicate fit pair.")
    require(all(r["n"] == len(r["membership"]) == (100 if r["parameters"] == "AB" else 200) for r in rows),
            "Changed fit population.")
    require(all(len(set(r["membership"])) == r["n"] and all(isinstance(p, int) and 0 <= p < 200 for p in r["membership"])
                and (r["parameters"] != "ABC" or r["membership"] == list(range(200))) for r in rows),
            "Invalid or duplicate sample membership.")
    return rows


def paired_results(rows, truths):
    groups = collections.defaultdict(dict)
    for row in rows:
        key = (row["state_id"], row["serial_id"], row["parameters"])
        require(row["prediction"] not in groups[key], "Duplicate prediction pair.")
        groups[key][row["prediction"]] = row
    pairs = []
    for (state, serial, fit), modes in sorted(groups.items()):
        require(set(modes) == {"analytic", "matched"}, "Incomplete prediction pair.")
        a, m = modes["analytic"], modes["matched"]
        require(a["membership"] == m["membership"] and a["input"] == m["input"] and a["identity"] == m["identity"],
                "Prediction pair uses different inputs or samples.")
        truth = truths[str(serial)]
        if fit == "AB":
            require(a["C"] == m["C"] == a["input"][2], "AB did not fix checkpoint C.")
        for i, p in enumerate(fit):
            for label, left, right in (("analytic-to-matched", a, m), ("input-to-analytic", None, a), ("input-to-matched", None, m)):
                x = a["input"][i] if left is None else left[p]
                y = right[p]
                references = [r for r in (left, right) if r is not None]
                available = all(finite(v) for v in (x, y, truth[i]))
                uncertainty_available = all(finite(r[p]) and finite(r["reference"][p]) for r in references)
                gain = abs(x-truth[i])-abs(y-truth[i]) if available else None
                uncertainty = sum(abs(r[p]-r["reference"][p]) for r in references) if uncertainty_available else None
                qualified = all(r["qualified"] for r in references)
                status = "unavailable" if not available or not uncertainty_available else "unqualified" if not qualified else (
                    "improved" if gain > uncertainty else "regressed" if gain < -uncertainty else "unresolved")
                pairs.append({"state_id": state, "serial_id": serial, "parameters": fit, "parameter": p,
                    "comparison": label, "left_error": x-truth[i] if finite(x) else None,
                    "right_error": y-truth[i] if finite(y) else None, "absolute_error_gain": gain,
                    "numerical_difference": uncertainty, "qualified": qualified, "status": status})
    return pairs


def population_stats(values, total):
    result = experiment.stats(v for v in values if finite(v))
    return {**result, "total": total, "missing": total-result["n"], "complete_population": result["n"] == total}


def fit_statistics(rows, truths):
    groups = collections.defaultdict(list)
    for row in rows:
        groups[(row["state_id"], row["parameters"], row["prediction"])].append(row)
    result = []
    for key, members in sorted(groups.items()):
        n = len(members)
        summary = dict(zip(("state_id", "parameters", "prediction"), key))
        summary.update(total=n, qualified=sum(r["qualified"] for r in members),
                       reasons=dict(collections.Counter(r["reason"] for r in members)))
        for i, p in enumerate(PARAMETERS):
            summary[p] = {"fixed": p == "C" and key[1] == "AB",
                "input": population_stats((r["input"][i]-truths[str(r["serial_id"])][i] for r in members), n),
                "estimate": population_stats((r[p]-truths[str(r["serial_id"])][i] if finite(r[p]) else None for r in members), n),
                "numerical_difference": population_stats((r[p]-r["reference"][p]
                    if finite(r[p]) and finite(r["reference"][p]) else None for r in members), n)}
        for field in ("condition", "stationarity", "input_loss", "loss"):
            summary[field] = population_stats((r[field] for r in members), n)
        for field in ("evaluations", "reference_evaluations", "linear_solves", "seconds"):
            summary[field] = sum(r[field] for r in members)
        result.append(summary)
    return result


def decomposition_error(values):
    require(all(finite(v) for v in values.values()), "Nonfinite subtraction evidence.")
    error = abs(values["total"]-values["estimation"]-values["operator"])
    require(error <= 64*sys.float_info.epsilon*max(1, *(abs(v) for v in values.values())),
            "Subtraction decomposition identity failed.")
    return error


def subtraction_statistics(directory, states):
    summaries, atom_summaries = [], []
    maximum_identity_error = 0.0
    for i, state in enumerate(states):
        groups, atoms, seen = collections.defaultdict(list), collections.defaultdict(list), set()
        with (directory / "samples" / f"{i}.csv").open() as stream:
            for row in csv.DictReader(stream):
                serial, sample = int(row["serial_id"]), int(row["sample"])
                require(row["state_id"] == state["id"] and (serial, sample) not in seen, "Duplicate or wrong sample state.")
                seen.add((serial, sample))
                distance = float(row["distance"])
                labels = ("all", "signal" if distance <= 1 else "tail",
                          "crosses-cutoff" if row["crosses_cutoff"] == "1" else "smooth-stencil",
                          "map-boundary" if row["boundary"] == "1" else "map-interior")
                for mode in ("analytic", "matched"):
                    values = {field: float(row[f"{mode}_{field}"]) for field in ("estimate", "truth", "estimation", "operator", "total")}
                    error = decomposition_error(values)
                    maximum_identity_error = max(maximum_identity_error, error)
                    for metric in ("estimation", "operator", "total"):
                        for label in labels:
                            groups[(mode, label, metric)].append(values[metric])
                        atoms[(serial, mode, metric)].append(values[metric])
        require(seen == {(serial, p) for serial in range(1, 169) for p in range(200)}, "Incomplete sample population.")
        for (mode, group, metric), values in sorted(groups.items()):
            summaries.append({"state_id": state["id"], "prediction": mode, "group": group, "metric": metric, **experiment.stats(values)})
        for (serial, mode, metric), values in sorted(atoms.items()):
            atom_summaries.append({"state_id": state["id"], "serial_id": serial, "prediction": mode,
                                   "metric": metric, **experiment.stats(values)})
    return summaries, atom_summaries, maximum_identity_error


def write_csv(path, rows):
    require(bool(rows), f"Empty CSV artifact: {path}")
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def summarize(directory):
    require(read(directory / "forward-status.json")["passed"], "Forward validation failed.")
    rows, truths = load_fits(directory), read(directory / "scoring-truth.json")
    states = read(directory / "state-index.json")["states"]
    pairs = paired_results(rows, truths)
    statistics = fit_statistics(rows, truths)
    neighbors, atom_neighbors, identity_error = subtraction_statistics(directory, states)
    counts = collections.defaultdict(collections.Counter)
    for p in pairs:
        counts[(p["state_id"], p["parameters"], p["parameter"], p["comparison"])][p["status"]] += 1
    comparisons = [{**dict(zip(("state_id", "parameters", "parameter", "comparison"), key)),
                    "total": sum(value.values()), "counts": dict(value)} for key, value in sorted(counts.items())]
    by_fit = {(r["state_id"], r["parameters"], r["prediction"]): r for r in statistics}
    for comparison in comparisons:
        key = (comparison["state_id"], comparison["parameters"])
        p = comparison["parameter"]
        left, right = comparison["comparison"].split("-to-")
        a = by_fit[(*key, "analytic" if left == "input" else left)][p]
        b = by_fit[(*key, right)][p]
        before, after = a["input" if left == "input" else "estimate"], b["estimate"]
        complete = before["complete_population"] and after["complete_population"]
        numerical = b["numerical_difference"]["complete_population"] and (left == "input" or a["numerical_difference"]["complete_population"])
        comparison.update(complete_population=complete, left_rmse=before.get("rmse"), right_rmse=after.get("rmse"),
                          rmse_ratio=after["rmse"]/before["rmse"] if complete and before.get("rmse") else None,
                          rmse_gain=before["rmse"]-after["rmse"] if complete else None,
                          numerical_difference_rms=((0 if left == "input" else a["numerical_difference"]["rmse"])
                          + b["numerical_difference"]["rmse"]) if numerical else None)
    worst = []
    for state in states:
        for parameter in PARAMETERS:
            candidates = [p for p in pairs if p["state_id"] == state["id"] and p["parameters"] == "ABC"
                          and p["parameter"] == parameter and p["comparison"] == "analytic-to-matched"
                          and finite(p["absolute_error_gain"]) and p["absolute_error_gain"] < 0]
            worst.extend(sorted(candidates, key=lambda p: (p["absolute_error_gain"], p["serial_id"]))[:5])
    summary = {"schema_version": 1, "experiment": "estimated-neighbor-sweep", "technical_complete": True,
        "states": states, "fit_count": len(rows), "failure_count": sum(not r["qualified"] for r in rows),
        "unavailable_fit_count": sum(any(not finite(r[p]) for p in r["parameters"]) for r in rows),
        "fit_statistics": statistics, "comparisons": comparisons, "neighbors": neighbors,
        "decomposition_max_error": identity_error, "worst_regressions": worst,
        "serial_100": [p for p in pairs if p["serial_id"] == 100],
        "production_quality_gate": "unchanged-uncalibrated", "convergence_assessed": False}
    experiment.write(directory / "results.json", summary)
    write_csv(directory / "paired-results.csv", pairs)
    write_csv(directory / "neighbor-atoms.csv", atom_neighbors)
    flat = []
    for row in rows:
        record = {k: row[k] for k in ("state_id", "serial_id", "parameters", "prediction", "n", "qualified", "reason",
                  "input_loss", "loss", "stationarity", "condition", "reference_difference", "evaluations", "reference_evaluations", "linear_solves", "seconds")}
        record.update({k: row["identity"][k] for k in IDENTITY if k != "position"})
        for i, p in enumerate(PARAMETERS):
            truth = truths[str(row["serial_id"])][i]
            record.update({p: row[p], f"input_{p}": row["input"][i], f"truth_{p}": truth,
                           f"reference_{p}": row["reference"][p], f"error_{p}": row[p]-truth if finite(row[p]) else None})
        flat.append(record)
    write_csv(directory / "fits.csv", flat)
    artifact_index(directory)
    print(json.dumps({k: summary[k] for k in ("technical_complete", "fit_count", "failure_count", "unavailable_fit_count")}, indent=2))


def compare(left, right, output):
    a, b = load_fits(left), load_fits(right)
    differences = [i for i, (x, y) in enumerate(zip(a, b))
                   if {k: v for k, v in x.items() if k != "seconds"} != {k: v for k, v in y.items() if k != "seconds"}]
    samples_equal = all((left / "samples" / f"{i}.csv").read_bytes() == (right / "samples" / f"{i}.csv").read_bytes() for i in range(8))
    statistics = [read(p / "results.json") for p in (left, right)]
    for value in statistics:
        for row in value["fit_statistics"]:
            row.pop("seconds")
    result = {"passed": not differences and samples_equal and statistics[0] == statistics[1]
              and read(left / "forward-status.json") == read(right / "forward-status.json"),
              "fit_differences": differences, "samples_equal": samples_equal, "statistics_equal": statistics[0] == statistics[1],
              "fit_count": len(a)}
    experiment.write(output, result)
    require(result["passed"], "j1/j4 sweep evidence differs.")


def plots(directory, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    result = read(directory / "results.json")
    output.mkdir(parents=True, exist_ok=True)
    states = [s["id"] for s in result["states"]]
    labels = [s.replace("baseline-", "base ").replace("failed-only-", "refined ") for s in states]
    stats = {(r["state_id"], r["parameters"], r["prediction"]): r for r in result["fit_statistics"]}
    colors = {"input": "#777777", "analytic": "#b56b45", "matched": "#246b86"}
    def save(fig, name):
        for extension in ("png", "pdf"):
            fig.savefig(output / f"{name}.{extension}", dpi=180)
        plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    for axis, p in zip(axes, PARAMETERS):
        for offset, mode in enumerate(colors):
            values = [stats[(s, "ABC", "analytic" if mode == "input" else mode)][p]["input" if mode == "input" else "estimate"].get("rmse", math.nan) for s in states]
            axis.bar([i+(offset-1)*.25 for i in range(8)], values, width=.25, label=mode, color=colors[mode])
        axis.set_xticks(range(8), labels, rotation=60, ha="right"); axis.set_yscale("log"); axis.set_title(f"{p} RMSE")
        axis.grid(axis="y", alpha=.2)
    axes[0].legend(); fig.suptitle("Fold-168 · frozen estimated neighbors · ABC one-sweep")
    save(fig, "parameter-errors")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    for axis, metric in zip(axes, ("estimation", "operator", "total")):
        for offset, mode in enumerate(("analytic", "matched")):
            values = [next(r["rmse"] for r in result["neighbors"] if r["state_id"] == s and r["prediction"] == mode
                           and r["metric"] == metric and r["group"] == "all") for s in states]
            axis.bar([i+(offset-.5)*.35 for i in range(8)], values, width=.35, label=mode, color=colors[mode])
        axis.set_xticks(range(8), labels, rotation=60, ha="right"); axis.set_title(f"Neighbor {metric} RMSE"); axis.grid(axis="y", alpha=.2)
    axes[0].legend(); fig.suptitle("Pointwise decomposition · RMSE components are not additive")
    save(fig, "neighbor-errors")
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), constrained_layout=True)
    for axis, p in zip(axes, PARAMETERS):
        for mode in colors:
            values = [stats[(f"failed-only-recovery-{i}", "ABC", "analytic" if mode == "input" else mode)][p]
                      ["input" if mode == "input" else "estimate"].get("rmse", math.nan) for i in range(10, 15)]
            axis.plot(range(10, 15), values, "o-", color=colors[mode], label=mode)
        axis.set_xticks(range(10, 15)); axis.set_xlabel("Recovery attempt"); axis.set_title(f"{p} RMSE"); axis.set_yscale("log"); axis.grid(alpha=.2)
    axes[0].legend(); fig.suptitle("Each sweep freezes its own recovery checkpoint")
    save(fig, "recovery-trajectory")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    runner = sub.add_parser("run")
    for name in ("executable", "model", "map", "manifest", "baseline-run", "refined-run", "output"):
        runner.add_argument("--"+name, type=Path, required=True)
    runner.add_argument("--source", type=Path, default=ROOT)
    runner.add_argument("--jobs", type=int, choices=(1, 4), default=4)
    summary = sub.add_parser("summarize"); summary.add_argument("directory", type=Path)
    comparison = sub.add_parser("compare")
    for name in ("left", "right", "output"):
        comparison.add_argument("--"+name, type=Path, required=True)
    plot = sub.add_parser("plots"); plot.add_argument("directory", type=Path); plot.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.operation == "run": run(args)
    elif args.operation == "summarize": summarize(args.directory)
    elif args.operation == "compare": compare(args.left, args.right, args.output)
    else: plots(args.directory, args.output)


if __name__ == "__main__":
    main()
