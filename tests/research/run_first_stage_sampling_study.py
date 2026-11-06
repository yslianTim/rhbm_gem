#!/usr/bin/env python3
"""Run the bounded, sequential study and save compact evidence using stdlib only."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess

POLICIES = ("F0", "F1", "Q0", "Q1")


def mean(values):
    return statistics.mean(values) if values else None


def summarize(data):
    atoms = data["atoms"]
    result = {key: data[key] for key in ("policy", "scenario", "bound_mode", "sampling",
              "interpolation", "grid_spacing", "calibration", "initialization_seconds",
              "initial_objective", "joint_seconds", "total_seconds", "peak_rss_bytes", "qscore")}
    for key in ("raw", "selected", "positive_raw", "positive_selected", "fitting_candidates", "retained"):
        result[key] = sum(a[key] for a in atoms)
    result["alpha_median"] = statistics.median(a["alpha"] for a in atoms)
    result["rank_counts"] = {str(rank): sum(a["rank"] == rank for a in atoms) for rank in (0, 1, 2)}
    result["distinct_fitting_radii_range"] = [min(a["distinct_fitting_radii"] for a in atoms),
                                             max(a["distinct_fitting_radii"] for a in atoms)]
    result["radial_counts"] = atoms[0]["radial_counts"]
    for estimator in ("ols", "mdpde"):
        for coordinate, index in (("amplitude", 0), ("width", 1)):
            errors = [a[estimator][index] - a["truth"][index] for a in atoms
                      if "truth" in a and a[estimator][index] is not None]
            if errors:
                result[f"{estimator}_{coordinate}_bias"] = mean(errors)
                result[f"{estimator}_{coordinate}_rmse"] = math.sqrt(mean([e * e for e in errors]))
        result[f"{estimator}_mean"] = [mean([a[estimator][k] for a in atoms if a[estimator][k] is not None])
                                       for k in (0, 1)]
    result["working_width_variance_median"] = statistics.median(
        a["working_information"]["amplitude_width_covariance"][3] for a in atoms
        if a["working_information"]) if any(a["working_information"] for a in atoms) else None
    result["working_information_min_eigenvalue_median"] = statistics.median(
        a["working_information"]["eigenvalues"][0] for a in atoms
        if a["working_information"]) if any(a["working_information"] for a in atoms) else None
    shells = [s for a in atoms for s in a.get("shells", [])]
    result["shells"] = {"successful": sum(s["successful"] for s in shells),
                        "failed": sum(not s["successful"] for s in shells),
                        "retries": sum(s["attempts"] - 1 for s in shells),
                        "requested": sum(s["requested"] for s in shells),
                        "final_candidates": sum(s["candidates"] for s in shells),
                        "used_accepted": sum(s["accepted"] for s in shells if s["successful"])}
    if shells:
        accepted = [s["accepted"] for s in shells if s["successful"]]
        result["shells"]["used_accepted_range"] = [min(accepted), max(accepted)] if accepted else []
    fit = data["joint"]
    result["joint"] = {key: fit.get(key) for key in ("objective", "runtime_qualification", "target_qualification")}
    result["joint"]["sweeps_max"] = max([c["sweeps"] for c in fit["components"]], default=0)
    result["joint"]["sweeps_sum"] = sum(c["sweeps"] for c in fit["components"])
    result["joint"]["stop_reasons"] = sorted(set(c["stop"] for c in fit["components"]))
    checks = fit["checks"] + [check for c in fit["components"] for check in c["checks"]]
    result["joint"]["failed_checks"] = [c for c in checks if c["status"] != 0]
    result["joint"]["replay_passed"] = all(c["status"] == 0 for c in checks if c["name"] == "independent-replay") and any(
        c["name"] == "independent-replay" for c in checks)
    result["joint"]["kkt_passed"] = all(c["status"] == 0 for c in checks if c["name"] in (
        "conditional-ac-optimality", "width-stationarity")) and any(c["name"] == "width-stationarity" for c in checks)
    result["joint"]["endpoint_passed"] = all(c["status"] == 0 for c in checks if c["name"] == "endpoint-verification") and any(
        c["name"] == "endpoint-verification" for c in checks)
    for key in ("common_bound", "common_projections", "downstream"):
        if key in data:
            result[key] = data[key]
    return result


def shared_equal(baseline, current):
    if isinstance(baseline, dict):
        return all(k in current and shared_equal(v, current[k]) for k, v in baseline.items())
    if isinstance(baseline, list):
        return len(baseline) == len(current) and all(shared_equal(a, b) for a, b in zip(baseline, current))
    return baseline == current


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--build", type=Path, default=Path("build/release"))
    parser.add_argument("--real-data", type=Path)
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    build = (root / args.build).resolve()
    runs = build / "first-stage-study"
    runs.mkdir(exist_ok=True)
    binary = build / "bin/first_stage_sampling_study"
    env = dict(os.environ, OMP_NUM_THREADS="1", OMP_THREAD_LIMIT="4", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    results = {}

    def run(name, scenario, policy, mode="production", paths=()):
        path = runs / f"{name}-{policy}-{mode}.json"
        if not args.summarize_only:
            print(path.name, flush=True)
            with path.with_suffix(".log").open("w") as log:
                subprocess.run([str(binary), str(path), scenario, policy, str(mode), *map(str, paths)],
                               cwd=root, env=env, stdout=log, stderr=log, check=True, timeout=600)
        data = json.loads(path.read_text())
        results[path.stem] = summarize(data)
        return data

    for scenario in ("fine", "coarse", "noisy", "crowded"):
        for replicate in range(3):
            name = f"{scenario}-r{replicate}"
            for policy in POLICIES:
                run(name, f"{scenario}:{replicate}", policy)
            if scenario in ("coarse", "crowded"):
                for policy in ("Q8", "S0", "S1"):
                    run(name, f"{scenario}:{replicate}", policy)
            if replicate == 0:
                bound = max(results[f"{name}-{p}-production"]["calibration"].get("bound", 0.) for p in POLICIES)
                if bound > 0.:
                    for policy in POLICIES:
                        run(name, f"{scenario}:{replicate}", policy, bound)
    provenance = {}
    if args.real_data:
        for model_name, map_name in (("6z6u.cif", "emd_11103_unsharpened.map"), ("8rqb.cif", "emd_19436_unsharpened.map")):
            paths = (args.real_data / model_name, args.real_data / map_name)
            if not all(p.is_file() for p in paths):
                provenance[model_name] = "unavailable"
                continue
            provenance[model_name] = {}
            for path in paths:
                digest = hashlib.sha256()
                with path.open("rb") as stream:
                    for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                        digest.update(chunk)
                provenance[model_name][path.name] = {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}
            name = model_name.split(".")[0]
            for policy in POLICIES:
                run(name, "real", policy, paths=paths)
            bound = max(results[f"{name}-{p}-production"]["calibration"].get("bound", 0.) for p in POLICIES)
            if bound > 0.:
                for policy in POLICIES:
                    run(name, "real", policy, bound, paths)
    baseline_dir = root / "docs/research/first-stage-sampling"
    parity = {}
    for old, new in (("baseline.json", "coarse-r0-F0-production.json"), ("baseline-6z6u.json", "6z6u-F0-production.json")):
        if (runs / new).exists():
            a, b = json.loads((baseline_dir / old).read_text()), json.loads((runs / new).read_text())
            parity[old] = all(shared_equal(a[key], b[key]) for key in ("atoms", "calibration", "joint", "qscore"))
            assert parity[old], f"F0 baseline changed: {old}"
    # One complete run per line keeps the checked-in numerical evidence reviewable.
    with (baseline_dir / "results.json").open("w") as stream:
        stream.write('{\n  "baseline_exact_parity": ' + json.dumps(parity) + ',\n  "real_data": ' + json.dumps(provenance))
        stream.write(',\n  "runs": {\n')
        stream.write(',\n'.join('    ' + json.dumps(name) + ': ' + json.dumps(value) for name, value in results.items()))
        stream.write('\n  }\n}\n')
    print(f"Saved {len(results)} study summaries; baseline parity: {parity}", flush=True)


if __name__ == "__main__":
    main()
