#!/usr/bin/env python3
"""Pool paired 6Z6U Halo-exclusion metric runs without hiding target qualification."""

import argparse
import csv
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

PARAMETERS = (
    ("A", "A_estimate_raw_units", "A_truth_raw_units", "A_error_raw_units"),
    ("B_angstrom", "B_estimate_angstrom", "B_truth_angstrom", "B_error_angstrom"),
    ("C", "C_estimate_raw_units", "C_truth_raw_units", "C_error_raw_units"),
)
STATUS_FIELDS = ("point_quality", "parameter_validity", "local_conditioning", "parameter_stability")


def error_stats(values):
    if not values:
        return {"count": 0, "bias": None, "variance_population": None,
                "variance_sample": None, "mse": None, "rmse": None}
    bias = statistics.fmean(values)
    variance_population = statistics.pvariance(values)
    variance_sample = statistics.variance(values) if len(values) > 1 else None
    mse = statistics.fmean(value * value for value in values)
    absolute = sorted(abs(value) for value in values)
    p95 = absolute[math.ceil(.95 * len(absolute)) - 1]
    return {"count": len(values), "bias": bias,
            "variance_population": variance_population, "variance_sample": variance_sample,
            "mse": mse, "mse_bias_plus_population_variance": bias * bias + variance_population,
            "rmse": math.sqrt(mse), "median_absolute_error": statistics.median(absolute),
            "p95_absolute_error": p95, "maximum_absolute_error": absolute[-1]}


def target_group_stats(target_groups, error_field, qualified_only):
    usable = []
    for atom, rows in target_groups:
        values = [row[error_field] for row in rows if row[error_field] is not None and
                  (not qualified_only or row["qualified_estimate"])]
        if values:
            usable.append((atom, rows, values))
    if not usable:
        return {"target_count": 0, "repetition_count": 0, "pooled_target_seed_error": error_stats([]),
                "mean_target_bias": None, "mean_target_squared_bias": None,
                "mean_target_population_variance": None, "mean_target_sample_variance": None,
                "mean_target_mse": None, "mean_target_bias_squared_plus_population_variance": None,
                "per_target": []}
    biases = [statistics.fmean(values) for _, _, values in usable]
    population_variances = [statistics.pvariance(values) for _, _, values in usable]
    sample_variances = [statistics.variance(values) for _, _, values in usable if len(values) > 1]
    mses = [statistics.fmean(value * value for value in values) for _, _, values in usable]
    per_target = []
    for index, (atom, rows, values) in enumerate(usable):
        per_target.append({"atom_id": atom, "coverage_ratio": rows[0]["coverage_ratio"],
            "repetition_count": len(values), "bias": biases[index],
            "population_variance": population_variances[index],
            "sample_variance": statistics.variance(values) if len(values) > 1 else None,
            "mse": mses[index]})
    return {"target_count": len(usable), "repetition_count": min(len(values) for _, _, values in usable),
        "pooled_target_seed_error": error_stats([value for _, _, values in usable for value in values]),
        "mean_target_bias": statistics.fmean(biases),
        "mean_target_squared_bias": statistics.fmean(value * value for value in biases),
        "mean_target_population_variance": statistics.fmean(population_variances),
        "mean_target_sample_variance": statistics.fmean(sample_variances) if sample_variances else None,
        "mean_target_mse": statistics.fmean(mses),
        "mean_target_bias_squared_plus_population_variance": statistics.fmean(
            bias * bias + variance for bias, variance in zip(biases, population_variances)),
        "median_absolute_target_bias": statistics.median(abs(value) for value in biases),
        "p95_absolute_target_bias": sorted(abs(value) for value in biases)[math.ceil(.95*len(biases))-1],
        "maximum_absolute_target_bias": max(abs(value) for value in biases),
        "per_target": per_target}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coverage_json", type=Path, help="P0 6Z6U mask audit JSON")
    parser.add_argument("output_json", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("metrics_json", type=Path, nargs="+")
    args = parser.parse_args()

    coverage = json.loads(args.coverage_json.read_text())
    low_coverage = {
        target["atom_id"] for target in coverage["local_identifiability"]["targets"]
        if target["retained_rows"] / target["original_rows"] < 0.5
    }
    coverage_by_target = {
        target["atom_id"]: target for target in coverage["local_identifiability"]["targets"]
    }
    metrics = [json.loads(path.read_text()) for path in args.metrics_json]
    seeds = [item["fit_noise"]["seed"] for item in metrics]
    if len(set(seeds)) != len(seeds) or len(seeds) < 2:
        raise ValueError("Metrics must contain at least two unique fixed-seed repetitions.")
    if any(item["fit_noise"]["sigma_normalized_map_units"] != .02 or
           not item["fit_noise"]["first_stage_uses_same_noisy_map"] for item in metrics):
        raise ValueError("This summary expects paired sigma=0.02 full-map noise runs.")

    records = []
    global_checks = defaultdict(Counter)
    for item in metrics:
        seed = item["fit_noise"]["seed"]
        row_counts = {run["method"]: run["observation_count"] for run in item["runs"]}
        for run in item["runs"]:
            global_checks[(run["method"], "runtime_convergence")][run["runtime_convergence"]] += 1
            global_checks[(run["method"], "target_runtime_convergence")][run["target_runtime_convergence"]] += 1
            for target in run["per_target_errors"]:
                target_coverage = coverage_by_target[target["atom_id"]]
                records.append({"method": run["method"], "seed": seed,
                    "noise_sha256": item["fit_noise"]["realization_sha256"],
                    "atom_id": target["atom_id"],
                    "coverage_ratio": target_coverage["retained_rows"] / target_coverage["original_rows"],
                    "coverage_group": "low_lt50pct" if target["atom_id"] in low_coverage else "other",
                    "distance_stratum": target["distance_stratum"],
                    "runtime_convergence": run["runtime_convergence"],
                    "target_runtime_convergence": run["target_runtime_convergence"],
                    "qualified_estimate": target["qualified_estimate"],
                    **{field: target[field] for field in STATUS_FIELDS},
                    "row_count": row_counts[run["method"]],
                    **{field: target[field] for _, estimate, truth, error in PARAMETERS
                       for field in (estimate, truth, error)}})

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    columns = list(records[0])
    with args.output_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(sorted(records, key=lambda row: (row["method"], row["atom_id"], row["seed"])))

    per_target = defaultdict(list)
    statuses = defaultdict(Counter)
    methods = sorted({record["method"] for record in records})
    for record in records:
        group = (record["method"], record["coverage_group"], record["atom_id"])
        per_target[group].append(record)
        for field in STATUS_FIELDS:
            statuses[(record["method"], record["coverage_group"], field)][record[field]] += 1

    summaries = {}
    for method in methods:
        summaries[method] = {}
        for group_name in ("all", "low_lt50pct", "other"):
            target_groups = [
                (atom, values) for (m, group, atom), values in per_target.items()
                if m == method and (group_name == "all" or group == group_name)
            ]
            parameter_summaries = {}
            for name, estimate, truth, error in PARAMETERS:
                parameter_summaries[name] = {
                    "finite_candidate": target_group_stats(target_groups, error, False),
                    "qualified_only": target_group_stats(target_groups, error, True)}
            status_summary = {field: dict(statuses[(method, group_name, field)])
                              for field in STATUS_FIELDS}
            if group_name == "all":
                status_summary = {field: dict(sum((statuses[(method, subgroup, field)]
                    for subgroup in ("low_lt50pct", "other")), Counter())) for field in STATUS_FIELDS}
            summaries[method][group_name] = {"target_count": len(target_groups),
                "repetition_count": len(seeds), "parameter_bias_variance_mse": parameter_summaries,
                "qualification_status_counts_across_target_repetitions": status_summary,
                "global_convergence_status_counts": {field: dict(global_checks[(method, field)])
                    for field in ("runtime_convergence", "target_runtime_convergence")}}

    record_by_key = {(row["method"], row["seed"], row["atom_id"]): row for row in records}
    paired_differences = {}
    for comparator in ("A", "R2"):
        comparison = f"B_minus_{comparator}_candidate_target_mse"
        paired_differences[comparison] = {}
        for group_name in ("all", "low_lt50pct", "other"):
            atom_ids = sorted({row["atom_id"] for row in records
                if group_name == "all" or row["coverage_group"] == group_name})
            differences = []
            for atom in atom_ids:
                deltas = {}
                paired_count = 0
                for _, _, _, error in PARAMETERS:
                    left, right = [], []
                    for seed in seeds:
                        b = record_by_key.get(("B", seed, atom))
                        base = record_by_key.get((comparator, seed, atom))
                        if b and base and b[error] is not None and base[error] is not None:
                            left.append(b[error]); right.append(base[error])
                    paired_count = min(paired_count or len(left), len(left)) if left else paired_count
                    deltas[error.replace("_error", "_mse_difference")] = (
                        statistics.fmean(value * value for value in left) -
                        statistics.fmean(value * value for value in right) if left and right else None)
                if paired_count:
                    differences.append({"atom_id": atom, "coverage_ratio": coverage_by_target[atom]["retained_rows"] /
                        coverage_by_target[atom]["original_rows"], "paired_repetition_count": paired_count, **deltas})
            parameter_deltas = {}
            worst_regressions = {}
            for _, _, _, error in PARAMETERS:
                field = error.replace("_error", "_mse_difference")
                values = [row[field] for row in differences if row[field] is not None]
                parameter_deltas[field] = error_stats(values)
                worst_regressions[field] = sorted(
                    (row for row in differences if row[field] is not None),
                    key=lambda row: row[field], reverse=True)[:20]
            paired_differences[comparison][group_name] = {
                "target_count_with_paired_values": len(differences),
                "mean_target_mse_difference": parameter_deltas,
                "largest_exclusion_regressions_by_parameter": worst_regressions}

    output = {"contract": "joint-halo-support-exclusion-6z6u-noise-statistics-v1",
        "noise_sigma_normalized_map_units": .02, "noise_layer": metrics[0]["fit_noise"]["layer"],
        "paired_seed_count": len(seeds), "seeds": sorted(seeds),
        "noise_realization_sha256_by_seed": {str(item["fit_noise"]["seed"]): item["fit_noise"]["realization_sha256"] for item in metrics},
        "low_coverage_definition": "exact-mask retained rows / original target rows < 0.5",
        "low_coverage_target_count": len(low_coverage),
        "target_coverage_source": args.coverage_json.name,
        "methods": summaries, "paired_target_mse_differences": paired_differences}
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
