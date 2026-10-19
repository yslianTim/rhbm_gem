#!/usr/bin/env python3
"""Run paired width-trust cases and audit policy scope and pre-intervention parity."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


SYNTHETIC_DATASETS = (
    "well-synthetic",
    "weak-synthetic",
    "moderate-synthetic",
    "noise-high-synthetic",
    "mixed-width-synthetic",
    "ac-correlated-synthetic",
)
STRATEGIES = ("A", "B1", "B2")
METRICS = (
    "source", "simulation_mode", "normalization", "truth_source", "truth_scale", "estimate_scale",
    "comparison_scale", "truth_available", "target_count", "observation_count",
    "policy_scope", "diagnostic_scope",
    "component_qualification", "global_qualification", "component_stop_reason",
    "search_completed", "independent_replay", "candidate_available", "candidate_objective",
    "candidate_ac_kkt", "candidate_width_gradient_inf_norm", "A_rmse", "B_rmse_angstrom",
    "C_rmse", "prediction_rmse", "maximum_accepted_width_A", "maximum_local_ac_condition",
    "accuracy_metric_status", "prediction_metric_status",
    "maximum_ac_cancellation_ratio", "minimum_accepted_projected_width_singular_value",
    "maximum_observed_width_A", "maximum_observed_local_ac_condition",
    "maximum_observed_ac_cancellation_ratio", "maximum_observed_abs_ac_basis_correlation",
    "minimum_observed_projected_width_singular_value", "extreme_width_incident",
    "policy_rejected_trials", "policy_statistics", "diagnostic_capture_limits", "search_sweeps",
    "profile_evaluations", "wall_seconds", "process_peak_rss_bytes",
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_run(binary, output_root, dataset, seed, strategy, model=None, map_path=None,
             execute=True):
    case = f"{dataset}/{seed if seed is not None else 'fixed-input'}"
    output = output_root / case / strategy
    output.parent.mkdir(parents=True, exist_ok=True)
    log_path = output.parent / f"{strategy}.log"
    command = [str(binary), dataset, str(output), strategy, str(seed or 0)]
    if model is not None:
        command.extend((str(model), str(map_path)))
    if execute:
        with log_path.open("w") as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
        if completed.returncode:
            tail = "\n".join(log_path.read_text().splitlines()[-40:])
            raise RuntimeError(f"Experiment failed ({completed.returncode}): {' '.join(command)}\n{tail}")
    elif not (output / "result.json").is_file():
        raise SystemExit(f"Experiment result is missing: {output / 'result.json'}")
    result = json.loads((output / "result.json").read_text())
    summary = {key: result.get(key) for key in METRICS}
    qualification = result.get("global_qualification")
    summary["accuracy_metric_status"] = result.get("accuracy_metric_status") or (
        "N/A — no truth" if not result.get("truth_available") else
        "qualified-estimate" if qualification == "Passed" else "unqualified-candidate-metrics"
    )
    summary["prediction_metric_status"] = result.get("prediction_metric_status") or (
        "Unavailable" if not result.get("candidate_available") else
        "qualified-fit" if qualification == "Passed" else "unqualified-candidate"
    )
    summary["seed"] = result.get("seed")
    summary["block_trace_count"] = len(result["block_update_trace"])
    return summary, result["block_update_trace"], command


def before_first_intervention(trace, statistics):
    sweep = statistics["first_rejection_sweep"]
    block = statistics["first_rejection_block"]
    component = statistics["first_rejection_component"]
    if sweep is None:
        return len(trace), None
    for index, event in enumerate(trace):
        if (event["sweep"], event["block"], event["component_index"]) == (sweep, block, component):
            return index, index
    raise RuntimeError("First policy rejection is absent from the block trace.")


def audit_pair(baseline, baseline_trace, candidate, candidate_trace):
    scope = candidate["policy_scope"]
    stats = candidate["policy_statistics"]
    prefix_count, intervention_index = before_first_intervention(candidate_trace, stats)
    candidate_prefix = candidate_trace[:prefix_count]
    baseline_prefix = baseline_trace[:prefix_count]
    target_blocks = {
        (block["component_index"], block["block_number_one_based"]) for block in scope
    }
    out_of_scope_prefix = lambda trace: [
        event for event in trace
        if (event["component_index"], event["block"]) not in target_blocks
    ]
    rejected_sweep = stats["first_rejection_sweep"]
    rejected_block = stats["first_rejection_block"]
    rejected_component = stats["first_rejection_component"]
    rejection_is_scoped = rejected_sweep is None or (
        rejected_component, rejected_block
    ) in target_blocks
    return {
        "policy_scope_matches_baseline": scope == baseline["policy_scope"],
        "out_of_scope_interventions": stats["out_of_scope_interventions"],
        "first_intervention_is_in_scope": rejection_is_scoped,
        "first_intervention_trace_index": intervention_index,
        "events_before_first_intervention": prefix_count,
        "out_of_scope_events_before_first_intervention": len(out_of_scope_prefix(candidate_prefix)),
        "pre_intervention_trace_matches_baseline": candidate_prefix == baseline_prefix,
        "out_of_scope_pre_intervention_trace_matches_baseline": (
            out_of_scope_prefix(candidate_prefix) == out_of_scope_prefix(baseline_prefix)
        ),
        "complete_trace_matches_baseline": candidate_trace == baseline_trace,
    }


def run_case(binary, output_root, dataset, seed, model=None, map_path=None, execute=True):
    runs = {}
    traces = {}
    commands = {}
    for strategy in STRATEGIES:
        runs[strategy], traces[strategy], commands[strategy] = load_run(
            binary, output_root, dataset, seed, strategy, model, map_path, execute
        )
    return {
        "dataset": dataset,
        "seed": seed,
        "strategies": runs,
        "scope_audit": {
            strategy: audit_pair(runs["A"], traces["A"], runs[strategy], traces[strategy])
            for strategy in ("B1", "B2")
        },
        "paired_inputs_match": all(
            runs["A"][key] == runs[strategy][key]
            for strategy in ("B1", "B2")
            for key in ("source", "simulation_mode", "normalization", "truth_available",
                        "target_count", "observation_count", "truth_source", "policy_scope")
        ),
        "commands": commands,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("model", type=Path)
    parser.add_argument("simulation_map", type=Path)
    parser.add_argument("unsharpened_map", type=Path)
    parser.add_argument("--seeds", type=int, nargs="+", default=[1001, 1002, 1003, 1004, 1005])
    parser.add_argument("--summarize-existing", action="store_true")
    args = parser.parse_args()
    if args.output.exists() and not args.summarize_existing:
        raise SystemExit(f"Choose a fresh output directory: {args.output}")
    if args.summarize_existing and not args.output.is_dir():
        raise SystemExit(f"Existing experiment directory is missing: {args.output}")
    for path in (args.binary, args.model, args.simulation_map, args.unsharpened_map):
        if not path.is_file():
            raise SystemExit(f"Input file is missing: {path}")
    if not args.summarize_existing:
        args.output.mkdir(parents=True)

    cases = []
    for dataset in SYNTHETIC_DATASETS:
        for seed in args.seeds:
            cases.append(run_case(args.binary, args.output, dataset, seed, execute=not args.summarize_existing))
    cases.append(run_case(args.binary, args.output, "6z6u-simulation", None,
                          args.model, args.simulation_map, execute=not args.summarize_existing))
    cases.append(run_case(args.binary, args.output, "6z6u-unsharpened", None,
                          args.model, args.unsharpened_map, execute=not args.summarize_existing))

    report = {
        "contract": "joint-width-trust-corrected-paired-results-v1",
        "seeds": args.seeds,
        "synthetic_datasets": list(SYNTHETIC_DATASETS),
        "strategies": list(STRATEGIES),
        "input_sha256": {
            str(path): sha256(path) for path in
            (args.model, args.simulation_map, Path(str(args.simulation_map) + ".simulation.json"),
             args.unsharpened_map)
        },
        "cases": cases,
    }
    report["all_scope_audits_pass"] = all(
        audit["out_of_scope_interventions"] == 0
        and audit["policy_scope_matches_baseline"]
        and audit["first_intervention_is_in_scope"]
        and audit["pre_intervention_trace_matches_baseline"]
        and audit["out_of_scope_pre_intervention_trace_matches_baseline"]
        for case in cases for audit in case["scope_audit"].values()
    )
    (args.output / "paired-summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "case_count": len(cases),
        "run_count": len(cases) * len(STRATEGIES),
        "summary": str(args.output / "paired-summary.json"),
        "all_scope_audits_pass": report["all_scope_audits_pass"],
    }, indent=2))


if __name__ == "__main__":
    main()
