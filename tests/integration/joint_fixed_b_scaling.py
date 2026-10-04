#!/usr/bin/env python3
"""Run and summarize the fixed-B chain/cube block-scaling campaign."""
import argparse
import csv
import json
import math
import subprocess
from pathlib import Path


TOPOLOGIES = ("chain", "cube")
SIZES = (128, 256, 512, 1024)
CORE_ATOMS = 128
GLOBAL_KKT_TOLERANCE = 1e-10
OBJECTIVE_TOLERANCE = 1e-12
AC_TOLERANCE = 1e-8


def classify_trend(rows):
    forward = {row["atoms"]: row for row in rows if row["order"] == "forward"}
    sizes = (256, 512, 1024)
    if any(size not in forward for size in sizes):
        return {"classification": "insufficient-evidence", "reason": "missing multi-block size"}
    selected = [forward[size] for size in sizes]
    if any(not row["global_kkt_passed"] for row in selected):
        return {"classification": "insufficient-evidence", "reason": "a multi-block case did not reach global KKT"}
    sweeps = [row["sweeps_to_global_kkt"] for row in selected]
    ratio = sweeps[-1] / sweeps[0] if sweeps[0] else math.inf
    if ratio >= 1.5 and sweeps[0] <= sweeps[1] <= sweeps[2]:
        classification = "growth-looking"
    elif ratio < 1.5 and max(sweeps) / min(sweeps) < 1.5:
        classification = "stable-looking"
    else:
        classification = "insufficient-evidence"
    return {
        "classification": classification,
        "multi_block_sizes": list(sizes),
        "sweeps_to_global_kkt": sweeps,
        "endpoint_sweep_ratio": ratio,
        "formal_growth_gate": "not evaluated; only three multi-block sizes",
    }


def compare_order(forward, reverse):
    f_converged = bool(forward["global_kkt_passed"])
    r_converged = bool(reverse["global_kkt_passed"])
    if not (f_converged and r_converged):
        return {
            "status": "non-convergence-asymmetry" if f_converged != r_converged else "not-converged",
            "compared_converged_endpoints": False,
            "passed": False,
        }
    parameter_difference = max(
        (abs(a - b) for a, b in zip(forward["block_scaled_parameters"], reverse["block_scaled_parameters"])),
        default=0.0,
    )
    objective_difference = abs(forward["final_objective"] - reverse["final_objective"])
    objective_tolerance = OBJECTIVE_TOLERANCE + 2e-12 * max(
        abs(forward["global_objective"]), abs(reverse["global_objective"])
    )
    passed = objective_difference <= objective_tolerance and parameter_difference <= AC_TOLERANCE
    return {
        "status": "converged-endpoint-comparison",
        "compared_converged_endpoints": True,
        "objective_difference": objective_difference,
        "scaled_ac_difference": parameter_difference,
        "sweep_count_difference": forward["sweeps"] - reverse["sweeps"],
        "objective_tolerance": objective_tolerance,
        "ac_tolerance": AC_TOLERANCE,
        "passed": passed,
    }


def analyze(rows):
    trends = {topology: classify_trend([row for row in rows if row["topology"] == topology])
              for topology in TOPOLOGIES}
    order_results = {}
    for topology in TOPOLOGIES:
        forward = next(row for row in rows if row["topology"] == topology and row["atoms"] == 512 and row["order"] == "forward")
        reverse = next(row for row in rows if row["topology"] == topology and row["atoms"] == 512 and row["order"] == "reverse")
        order_results[topology] = compare_order(forward, reverse)
    cases_passed = all(row["success"] and row["objective_parity"] and row["global_kkt_passed"] for row in rows)
    stable = all(value["classification"] == "stable-looking" for value in trends.values())
    order_parity = all(value["passed"] for value in order_results.values())
    return {
        "fixed_b_scaling": "passed" if cases_passed else "failed",
        "all_cases_converged_and_objective_parity": cases_passed,
        "sweep_trends": trends,
        "forward_reverse_512": order_results,
        "no_severe_sweep_growth": stable,
        "no_material_converged_order_dependence": order_parity,
        "resume_nonlinear_f2": cases_passed and stable and order_parity,
    }


def _run(command):
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    if result.stdout.strip():
        print(result.stdout.strip())
    if result.stderr.strip():
        print(result.stderr.strip())


def run_campaign(executable, output_dir):
    executable = Path(executable).resolve()
    root = Path(output_dir)
    references_dir = root / "global-references"
    individuals_dir = root / "individual-results"
    references_dir.mkdir(parents=True, exist_ok=True)
    individuals_dir.mkdir(parents=True, exist_ok=True)

    references = {}
    for topology in TOPOLOGIES:
        for atoms in SIZES:
            path = references_dir / f"{topology}-{atoms}.json"
            print(f"global reference: {topology}-{atoms}")
            _run([str(executable), "--global-reference", str(path), topology, str(atoms)])
            references[(topology, atoms)] = json.loads(path.read_text())

    cases = [(topology, atoms, "forward") for topology in TOPOLOGIES for atoms in SIZES]
    cases.extend((("chain", 512, "reverse"), ("cube", 512, "reverse")))
    rows = []
    for topology, atoms, order in cases:
        path = individuals_dir / f"{topology}-{atoms}-{order}.json"
        reference = references[(topology, atoms)]
        print(f"fixed-B block search: {topology}-{atoms}-{order}")
        _run([str(executable), "--scaling-case", str(path), topology, str(atoms), order,
              format(reference["objective"], ".17g")])
        row = json.loads(path.read_text())
        row["reference_peak_rss_mb"] = reference["peak_rss_mb"]
        rows.append(row)

    root.mkdir(parents=True, exist_ok=True)
    result = analyze(rows)
    summary_rows = [{key: value for key, value in row.items()
                     if key not in ("block_scaled_parameters", "sweep_objective", "sweep_global_kkt")}
                    for row in rows]
    (root / "summary.json").write_text(json.dumps({"cases": summary_rows}, indent=2) + "\n")
    (root / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    (root / "campaign-manifest.json").write_text(json.dumps({
        "phase": "F1.5 fixed-B block scaling",
        "topologies": list(TOPOLOGIES),
        "atom_counts": list(SIZES),
        "reverse_controls": ["chain-512", "cube-512"],
        "core_atoms": CORE_ATOMS,
        "block_order": "forward serial Gauss-Seidel; reverse only at 512 controls",
        "local_solver": "existing constrained SolveLinear",
        "maximum_sweeps": 200,
        "global_kkt_tolerance": GLOBAL_KKT_TOLERANCE,
        "objective_parity_tolerance": "1e-12 + 2e-12 * abs(global_objective)",
        "formal_growth_gate": "not evaluated; three multi-block sizes, trend classification only",
    }, indent=2) + "\n")
    fields = ("topology", "atoms", "order", "blocks", "sweeps", "sweeps_to_objective_parity",
              "sweeps_to_global_kkt", "block_solves", "linear_solves", "factorizations", "factor_seconds",
              "search_seconds", "peak_rss_mb", "reference_peak_rss_mb", "initial_objective", "final_objective",
              "final_global_kkt", "maximum_block_rows", "maximum_block_columns", "boundary_row_fraction",
              "success", "objective_parity", "global_kkt_passed")
    with (root / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row[field] for field in fields} for row in rows)
    print(f"fixed-B scaling={result['fixed_b_scaling']} F2 gate={result['resume_nonlinear_f2']}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    result = run_campaign(args.executable, args.output_dir)
    raise SystemExit(0 if result["fixed_b_scaling"] == "passed" else 1)


if __name__ == "__main__":
    main()
