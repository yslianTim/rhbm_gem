"""Compare only confirmed FixedNeighbor endpoints across block orders."""

import argparse
import csv
import json
from pathlib import Path


EXPECTED_CASES = {(topology, atoms) for topology, atoms in (
    ("chain", 256), ("chain", 512), ("cube", 256), ("cube", 512), ("cube", 1024))}
ASSESSED_CASES = {case for case in EXPECTED_CASES if case[1] <= 512}
ETA_DIFFERENCE_THRESHOLD = 1e-10
AC_DIFFERENCE_THRESHOLD = 1e-10


def _inf_difference(lhs, rhs):
    if not isinstance(lhs, list) or not isinstance(rhs, list) or len(lhs) != len(rhs):
        return None
    return max((abs(a - b) for a, b in zip(lhs, rhs)), default=0.0)


def _direction(fixed, order):
    global_kkt = fixed.get("global_kkt")
    return {
        "block_order": order,
        "search_converged": fixed.get("search_converged") is True,
        "search_reason": fixed.get("search_reason"),
        "sweeps": fixed.get("sweeps"),
        "cheap_sweep": fixed.get("first_order_stationarity_sweep"),
        "confirmed_sweep": fixed.get("confirmed_stationarity_sweep"),
        "objective": fixed.get("objective"),
        "global_kkt": global_kkt.get("value") if isinstance(global_kkt, dict) else fixed.get("final_global_ac_kkt"),
        "width_gradient_inf_norm": fixed.get("width_gradient_inf_norm", fixed.get("final_raw_width_gradient_inf_norm")),
        "endpoint_certified": fixed.get("endpoint_certified"),
        "runtime_convergence": fixed.get("runtime_convergence"),
    }


def analyze(cases, expected=EXPECTED_CASES):
    grouped = {}
    for case in cases:
        key = (case["topology"], case["atoms"])
        order = case["fixed_neighbor"].get("block_order", case.get("block_order"))
        if order in ("forward", "reverse"):
            grouped.setdefault(key, {})[order] = case["fixed_neighbor"]

    rows = []
    failures = []
    for key in sorted(expected):
        topology, atoms = key
        directions = grouped.get(key, {})
        missing_orders = [order for order in ("forward", "reverse") if order not in directions]
        if missing_orders:
            failures.append({"topology": topology, "atoms": atoms, "reason": "missing-order"})
            rows.append({"topology": topology, "atoms": atoms, "missing_orders": missing_orders,
                "comparison_status": "not-comparable", "convergence_asymmetry": None})
            continue

        forward, reverse = directions["forward"], directions["reverse"]
        fwd, rev = _direction(forward, "forward"), _direction(reverse, "reverse")
        fwd_confirmed = (fwd["search_converged"] and isinstance(fwd["confirmed_sweep"], int) and
            0 < fwd["confirmed_sweep"] <= fwd["sweeps"])
        rev_confirmed = (rev["search_converged"] and isinstance(rev["confirmed_sweep"], int) and
            0 < rev["confirmed_sweep"] <= rev["sweeps"])
        asymmetry = fwd_confirmed != rev_confirmed
        if asymmetry:
            failures.append({"topology": topology, "atoms": atoms, "reason": "convergence-asymmetry"})
        if not fwd_confirmed or not rev_confirmed:
            if not asymmetry:
                failures.append({"topology": topology, "atoms": atoms, "reason": "both-orders-unconfirmed"})
            rows.append({"topology": topology, "atoms": atoms, "directions": {"forward": fwd, "reverse": rev},
                "comparison_status": "convergence-asymmetry" if asymmetry else "unconfirmed",
                "convergence_asymmetry": asymmetry})
            continue

        eta_difference = _inf_difference(forward.get("final_eta"), reverse.get("final_eta"))
        beta_difference = _inf_difference(forward.get("final_beta"), reverse.get("final_beta"))
        forward_weights = forward.get("final_ac_scaling_weights")
        reverse_weights = reverse.get("final_ac_scaling_weights")
        if (beta_difference is None or not isinstance(forward_weights, list) or
                not isinstance(reverse_weights, list) or len(forward_weights) != len(reverse_weights) or
                len(forward_weights) != len(forward.get("final_beta", []))):
            ac_difference = None
        else:
            ac_difference = max((abs(a - b) * max(wf, wr) for a, b, wf, wr in zip(
                forward["final_beta"], reverse["final_beta"], forward_weights, reverse_weights)), default=0.0)
        objective_difference = (abs(fwd["objective"] - rev["objective"])
            if isinstance(fwd["objective"], (int, float)) and isinstance(rev["objective"], (int, float)) else None)
        if eta_difference is None or ac_difference is None or objective_difference is None:
            failures.append({"topology": topology, "atoms": atoms, "reason": "missing-final-state-metrics"})
        elif eta_difference > ETA_DIFFERENCE_THRESHOLD or ac_difference > AC_DIFFERENCE_THRESHOLD:
            failures.append({"topology": topology, "atoms": atoms, "reason": "material-parameter-divergence"})

        assessed = key in ASSESSED_CASES
        endpoint_passed = None
        endpoint_parity = None
        runtime_parity = None
        if assessed:
            endpoint_parity = fwd["endpoint_certified"] == rev["endpoint_certified"]
            runtime_parity = fwd["runtime_convergence"] == rev["runtime_convergence"]
            endpoint_passed = all(direction["endpoint_certified"] is True and
                direction["runtime_convergence"] == "Passed" for direction in (fwd, rev))
            if not endpoint_passed:
                failures.append({"topology": topology, "atoms": atoms, "reason": "endpoint-or-runtime-failed"})

        rows.append({
            "topology": topology,
            "atoms": atoms,
            "measurement_scope": "endpoint-assessed" if assessed else "fixed-neighbor-search-only",
            "directions": {"forward": fwd, "reverse": rev},
            "comparison_status": "compared",
            "convergence_asymmetry": False,
            "objective_abs_difference": objective_difference,
            "eta_inf_difference": eta_difference,
            "scaled_ac_inf_difference": ac_difference,
            "cheap_sweep_difference": fwd["cheap_sweep"] - rev["cheap_sweep"],
            "confirmed_sweep_difference": fwd["confirmed_sweep"] - rev["confirmed_sweep"],
            "endpoint_certification_parity": endpoint_parity,
            "runtime_convergence_parity": runtime_parity,
            "endpoint_and_runtime_passed_both": endpoint_passed,
            "no_material_parameter_divergence": (eta_difference is not None and ac_difference is not None and
                eta_difference <= ETA_DIFFERENCE_THRESHOLD and ac_difference <= AC_DIFFERENCE_THRESHOLD),
        })

    missing_cases = sorted(expected - set(grouped))
    if missing_cases:
        for topology, atoms in missing_cases:
            failures.append({"topology": topology, "atoms": atoms, "reason": "missing-case"})
    return {
        "order_sensitivity_gate": "passed" if not failures else "failed",
        "confirmation_rule": "A/C KKT <=1e-10 and width-gradient <=1e-12 plus eta-change <=1e-10",
        "parameter_difference_thresholds": {
            "eta_inf": ETA_DIFFERENCE_THRESHOLD,
            "scaled_ac_inf": AC_DIFFERENCE_THRESHOLD,
            "scaled_ac_definition": "max absolute beta difference times the larger endpoint scale weight",
        },
        "expected_cases": [{"topology": t, "atoms": n} for t, n in sorted(expected)],
        "missing_cases": [{"topology": t, "atoms": n} for t, n in missing_cases],
        "failures": failures,
        "cases": rows,
    }


def write_report(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "order-analysis.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = ["topology", "atoms", "measurement_scope", "comparison_status", "objective_abs_difference",
        "eta_inf_difference", "scaled_ac_inf_difference", "cheap_sweep_difference",
        "confirmed_sweep_difference", "endpoint_certification_parity", "runtime_convergence_parity",
        "endpoint_and_runtime_passed_both", "no_material_parameter_divergence"]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(report["cases"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, help="order-sensitivity artifact directory")
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args()
    cases = [json.loads(Path(path).read_text()) for path in args.cases]
    report = analyze(cases)
    write_report(report, Path(args.output))
    print(report["order_sensitivity_gate"])
    return 0 if report["order_sensitivity_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
