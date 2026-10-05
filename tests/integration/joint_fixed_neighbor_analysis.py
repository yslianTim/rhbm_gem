import argparse
import json
from pathlib import Path


EXPECTED_CASES = {(topology, atoms) for topology in ("chain", "cube") for atoms in (32, 128, 256)}


def analyze(cases, expected=EXPECTED_CASES):
    rows = []
    seen = set()
    for case in cases:
        key = (case["topology"], case["atoms"])
        seen.add(key)
        legacy = case["global_legacy_compact"]
        operator = case["global_operator_pcg"]
        fixed = case["fixed_neighbor"]
        comparison = next((item for item in case["comparisons"] if item["method"] == "FixedNeighbor"), None)
        checks = [legacy, operator, fixed]
        reasons = []
        if any(result["runtime_convergence"] != "Passed" for result in checks):
            reasons.append("runtime-convergence")
        if not fixed["search_converged"] or not fixed["endpoint_certified"]:
            reasons.append("search-or-endpoint-certification")
        first_order_sweep = fixed.get("first_order_stationarity_sweep", 0)
        confirmed_sweep = fixed.get("confirmed_stationarity_sweep", 0)
        if (first_order_sweep == 0 or confirmed_sweep == 0 or confirmed_sweep < first_order_sweep or
                confirmed_sweep > fixed["sweeps"]):
            reasons.append("stationarity-unconfirmed")
        for name in ("assessment_inner", "assessment_gradient", "assessment_local", "assessment_identified"):
            if fixed[name] is None or fixed[name]["status"] != "Passed":
                reasons.append(name)
        if comparison is None:
            reasons.append("missing-legacy-comparison")
        else:
            tolerance = 1e-12 + 2e-12 * abs(legacy["objective"])
            if comparison["objective_difference"] > tolerance:
                reasons.append("objective-parity")
            if comparison["scaled_ac_inf_difference"] > 1e-8:
                reasons.append("scaled-ac-parity")
        rows.append({"topology": key[0], "atoms": key[1], "passed": not reasons, "reasons": reasons,
            "sweeps": fixed["sweeps"], "objective": fixed["objective"],
            "first_order_stationarity_sweep": first_order_sweep,
            "confirmed_stationarity_sweep": confirmed_sweep,
            "confirmation_extra_sweeps": confirmed_sweep-first_order_sweep if confirmed_sweep else None,
            "global_kkt": fixed["global_kkt"]["value"], "width_gradient_inf_norm": fixed["width_gradient_inf_norm"],
            "search_seconds": fixed["search_seconds"], "peak_rss_mb": case["peak_rss_mb"]})
    missing = sorted(expected - seen)
    return {"fixed_neighbor_f2_gate": "passed" if not missing and all(row["passed"] for row in rows) else "failed",
        "expected_cases": [{"topology": topology, "atoms": atoms} for topology, atoms in sorted(expected)],
        "missing_cases": [{"topology": topology, "atoms": atoms} for topology, atoms in missing],
        "cases": rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args()
    cases = [json.loads(Path(path).read_text()) for path in args.cases]
    result = analyze(cases)
    Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["fixed_neighbor_f2_gate"])
    return 0 if result["fixed_neighbor_f2_gate"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
