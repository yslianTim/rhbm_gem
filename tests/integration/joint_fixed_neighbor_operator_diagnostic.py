"""Analyze per-sweep endpoint evidence for the local OperatorPcg candidate."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


LOCAL_TOLERANCE = 1e-10
WIDTH_TOLERANCE = 1e-12


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _fixed(case):
    return case["fixed_neighbor"]


def _last_assessment(case):
    fixed = _fixed(case)
    snapshots = case.get("endpoint_assessment_by_sweep", [])
    if snapshots:
        return snapshots[-1].get("assessment", {})
    return fixed.get("endpoint_assessment", {})


def _correction(assessment):
    correction = assessment.get("correction", {})
    vector = correction.get("vector")
    if vector is None:
        return []
    return vector


def _trajectory(case):
    fixed = _fixed(case)
    assessments = {item.get("sweep"): item
                   for item in case.get("endpoint_assessment_by_sweep", [])}
    rows = []
    for sweep in fixed.get("sweep_telemetry", []):
        number = sweep.get("sweep")
        endpoint = assessments.get(number, {})
        assessment = endpoint.get("assessment", {})
        correction = assessment.get("correction", {})
        projected = assessment.get("projected_width") or {}
        jacobian = assessment.get("corrected_jacobian") or {}
        rows.append({
            "sweep": number,
            "objective": sweep.get("objective_after", sweep.get("objective")),
            "global_ac_kkt": sweep.get("global_ac_kkt"),
            "global_width_gradient_inf_norm": sweep.get("global_width_gradient_inf_norm"),
            "eta_change_inf": sweep.get("eta_change_inf"),
            "beta_scaled_change": sweep.get("beta_scaled_change"),
            "search_stop_reason": endpoint.get("search_stop_reason",
                                               fixed.get("search_reason")),
            "assessment_inner": assessment.get("inner"),
            "assessment_gradient": assessment.get("gradient"),
            "assessment_local": assessment.get("local"),
            "assessment_identified": assessment.get("identified"),
            "endpoint_trust": (assessment.get("endpoint_trust") or {}).get("passed"),
            "assessment_local_correction_inf_norm": correction.get("inf_norm",
                                                                      assessment.get("local_correction_inf_norm")),
            "reference_coefficient_difference": assessment.get("coefficient_difference"),
            "primary_gradient_inf_norm": assessment.get("primary_gradient_inf_norm"),
            "reference_gradient_inf_norm": assessment.get("reference_gradient_inf_norm"),
            "width_rank": projected.get("rank"),
            "width_minimum_singular_value": projected.get("minimum_singular_value"),
            "corrected_jacobian_rank": jacobian.get("rank"),
            "corrected_jacobian_minimum_singular_value": jacobian.get("minimum_singular_value"),
            "correction_max_coordinate": correction.get("max_coordinate"),
        })
    return rows


def classify_correction(trajectory):
    values = [row["assessment_local_correction_inf_norm"] for row in trajectory
              if _finite(row.get("assessment_local_correction_inf_norm"))]
    ranks = [(row.get("width_rank"), row.get("corrected_jacobian_rank"))
             for row in trajectory if row.get("width_rank") is not None or
             row.get("corrected_jacobian_rank") is not None]
    if len(set(ranks)) > 1:
        return "rank/active-face driven"
    if len(values) < 2:
        return "unclear"
    differences = [right - left for left, right in zip(values, values[1:])]
    scale = max(max(abs(value) for value in values), 1e-300)
    noise = max(1e-13, 0.01 * scale)
    tail = values[-min(3, len(values)):]
    if max(tail) - min(tail) <= noise:
        return "plateaued"
    if all(value < -noise for value in differences):
        return "still decreasing"
    if any(value > noise for value in differences) and any(value < -noise for value in differences):
        return "oscillatory"
    return "unclear"


def _inf_difference(lhs, rhs):
    if len(lhs) != len(rhs) or not lhs:
        return None
    return max(abs(left - right) for left, right in zip(lhs, rhs))


def _scaled_ac_difference(operator, legacy):
    beta = operator.get("beta", [])
    reference = legacy.get("beta", [])
    operator_weights = operator.get("ac_scaling_weights", [])
    legacy_weights = legacy.get("ac_scaling_weights", [])
    if len(beta) != len(reference) or len(beta) != len(operator_weights) or len(beta) != len(legacy_weights):
        return None
    return max(abs(left - right) * max(abs(left_weight), abs(right_weight))
               for left, right, left_weight, right_weight in
               zip(beta, reference, operator_weights, legacy_weights))


def _state(case):
    fixed = _fixed(case)
    assessment = _last_assessment(case)
    correction = assessment.get("correction", {})
    return {
        "eta": fixed.get("final_eta", []),
        "beta": fixed.get("final_beta", []),
        "objective": fixed.get("objective"),
        "ac_scaling_weights": fixed.get("final_ac_scaling_weights", []),
        "assessment": assessment,
        "correction_vector": correction.get("vector", []),
        "runtime_convergence": fixed.get("runtime_convergence"),
        "endpoint_certified": fixed.get("endpoint_certified"),
    }


def _largest_corrections(state, limit=5):
    vector = state.get("correction_vector", [])
    eta = state.get("eta", [])
    entries = []
    for atom, value in enumerate(vector):
        if not _finite(value):
            continue
        width = math.exp(eta[atom]) if atom < len(eta) and _finite(eta[atom]) else None
        entries.append({"atom": atom, "value": value, "width": width,
                        "absolute_value": abs(value)})
    entries.sort(key=lambda item: (-item["absolute_value"], item["atom"]))
    for item in entries:
        item.pop("absolute_value")
    return entries[:limit]


def _state_summary(state):
    assessment = state["assessment"]
    correction = assessment.get("correction", {})
    return {
        "eta": state["eta"],
        "beta": state["beta"],
        "objective": state["objective"],
        "reference_coefficient_difference": assessment.get("coefficient_difference"),
        "primary_gradient_inf_norm": assessment.get("primary_gradient_inf_norm"),
        "reference_gradient_inf_norm": assessment.get("reference_gradient_inf_norm"),
        "correction_inf_norm": correction.get("inf_norm", assessment.get("local_correction_inf_norm")),
        "correction_max_coordinate": correction.get("max_coordinate"),
        "largest_corrections": _largest_corrections(state),
        "runtime_convergence": state["runtime_convergence"],
        "endpoint_certified": state["endpoint_certified"],
    }


def analyze_case(operator_case, legacy_case):
    operator_fixed = _fixed(operator_case)
    legacy_fixed = _fixed(legacy_case)
    if (operator_case.get("topology"), operator_case.get("atoms")) != \
            (legacy_case.get("topology"), legacy_case.get("atoms")):
        raise ValueError("operator and legacy cases must have the same topology and atom count")
    operator = _state(operator_case)
    legacy = _state(legacy_case)
    operator_vector = _correction(operator["assessment"])
    legacy_vector = _correction(legacy["assessment"])
    comparison = {
        "eta_inf_difference": _inf_difference(operator["eta"], legacy["eta"]),
        "scaled_ac_inf_difference": _scaled_ac_difference(
            {"beta": operator["beta"], "ac_scaling_weights": operator["ac_scaling_weights"]},
            {"beta": legacy["beta"], "ac_scaling_weights": legacy["ac_scaling_weights"]}),
        "objective_difference": (abs(operator["objective"] - legacy["objective"])
                                 if _finite(operator["objective"]) and _finite(legacy["objective"])
                                 else None),
        "correction_vector_inf_difference": _inf_difference(operator_vector, legacy_vector),
        "reference_coefficient_difference": {
            "operator": operator["assessment"].get("coefficient_difference"),
            "legacy": legacy["assessment"].get("coefficient_difference"),
        },
    }
    trajectory = _trajectory(operator_case)
    return {
        "topology": operator_case.get("topology"),
        "atoms": operator_case.get("atoms"),
        "policy": {
            "outer_search": "FixedNeighbor",
            "outer_core_atoms": operator_fixed.get("outer_core_atoms",
                                                   operator_fixed.get("core_atoms",
                                                                      operator_case.get("core_atoms"))),
            "block_order": operator_fixed.get("block_order"),
            "local_work": operator_fixed.get("local_work_policy"),
            "operator_local_search": operator_fixed.get("local_search_method"),
            "operator_preconditioner": operator_fixed.get("local_preconditioner"),
            "local_schwarz_core_atoms": operator_fixed.get("local_schwarz_core_atoms"),
            "local_schwarz_overlap_hops": operator_fixed.get("local_schwarz_overlap_hops"),
            "local_schwarz_max_block_atoms": operator_fixed.get("local_schwarz_max_block_atoms"),
        },
        "operator": _state_summary(operator),
        "legacy": _state_summary(legacy),
        "comparison": comparison,
        "correction_diagnosis": {
            "classification": classify_correction(trajectory),
            "trajectory": trajectory,
            "thresholds": {"local_correction": LOCAL_TOLERANCE,
                           "width_gradient": WIDTH_TOLERANCE},
        },
    }


def write_report(report, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    summary = output.with_suffix(".csv")
    columns = ["topology", "atoms", "classification", "operator_correction", "legacy_correction",
               "eta_inf_difference", "scaled_ac_inf_difference", "objective_difference",
               "correction_vector_inf_difference"]
    with summary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            writer.writerow({
                "topology": case["topology"], "atoms": case["atoms"],
                "classification": case["correction_diagnosis"]["classification"],
                "operator_correction": case["operator"]["correction_inf_norm"],
                "legacy_correction": case["legacy"]["correction_inf_norm"],
                "eta_inf_difference": case["comparison"]["eta_inf_difference"],
                "scaled_ac_inf_difference": case["comparison"]["scaled_ac_inf_difference"],
                "objective_difference": case["comparison"]["objective_difference"],
                "correction_vector_inf_difference": case["comparison"]["correction_vector_inf_difference"],
            })


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--case", action="append", required=True,
                        help="topology:operator-json:legacy-json; may be repeated")
    args = parser.parse_args(argv)
    cases = []
    for specification in args.case:
        topology, operator_path, legacy_path = specification.split(":", 2)
        operator_case = json.loads(Path(operator_path).read_text())
        legacy_case = json.loads(Path(legacy_path).read_text())
        if operator_case.get("topology") != topology or legacy_case.get("topology") != topology:
            parser.error(f"case topology does not match {specification}")
        cases.append(analyze_case(operator_case, legacy_case))
    report = {"phase": "P1 fixed-neighbor operator endpoint correction attribution",
              "cases": sorted(cases, key=lambda item: (item["topology"], item["atoms"]))}
    write_report(report, args.output)
    print(json.dumps({"output": str(args.output),
                      "cases": len(report["cases"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
