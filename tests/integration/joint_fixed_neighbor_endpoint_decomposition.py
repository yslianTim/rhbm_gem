"""Classify same-width FixedNeighbor endpoint diagnostics using fixed thresholds."""
import argparse
import json
from pathlib import Path


INNER_TOLERANCE = 1e-10
WIDTH_TOLERANCE = 1e-12
LOCAL_TOLERANCE = 1e-10
RANK_BOUNDARY_FACTOR = 2.0


def _passed(state):
    assessment = state["assessment"]
    trust = assessment["endpoint_trust"]
    return (assessment["inner"] and assessment["gradient"] and assessment["local"] and
            assessment["identified"] and trust["passed"])


def _near_rank_boundary(state):
    assessment = state["assessment"]
    if assessment["identified"]:
        return False
    for key in ("projected_width", "corrected_jacobian", "normalized_width"):
        spectrum = assessment[key]
        if spectrum is None:
            continue
        minimum = spectrum["minimum_singular_value"]
        threshold = spectrum["rank_threshold"]
        if minimum is not None and threshold is not None and threshold > 0:
            if threshold / RANK_BOUNDARY_FACTOR <= minimum <= threshold * RANK_BOUNDARY_FACTOR:
                return True
    return False


def _local_correction(state):
    return state["assessment"]["local_correction_inf_norm"]


def _width_gradient(state):
    return state["width_gradient_inf_norm"]


def _materially_improved(raw, profiled):
    for metric in (_local_correction, _width_gradient):
        before, after = metric(raw), metric(profiled)
        if before is not None and after is not None and before > 0 and after <= 0.5 * before:
            return True
    return False


def classify(snapshots):
    rank_boundary = False
    width_gap = False
    mixed = False
    ac_profile_gap = False
    evidence = []
    for snapshot in snapshots:
        raw = snapshot["raw"]
        primary = snapshot["same_eta_primary"]
        reference = snapshot["same_eta_reference"]
        if any(_near_rank_boundary(state) for state in (raw, primary, reference)):
            rank_boundary = True

        raw_failed = not _passed(raw)
        primary_passed = _passed(primary)
        reference_passed = _passed(reference)
        block_primary = snapshot["coefficient_difference_block_primary"]
        primary_reference = snapshot["coefficient_difference_primary_reference"]
        profile_inner = primary["assessment"]["inner"] and reference["assessment"]["inner"]
        profile_width_or_local_failed = not (
            primary["assessment"]["gradient"] and primary["assessment"]["local"] and
            reference["assessment"]["gradient"] and reference["assessment"]["local"])

        if raw_failed and (primary_passed or reference_passed):
            ac_profile_gap = True
            evidence.append({"sweep": snapshot["sweep"], "rule": "raw-fails-profile-passes"})
        elif (block_primary > INNER_TOLERANCE and primary_reference <= INNER_TOLERANCE and
              _materially_improved(raw, primary) and primary_passed and reference_passed):
            ac_profile_gap = True
            evidence.append({"sweep": snapshot["sweep"], "rule": "same-eta-reprofile-materially-improves"})

        if (block_primary > INNER_TOLERANCE and profile_inner and profile_width_or_local_failed):
            mixed = True
            evidence.append({"sweep": snapshot["sweep"], "rule": "beta-gap-and-profile-width-or-correction-fails"})
        elif (profile_inner and profile_width_or_local_failed and
              max(_width_gradient(primary) or 0.0, _width_gradient(reference) or 0.0) > WIDTH_TOLERANCE):
            width_gap = True
            evidence.append({"sweep": snapshot["sweep"], "rule": "profiled-width-gradient-exceeds-threshold"})

    if rank_boundary and not (ac_profile_gap or mixed or width_gap):
        classification = "rank-boundary-related"
    elif ac_profile_gap and not mixed:
        classification = "ac-profile-manifold-gap"
    elif mixed:
        classification = "mixed"
    elif width_gap:
        classification = "width-stationarity-gap"
    else:
        classification = "unidentified"
    return {"classification": classification, "evidence": evidence,
        "thresholds": {"inner_coefficient_difference": INNER_TOLERANCE,
            "width_stationarity": WIDTH_TOLERANCE, "local_correction": LOCAL_TOLERANCE}}


def analyze(cases):
    selected = {}
    for case in cases:
        selected[case["topology"] + "-" + str(case["atoms"])] = case
    control = selected.get("chain-256")
    target = selected.get("cube-256")
    if target is None:
        return {"classification": "unidentified", "evidence": [], "missing": ["cube-256"]}
    result = classify(target.get("endpoint_decomposition", []))
    result["control_case"] = None if control is None else {
        "topology": "chain", "atoms": 256,
        "snapshots": len(control.get("endpoint_decomposition", []))}
    result["target_case"] = {"topology": "cube", "atoms": 256,
        "snapshots": len(target.get("endpoint_decomposition", []))}
    if control is None:
        result["missing_control"] = "chain-256"
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args()
    cases = [json.loads(Path(path).read_text()) for path in args.cases]
    result = analyze(cases)
    Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["classification"])
    return 0 if result["classification"] != "unidentified" else 1


if __name__ == "__main__":
    raise SystemExit(main())
