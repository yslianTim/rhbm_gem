"""Check paired research audits reject infeasible or uncertified bounds."""

import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tests/experiments/run_joint_global_width_bound_paired.py"
SPEC = importlib.util.spec_from_file_location("global_width_bound_runner", SCRIPT)
RUNNER = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
SPEC.loader.exec_module(RUNNER)


def group(label: str, bound: float | None) -> dict:
    calibration = {"status": "Available", "p95": 1.0, "bmax2": 2.0, "bmax4": 4.0}
    bounded = label in ("B2", "B4")
    return {
        "group": label,
        "strategy": "active-set-bounded-lm" if bounded else "production-profile-lm",
        "input_fingerprint": "input",
        "model_sha256": "model",
        "map_sha256": "map",
        "target_count": 10,
        "observation_count": 100,
        "normalization_divisor": 1.0,
        "chain_selection": "Chain A non-hydrogen",
        "pure_halo_support_exclusion": True,
        "first_stage_initialization_fingerprint": "first-stage",
        "effective_initialization_fingerprint": "first-stage" if label == "A" else f"{label}-initial",
        "calibration": calibration,
        "problem_objective": "same objective",
        "candidate_available": True,
        "global_width_upper_bound": bound,
        "global_eta_upper_bound": math.log(bound) if bound else None,
        "maximum_bound_violation_log_width": 0.0,
        "independent_replay_check": {"status": "Passed"} if bounded else None,
        "width_stationarity": {"status": "Passed"} if bounded else None,
        "active_set_subproblem_failures": 0 if bounded else None,
        "truth_leakage": False,
    }


def main() -> None:
    groups = {"A": group("A", None), "B2": group("B2", 2.0), "B4": group("B4", 4.0)}
    assert RUNNER.audit_pair(groups)["all_audits_pass"]

    infeasible = {name: dict(value) for name, value in groups.items()}
    infeasible["B2"]["maximum_bound_violation_log_width"] = 1e-8
    assert not RUNNER.audit_pair(infeasible)["all_audits_pass"]

    replay_failure = {name: dict(value) for name, value in groups.items()}
    replay_failure["B4"]["independent_replay_check"] = {"status": "Failed"}
    assert not RUNNER.audit_pair(replay_failure)["all_audits_pass"]
    print("global width bound paired-audit contract passed")


if __name__ == "__main__":
    main()
