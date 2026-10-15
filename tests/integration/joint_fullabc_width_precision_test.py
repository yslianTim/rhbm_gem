"""Run the bounded Block 21 FullABC precision and degeneracy regression."""
import argparse
from decimal import Decimal
import json
import subprocess
import tempfile
from pathlib import Path


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def nearest(rows, width):
    return min(rows, key=lambda row: abs(row["width"] - width))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--fixture", required=True, type=Path)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="joint-fullabc-width-") as temp:
        output = Path(temp) / "study"
        subprocess.run([str(args.executable), str(args.fixture), str(output)], check=True)
        summary = json.loads((output / "summary.json").read_text())

    require(summary["local_rows"] == 930 and summary["local_columns"] == 24,
            "Block 21 fixture shape changed")
    require(summary["profile_reason"] == "kkt-failed" and not summary["profile_valid"],
            "Block 21's unqualified production profile must remain rejected")
    require(summary["production_scaled_kkt"] > 1e-10,
            "Block 21 regression must continue to fail the unchanged KKT gate")
    require(summary["full_normalized_design"]["rank"] == 24,
            "Full design must retain its observed full-rank classification")
    require(summary["free_normalized_design"]["rank"] == 23,
            "Active A face must retain its free-face rank classification")

    widths_11721 = summary["fixed_ac_profile_width_11721"]
    widths_11650 = summary["fixed_ac_profile_width_11650"]
    require(len(widths_11721) == 39 and len(widths_11650) == 24,
            "The bounded profile grids changed")
    moderate = nearest(widths_11721, 0.4942385633)
    require(moderate["profile_valid"] and moderate["inner_kkt"] <= 1e-10
            and moderate["local_full_condition"] < 20,
            "The moderate-width identifiable control must remain qualified")
    require(Decimal(moderate["fixed_double_design_hp_kkt"]) < Decimal("1e-50")
            and abs(Decimal(moderate["fixed_double_design_hp_loss"])
                    - Decimal(str(moderate["profile_loss"]))) < Decimal("1e-18"),
            "Well-conditioned binary64 profile must agree with the constrained high-precision reference")

    captured = nearest(widths_11721, 85164.2929583)
    pair = captured["pair"]
    require(pair["normalized_correlation"] == 1.0
            and pair["normalized_pair_sigma_min"] > 0.0
            and pair["normalized_pair_sigma_min"] < 4e-11,
            "Large-width G/K dependence must remain distinct from rank deficiency")
    require(captured["local_full_rank"] == 24 and captured["inner_kkt"] > 1e-10,
            "The captured full-rank A/C candidate must remain KKT-unqualified")
    require(abs(captured["normalized_A"] + captured["normalized_C"]) < 2.0
            and abs(captured["normalized_A"]) > 9e8,
            "Severe A/C cancellation regression changed")
    far_tail = nearest(widths_11721, 888611.05)
    require(abs(far_tail["profile_loss"] - captured["profile_loss"]) < 1e-8
            and abs(far_tail["normalized_A"]) > 1e11,
            "Large-width profile flatness and coefficient growth changed")

    small_captured = nearest(widths_11650, 0.123236)
    gnorm = small_captured["pair"]["gaussian_norm"]
    require(1e-39 < gnorm < 1e-37
            and summary["atom_11650_gaussian_norm_over_observation_scale"] < 1e-39,
            "Captured small-width Gaussian information changed")
    underflow = nearest(widths_11650, 0.01)
    require(underflow["pair"]["gaussian_norm"] == 0.0
            and underflow["pair"]["gaussian_nonzero_entries"] == 0
            and underflow["profile_reason"] == "rank-deficient"
            and Decimal(underflow["functional_kernel_hp_width_gradient"]) != 0,
            "Binary64 kernel underflow must remain distinct from functional response")

    reference = summary["high_precision_constrained_reference"]
    fixed = reference["fixed_rounded_double_design"]
    functional = reference["functional_kernel_design"]
    require(fixed["active_a_atom_indices"] == [5]
            and functional["active_a_atom_indices"] == [5],
            "High-precision references must retain the production A>=0 active face")
    require(Decimal(fixed["high_precision_kkt"]) < Decimal("1e-50")
            and Decimal(functional["high_precision_kkt"]) < Decimal("1e-50"),
            "Constrained Householder QR references did not qualify")
    require(Decimal(fixed["production_beta_high_precision_kkt"]) > Decimal("1e-10"),
            "High precision must continue to reject the production coefficients")
    require(fixed["rounded_hp_solution_double_beta_binary64_kkt"] < 1e-10
            and Decimal(fixed["rounded_hp_solution_double_beta_high_precision_kkt"]) > Decimal("1e-10"),
            "Near-threshold binary64/high-precision KKT discrepancy changed")

    print("Block 21 FullABC precision regression passed: full rank, unresolved double KKT, constrained QR references, large/small-width controls.")


if __name__ == "__main__":
    main()
