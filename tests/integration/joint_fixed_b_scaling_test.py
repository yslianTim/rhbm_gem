import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from joint_fixed_b_scaling import classify_trend, compare_order


class FixedBScalingAnalysisTest(unittest.TestCase):
    def test_classifies_bounded_three_point_sweep_trend(self):
        rows = [
            {"atoms": atoms, "order": "forward", "global_kkt_passed": True,
             "sweeps_to_global_kkt": sweeps}
            for atoms, sweeps in ((256, 8), (512, 9), (1024, 10))
        ]
        self.assertEqual(classify_trend(rows)["classification"], "stable-looking")

    def test_classifies_growth_looking_trend(self):
        rows = [
            {"atoms": atoms, "order": "forward", "global_kkt_passed": True,
             "sweeps_to_global_kkt": sweeps}
            for atoms, sweeps in ((256, 8), (512, 10), (1024, 13))
        ]
        self.assertEqual(classify_trend(rows)["classification"], "growth-looking")

    def test_marks_missing_or_unconverged_points_insufficient(self):
        missing = [{"atoms": 256, "order": "forward", "global_kkt_passed": True,
                    "sweeps_to_global_kkt": 8}]
        self.assertEqual(classify_trend(missing)["classification"], "insufficient-evidence")
        unconverged = [
            {"atoms": atoms, "order": "forward", "global_kkt_passed": atoms != 512,
             "sweeps_to_global_kkt": 8}
            for atoms in (256, 512, 1024)
        ]
        self.assertEqual(classify_trend(unconverged)["classification"], "insufficient-evidence")

    def test_compares_only_converged_order_endpoints(self):
        forward = {
            "global_kkt_passed": True, "final_objective": 1.0e-10, "global_objective": 1.0e-10,
            "sweeps": 8, "block_scaled_parameters": [1.0, 2.0],
        }
        reverse = {
            "global_kkt_passed": True, "final_objective": 1.0e-10 + 1.0e-15,
            "global_objective": 1.0e-10, "sweeps": 9, "block_scaled_parameters": [1.0 + 1.0e-9, 2.0],
        }
        result = compare_order(forward, reverse)
        self.assertTrue(result["compared_converged_endpoints"])
        self.assertTrue(result["passed"])
        reverse["global_kkt_passed"] = False
        failed = compare_order(forward, reverse)
        self.assertFalse(failed["compared_converged_endpoints"])
        self.assertEqual(failed["status"], "non-convergence-asymmetry")


if __name__ == "__main__":
    unittest.main()
