import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import joint_abc_coverage as coverage
import numpy as np


class JointABCCoverageRunnerTest(unittest.TestCase):
    def initialization(self):
        atoms = [{"serial_id": 8}, {"serial_id": 3}]
        value = {"source": "first-stage-float32-map", "sampling_method": "FibonacciDeterministic", "jobs": 1,
                 "seed_abc": [0, 1, 0], "uses_peeling": False, "uses_truth": False, "valid": True,
                 "reason": "valid-widths", "b0": [.6, .7], "unqualified_local_count": 1,
                 "atoms": [{"identity": a, "mdpde": [2, b, 0], "ols": [2.1, b+.01, 0], "alpha": .2,
                            "sample_count": 10, "raw_sample_count": 20, "diagnostics": {"qualification": k}}
                           for k, (a, b) in enumerate(zip(atoms, [.6, .7]))]}
        return atoms, value

    def test_first_stage_initialization_keeps_unqualified_but_legal_widths(self):
        atoms, value = self.initialization()
        np.testing.assert_array_equal(coverage.validate_initialization(value, atoms), [.6, .7])
        np.testing.assert_allclose(coverage.initial_widths(value["b0"], atoms, "mixed"), [.72, .56])
        with self.assertRaisesRegex(RuntimeError, "Invalid first-stage B"):
            coverage.initial_widths([0, .7], atoms, "first-stage")

    def test_wrong_identity_or_truth_start_is_rejected(self):
        atoms, value = self.initialization(); wrong = copy.deepcopy(value)
        wrong["atoms"][0]["identity"]["serial_id"] = 1
        with self.assertRaisesRegex(RuntimeError, "identity mismatch"): coverage.validate_initialization(wrong, atoms)
        value["uses_truth"] = True
        with self.assertRaisesRegex(RuntimeError, "contract"): coverage.validate_initialization(value, atoms)

    def test_initial_width_must_match_mdpde_and_illegal_width_cannot_pass(self):
        atoms, value = self.initialization(); value["b0"][0] = .5
        with self.assertRaisesRegex(RuntimeError, "does not come"): coverage.validate_initialization(value, atoms)
        atoms, value = self.initialization(); value["b0"][0] = 0
        with self.assertRaisesRegex(RuntimeError, "validity"): coverage.validate_initialization(value, atoms)
        value["valid"] = False; value["atoms"][0]["mdpde"][1] = 0
        self.assertIsNone(coverage.validate_initialization(value, atoms))

    def test_wrong_width_contract_rejected_even_with_matching_hash_adapter(self):
        fixture = {"input_hashes": {}, "generation_contract": {}, "width_contract": {"b": [.5]*168}}
        manifest = {"atoms": [{"serial_id": k+1, "preparation_index": k, "element": 8} for k in range(168)]}
        with patch.object(coverage, "read", return_value=fixture), \
                patch.object(coverage.fixed.fold, "validate_input_hashes", return_value={}), \
                patch.object(coverage.fixed.fold, "load_simulation_manifest", return_value=manifest):
            with self.assertRaisesRegex(RuntimeError, "width rule"):
                coverage.validate_inputs({"manifest": Path("unused")}, fixture)

    def test_missing_case_rejected_before_dataset_scoring(self):
        with patch.object(coverage, "read", return_value={"complete": True, "cases": coverage.CASES[:-1]}):
            with self.assertRaisesRegex(RuntimeError, "Missing completed coverage case"):
                coverage.summarize_dataset(Path("unused"), "baseline", {}, {}, {})

    def test_independent_svd_exposes_zero_signal_and_duplicate_rank(self):
        axis = np.arange(-2.4, 4.9, .6)
        points = np.array(list(__import__("itertools").product(axis, axis, axis)))
        for name in ("baseline", "zero-signal", "duplicate"):
            atoms, truth = coverage.synthetic_truth(name)
            y = coverage.fixed.predict(coverage.fixed.columns(points, atoms, truth[:, 1]), truth[:, [0, 2]].ravel(), len(points))
            evidence = coverage.truth_identifiability(points, atoms, truth, y)
            if name == "baseline":
                self.assertEqual(evidence["design_rank"], 24); self.assertEqual(evidence["width_rank"], 12)
            elif name == "zero-signal":
                self.assertEqual(evidence["width_column_norms"][5], 0); self.assertEqual(evidence["width_rank"], 11)
            else: self.assertLess(evidence["design_rank"], 24)

    def test_fused_coordinates_resolve_hard_support_boundary_without_radius_slack(self):
        atoms, _ = coverage.synthetic_truth("near-0.10")
        dataset = {"atoms": atoms, "grid_size": [33, 29, 29],
                   "generation_origin": [-3.6]*3, "generation_spacing": [.3]*3}
        ids, counts, _ = coverage.sphere_membership(dataset)
        self.assertEqual(len(ids), 7243); self.assertEqual(int(counts.sum()), 29595)
        axes = coverage.grid_axes(dataset)
        indices = np.array([[7, 16, 7], [19, 9, 8]])
        points = np.column_stack([axes[k][indices[:, k]] for k in range(3)])
        basis = coverage.columns(points, [atoms[3]], [.5])
        np.testing.assert_array_equal(basis[0][0], [1])
        # An outward rounded boundary point is excluded, never tolerated inward.
        self.assertGreater(coverage.squared_distances(points-atoms[3]["position"])[0], 6.25)

    def test_failed_branches_stay_in_comparisons_but_cannot_be_representative(self):
        fits = {c: {"joint_qualified": True, "primary": {"rss": k+1}} for k, c in enumerate(coverage.CASES)}
        estimates = {c: np.array([[1., .5, -.2]]) for c in coverage.CASES}
        fits[coverage.CASES[0]]["joint_qualified"] = False
        ok, pairs, representatives = coverage.compare_starts(fits, estimates)
        self.assertFalse(ok); self.assertEqual(len(pairs), 12)
        self.assertEqual(representatives["double"], "narrower-double")
        self.assertFalse(pairs[0]["both_qualified"])

    def test_fresh_output_required(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(RuntimeError, "fresh output"):
                coverage.run(type("Args", (), {"output": Path(directory)})())

    def test_numerical_replay_failure_is_retained_without_promoting_endpoint(self):
        failures = []
        coverage.numerical_check(False, failures, False, "Independent A/C KKT differs.")
        self.assertEqual(failures, ["Independent A/C KKT differs."])
        with self.assertRaisesRegex(RuntimeError, "KKT differs"):
            coverage.numerical_check(True, [], False, "Independent A/C KKT differs.")


if __name__ == "__main__": unittest.main()
