import copy
import csv
import tempfile
import unittest
from pathlib import Path

import atom_centered_voxel_union as union
import unique_stencil_grid_test as previous
import mdpde_experiment as experiment


class AtomCenteredVoxelUnionRunnerTest(unittest.TestCase):
    def test_scope_is_one_state_and_explicit_alpha(self):
        index = {"selected_state_ids": ["baseline-best-28"], "alphas": [0, .1, .5, 1]}
        self.assertEqual(len(union.expected_cases(index)), 4)
        self.assertEqual(union.expected_cases(dict(index, alphas=[0])), ["baseline-best-28-alpha-0"])
        for invalid in ([0, 0], [.1, 0], [], [.2]):
            with self.assertRaises(RuntimeError): union.expected_cases(dict(index, alphas=invalid))
        with self.assertRaises(RuntimeError): union.expected_cases(dict(index, selected_state_ids=["a", "b"]))

    def test_geometry_uses_generation_not_header_and_deduplicates(self):
        data = {"grid_size": [9, 9, 9], "generation_origin": [-2., -2., -2.], "generation_spacing": [.5]*3,
                "sampling_origin": [-2.1]*3, "sampling_spacing": [.7]*3, "radius": 2.5,
                "atoms": [{"position": [0, 0, 0]}, {"position": [0, 0, 0]}]}
        indices, count, nearest = union.sphere_membership(data)
        self.assertEqual(len(indices), len(set(indices)))
        self.assertTrue(all(count == 2)); self.assertEqual(min(nearest), 0)
        for i, distance in zip(indices, nearest):
            xyz = [i % 9, i//9 % 9, i//81]
            expected = sum((-2+.5*x)**2 for x in xyz)**.5
            self.assertAlmostEqual(distance, expected)
            self.assertLessEqual(distance, 2.5)

    def test_duplicate_and_missing_voxels_are_rejected(self):
        data = {"experiment": union.NAME, "radius": 2.5, "membership_geometry": "generation",
                "grid_size": [3]*3, "generation_origin": [-1.]*3, "generation_spacing": [1.]*3,
                "sampling_origin": [-1.]*3, "sampling_spacing": [1.]*3,
                "atoms": [{"position": [0, 0, 0]}], "row_count": 27, "alphas": union.ALPHAS}
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            (directory/"slots.csv").write_text("")
            (directory/"voxels.csv").write_text("row,index,x,y,z,observed,multiplicity,reference_double,quantization_bound,nearest_distance,in_stencil\n0,0,-1,-1,-1,0,1,0,0,1.7320508075688772,0\n1,0,-1,-1,-1,0,1,0,0,1.7320508075688772,0\n")
            data.update(voxel_table_sha256=union.fold.sha256_file(directory/"voxels.csv"), slot_table_sha256=union.fold.sha256_file(directory/"slots.csv"))
            with self.assertRaisesRegex(RuntimeError, "voxel index"): union.validate_dataset(directory, data, 27)

    def test_b_and_qualification_evidence_are_required(self):
        fit, context = previous.UniqueStencilGridRunnerTest().fixture()
        fit.update(experiment=union.NAME, svd_preconditioner="independent-tsqr-8192")
        union.validate_fit(fit, context, .1, 30)
        changed = copy.deepcopy(fit); changed["abc"][0][1] += .01
        with self.assertRaisesRegex(RuntimeError, "Width changed"): union.validate_fit(changed, context, .1, 30)
        changed = copy.deepcopy(fit); changed["qualified"] = True
        with self.assertRaisesRegex(RuntimeError, "Unverified"): union.validate_fit(changed, context, .1, 30)
        changed = dict(fit, svd_preconditioner="primary-R")
        with self.assertRaisesRegex(RuntimeError, "independent SVD"): union.validate_fit(changed, context, .1, 30)

    def test_missing_completed_case_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            experiment.write(p/"state-index.json", {"selected_state_ids": ["baseline-best-28"], "alphas": [0]})
            experiment.write(p/"completion.json", {"complete": False, "cases": []})
            experiment.write(p/"fit-index.json", {"cases": [], "alphas": [0]})
            with self.assertRaisesRegex(RuntimeError, "Missing requested cases"): union.summarize(p)

    def test_replay_exclusions_do_not_hide_scientific_changes(self):
        data = {"seconds": 1, "design_build_seconds": 2, "independent_validation_seconds": 3,
                "peak_child_rss_bytes": 4, "iterations": 100, "variance": .01, "beta": [2]}
        self.assertEqual(union.scientific(data), {"iterations": 100, "variance": .01, "beta": [2]})
        self.assertNotEqual(union.scientific(data), union.scientific(dict(data, variance=.02)))


if __name__ == "__main__": unittest.main()
