import copy
import tempfile
import unittest
from pathlib import Path

import matched_joint_ac as ac
import mdpde_experiment as experiment


class MatchedJointACRunnerTest(unittest.TestCase):
    def test_pairing_respects_qualification_and_numerical_uncertainty(self):
        self.assertEqual(ac.pair_status(.2, .1, 0, 0, .001, True), "improved")
        self.assertEqual(ac.pair_status(.2, .1, 0, 0, .001, False), "unqualified")
        self.assertEqual(ac.pair_status(.2, .1, 0, .1, .1, True), "unresolved")
        self.assertEqual(ac.pair_status(.2, None, 0, 0, None, False), "unavailable")
        self.assertEqual(ac.pair_status(.1, .2, 0, 0, .001, True), "regressed")

    def test_case_requires_same_full_domain_fixed_b_and_exact_assembly(self):
        dataset = {"component_atoms": [[0]], "component_rows": [[0, 1, 2]], "rows": [{"owner": 0}]*3}
        context = {"state": [[6, .51, -.2]], "atoms": [{"alpha": .1}]}
        blocks = [{"owner": 0, "alpha": .1, "rows": 3, "lambda": 1.0}]
        assembled = [{"mode": m, "abc": [[6 if m == "input" else 6.1, .51, -.2]]}
                     for m in ("input", "frozen", "joint")]
        case = {"component": 0, "members": [0], "row_count": 3, "blocks": blocks, "state_id": "s", "assembled": assembled}
        base = {"schema_version": 2, "variances": [.01], "columns": 2, "rows": 3, "blocks": blocks, "state_id": "s", "component": 0,
                "beta": [6.1, -.2], "qualified": True, "uncertainty": [1e-8, 1e-8]}
        fits = [dict(base, mode="joint"), dict(base, mode="frozen", atom_index=0)]
        ac.validate_case(case, dataset, fits, context)
        changed = copy.deepcopy(case); changed["assembled"][2]["abc"][0][1] = .5
        with self.assertRaisesRegex(RuntimeError, "Width changed"):
            ac.validate_case(changed, dataset, fits, context)
        changed = copy.deepcopy(fits); changed[1]["rows"] = 2
        with self.assertRaisesRegex(RuntimeError, "domain"):
            ac.validate_case(case, dataset, changed, context)
        changed = copy.deepcopy(case); changed["assembled"][2]["abc"][0][0] = 6.2
        with self.assertRaisesRegex(RuntimeError, "Assembly"):
            ac.validate_case(changed, dataset, fits, context)

    def test_alpha_mapping_and_unexecuted_comparisons_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            experiment.write(p / "fit-index.json", {"selected_state_ids": ["executed"]})
            with self.assertRaisesRegex(RuntimeError, "not executed"):
                ac.compare(p, p, p / "comparison.json", "missing")

    def test_population_keeps_missing_and_unqualified_atoms(self):
        stats = ac.sweep.population_stats([.1, None], 168)
        self.assertEqual(stats["total"], 168)
        self.assertEqual(stats["missing"], 167)
        self.assertFalse(stats["complete_population"])

    def test_truth_scoring_is_read_only(self):
        values = [.2, .1, 0, 1e-9, 1e-9, True]
        saved = copy.deepcopy(values)
        ac.pair_status(*values)
        self.assertEqual(values, saved)
        self.assertEqual(ac.pair_status(.2, .1, .3, 1e-9, 1e-9, True), "regressed")

    def test_thread_comparison_ignores_only_elapsed_times(self):
        value = {"seconds": 1, "beta": [1, 2], "reference": {"seconds": 3, "variance": .1}, "worker_seconds": 4}
        self.assertEqual(ac.scientific(value), {"beta": [1, 2], "reference": {"variance": .1}})
        self.assertNotEqual(ac.scientific(value), ac.scientific(dict(value, beta=[1, 3])))

    def test_incomplete_matrix_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            experiment.write(p / "fit-index.json", {"states": []})
            experiment.write(p / "dataset.json", {})
            experiment.write(p / "state-index.json", {"states": []})
            with self.assertRaisesRegex(RuntimeError, "State index"):
                ac.summarize(p)


if __name__ == "__main__":
    unittest.main()
