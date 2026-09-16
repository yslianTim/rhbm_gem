import copy
import tempfile
import unittest
from pathlib import Path

import mdpde_experiment as experiment
import observation_matching as matched


def row(prediction="analytic", qualified=True):
    return {"case": "fold", "h": .1, "serial_id": 1, "layer": "file", "prediction": prediction,
            "parameters": "ABC", "truth": [6, .5, 0], "A": 6.1, "B": .51, "C": -.1,
            "reference": {"A": 6.100001, "B": .510001, "C": -.100001}, "membership": [0, 1],
            "qualified": qualified, "reason": "qualified" if qualified else "nonstationary",
            "condition": 3, "evaluations": 150, "reference_evaluations": 280, "linear_solves": 430}


class ObservationMatchingRunnerTest(unittest.TestCase):
    def test_self_consistency_cannot_pass_with_wrong_parameters(self):
        value = row()
        value["layer"] = "analytic"
        self.assertFalse(matched.control_checks([value])["passed"])
        value.update(A=6, B=.5, C=0)
        self.assertTrue(matched.control_checks([value])["passed"])

    def test_unqualified_fits_remain_in_statistics(self):
        values = [row(), row(qualified=False)]
        values[1]["A"] = 7
        summary = matched.fit_statistics(values)[0]
        self.assertEqual(summary["total"], 2)
        self.assertEqual(summary["qualified"], 1)
        self.assertEqual(summary["A"]["n"], 2)
        self.assertGreater(summary["A"]["rmse"], .7)

    def test_missing_fit_population_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            experiment.write(root / "fit-index.json", {"targets": 0, "fits": 0})
            with self.assertRaisesRegex(RuntimeError, "Incomplete experiment matrix"):
                matched.load_fits(root)

    def test_pair_cannot_change_membership(self):
        a, b = row(), row("matched")
        b["membership"] = [0]
        with self.assertRaisesRegex(RuntimeError, "different samples"):
            matched.paired_results([a, b])

    def test_zero_charge_uses_absolute_error_and_keeps_failure(self):
        a, b = row(), row("matched", qualified=False)
        b["C"] = -.0001
        b["reference"]["C"] = -.00010001
        pair = matched.paired_results([a, b])[0]
        self.assertEqual(pair["C"]["status"], "improved")
        self.assertFalse(pair["both_qualified"])
        self.assertNotIn("relative_error", pair["C"])

    def test_better_accuracy_does_not_override_qualification_or_numerical_uncertainty(self):
        rows = []
        for parameters in ("AB", "ABC"):
            for prediction in ("analytic", "matched"):
                for i in range(168):
                    r = row(prediction)
                    r["parameters"], r["serial_id"] = parameters, i + 1
                    if prediction == "matched":
                        r.update(A=6.0001, B=.5001, C=.0001)
                        r["reference"] = {p: r[p] + 1e-8 for p in matched.PARAMETERS}
                    rows.append(r)
        self.assertTrue(matched.integration_decision(matched.fit_statistics(rows))["recommend_peeling_experiment"])
        rows[-1]["qualified"] = False
        self.assertFalse(matched.integration_decision(matched.fit_statistics(rows))["recommend_peeling_experiment"])
        rows[-1]["qualified"] = True
        uncertain = copy.deepcopy(rows)
        for r in uncertain:
            r["reference"]["A"] = r["A"] + 1
        self.assertFalse(matched.integration_decision(matched.fit_statistics(uncertain))["recommend_peeling_experiment"])


if __name__ == "__main__":
    unittest.main()
