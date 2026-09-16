import copy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import estimated_neighbor_sweep as sweep
import mdpde_experiment as experiment


def row(mode="analytic", fit="ABC"):
    return {"state_id": "baseline-best-28", "serial_id": 1, "identity": {"serial_id": 1},
            "parameters": fit, "prediction": mode, "membership": [0, 1], "input": [6.02, .501, .01],
            "A": 6.1, "B": .51, "C": .01 if fit == "AB" else .1,
            "reference": {"A": 6.10000001, "B": .51000001, "C": .01 if fit == "AB" else .10000001},
            "qualified": True, "reason": "qualified", "condition": 3, "stationarity": 1e-10,
            "input_loss": .1, "loss": .01, "evaluations": 129, "reference_evaluations": 257,
            "linear_solves": 386, "seconds": 1.0}


def context_fixture():
    atoms, manifest_atoms = [], []
    for i in range(168):
        identity = {"serial_id": i+1, "chain_id": "A", "sequence_id": i+1, "component_id": "ALA",
                    "atom_id": "CA", "alternate_indicator": ".", "position": [float(i), 0, 0]}
        samples = [{"distance": .5 if p < 100 else 1.5, "position": [float(i)+p*.001, 0, 0],
                    "selected": p < 100, "response": 1.0} for p in range(200)]
        atoms.append({"index": i, "identity": identity, "samples": samples})
        manifest_atoms.append(copy.deepcopy(identity))
    return {"schema_version": 1, "atoms": atoms, "state": [[6, .5, .1] for _ in atoms]}, {"atoms": manifest_atoms}


class EstimatedNeighborSweepRunnerTest(unittest.TestCase):
    def test_population_and_identity_checks(self):
        context, manifest = context_fixture()
        sweep.validate_context(context, manifest)
        reference = copy.deepcopy(context)
        context["atoms"][0]["identity"]["atom_id"] = "CB"
        with self.assertRaisesRegex(RuntimeError, "identity mismatch"):
            sweep.validate_context(context, manifest)
        context = copy.deepcopy(reference)
        context["atoms"][0]["samples"][0]["response"] += .01
        with self.assertRaisesRegex(RuntimeError, "samples or identities differ"):
            sweep.validate_context(context, manifest, reference)
        context = copy.deepcopy(reference)
        context["state"][0][2] = float("nan")
        with self.assertRaisesRegex(RuntimeError, "Nonfinite"):
            sweep.validate_context(context, manifest)
        context = copy.deepcopy(reference)
        context["atoms"].pop()
        with self.assertRaisesRegex(RuntimeError, "population"):
            sweep.validate_context(context, manifest)

    def test_missing_checkpoint_cannot_be_replaced(self):
        with tempfile.TemporaryDirectory() as temp:
            args = SimpleNamespace(baseline_run=Path(temp), refined_run=Path(temp))
            with self.assertRaisesRegex(RuntimeError, "Missing checkpoint"):
                sweep.prepare_index(args, {}, {})

    def test_incomplete_matrix_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            experiment.write(directory / "fit-index.json", {"targets": 1344, "fits": 5376, "states": []})
            with self.assertRaisesRegex(RuntimeError, "Incomplete sweep matrix"):
                sweep.load_fits(directory)

    def test_pair_inputs_and_fixed_charge_are_checked(self):
        truth = {"1": [6, .5, 0]}
        a, b = row(), row("matched")
        b["membership"] = [0]
        with self.assertRaisesRegex(RuntimeError, "different inputs or samples"):
            sweep.paired_results([a, b], truth)
        with self.assertRaisesRegex(RuntimeError, "Incomplete prediction pair"):
            sweep.paired_results([a], truth)
        a, b = row(fit="AB"), row("matched", "AB")
        self.assertEqual(len(sweep.paired_results([a, b], truth)), 6)
        b["C"] = 0
        with self.assertRaisesRegex(RuntimeError, "fix checkpoint C"):
            sweep.paired_results([a, b], truth)

    def test_scoring_truth_cannot_mutate_fits_and_unqualified_is_not_a_win(self):
        a, b = row(), row("matched")
        b.update(A=6.001, B=.50001, C=.00001, qualified=False, reason="nonstationary")
        b["reference"] = {p: b[p]+1e-9 for p in sweep.PARAMETERS}
        rows = [a, b]; before = copy.deepcopy(rows)
        pairs = sweep.paired_results(rows, {"1": [6, .5, 0]})
        matched = [p for p in pairs if p["comparison"] == "analytic-to-matched"]
        self.assertTrue(all(p["absolute_error_gain"] > 0 and p["status"] == "unqualified" for p in matched))
        sweep.paired_results(rows, {"1": [5, .45, -.5]})
        self.assertEqual(rows, before)
        stats = sweep.fit_statistics(rows, {"1": [6, .5, 0]})
        self.assertEqual(stats[1]["qualified"], 0)
        self.assertEqual(stats[1]["C"]["estimate"]["n"], 1)

    def test_missing_parameters_keep_denominator_and_numerical_ties_are_unresolved(self):
        a, b = row(), row("matched")
        b["A"] = None
        pairs = sweep.paired_results([a, b], {"1": [6, .5, 0]})
        self.assertEqual(pairs[0]["status"], "unavailable")
        stats = sweep.fit_statistics([a, b], {"1": [6, .5, 0]})[1]["A"]["estimate"]
        self.assertEqual((stats["total"], stats["n"], stats["missing"]), (1, 0, 1))
        self.assertFalse(stats["complete_population"])
        b = row("matched")
        pairs = sweep.paired_results([a, b], {"1": [6, .5, 0]})
        self.assertTrue(all(p["status"] == "unresolved" for p in pairs if p["comparison"] == "analytic-to-matched"))

    def test_matched_advantage_does_not_imply_improvement_over_input(self):
        a, b = row(), row("matched")
        b.update(A=6.05, B=.505, C=.05)
        b["reference"] = {p: b[p]+1e-9 for p in sweep.PARAMETERS}
        pairs = sweep.paired_results([a, b], {"1": [6, .5, 0]})
        self.assertTrue(all(p["status"] == "improved" for p in pairs if p["comparison"] == "analytic-to-matched"))
        self.assertTrue(all(p["status"] == "regressed" for p in pairs if p["comparison"] == "input-to-matched"))

    def test_neighbor_components_are_added_pointwise_not_by_rmse(self):
        values = {"estimate": 2., "truth": 1., "estimation": 1., "operator": -.75, "total": .25}
        self.assertEqual(sweep.decomposition_error(values), 0)
        values["total"] = 1.75
        with self.assertRaisesRegex(RuntimeError, "decomposition identity"):
            sweep.decomposition_error(values)


if __name__ == "__main__":
    unittest.main()
