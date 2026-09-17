import copy
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import fixed_b_oracle as oracle


class FixedBOracleRunnerTest(unittest.TestCase):
    def context(self):
        atoms = [{**{k: "A" for k in oracle.IDENTITY}, "serial_id": k+1, "position": [k, 0, 0]} for k in range(168)]
        manifest = {"atoms": atoms}
        context = {"schema_version": 1, "atoms": [{"identity": a, "index": k} for k, a in enumerate(atoms)],
                   "state": [[6, .5+k*1e-6, -.2] for k in range(168)]}
        return manifest, context

    def test_checkpoint_pairs_by_identity_and_only_reads_b(self):
        manifest, context = self.context()
        expected = oracle.checkpoint_widths(manifest, context)
        context["atoms"].reverse(); context["state"].reverse()
        for k, atom in enumerate(context["atoms"]):
            atom["index"] = k; context["state"][k][0] = "unused"; context["state"][k][2] = "unused"
        self.assertEqual(oracle.checkpoint_widths(manifest, context), expected)

    def test_duplicate_identity_changed_position_and_invalid_width_fail(self):
        manifest, context = self.context()
        for change in ("duplicate", "position", "width"):
            c = copy.deepcopy(context)
            if change == "duplicate": c["atoms"][0]["identity"] = c["atoms"][1]["identity"]
            elif change == "position": c["atoms"][0]["identity"]["position"][0] += 1
            else: c["state"][0][1] = 0
            with self.assertRaises(RuntimeError): oracle.checkpoint_widths(manifest, c)

    def test_wrong_hash_rejected_before_manifest_load(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)/"model"; p.write_text("bad")
            with self.assertRaisesRegex(RuntimeError, "SHA-256"):
                oracle.validate_inputs({"model": p}, {"input_hashes": {"model": "0"*64}})

    def test_existing_output_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(RuntimeError, "fresh output"):
                oracle.run(SimpleNamespace(output=Path(temp)))

    def test_missing_case_is_not_a_numerical_failure(self):
        fixture = {"input_hashes": {}}
        data = {"fixture.json": fixture, "inputs.json": {}, "input-hashes.json": {},
                "completion.json": {"complete": True, "cases": oracle.CASES[:-1]}}
        with patch.object(oracle, "read", side_effect=lambda p: fixture if p == oracle.FIXTURE else data[p.name]), \
                patch.object(oracle, "validate_inputs", return_value=({}, [])):
            with self.assertRaisesRegex(RuntimeError, "Missing completed case"):
                oracle.summarize(Path("unused"))

    def test_basis_cutoff_center_and_negative_charge(self):
        points = np.array([[0., 0, 0], [2.5, 0, 0], [2.50001, 0, 0]])
        basis = oracle.columns(points, [{"position": [0, 0, 0]}], [.5])
        self.assertEqual(basis[0][0].tolist(), [0, 1])
        y = oracle.predict(basis, [0, -1], len(points))
        self.assertLess(y[0], 0); self.assertLess(y[1], 0); self.assertEqual(y[2], 0)

    def test_quantization_sign_and_mismatch(self):
        y64 = np.array([1.00000007, 1.00000003, -1.00000007, 0.])
        y32 = y64.astype(np.float32).astype(float); delta = y32-y64
        oracle.validate_quantization(y64, y32, delta)
        self.assertTrue(np.any(delta < 0) and np.any(delta > 0))
        for observed, q in ((y32+1e-12, delta), (y32, np.abs(delta))):
            with self.assertRaisesRegex(RuntimeError, "Quantization"):
                oracle.validate_quantization(y64, observed, q)

    def test_fit_rejects_changed_b_before_scoring(self):
        fit = {"case": oracle.CASES[0], "experiment": "fixed-b-oracle", "alpha": 0, "b": [.51]*168, "row_count": 1}
        with self.assertRaisesRegex(RuntimeError, "Changed fixed B"):
            oracle.validate_fit(Path("unused"), oracle.CASES[0], fit, [], np.ones(1), [.5]*168)

    def test_repeat_excludes_only_timing_from_scientific_files(self):
        data = {"seconds": 1, "primary": {"seconds": 2, "beta": [1]}, "qualified": False, "oracle_recovery": {"passed": False}}
        self.assertEqual(oracle.scientific(data), {"primary": {"beta": [1]}, "qualified": False, "oracle_recovery": {"passed": False}})
        self.assertNotEqual(oracle.scientific(data), oracle.scientific(dict(data, qualified=True)))


if __name__ == "__main__": unittest.main()
