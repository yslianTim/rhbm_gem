import copy
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import joint_abc_profile as joint
import numpy as np


class JointABCProfileRunnerTest(unittest.TestCase):
    def test_starts_use_serial_identity_without_truth_or_common_width(self):
        b = np.array([.51, .49, .53]); atoms = [{"serial_id": n} for n in (8, 3, 5)]
        for start, factors in (("checkpoint", [1, 1, 1]), ("narrower", [.8]*3),
                               ("wider", [1.2]*3), ("mixed", [1.2, .8, .8])):
            np.testing.assert_array_equal(joint.initial_widths(b, atoms, start), b*factors)
        np.testing.assert_array_equal(b, [.51, .49, .53])

    def test_raw_gradient_matches_objective_with_fixed_ac(self):
        points = np.array([[0., 0, 0], [1e-6, 0, 0], [.3, .2, 0], [1.1, 0, 0], [2.5, 0, 0], [3, 0, 0]])
        atoms = [{"position": [0, 0, 0]}]; endpoint = {"beta": [1.2, -.3], "b": [.5]}; y = np.linspace(-.2, .3, len(points))
        computed = joint.raw_certificate(points, atoms, y, endpoint)
        h = 1e-5; plus = dict(endpoint, b=[.5*np.exp(h)]); minus = dict(endpoint, b=[.5*np.exp(-h)])
        rp = joint.raw_certificate(points, atoms, y, plus)["residual"]
        rm = joint.raw_certificate(points, atoms, y, minus)["residual"]
        self.assertAlmostEqual(computed["b_gradient"][0], (rp@rp-rm@rm)/(4*h*max(1, np.linalg.norm(y))**2), places=7)
        self.assertEqual(computed["prediction"][-1], 0)

    def test_fresh_output_and_wrong_hash_still_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(RuntimeError, "fresh output"):
                joint.fixed.run(SimpleNamespace(output=Path(tmp)), "joint-abc-profile", joint.summarize)
            path = Path(tmp)/"model"; path.write_text("bad")
            with self.assertRaisesRegex(RuntimeError, "SHA-256"):
                joint.fixed.validate_inputs({"model": path}, {"input_hashes": {"model": "0"*64}})

    def test_missing_case_rejected_before_scientific_scoring(self):
        fixture = {}; data = {"fixture.json": fixture, "inputs.json": {}, "input-hashes.json": {},
                              "completion.json": {"complete": True, "execution_complete": True, "cases": joint.CASES[:-1]}}
        with patch.object(joint, "read", side_effect=lambda p: fixture if p == joint.fixed.FIXTURE else data[p.name]), \
                patch.object(joint.fixed, "validate_inputs", return_value=({}, [])):
            with self.assertRaisesRegex(RuntimeError, "Missing completed joint case"):
                joint.summarize(Path("unused"))

    def candidate(self):
        endpoint = {"valid": True, "feasible": True, "kkt_passed": True, "projected_kkt": 0, "b_gradient_inf": 0,
                    "beta": [1, -.1, 2, .2], "b": [.5, .5]}
        check = {"passed": True, "same_active_face": True, "plus_valid": True, "minus_valid": True, "relative_l2_difference": 1e-8}
        return {"primary": endpoint, "reference": copy.deepcopy(endpoint), "scaled_reference_difference": 0,
                "design_spectrum": {"rank": 4}, "local_correction_inf": 0, "width_spectrum": {"rank": 2},
                "profile_jacobian_spectrum": {"rank": 2}, "derivative_checks": [copy.deepcopy(check) for _ in range(6)],
                "derivative_verified": True, "qualification_checks": {k: True for k in ("inner", "b_gradient", "local_correction", "identified", "derivative")}}

    def test_native_lm_success_cannot_override_failed_joint_checks(self):
        fit = self.candidate(); fit["lm_status"] = 1; self.assertTrue(joint.qualification(fit))
        fit["primary"]["valid"] = False; self.assertFalse(joint.qualification(fit))
        fit = self.candidate(); fit["local_correction_inf"] = 1e-8; fit["qualification_checks"]["local_correction"] = False
        self.assertFalse(joint.qualification(fit))
        fit = self.candidate(); fit["width_spectrum"]["rank"] = 1; fit["qualification_checks"]["identified"] = False
        self.assertFalse(joint.qualification(fit))
        fit = self.candidate(); fit["derivative_checks"][0]["same_active_face"] = False
        fit["derivative_checks"][0]["passed"] = False; fit["derivative_verified"] = False; fit["qualification_checks"]["derivative"] = False
        self.assertFalse(joint.qualification(fit))

    def test_representative_uses_rss_and_all_starts_must_qualify(self):
        fits = {name: {"joint_qualified": True, "primary": {"rss": k+1.}} for k, name in enumerate(joint.CASES)}
        abc = {name: np.array([[1., .5, -.1]]) for name in joint.CASES}
        ok, pairs, representatives = joint.compare_starts(fits, abc)
        self.assertTrue(ok); self.assertEqual(len(pairs), 12)
        self.assertEqual(representatives, {"double": "checkpoint-double", "float32": "checkpoint-float32"})
        abc["mixed-double"][0, 1] += 1e-7
        self.assertFalse(joint.compare_starts(fits, abc)[0])
        fits["mixed-double"]["joint_qualified"] = False
        self.assertFalse(joint.compare_starts(fits, abc)[0])

    def test_exact_repeat_only_excludes_resources_and_time(self):
        data = {"seconds": 1, "process_peak_rss_bytes": 123, "primary": {"beta": [1], "seconds": 2}, "lm_status": 5,
                "execution_complete": True, "joint_qualified": False}
        expected = {"primary": {"beta": [1]}, "lm_status": 5, "execution_complete": True, "joint_qualified": False}
        self.assertEqual(joint.scientific(data), expected)
        self.assertNotEqual(joint.scientific(dict(data, joint_qualified=True)), expected)


if __name__ == "__main__": unittest.main()
