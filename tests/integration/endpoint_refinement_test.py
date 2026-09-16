import unittest

import endpoint_refinement as experiment


class EndpointRunnerTest(unittest.TestCase):
    def test_budget_is_frozen_from_full_cost_with_twofold_margin(self):
        self.assertEqual(experiment.choose_budget(16),32)
        self.assertEqual(experiment.choose_budget(35),128)
        self.assertEqual(experiment.choose_budget(128),256)
        self.assertIsNone(experiment.choose_budget(129))

    def test_native_success_cannot_hide_rejection(self):
        record={"triggered":True,"effective_status":0,
                "refinement":{"accepted":False,"reason":"branch-mismatch","candidate_equation_evaluations":20}}
        with self.assertRaisesRegex(RuntimeError,"inherited SUCCESS"):
            experiment.validate_refinements([record],128)
        record["effective_status"]=4
        self.assertEqual(experiment.validate_refinements([record],128)["rejections"],{"branch-mismatch":1})

    def test_acceptance_requires_all_independent_evidence(self):
        r={"candidate_equation_evaluations":20,"accepted":True,"branch":{"pass":True},
           "root":{"equation_pass":True},"reference":{"reference_pass":True}}
        record={"triggered":True,"effective_status":0,"refinement":r}
        self.assertEqual(experiment.validate_refinements([record],128)["accepted"],1)
        r["reference"]["reference_pass"]=False
        with self.assertRaisesRegex(RuntimeError,"incomplete"):
            experiment.validate_refinements([record],128)
        r["reference"]["reference_pass"]=True
        r["candidate_equation_evaluations"]=129
        with self.assertRaisesRegex(RuntimeError,"budget"):
            experiment.validate_refinements([record],128)


if __name__ == "__main__":
    unittest.main()
