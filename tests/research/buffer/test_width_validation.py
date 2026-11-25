import unittest

from confirm_truth import shared_panels, study_case
from truth import fitted_points, stratify_scenario, width_bin


def fixture():
    return {"spacing":.8,"upper":3.,"dimensions":[2,2,2],"origin":[0.,0.,0.],
            "positions":[[0.,0.,0.]]*4,"truth":[[150.,.24,-3.],[165.,.72,-2.6],[250.,.79,-3.6],[170.,.64,-1.8]],
            "noise":[[i,0.,1.] for i in range(8)],"noise_model":{"rms":1.}}


class WidthValidationTest(unittest.TestCase):
    def test_sweep_only_changes_requested_truth(self):
        base=fixture()
        case=study_case(base,"2-A",{1:{1:2.}})
        self.assertEqual(case["truth"][0],[150.,2.,-3.])
        self.assertEqual(base["truth"][0],[150.,.24,-3.])
        self.assertEqual(case["truth"][1:],base["truth"][1:])
        for key in ("positions","spacing","origin","dimensions","noise","upper"):
            self.assertEqual(case[key],base[key])

    def test_feasible_range_uses_b1_not_the_reporting_boundary(self):
        self.assertEqual(study_case(fixture(),"broad",{1:{1:2.85}})["truth"][0][1],2.85)
        for width in (.1722,3.1):
            with self.assertRaisesRegex(ValueError,"B1 interval"):
                study_case(fixture(),"out of model",{1:{1:width}})

    def test_report_bins_include_exact_boundary(self):
        self.assertEqual([width_bin(t) for t in (.5,.75,1.,1.25,1.75,2.,2.0001)],
                         ["≤0.5","0.5–1.0","0.5–1.0","1.0–1.5","1.5–2.0","1.5–2.0",">2.0"])

    def test_shared_noise_is_reproducible_and_validation_is_independent(self):
        inputs={"nominal":fixture(),"coarse":fixture()}
        first=shared_panels(inputs,1,2,6)
        self.assertEqual(first,shared_panels(inputs,1,2,6))
        self.assertEqual(first["nominal"]["noise"],first["coarse"]["noise"])
        self.assertNotEqual(first["nominal"]["noise"],first["nominal"]["validation_noise"])
        self.assertEqual([d[0] for d in first["nominal"]["validation_noise"]],list(range(8)))

    def test_mixed_width_classifies_truth_and_preserves_component_failure(self):
        true=[[8.,2.,-1.],[10.,.7,0.],[12.,2.85,-1.],[4.,.6,0.]]
        fits=[{"replicate":0,"strategy":"G0","qualified":True,"points":[[8.,2.7,-1.],true[1]]},
              {"replicate":0,"strategy":"single","qualified":True,"points":[[8.,2.8,-1.],true[1],true[2]]},
              {"replicate":1,"strategy":"G0","qualified":True,"points":true[:2]},
              {"replicate":1,"strategy":"single","qualified":False,"points":true[:3]}]
        result=stratify_scenario({"fits":fits,"truth":true})
        self.assertEqual(result["population"],"B")
        self.assertTrue(result["atoms"]["1"]["in_interest"])
        self.assertAlmostEqual(result["atoms"]["1"]["paired"]["single"]["rmse"][1],.8)
        self.assertEqual(result["transitions"]["n10"],[1])
        self.assertEqual(result["atoms"]["3"]["coverage"],{"G0":0,"single":1})
        self.assertEqual(result["bins"]["1.0–1.5"],0)
        self.assertEqual(fitted_points(fits[-1]),{})


if __name__ == "__main__":
    unittest.main()
