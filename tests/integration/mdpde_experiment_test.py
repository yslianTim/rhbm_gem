import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import mdpde_experiment as experiment


class ExperimentRunnerTest(unittest.TestCase):
    def test_empty_solver_results_are_not_a_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaisesRegex(RuntimeError,"No solver results"):
                experiment.summarize(SimpleNamespace(output=root/"out",solver=root,forward=root,plots=False))

    def test_statistics_preserve_bias_and_absolute_tail(self):
        result = experiment.stats([-2.0,1.0,1.0,float("nan")])
        self.assertEqual(result["n"],3)
        self.assertEqual(result["bias"],0.0)
        self.assertEqual(result["p99"],2.0)
        self.assertAlmostEqual(result["rmse"],2**.5)


if __name__ == "__main__":
    unittest.main()
