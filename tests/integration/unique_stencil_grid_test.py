import copy
import tempfile
import unittest
from pathlib import Path

import unique_stencil_grid as grid
import mdpde_experiment as experiment


class UniqueStencilGridRunnerTest(unittest.TestCase):
    def fixture(self):
        context = {"state": [[6, .51, -.2]]}
        fit = {"schema_version": 2, "experiment": "unique-stencil-grid", "rows": 30, "columns": 2,
            "linear_solver": "sparse-qr", "sparse_row_reduction": "householder-qr-1024", "design_nonzeros": 60,
            "alpha": .1, "iteration_budget": 100, "refinement_budget": 100,
            "blocks": [{"scope": "global", "alpha": .1, "rows": 30, "lambda": 1}],
            "initial_beta": [6, -.2], "abc": [[6.1, .51, -.21]], "beta": [6.1, -.21],
            "qualified": False, "reason": "budget-exhausted", "branches": []}
        return fit, context

    def test_fixed_width_and_full_population(self):
        fit, context = self.fixture(); grid.validate_fit(fit, context, .1, 30)
        for key, value, message in (("rows", 29, "population"), ("linear_solver", "dense-qr", "sparse"),
                                    ("design_nonzeros", 0, "sparse"), ("blocks", [], "global block"),
                                    ("iteration_budget", 101, "budget"), ("refinement_budget", 4000, "budget")):
            changed = copy.deepcopy(fit); changed[key] = value
            with self.assertRaisesRegex(RuntimeError, message): grid.validate_fit(changed, context, .1, 30)
        changed = copy.deepcopy(fit); changed["abc"][0][1] = .5
        with self.assertRaisesRegex(RuntimeError, "Width changed"): grid.validate_fit(changed, context, .1, 30)

    def test_qualification_requires_both_phases_and_svd(self):
        fit, context = self.fixture(); fit["qualified"] = True
        with self.assertRaisesRegex(RuntimeError, "Unverified"): grid.validate_fit(fit, context, .1, 30)
        fit["qualified"] = False
        phase = {"iterations": 100, "stop": "budget-exhausted", "stationarity": 1e-5,
                 "beta": fit["beta"], "variances": [.01]}
        branch = {"primary": phase, "reference": phase, "trace": [], "reference_trace": [], "qualified": False}
        fit["branches"] = [dict(copy.deepcopy(branch), seed=seed) for seed in ("checkpoint", "constrained-ls")]
        fit["selected_seed"] = 0; fit["variances"] = [.01]
        grid.validate_fit(fit, context, .1, 30)
        fit["branches"][0]["primary"]["iterations"] = 101
        with self.assertRaisesRegex(RuntimeError, "budget exceeded"): grid.validate_fit(fit, context, .1, 30)

    def test_duplicate_voxel_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            (p/"voxels.csv").write_text("row,index,x,y,z,observed,multiplicity,reference_double,quantization_bound\n0,0,0,0,0,0,1,0,0\n1,0,0,0,0,0,1,0,0\n")
            (p/"slots.csv").write_text("")
            dataset = {"schema_version": 1, "experiment": "unique-stencil-grid", "atoms": [None]*168,
                "samples": [None]*33600, "row_count": grid.ROWS, "slot_count": 33600*64, "alphas": grid.ALPHAS,
                "grid_size": [10]*3, "generation_origin": [0]*3, "generation_spacing": [1]*3,
                "voxel_table_sha256": grid.fold.sha256_file(p/"voxels.csv"), "slot_table_sha256": grid.fold.sha256_file(p/"slots.csv")}
            with self.assertRaisesRegex(RuntimeError, "Duplicate"): grid.validate_dataset(p, dataset)

    def test_missing_cases_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            experiment.write(p/"fit-index.json", {"schema_version": 1, "states": []})
            experiment.write(p/"state-index.json", {"states": [], "selected_state_ids": []})
            experiment.write(p/"dataset.json", {})
            with self.assertRaisesRegex(RuntimeError, "State index"): grid.summarize(p)

    def test_unexecuted_comparison_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp); experiment.write(p/"fit-index.json", {"selected_state_ids": []})
            with self.assertRaisesRegex(RuntimeError, "not executed"): grid.compare(p, p, p/"compare.json", "missing")

    def test_compare_excludes_only_resources(self):
        a = {"seconds": 1, "peak_child_rss_bytes": 42, "beta": [1], "variance": .1, "iterations": 100}
        self.assertEqual(grid.scientific(a), {"beta": [1], "variance": .1, "iterations": 100})
        self.assertNotEqual(grid.scientific(a), grid.scientific(dict(a, iterations=101)))

    def test_weight_replay_allows_residual_rounding_but_rejects_changed_weights(self):
        prediction, residual = 2.8691596224482487, -0.00016634013013350213
        observed = prediction+residual
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp); (p/"residuals").mkdir()
            (p/"residuals/case-samples.csv").write_text(f"sample,prediction,residual\n0,{prediction},{residual}\n")
            expected = {"grid_rmse": abs(residual), "matched_sample_rmse": abs(residual)}
            fit = {"alpha": 1., "variances": [2.6476431666842024e-8]}
            for weight, valid in [(0.59302302449303623, True), (0.59303302449303623, False)]:
                (p/"residuals/case.csv").write_text(f"row,prediction,residual,weight\n0,{prediction},{residual},{weight}\n")
                if valid: grid.residual_statistics(p, "case", [observed], [{"response": observed}], expected, fit)
                else:
                    with self.assertRaisesRegex(RuntimeError, "weight replay"):
                        grid.residual_statistics(p, "case", [observed], [{"response": observed}], expected, fit)


if __name__ == "__main__": unittest.main()
