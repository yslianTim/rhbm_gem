import copy
from pathlib import Path
import tempfile
import unittest

import numpy as np
import joint_abc_component_records as records
import joint_abc_components as runner


class ComponentRecordsTest(unittest.TestCase):
    def test_census_reconstructs_structure_and_rejects_lost_rows(self):
        data = {"ids": ["a", "b", "c", "d"], "hash": "snapshot", "y64": np.ones(4),
                "table": np.array([(0, 0, 0.), (1, 0, 6.25), (1, 1, 0.), (2, 2, 1.)],
                                  dtype=[("row", int), ("atom", int), ("square", float)])}
        census = {"component_count": 3, "snapshot_sha256": "snapshot",
                  "components": [{"id": "a", "atoms": [0, 1], "rows": [0, 1], "memberships": 3},
                                 {"id": "c", "atoms": [2], "rows": [2], "memberships": 1},
                                 {"id": "d", "atoms": [3], "rows": [], "memberships": 0}],
                  "atom_component": [0, 0, 1, 2], "row_component": [0, 0, 1, -1],
                  "constant_rows": [3], "unobserved_atoms": [3]}
        for c in census["components"]:
            c["parent_atom_to_local"] = [c["atoms"].index(a) if a in c["atoms"] else -1 for a in range(4)]
            c["parent_row_to_local"] = [c["rows"].index(r) if r in c["rows"] else -1 for r in range(4)]
        self.assertTrue(records.verify_census(data, census))
        bad = copy.deepcopy(census); bad["components"][0]["parent_row_to_local"][1] = 0
        with self.assertRaisesRegex(RuntimeError, "local mapping"):
            records.verify_census(data, bad)
        bad = copy.deepcopy(census); bad["components"][0]["rows"] = [0]
        with self.assertRaisesRegex(RuntimeError, "membership"):
            records.verify_census(data, bad)

    def test_rehashed_recipe_cannot_hide_changed_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); source = root/"source"; source.mkdir()
            (source/"voxels.csv").write_text("changed observations")
            runner.write(source/"snapshot.json", {"voxels_sha256": "0"*64})
            runner.write(root/"composition.json", {
                "kind": "disjoint-frozen-snapshot-composition", "sources": [
                    {"path": "source", "identity_prefix": "x:", "snapshot_sha256": records.digest(source/"snapshot.json")} ]})
            with self.assertRaisesRegex(RuntimeError, "hash mismatch"):
                records.load(root)

    def test_mixed_start_uses_recorded_vector(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            runner.write(root/"snapshot.json", {"identity": "unchanged"})
            runner.write(root/"guarded/fits/mixed-double.json", {
                "observation_snapshot_sha256": records.digest(root/"snapshot.json"), "initial_b": [.12, .83]})
            data = {"sources": [root, root], "ids": ["remapped3", "remapped2", "remapped1", "remapped0"]}
            np.testing.assert_array_equal(records.frozen_start(data, "mixed-double"), [.12, .83, .12, .83])

    def test_regular_pair_cannot_pass_without_both_certificates_and_states(self):
        certificates = [{"regular_qualified": True}, {"regular_qualified": False}]
        row = runner.endpoint_parity({}, np.ones(1), {}, {}, certificates, True)
        self.assertFalse(row["passed"])
        certificates[1]["regular_qualified"] = True
        row = runner.endpoint_parity({}, np.ones(1), {}, {}, certificates, False)
        self.assertFalse(row["passed"])

    def test_final_scope_rejects_initial_direction_and_renormalized_child(self):
        fit = {"width_spectrum": {"weak_directions": [[0., 1.]]}}
        context = {"audit": {"directions": [[.5, .5], [.5, -.5], [1., 0.]]}}
        with self.assertRaisesRegex(RuntimeError, "endpoint weakest"):
            runner.validate_endpoint_scope(fit, context)
        context["audit"]["directions"][2] = [0., 1.]
        self.assertTrue(runner.validate_endpoint_scope(fit, context))
        child = {"audit": {"directions": [[.5], [-.5], [1.]]}}
        self.assertTrue(runner.validate_endpoint_scope({}, child, context, [1]))
        child["audit"]["directions"][0] = [1.]
        with self.assertRaisesRegex(RuntimeError, "renormalized"):
            runner.validate_endpoint_scope({}, child, context, [1])

    def test_raw_component_replay_uses_parent_scale(self):
        data = {"ids": ["a"], "table": np.array([(0, 0, 0.), (1, 0, 1.)],
                dtype=[("row", int), ("atom", int), ("square", float)])}
        y = np.array([1., .5]); endpoint = {"beta": [1., -.2], "b": [.5]}
        local = runner.replay(data, y, endpoint, np.linalg.norm(y))
        parent = runner.replay(data, y, endpoint, 100.)
        np.testing.assert_array_equal(local["prediction"], parent["prediction"])
        self.assertAlmostEqual(parent["projected_kkt"], local["projected_kkt"]*np.linalg.norm(y)/100.)
        np.testing.assert_allclose(parent["b_gradient"], local["b_gradient"]*(np.linalg.norm(y)/100.)**2)


if __name__ == "__main__":
    unittest.main()
