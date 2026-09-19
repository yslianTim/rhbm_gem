import copy
import tempfile
import unittest
from pathlib import Path
from joint_runtime_support import differences, scientific, unpack, read
from joint_component_runtime import CATALOG, compare
import joint_fixture_records as records
import joint_offline_support as runner
import numpy as np


class JointRuntimeSupportTest(unittest.TestCase):
    def test_missing_reference_cannot_pass_empty_comparison(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(RuntimeError, 'Missing comparison'):
                compare(Path(tmp)/'missing', Path(tmp)/'also-missing')

    def test_decisions_and_missing_states_are_exact(self):
        self.assertTrue(differences({'accepted': False}, {'accepted': True}))
        self.assertTrue(differences({'state': None}, {'state': []}))
        self.assertTrue(differences({'rank': 1}, {'rank': 2}))
        self.assertTrue(differences({'width': .5}, {'width': .6}))
        self.assertFalse(differences({'width': .5}, {'width': .5+1e-13}))
        self.assertTrue(differences({'width': .5}, {'width': float('nan')}))

    def test_only_costs_are_removed(self):
        self.assertEqual(scientific({'search_seconds': 2, 'process_peak_rss_bytes': 1,
                                     'state': {'objective': .2}, 'regular': False}),
                         {'state': {'objective': .2}, 'regular': False})

    def test_fixture_is_self_contained_and_restores_modified_member(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = unpack(CATALOG, 'baseline', tmp)
            self.assertEqual(len(read(root/'cases.json')), 8)
            expected = (root/'snapshot.json').read_bytes()
            (root/'snapshot.json').write_text('{}')
            unpack(CATALOG, 'baseline', tmp)
            self.assertEqual((root/'snapshot.json').read_bytes(), expected)


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

    def test_raw_component_replay_uses_parent_scale(self):
        data = {"ids": ["a"], "table": np.array([(0, 0, 0.), (1, 0, 1.)],
                dtype=[("row", int), ("atom", int), ("square", float)])}
        y = np.array([1., .5]); endpoint = {"beta": [1., -.2], "b": [.5]}
        local = runner.replay(data, y, endpoint, np.linalg.norm(y))
        parent = runner.replay(data, y, endpoint, 100.)
        np.testing.assert_array_equal(local["prediction"], parent["prediction"])
        self.assertAlmostEqual(parent["projected_kkt"], local["projected_kkt"]*np.linalg.norm(y)/100.)
        np.testing.assert_allclose(parent["b_gradient"], local["b_gradient"]*(np.linalg.norm(y)/100.)**2)

if __name__ == '__main__': unittest.main()
