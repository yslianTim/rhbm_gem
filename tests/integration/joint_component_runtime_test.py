import copy
import tempfile
import unittest
from pathlib import Path
from joint_runtime_support import differences, scientific, unpack, read
from joint_component_runtime import CATALOG, compare, runtime_expected, backend_differences
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

    def test_runtime_projection_preserves_search_and_separates_derivative_audit(self):
        expected = {'usable_state': True, 'search_success': False,
                    'qualification_checks': dict(inner=True, b_gradient=True, local_correction=True,
                                                 identified=True, derivative=False),
                    'qualification_failure': 'derivative-unverified'}
        original = copy.deepcopy(expected)
        runtime = runtime_expected(expected)
        self.assertEqual(expected, original)
        self.assertEqual(runtime['runtime_convergence'], 'passed')
        self.assertFalse(runtime['search_success'])
        self.assertNotIn('derivative', runtime['runtime_checks'])
        expected['qualification_checks']['local_correction'] = False
        self.assertEqual(runtime_expected(expected)['runtime_convergence'], 'failed')
        expected['usable_state'] = False
        self.assertEqual(runtime_expected(expected)['runtime_convergence'], 'unavailable')
        expected = copy.deepcopy(original)
        del expected['qualification_checks']['local_correction']
        self.assertEqual(runtime_expected(expected)['runtime_convergence'], 'unavailable')
        expected['qualification_checks'] = {}
        self.assertEqual(runtime_expected(expected)['runtime_convergence'], 'unavailable')

    def test_backend_contract_allows_trace_changes_but_not_lost_evidence(self):
        record = dict(usable_state=True, runtime_convergence='passed',
                      runtime_checks=dict(inner=True, b_gradient=True, local_correction=True, identified=True),
                      search_success=True, stop_reason='native-lm-stop', accepted_updates=3, profile_evaluations=4,
                      trials=[dict(accepted=True, valid=True, trust_passed=True)],
                      last_trusted_state=dict(valid=True, feasible=True, kkt_passed=True, free_rank=2,
                                              objective=1e-15, beta=[2., .1], eta=[-.7], b=[.5]))
        actual = copy.deepcopy(record)
        actual.update(accepted_updates=4, profile_evaluations=5, stop_reason='unrepresentable-step', search_success=False)
        self.assertFalse(backend_differences(record, actual, 1.))
        actual['last_trusted_state']['b'][0] += 1e-6
        self.assertTrue(backend_differences(record, actual, 1.))
        actual = copy.deepcopy(record); del actual['runtime_checks']['inner']
        self.assertTrue(backend_differences(record, actual, 1.))
        actual = copy.deepcopy(record); actual['trials'][0]['trust_passed'] = False
        self.assertTrue(backend_differences(record, actual, 1.))
        actual = copy.deepcopy(record); actual['profile_evaluations'] = 201
        self.assertTrue(backend_differences(record, actual, 1.))
        actual = copy.deepcopy(record); actual['last_trusted_state']['free_rank'] = 1
        self.assertTrue(backend_differences(record, actual, 1.))

    def test_backend_nonconverged_endpoint_must_not_worsen_objective(self):
        expected = dict(usable_state=True, runtime_convergence='failed', runtime_checks={'identified': False},
                        search_success=False, stop_reason='profile-budget', accepted_updates=1, profile_evaluations=200,
                        last_trusted_state=dict(valid=True, feasible=True, kkt_passed=True, free_rank=2, objective=.1))
        actual = copy.deepcopy(expected)
        actual['last_trusted_state']['objective'] = .09
        self.assertFalse(backend_differences(expected, actual, 10.))
        expected['last_trusted_state']['active_atoms'] = []
        actual['last_trusted_state'].update(active_atoms=[0], free_rank=1)
        self.assertFalse(backend_differences(expected, actual, 10.))
        actual['last_trusted_state']['objective'] = .10001
        self.assertTrue(backend_differences(expected, actual, 10.))
        actual['last_trusted_state'] = None
        self.assertTrue(backend_differences(expected, actual, 10.))

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
