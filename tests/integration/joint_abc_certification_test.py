import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import joint_abc_certification as cert
import simulation_contract as contract


class CertificationTest(unittest.TestCase):
    def manifest(self, directory, model='single_gaus'):
        paths = {k: directory/k for k in ('model', 'map', 'manifest')}
        for k in ('model', 'map'): paths[k].write_text(k)
        widths = [.5, .45, .4]
        value = {'schema_version': 2, 'source': {'model_sha256': cert.fixed.fold.sha256_file(paths['model'])},
                 'output': {'map_sha256': cert.fixed.fold.sha256_file(paths['map'])},
                 'generator': dict.fromkeys(('source_sha256', 'configuration_sha256', 'build_sha256'), 'a'*64),
                 'settings': {'potential_model': model, 'blurring_width': .5, 'cutoff_distance': 2.5},
                 'kernel': {'version': model+'-v1', 'width_policy': contract.POLICY if model == 'single_gaus' else {'id': 'model-specific-v1'},
                            'near_zero_distance': None if model == 'single_gaus_user' else 1e-5,
                            'charge_term_cutoff': {'single_gaus': 2.5, 'five_gaus_charge': 3., 'single_gaus_user': None}[model],
                            'minimum_charge_width': 1e-5 if model == 'five_gaus_charge' else None},
                 'support': dict(contract.SUPPORT, outer_cutoff=2.5),
                 'execution': {'accumulation': 'z_planes_in_preparation_order'}, 'atom_count': 3, 'atoms': []}
        for k in range(3):
            value['atoms'].append({'serial_id': k+1, 'sequence_id': k+1, 'chain_id': 'A', 'component_id': 'ALA',
                'atom_id': 'CNO'[k], 'alternate_indicator': '.', 'structure': 0, 'position': [k, 0., 0.],
                'element': 6+k, 'preparation_index': k,
                'effective_gaussian_width': widths[k] if model == 'single_gaus' else None,
                'effective_charge_width': widths[k] if model == 'single_gaus' else .5 if model == 'five_gaus_charge' else None})
        cert.write(paths['manifest'], value)
        return paths, value

    def test_three_models_have_explicit_applicability(self):
        with tempfile.TemporaryDirectory() as d:
            for model in ('single_gaus', 'five_gaus_charge', 'single_gaus_user'):
                paths, _ = self.manifest(Path(d), model)
                _, result = contract.resolve(paths['manifest'], paths)
                self.assertEqual(result['widths'], [.5, .45, .4] if model == 'single_gaus' else [None]*3)

    def test_wrong_width_policy_hash_identity_and_version_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            paths, value = self.manifest(Path(d))
            mutations = [lambda v: v['atoms'][1].update(effective_charge_width=.5),
                lambda v: v['kernel']['width_policy'].update(nitrogen=1),
                lambda v: v['source'].update(model_sha256='b'*64),
                lambda v: v.update(schema_version=3), lambda v: v.update(schema_version=1),
                lambda v: v['atoms'][1].update(preparation_index=0),
                lambda v: v['atoms'][1].update(**{k: v['atoms'][0][k] for k in ('serial_id', 'chain_id', 'sequence_id', 'component_id', 'atom_id', 'alternate_indicator')})]
            for mutate in mutations:
                bad = copy.deepcopy(value); mutate(bad); cert.write(paths['manifest'], bad)
                with self.assertRaises(RuntimeError): contract.resolve(paths['manifest'], paths)

    def test_failure_matrix_keeps_multiple_failures_and_unavailable(self):
        fit = {'execution_complete': True, 'qualification_checks': {'inner': True, 'b_gradient': True, 'local_correction': False, 'identified': False}}
        audit = {'derivative_status': 'unavailable', 'derivative_verified': False, 'reason': 'no width spectrum'}
        result = cert.certificate(fit, audit, False)
        self.assertEqual(result['failures'], ['local_correction', 'identified', 'endpoint_replay'])
        self.assertEqual(result['unavailable'], ['derivative', 'endpoint_trust', 'regular_active_face'])
        self.assertFalse(result['regular_qualified'])
        self.assertFalse(result['truth_used'])

    def ladder(self):
        ladder = {'direction': 0, 'samples': [{'h': .01*2**-k} for k in range(17)],
                  'candidates': [{'index': 0, 'estimated_relative_error': 1e-9}, {'index': 1, 'estimated_relative_error': 1e-8}],
                  'selected': 0, 'first_relative_error': 1e-8, 'second_relative_error': 1e-8, 'passed': True}
        audit = {'schema_version': 2, 'case': 'x', 'dataset': 'y', 'legacy_joint_qualified': False,
                 'ladders': [copy.deepcopy(ladder) for _ in range(3)], 'derivative_verified': False}
        return audit, {'case': 'x', 'dataset': 'y', 'joint_qualified': False}

    def test_step_selection_cannot_choose_analytic_match_or_truth(self):
        audit, fit = self.ladder(); cert.validate_audit(audit, fit)
        audit['truth_error'] = 0; audit['ladders'][0]['selected'] = 1
        with self.assertRaisesRegex(RuntimeError, 'selection'): cert.validate_audit(audit, fit)
        audit['ladders'][0]['selected'] = 0; audit['derivative_verified'] = True
        with self.assertRaisesRegex(RuntimeError, 'False derivative'): cert.validate_audit(audit, fit)

    def test_false_precision_agreement_is_rejected(self):
        audit, fit = self.ladder()
        audit['precision'] = {'valid50': True, 'valid100': True, 'agreement_passed': True, 'derivative_passed': False,
            'maximum_scaled_precision_difference': '1e-15'}
        with self.assertRaisesRegex(RuntimeError, 'precision agreement'): cert.validate_audit(audit, fit)

    def test_legacy_derivative_success_cannot_skip_precision_target_ladder(self):
        audit, fit = self.ladder()
        fit['derivative_verified'] = True; audit['derivative_verified'] = True
        audit['ladders'] = []; audit['precision'] = {'agreement_passed': True, 'derivative_passed': True}
        with self.assertRaisesRegex(RuntimeError, 'full derivative ladder'): cert.validate_audit(audit, fit)

    def test_modified_snapshot_rejected_before_membership_reconstruction(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); (root/'voxels.csv').write_text('changed')
            cert.write(root/'snapshot.json', {'support_policy': 'sphere-fma-v1', 'rows': 1, 'atoms': 1, 'voxels_sha256': 'a'*64})
            with self.assertRaisesRegex(RuntimeError, 'hash mismatch'): cert.snapshot(root, np.zeros((1,3)), [{}])

    def test_snapshot_replay_keeps_zero_negative_observations(self):
        table = np.array([(0,0,0.), (1,0,1.), (2,0,6.25)], dtype=[('row',int), ('atom',int), ('square',float)])
        y = np.array([0., -1., 2.])
        result = cert.snapshot_replay(table, 1, y, {'b': [.5], 'beta': [0., -.2]})
        np.testing.assert_array_equal(result['residual'], result['prediction']-y)
        self.assertEqual(len(result['residual']), 3)

    def test_failed_initialization_has_no_accepted_trials(self):
        self.assertEqual(cert.validate_accepted({"variant": "guarded", "execution_complete": False}, None, 0, None), [])

    def test_untrusted_accepted_state_is_never_success(self):
        with self.assertRaisesRegex(RuntimeError, 'untrusted'):
            cert.validate_accepted({'variant': 'guarded', 'trials': [{'accepted': True, 'trust': {'passed': False}}]}, None, 1, None)

    def test_missing_matrix_and_existing_directory_fail(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(RuntimeError, 'fresh output'): cert.run(type('Args', (), {'output': Path(d)})())
            with patch.object(cert, 'paths_and_manifest', return_value=({}, {}, {}, {})), patch.object(cert, 'read', return_value={'complete': False}):
                with self.assertRaisesRegex(RuntimeError, 'Incomplete'): cert.summarize(Path(d))


if __name__ == '__main__': unittest.main()
