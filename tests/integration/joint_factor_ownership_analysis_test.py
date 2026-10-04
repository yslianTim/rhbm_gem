"""Tests for search factor ownership evidence and trajectory parity."""
import unittest

import joint_factor_ownership_analysis as analysis


def report(owner, *, eta=None, status='completed'):
    return {
        'case': 'cube-512',
        'metadata': {'solver_policy': {'operator_factor_ownership': owner}},
        'execution': {'status': status, 'runs': [{'process_wall_seconds': 12.0,
            'stages': [{'sampled_tree_peak_rss_bytes': 1024}]}]},
        'result': {'details': {
            'operator_factor_ownership': owner,
            'search_seconds': 10.0,
            'stop_reason': 'operator-gradient-stop',
            'profile_evaluations': 5,
            'accepted_updates': 4,
            'accepted_objective': 0.25,
            'accepted_gradient_inf_norm': 1e-12,
            'returned_search_state': eta if eta is not None else [0.1, 0.2],
            'search_work': {
                'factor_residency': {'maximum_concurrent_factor_count': 2,
                    'maximum_concurrent_owned_bytes_estimate': 4096,
                    'maximum_concurrent_factor_stage': 'trial-evaluation'},
                'pcg_iteration_counts': [4, 5],
                'operator_rank_status': 'full-rank',
                'operator_rank_certificate': 'local-support',
                'operator_rank_reason': 'rank-verified-full',
                'accepted_factor_reuse_attempts': 5 if owner.startswith('reuse') else 0,
                'accepted_factor_reuse_accepted': 5 if owner.startswith('reuse') else 0,
                'accepted_factor_reuse_fallbacks': 0,
            },
            'spqr_factorization': {'symbolic': {'calls': 3, 'seconds': .2},
                                   'numeric': {'calls': 6, 'seconds': .4}},
        }},
    }


class FactorOwnershipAnalysisTest(unittest.TestCase):
    def test_matching_search_trajectory_is_reported(self):
        current = analysis._run(report('dedicated-fixed'))
        candidate = analysis._run(report('reuse-accepted-handoff'))
        parity = analysis.trajectory_parity(current, candidate)
        self.assertTrue(parity['passed'], parity)
        self.assertEqual(candidate['maximum_concurrent_factor_count'], 2)
        self.assertEqual(candidate['peak_rss_bytes'], 1024)
        self.assertEqual(candidate['reuse_accepted'], 5)
        self.assertEqual(candidate['derived_symbolic_reuses'], 3)

    def test_changed_trace_fails_parity(self):
        current = analysis._run(report('dedicated-fixed'))
        candidate = analysis._run(report('reuse-accepted-copy-on-write', eta=[0.1, 0.3]))
        parity = analysis.trajectory_parity(current, candidate)
        self.assertFalse(parity['passed'])
        self.assertFalse(parity['checks']['returned_search_state'])

    def test_incomplete_or_failed_campaign_row_keeps_status(self):
        row = analysis._run(report('reuse-accepted-handoff', status='timeout'))
        self.assertEqual(row['status'], 'timeout')
        self.assertEqual(row['maximum_concurrent_factor_count'], 2)


if __name__ == '__main__':
    unittest.main()
