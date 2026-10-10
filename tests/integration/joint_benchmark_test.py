"""Contract tests for the current FixedNeighbor Joint benchmark surface."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import joint_benchmark as benchmark


class JointBenchmarkContractTest(unittest.TestCase):
    def test_surface_has_only_current_profiles_and_options(self):
        parser = benchmark.build_parser()
        profile_action = next(action for action in parser._actions if action.dest == 'profile')
        self.assertEqual(set(profile_action.choices),
                         {'search', 'solve', 'workflow', 'postprocess', 'command'})
        options = {option for action in parser._actions for option in action.option_strings}
        self.assertEqual(options, {
            '-h', '--help', '--profile', '--case', '--build-dir', '--output', '--repeat',
            '--warmup', '--timeout', '--rss-limit', '--cli', '--model', '--map',
        })

    def test_default_solver_metadata_locks_fixed_neighbor_production_policy(self):
        args = benchmark.build_parser().parse_args([
            '--profile', 'solve', '--case', 'chain-8', '--build-dir', 'build/debug-tests',
            '--output', 'result.json'])
        self.assertEqual(benchmark.solver_policy_metadata(), {
            'search_method': 'FixedNeighbor',
            'sparse_backend': 'EIGEN',
            'fixed_neighbor_core_atoms': 12,
            'fixed_neighbor_block_order': 'Forward',
            'fixed_neighbor_maximum_sweeps': 30,
            'fixed_neighbor_local_search': 'LegacyCompact',
            'fixed_neighbor_update_policy': 'OneAccepted',
        })

    def test_profile_commands_use_fixed_neighbor_driver(self):
        parser = benchmark.build_parser()
        with tempfile.TemporaryDirectory() as temporary:
            build = Path(temporary)
            (build / 'bin').mkdir()
            driver = build / 'bin/joint_fixed_neighbor_experiment'
            driver.touch()
            for profile, mode in (('search', '--scaling-only'), ('solve', '--case')):
                args = parser.parse_args([
                    '--profile', profile, '--case', 'chain-8', '--build-dir', str(build),
                    '--output', 'result.json'])
                command = benchmark.command_for_profile(args, Path('result.json'), build)
                self.assertEqual(command[0], str(driver))
                self.assertEqual(command[1], mode)
                self.assertEqual(command[-2:], ['chain', '8'])
                self.assertEqual(len(command), 5)

    def test_custom_core_size_option_is_not_exposed(self):
        with self.assertRaises(SystemExit):
            benchmark.build_parser().parse_args([
                '--profile', 'search', '--case', 'chain-8', '--build-dir', 'build',
                '--output', 'result.json', '--fixed-core-atoms', '16'])

    def test_normalization_preserves_search_and_endpoint_contracts(self):
        search = benchmark.normalize_result('search', {'fixed_neighbor': {
            'search_converged': True, 'search_reason': 'block-stationary', 'sweeps': 2,
            'block_solves': 8, 'profile_evaluations': 16, 'accepted_blocks': 4,
            'objective': .2, 'final_global_ac_kkt': 1e-12,
            'final_raw_width_gradient_inf_norm': 1e-13,
        }})
        self.assertIsNone(search['qualified'])
        self.assertEqual(search['scientific_status'], 'search-only')
        self.assertTrue(search['details']['search_completed'])
        solve = benchmark.normalize_result('solve', {'fixed_neighbor': {
            'search_converged': True, 'search_reason': 'block-stationary',
            'endpoint_certified': True, 'runtime_convergence': 'Passed',
            'objective': .1, 'sweep_telemetry': [],
        }})
        self.assertTrue(solve['qualified'])
        self.assertNotIn('global_legacy_compact', solve['details'])


def smoke(build):
    script = Path(__file__).with_name('joint_benchmark.py')
    with tempfile.TemporaryDirectory(prefix='joint-benchmark-smoke-') as temporary:
        root = Path(temporary)
        for profile in ('search', 'solve'):
            output = root / f'{profile}.json'
            subprocess.run([
                sys.executable, str(script), '--profile', profile, '--case', 'chain-8',
                '--build-dir', str(build), '--output', str(output), '--timeout', '60',
            ], check=True, capture_output=True, text=True, timeout=120)
            report = json.loads(output.read_text())
            if report['execution']['status'] != 'completed':
                raise AssertionError(f'{profile} did not complete: {report["execution"]}')
            policy = report['metadata']['solver_policy']
            if (policy['search_method'] != 'FixedNeighbor' or
                    policy['fixed_neighbor_core_atoms'] != 12 or
                    policy['fixed_neighbor_maximum_sweeps'] != 30 or
                    policy['fixed_neighbor_block_order'] != 'Forward' or
                    policy['fixed_neighbor_local_search'] != 'LegacyCompact'):
                raise AssertionError(f'FixedNeighbor policy metadata is incomplete: {policy}')
            if profile == 'search' and report['result']['qualified'] is not None:
                raise AssertionError('search-only benchmark must not claim endpoint qualification')
    return 0


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--smoke':
        raise SystemExit(smoke(Path(sys.argv[2]).resolve()))
    unittest.main()
