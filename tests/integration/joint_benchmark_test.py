"""Contract checks and deterministic smoke runs for the unified Joint benchmark."""
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import joint_benchmark as benchmark


class JointBenchmarkContract(unittest.TestCase):
    def test_profiles_have_distinct_measurement_scopes(self):
        expected = {'prepare', 'fixed', 'solve', 'rank', 'workflow', 'postprocess', 'command'}
        self.assertEqual(set(benchmark.SCOPES), expected)
        self.assertEqual(len(set(benchmark.SCOPES.values())), len(expected))

    def test_process_status_keeps_limits_and_incomplete_work_explicit(self):
        for source, result in (('completed', 'completed'), ('process-failure', 'process_error'),
                               ('time-limit', 'timeout'), ('rss-limit', 'rss_limit'),
                               ('not-run-budget', 'not_run')):
            self.assertEqual(benchmark.process_status({'status': source}), result)

    def test_profile_timers_use_their_named_work(self):
        self.assertEqual(benchmark.raw_elapsed('prepare',
            {'construction_seconds': .2, 'basis_seconds': .3}, 99), .5)
        self.assertEqual(benchmark.raw_elapsed('fixed',
            {'steps': [{'kind': 'schwarz', 'fixed_step_wall_seconds': .4}]}, 99), .4)
        self.assertAlmostEqual(benchmark.raw_elapsed('solve',
            {'search_seconds': .7, 'assessment_seconds': .2}, 99), .9)
        self.assertAlmostEqual(benchmark.raw_elapsed('postprocess',
            {'phases': {'uncertainty_seconds': .1, 'save_seconds': .2}}, 99), .3)

    def test_process_completion_does_not_imply_numerical_qualification(self):
        result = benchmark.normalize_result('solve', {
            'returned_assessment': {'runtime_convergence': 'failed'},
        })
        self.assertFalse(result['qualified'])
        self.assertEqual(result['scientific_status'], 'failed')

    def test_operator_and_schwarz_defaults_are_stable(self):
        parser = benchmark.build_parser()
        args = parser.parse_args(['--profile', 'solve', '--case', 'chain-8',
                                  '--build-dir', 'build/debug', '--output', 'result.json'])
        benchmark.validate_args(parser, args)
        self.assertEqual((args.operator_rank, args.schwarz_core_atoms, args.schwarz_overlap_hops,
                          args.schwarz_max_block_atoms, args.schwarz_storage_mib,
                          args.schwarz_scratch_mib), ('auto', 128, 1, 512, 512, 256))
        policy = benchmark.solver_policy_metadata(args, 'EIGEN')
        self.assertEqual(policy['search_method'], 'OperatorPcg')
        self.assertEqual(policy['resolved_rank_backend'], 'Dense')
        self.assertEqual(policy['schwarz_overlap_hops'], 1)

    def test_operator_policy_options_reach_only_sparse_driver_routes(self):
        parser = benchmark.build_parser()
        args = parser.parse_args(['--profile', 'solve', '--case', 'chain-8',
                                  '--build-dir', 'build/debug', '--output', 'result.json',
                                  '--operator-rank', 'spqr-bounds', '--schwarz-core-atoms', '2',
                                  '--schwarz-overlap-hops', '0', '--schwarz-max-block-atoms', '8',
                                  '--schwarz-storage-mib', '64', '--schwarz-scratch-mib', '32'])
        benchmark.validate_args(parser, args)
        with tempfile.TemporaryDirectory() as temporary:
            build = Path(temporary)
            driver = build / 'bin/joint_sparse_benchmark'
            driver.parent.mkdir()
            driver.touch()
            command = benchmark.command_for_profile(args, {'kind': 'synthetic', 'topology': 'chain', 'atoms': 8},
                                                    Path('result.json'), build)
            for option, value in (('--operator-rank', 'spqr-bounds'), ('--schwarz-core-atoms', '2'),
                                  ('--schwarz-overlap-hops', '0'), ('--schwarz-max-block-atoms', '8'),
                                  ('--schwarz-storage-mib', '64'), ('--schwarz-scratch-mib', '32')):
                self.assertEqual(command[command.index(option) + 1], value)
            self.assertEqual(benchmark.solver_policy_metadata(args, 'SPQR')['resolved_rank_backend'], 'SpqrBounds')

    def test_non_sparse_profiles_do_not_receive_operator_policy_options(self):
        parser = benchmark.build_parser()
        args = parser.parse_args(['--profile', 'workflow', '--case', 'full',
                                  '--build-dir', 'build/debug', '--output', 'result.json'])
        with tempfile.TemporaryDirectory() as temporary:
            build = Path(temporary)
            driver = build / 'bin/joint_postprocessing_benchmark'
            driver.parent.mkdir()
            driver.touch()
            command = benchmark.command_for_profile(args, {}, Path('result.sqlite'), build)
        self.assertEqual(command, [str(driver), 'full', 'workflow', 'result.sqlite'])

    def test_invalid_schwarz_policy_values_are_rejected(self):
        base = ['--profile', 'solve', '--case', 'chain-8', '--build-dir', 'build/debug',
                '--output', 'result.json']
        for option, value in (('--schwarz-core-atoms', '0'), ('--schwarz-overlap-hops', '-1'),
                              ('--schwarz-max-block-atoms', '127'), ('--schwarz-storage-mib', '0'),
                              ('--schwarz-scratch-mib', '-1')):
            parser = benchmark.build_parser()
            args = parser.parse_args(base + [option, value])
            with self.assertRaises(SystemExit):
                benchmark.validate_args(parser, args)

    def test_missing_telemetry_does_not_change_result_normalization(self):
        raw = {'search': {'stop_reason': 'converged', 'profile_evaluations': 3, 'accepted_updates': 2},
               'returned_assessment': {'runtime_convergence': 'passed', 'primary': {'valid': True}}}
        details = benchmark.normalize_result('solve', raw)['details']
        self.assertEqual(set(details), {'search_completed', 'stop_reason', 'profile_evaluations',
                                        'accepted_updates', 'endpoint_valid'})

def smoke(build):
    script = Path(__file__).with_name('joint_benchmark.py')
    cli = build / 'bin/RHBM-GEM'
    source_model = Path(__file__).resolve().parents[1] / 'fixtures/test_model.cif'
    with tempfile.TemporaryDirectory(prefix='joint-benchmark-smoke-') as temporary:
        root = Path(temporary)
        model = root / 'input.cif'
        shutil.copyfile(source_model, model)
        model.write_text(model.read_text().rsplit('#', 1)[0] +
            'ATOM 2 C CB . ALA A 1 1.2 0.0 0.0 1.0 0.0 1\n'
            'ATOM 3 C CB . ALA A 2 12.0 0.0 0.0 1.0 0.0 1\n#\n')
        subprocess.run([str(cli), 'map_simulation', '-a', str(model), '-o', str(root),
                        '--potential-model', 'single', '--blurring-width', '.5', '-g', '.3', '-v', '0'],
                       check=True, capture_output=True, text=True, timeout=60)
        map_path = next(root.glob('*.map'))
        profiles = ('prepare', 'fixed', 'solve', 'rank', 'workflow', 'postprocess')
        for profile in profiles:
            output = root / f'{profile}.json'
            command = [sys.executable, str(script), '--profile', profile, '--case',
                            'full' if profile in ('workflow', 'postprocess') else 'chain-8',
                            '--build-dir', str(build), '--output', str(output), '--timeout', '60']
            if profile == 'fixed':
                command.extend(['--svd-mode', 'auto'])
            subprocess.run(command, check=True, capture_output=True, text=True, timeout=120)
            report = json.loads(output.read_text())
            if report['execution']['status'] != 'completed':
                raise AssertionError(f'{profile} did not complete: {report["execution"]}')
            if report['profile'] != profile or not report['measurement_scope']:
                raise AssertionError(f'{profile} report is missing its scope')
            direct = root / f'direct-{profile}.json'
            if profile == 'prepare':
                direct_command = [str(build / 'bin/joint_sparse_benchmark'), 'synthetic', 'chain', '8',
                                  'prepare', str(direct), '--resources']
            elif profile == 'fixed':
                driver = build / 'bin/joint_sparse_benchmark'
                frozen = root / 'direct-frozen.json'
                subprocess.run([str(driver), 'synthetic', 'chain', '8', 'fixed', str(frozen),
                                '--fixed', 'freeze', '--svd-mode', 'auto', '--resources'],
                               check=True, capture_output=True, text=True, timeout=60)
                direct_command = [str(driver), 'synthetic', 'chain', '8', 'fixed', str(direct),
                                  '--fixed', 'normal', '--state', str(frozen),
                                  '--fixed-preconditioner', 'schwarz', '--svd-mode', 'auto', '--resources']
            elif profile == 'solve':
                direct_command = [str(build / 'bin/joint_sparse_benchmark'), 'synthetic', 'chain', '8',
                                  'fixed', str(direct), '--search', 'schwarz', '--resources']
            elif profile == 'rank':
                direct_command = [str(build / 'bin/joint_sparse_benchmark'), 'synthetic', 'chain', '8',
                                  'rank', str(direct), '--resources']
            else:
                direct_command = [str(build / 'bin/joint_postprocessing_benchmark'), 'full',
                                  'workflow' if profile == 'workflow' else 'post', str(root / f'{profile}.sqlite')]
            direct_run = subprocess.run(direct_command, check=True, capture_output=True, text=True, timeout=60)
            if profile in ('workflow', 'postprocess'):
                previous = json.loads(direct_run.stdout.splitlines()[-1])
            else:
                previous = json.loads(direct.read_text())
            actual_result = json.loads(json.dumps(report['result']))
            expected_result = benchmark.normalize_result(profile, previous)
            if profile == 'fixed':
                timing_fields = ('preparation_seconds', 'build_seconds', 'solve_seconds', 'fixed_step_wall_seconds')
                for result in (actual_result, expected_result):
                    for step in result['details']['steps']:
                        for field in timing_fields:
                            step.pop(field, None)
            if report['problem'] != benchmark.problem_result(previous) or \
                    report['numerics'] != benchmark.numerics_result(previous) or \
                    actual_result != expected_result:
                raise AssertionError(f'{profile} changed deterministic driver results')
            if profile == 'solve':
                assessment = previous['returned_assessment']
                state = previous['returned_state']
                expected = (assessment['runtime_convergence'], assessment['design_spectrum']['rank'],
                            state['objective'], previous['search']['profile_evaluations'],
                            previous['search']['accepted_updates'])
                actual = (report['numerics']['runtime_convergence'], report['numerics']['rank'],
                          report['numerics']['objective'], report['result']['details']['profile_evaluations'],
                          report['result']['details']['accepted_updates'])
                if actual != expected:
                    raise AssertionError(f'solve profile changed deterministic driver results: {actual} != {expected}')

        output = root / 'command.json'
        subprocess.run([sys.executable, str(script), '--profile', 'command', '--case', 'command-smoke',
                        '--build-dir', str(build), '--model', str(model), '--map', str(map_path),
                        '--output', str(output), '--timeout', '60'],
                       check=True, capture_output=True, text=True, timeout=180)
        report = json.loads(output.read_text())
        if report['execution']['status'] != 'completed' or not report['result']['persistence'] or not report['result']['export']:
            raise AssertionError(f'command profile did not persist and export: {report["execution"]}')

    return 0


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--smoke':
        raise SystemExit(smoke(Path(sys.argv[2]).resolve()))
    unittest.main()
