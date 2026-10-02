"""Contract checks and deterministic smoke runs for the unified Joint benchmark."""
import json
import io
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path

import joint_benchmark as benchmark
import joint_schwarz_sweep as schwarz_sweep


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

            args.profile = 'fixed'
            args.state_path = Path('state.json')
            fixed = benchmark.command_for_profile(args, {'kind': 'synthetic', 'topology': 'chain', 'atoms': 8},
                                                  Path('fixed.json'), build)
            self.assertEqual(fixed[fixed.index('--schwarz-overlap-hops') + 1], '0')

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
            with redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    benchmark.validate_args(parser, args)

    def test_missing_telemetry_does_not_change_result_normalization(self):
        raw = {'search': {'stop_reason': 'converged', 'profile_evaluations': 3, 'accepted_updates': 2},
               'returned_assessment': {'runtime_convergence': 'passed', 'primary': {'valid': True}}}
        details = benchmark.normalize_result('solve', raw)['details']
        self.assertEqual(set(details), {'search_completed', 'stop_reason', 'profile_evaluations',
                                        'accepted_updates', 'endpoint_valid'})

    def test_solve_result_preserves_pcg_distribution_and_operator_counts(self):
        work = {'pcg_solves': 2, 'pcg_iterations': 7, 'pcg_iteration_counts': [3, 4],
                'pcg_iterations_min': 3, 'pcg_iterations_median': 3.5,
                'pcg_iterations_max': 4, 'pcg_iterations_mean': 3.5,
                'operator_normals': 9, 'operator_applications': 2, 'operator_adjoints': 3}
        details = benchmark.normalize_result('solve', {'search_work': work})['details']
        self.assertEqual(details['search_work'], work)

    def test_sweep_does_not_expand_identity_or_diagonal_over_schwarz_dimensions(self):
        parser = schwarz_sweep.build_parser()
        args = parser.parse_args(['--build-dir', 'build/debug', '--output', 'sweep.json',
                                  '--atoms', '8', '--cores', '2', '4', '--overlaps', '0', '1',
                                  '--preconditioners', 'identity', 'diagonal', 'schwarz'])
        schwarz_sweep.validate_args(parser, args)
        configurations = schwarz_sweep.build_configurations(args)
        self.assertEqual(len(configurations), 6)
        controls = [item for item in configurations if item['preconditioner'] != 'schwarz']
        self.assertEqual(len(controls), 2)
        self.assertTrue(all(item['core_atoms'] is None and item['overlap_hops'] is None for item in controls))
        self.assertEqual(len({schwarz_sweep.configuration_key(item) for item in configurations}), 6)

    def test_sweep_command_passes_overlap_only_to_schwarz(self):
        args = schwarz_sweep.build_parser().parse_args(['--build-dir', 'build/debug', '--output', 'sweep.json'])
        config = {'topology': 'chain', 'atoms': 8, 'preconditioner': 'identity',
                  'core_atoms': None, 'overlap_hops': None, 'max_block_atoms': None,
                  'storage_mib': None, 'scratch_mib': None, 'operator_rank': 'auto', 'repeat': 1}
        command = schwarz_sweep.benchmark_command(Path('joint_benchmark.py'), config, args, Path('run.json'))
        self.assertNotIn('--schwarz-core-atoms', command)
        self.assertNotIn('--schwarz-overlap-hops', command)
        self.assertEqual(command[command.index('--warmup') + 1], '0')

    def test_sweep_aggregates_completed_repetitions_and_excludes_warmup(self):
        def run(kind, index, counts, scale, search_seconds, wall_seconds, rss_bytes):
            work = {
                'pcg_solves': 2 * scale,
                'pcg_iterations': sum(counts),
                'pcg_iteration_counts': counts,
                'operator_normals': 7 + 2 * scale,
                'operator_applications': 2 * scale,
                'operator_adjoints': 3 * scale,
                'linearizations': scale,
                'damping_trials': 2 * scale,
                'local_builds': scale,
                'inverse_actions': 4 * scale,
                'partition_seconds': .1 * scale,
                'operator_prepare_seconds': .2 * scale,
                'metric_seconds': .3 * scale,
                'local_seconds': .4 * scale,
                'factor_seconds': .5 * scale,
                'inverse_seconds': .05 * scale,
                'pcg_seconds': .6 * scale,
                'operator_normal_seconds': .15 * scale,
                'operator_apply_seconds': .02 * scale,
                'operator_adjoint_seconds': .03 * scale,
                'operator_rank_seconds': .04 * scale,
                'factor_builds': 1,
            }
            details = {'search_work': work, 'search_seconds': search_seconds,
                       'assessment_seconds': 1., 'accepted_updates': scale,
                       'profile_evaluations': 3 * scale, 'stop_reason': 'converged'}
            return {'kind': kind, 'index': index, 'status': 'completed',
                    'process_wall_seconds': wall_seconds, 'peak_rss_bytes': rss_bytes,
                    'result': {'details': details}}

        config = {'topology': 'chain', 'atoms': 8, 'preconditioner': 'schwarz',
                  'core_atoms': 2, 'overlap_hops': 1, 'max_block_atoms': 8,
                  'operator_rank': 'auto', 'repeat': 3, 'warmup': 1}
        report = {
            'execution': {'status': 'completed', 'runs': [
                run('warmup', 1, [1000], 100, 1000., 1000., 1000),
                run('measurement', 1, [2, 4], 1, 4., 5., 100),
                run('measurement', 2, [4, 6], 2, 8., 10., 400),
                run('measurement', 3, [6, 10], 3, 12., 15., 300),
            ]},
            'metadata': {'solver_policy': {'resolved_rank_backend': 'Dense'}},
            'problem': {'atoms': 8, 'voxels': 120},
            'numerics': {'objective': 1.2, 'runtime_convergence': 'passed'},
        }
        row = schwarz_sweep.summarize(config, report, 'run.json')
        self.assertEqual(row['completed_repetitions'], 3)
        self.assertEqual(row['failed_repetitions'], 0)
        self.assertTrue(row['measurements_complete'])
        self.assertEqual(row['measurements']['search_seconds'], {'min': 4., 'median': 8., 'max': 12.})
        self.assertEqual(row['measurements']['partition_seconds'], {'min': .1, 'median': .2, 'max': .3})
        self.assertEqual(row['measurements']['peak_rss_bytes'], {'min': 100, 'median': 300, 'max': 400})
        self.assertEqual(row['peak_rss_bytes'], 400)
        self.assertEqual(row['measurements']['operator_normals'], {'min': 9, 'median': 11, 'max': 13})
        self.assertEqual(row['measurements']['pcg_solves'], {'min': 2, 'median': 4, 'max': 6})
        self.assertEqual(row['measurements']['factor_builds'], {'min': 1, 'median': 1, 'max': 1})
        self.assertEqual(row['pcg_iterations_per_solve'], 5.)
        self.assertEqual(row['pcg_iterations_per_solve_pooled_median'], 5.)
        self.assertEqual(row['pcg_iteration_counts_by_repetition'], [[2, 4], [4, 6], [6, 10]])

    def test_sweep_marks_failed_measurement_repetition_incomplete(self):
        config = {'topology': 'chain', 'atoms': 8, 'preconditioner': 'schwarz',
                  'core_atoms': 2, 'overlap_hops': 1, 'max_block_atoms': 8,
                  'operator_rank': 'auto', 'repeat': 3, 'warmup': 0}
        details = {'search_work': {'pcg_solves': 1, 'pcg_iterations': 2,
                                   'pcg_iteration_counts': [2]},
                   'search_seconds': 1., 'assessment_seconds': .1}
        runs = [
            {'kind': 'measurement', 'status': 'completed', 'process_wall_seconds': 1.,
             'peak_rss_bytes': 100, 'result': {'details': details}},
            {'kind': 'measurement', 'status': 'timeout'},
            {'kind': 'measurement', 'status': 'completed', 'process_wall_seconds': 2.,
             'peak_rss_bytes': 200, 'result': {'details': details}},
        ]
        row = schwarz_sweep.summarize(config, {'execution': {'status': 'timeout', 'runs': runs}}, 'run.json')
        self.assertEqual(row['completed_repetitions'], 2)
        self.assertEqual(row['failed_repetitions'], 1)
        self.assertFalse(row['measurements_complete'])
        self.assertEqual(row['measurements']['pcg_solves'], {'min': 1, 'median': 1., 'max': 1})

    def test_sweep_keeps_failed_and_unavailable_statuses(self):
        config = {'topology': 'chain', 'atoms': 8, 'preconditioner': 'schwarz',
                  'core_atoms': 2, 'overlap_hops': 0, 'max_block_atoms': 8,
                  'storage_mib': 512, 'scratch_mib': 256, 'operator_rank': 'auto', 'repeat': 1}
        unavailable = schwarz_sweep.summarize(config, {'execution': {'status': 'unavailable'}}, 'run.json')
        failed = schwarz_sweep.summarize(config, {'execution': {'status': 'timeout'}}, 'run.json')
        self.assertEqual(unavailable['status'], 'unavailable')
        self.assertEqual(failed['status'], 'timeout')


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
                if not report['result']['details'].get('rank_diagnostics') or \
                        not report['result']['details'].get('preparation_work'):
                    raise AssertionError('fixed profile omitted rank or preparation work telemetry')

            def strip_timings(value):
                if isinstance(value, dict):
                    for field in list(value):
                        if field.endswith('_seconds'):
                            value.pop(field)
                        else:
                            strip_timings(value[field])
                elif isinstance(value, list):
                    for item in value:
                        strip_timings(item)

            strip_timings(actual_result)
            strip_timings(expected_result)
            if report['problem'] != benchmark.problem_result(previous) or \
                    report['numerics'] != benchmark.numerics_result(previous) or \
                    actual_result != expected_result:
                raise AssertionError(f'{profile} changed deterministic driver results: '
                                     f'{actual_result.get("details")} != {expected_result.get("details")}')
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
                policy = report['metadata']['solver_policy']
                if policy['resolved_rank_backend'] == 'SpqrBounds':
                    work = report['result']['details']['search_work']
                    if (work['operator_rank_status'] != 'full-rank' or work['operator_rank_checks'] <= 0 or
                            work['operator_rank_compact_extractions'] != 0 or
                            work['operator_rank_free_design_svds'] != 0):
                        raise AssertionError(f'SPQR bounded rank used an unexpected path: {work}')

        output = root / 'command.json'
        subprocess.run([sys.executable, str(script), '--profile', 'command', '--case', 'command-smoke',
                        '--build-dir', str(build), '--model', str(model), '--map', str(map_path),
                        '--output', str(output), '--timeout', '60'],
                       check=True, capture_output=True, text=True, timeout=180)
        report = json.loads(output.read_text())
        if report['execution']['status'] != 'completed' or not report['result']['persistence'] or not report['result']['export']:
            raise AssertionError(f'command profile did not persist and export: {report["execution"]}')

        sweep_output = root / 'sweep.json'
        sweep_csv = root / 'sweep.csv'
        subprocess.run([sys.executable, str(Path(__file__).with_name('joint_schwarz_sweep.py')),
                        '--build-dir', str(build), '--output', str(sweep_output), '--csv', str(sweep_csv),
                        '--topologies', 'chain', '--atoms', '8', '--cores', '2', '--overlaps', '0', '1',
                        '--preconditioners', 'schwarz', '--timeout', '60'],
                       check=True, capture_output=True, text=True, timeout=180)
        aggregate = json.loads(sweep_output.read_text())
        configurations = aggregate['configurations']
        if len(configurations) != 2 or len({row['configuration']['overlap_hops'] for row in configurations}) != 2:
            raise AssertionError(f'sweep contains duplicate or missing configurations: {configurations}')
        if not sweep_csv.is_file():
            raise AssertionError('sweep CSV summary was not written')
        for row in configurations:
            individual = root / row['individual_json']
            if not individual.is_file() or row['status'] != 'completed':
                raise AssertionError(f'sweep run was not completed: {row}')
            item = json.loads(individual.read_text())
            if row['objective'] != item['numerics']['objective'] or row['runtime_convergence'] != item['numerics']['runtime_convergence']:
                raise AssertionError(f'sweep summary does not match individual result: {row}')
            policy = item['metadata']['solver_policy']
            if policy['schwarz_core_atoms'] != 2 or policy['schwarz_overlap_hops'] != row['configuration']['overlap_hops']:
                raise AssertionError(f'sweep policy was not propagated: {policy}')
    return 0


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--smoke':
        raise SystemExit(smoke(Path(sys.argv[2]).resolve()))
    unittest.main()
