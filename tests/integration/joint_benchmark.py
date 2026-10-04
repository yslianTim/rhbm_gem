"""Run explicit Joint benchmark profiles with shared process and provenance support."""
from __future__ import annotations

import argparse
import json
import math
import re
import statistics
import subprocess
import tempfile
import time
from pathlib import Path

from experiment_io import ROOT, digest, sha, write
from experiment_process import ENV, monitored
from experiment_provenance import source_hash
from joint_fixture_records import CATALOG
from joint_runtime_support import read, unpack


SCOPES = {
    'prepare': 'joint_problem_construction_and_basis_preparation',
    'fixed': 'fixed_state_operator_preparation_and_selected_step',
    'search': 'search-only',
    'solve': 'joint_search_and_returned_state_assessment',
    'rank': 'bounded_free_design_rank_evaluation',
    'workflow': 'prepare_estimate_postprocess_and_persist',
    'postprocess': 'uncertainty_peeling_summary_and_persist',
    'command': 'cli_analysis_sqlite_reload_and_joint_export',
}


def cache_value(cache: str, name: str, default=None):
    prefix = name + ':'
    for line in cache.splitlines():
        if line.startswith(prefix) and '=' in line:
            return line.split('=', 1)[1]
    return default


def build_metadata(build: Path):
    cache_path = build / 'CMakeCache.txt'
    if not cache_path.is_file():
        raise ValueError(f'No CMake build found at {build}')
    cache = cache_path.read_text()
    configuration = build / 'generated' / cache_value(cache, 'CMAKE_BUILD_TYPE', 'Debug') / 'SimulationConfiguration-CXX.txt'
    if not configuration.is_file():
        configurations = sorted((build / 'generated').glob('*/SimulationConfiguration-CXX.txt'))
        configuration = configurations[0] if configurations else configuration
    values = {}
    if configuration.is_file():
        for line in configuration.read_text().splitlines():
            if '=' in line:
                key, value = line.split('=', 1)
                values[key] = value
    compiler = values.get('compiler', '').split(';')
    return {
        'backend': cache_value(cache, 'RHBM_GEM_JOINT_SPARSE_BACKEND', 'unknown'),
        'build_type': cache_value(cache, 'CMAKE_BUILD_TYPE', 'multi-config'),
        'compiler': ' '.join(compiler[1:]) if len(compiler) > 1 else None,
    }


def process_status(record):
    status = record.get('status')
    if status == 'completed':
        return 'completed'
    if status in ('time-limit', 'rss-limit', 'rss-limit-observed-after-exit'):
        return 'timeout' if status == 'time-limit' else 'rss_limit'
    if status == 'not-run-budget':
        return 'not_run'
    return 'process_error'


def process_record(record):
    result = {key: record.get(key) for key in (
        'status', 'exit_code', 'wall_seconds', 'sampled_tree_peak_rss_bytes',
        'os_process_peak_rss_bytes', 'samples', 'maximum_sampling_gap_seconds')}
    result['raw_status'] = result['status']
    result['status'] = process_status(record)
    return result


def spqr_ordering_unavailable_reason(stderr):
    return next((line.strip() for line in stderr.splitlines()
                 if 'spqr-ordering-unavailable:' in line), None)


def synthetic_case(case):
    match = re.fullmatch(r'(chain|cube)-(\d+)', case)
    if not match:
        return None
    return match.group(1), int(match.group(2))


def case_input(args, work):
    synthetic = synthetic_case(args.case)
    if synthetic:
        topology, atoms = synthetic
        return {
            'kind': 'synthetic', 'topology': topology, 'atoms': atoms,
            'input_sha256': digest({
                'generator': sha(ROOT / 'tests/support/JointOperatorWorkload.cpp'),
                'topology': topology, 'atoms': atoms,
            }),
        }
    if args.profile == 'prepare' or ':' not in args.case:
        raise ValueError('Use chain-N or cube-N; fixture cases use DATASET:CASE.')
    dataset, fixture_case = args.case.split(':', 1)
    catalog = read(CATALOG)
    entry = catalog['datasets'].get(dataset)
    if not entry:
        raise ValueError(f'Unknown Joint fixture: {dataset}')
    fixture = unpack(CATALOG, dataset, work)
    cases = read(fixture / 'cases.json')
    if fixture_case not in cases:
        raise ValueError(f'Unknown Joint fixture case: {args.case}')
    files = {path.name: sha(path) for path in sorted(fixture.iterdir()) if path.is_file()}
    return {
        'kind': 'fixture', 'directory': fixture, 'fixture_case': fixture_case,
        'input_sha256': digest(files), 'fixture_sha256': entry['sha256'],
    }


def sparse_command(driver, case, phase, output, options=(), svd_mode=None):
    if case['kind'] == 'synthetic':
        command = [driver, 'synthetic', case['topology'], str(case['atoms']), phase, output]
    else:
        command = [driver, 'fixture', case['directory'], case['fixture_case'], output]
    mode_options = ('--svd-mode', svd_mode) if svd_mode else ()
    return [*map(str, command), *options, *mode_options, '--resources']


def run_process(command, directory, deadline, rss_limit, seconds):
    result = monitored(command, directory, deadline, rss_limit=rss_limit, seconds=seconds)
    return process_record(result), result


def driver_json(directory, output=None):
    if output and output.is_file():
        return read(output)
    lines = (directory / 'stdout.txt').read_text(errors='replace').splitlines()
    for line in reversed(lines):
        try:
            return json.loads(line)
        except json.JSONDecodeError:
            continue
    raise ValueError('Benchmark driver did not produce a JSON result')


def partial_driver_json(profile, directory, output=None):
    if profile not in ('workflow', 'postprocess'):
        try:
            return driver_json(directory, output)
        except (OSError, ValueError, json.JSONDecodeError):
            return None
    return None


def command_for_profile(args, case, output, build):
    sparse = build / 'bin/joint_sparse_benchmark'
    workflow = build / 'bin/joint_postprocessing_benchmark'
    if args.profile in ('prepare', 'fixed', 'search', 'solve', 'rank'):
        if not sparse.is_file():
            raise ValueError('Build joint_sparse_benchmark with RHBM_GEM_BUILD_BENCHMARKS=ON')
        if args.profile == 'prepare':
            return sparse_command(sparse, case, 'prepare', output, svd_mode=args.svd_mode)
        if args.profile == 'fixed':
            return sparse_command(sparse, case, 'fixed', output, (
                '--fixed', args.fixed_action, '--state', args.state_path,
                '--fixed-preconditioner', args.preconditioner,
                *operator_policy_options(args)), args.svd_mode)
        if args.profile == 'search':
            ordering = ('--spqr-ordering', args.spqr_ordering) if args.spqr_ordering else ()
            representation = (('--operator-factor-representation', args.operator_factor_representation)
                              if args.operator_factor_representation != 'exported-fixed' else ())
            ownership_name = args.operator_factor_ownership or 'reuse-accepted-copy-on-write'
            ownership = (('--operator-factor-ownership', ownership_name)
                         if ownership_name != 'reuse-accepted-copy-on-write' else ())
            return sparse_command(sparse, case, 'fixed', output,
                                  ('--search', args.preconditioner, '--search-only',
                                   *operator_policy_options(args), *ordering, *representation, *ownership), args.svd_mode)
        if args.profile == 'solve':
            return sparse_command(sparse, case, 'fixed', output,
                                  ('--search', args.preconditioner, *operator_policy_options(args),
                                   '--assessment-reduction', args.assessment_reduction,
                                   '--projected-reduction', args.projected_reduction), args.svd_mode)
        return sparse_command(sparse, case,
                              'rank' if args.rank_mode == 'prototype' else 'rank-oracle', output,
                              rank_budget_options(args), args.svd_mode)
    if args.profile in ('workflow', 'postprocess'):
        if not workflow.is_file():
            raise ValueError('Build joint_postprocessing_benchmark with RHBM_GEM_BUILD_BENCHMARKS=ON')
        return [str(workflow), args.case, 'workflow' if args.profile == 'workflow' else 'post', str(output)]
    cli = args.cli or build / 'bin/RHBM-GEM'
    if not cli.is_file():
        raise ValueError('Build the RHBM-GEM CLI before using the command profile')
    if not args.model or not args.map_path:
        return None
    database = output / 'joint.sqlite'
    export = output / 'export'
    analysis = [str(cli), 'potential_analysis', '--estimator', 'joint-components',
                '--only-backbone', 'true', '-a', str(args.model), '-m', str(args.map_path),
                '-d', str(database), '-k', args.case, '--map-normalization', 'false', '-v', '0']
    dump = [str(cli), 'result_dump', '--printer', 'joint', '-d', str(database),
            '-k', args.case, '-o', str(export)]
    return analysis, dump, database, export


def rank_budget_options(args):
    return (
        '--operator-rank-seconds', str(args.operator_rank_seconds),
        '--operator-rank-work-entries', str(args.operator_rank_work_entries),
        '--operator-rank-workspace-mib', str(args.operator_rank_workspace_mib),
    )


def operator_policy_options(args):
    return (
        '--operator-rank', args.operator_rank,
        *rank_budget_options(args),
        '--schwarz-core-atoms', str(args.schwarz_core_atoms),
        '--schwarz-overlap-hops', str(args.schwarz_overlap_hops),
        '--schwarz-max-block-atoms', str(args.schwarz_max_block_atoms),
        '--schwarz-storage-mib', str(args.schwarz_storage_mib),
        '--schwarz-scratch-mib', str(args.schwarz_scratch_mib),
    )


def solver_policy_metadata(args, backend):
    rank_profile = args.profile == 'rank'
    active = (args.profile == 'fixed' or
              (args.profile in ('search', 'solve') and args.preconditioner != 'legacy') or rank_profile)
    resolved = None
    backend_name = backend.upper()
    if rank_profile:
        resolved = 'Dense' if args.rank_mode == 'oracle' else 'SpqrBounds'
        rank_mode = 'dense' if args.rank_mode == 'oracle' else 'spqr-bounds'
    elif active:
        if args.operator_rank == 'dense' or (backend_name == 'EIGEN' and args.operator_rank == 'auto'):
            resolved = 'Dense'
        elif args.operator_rank in ('auto', 'spqr-bounds') and backend_name == 'SPQR':
            resolved = 'SpqrBounds'
        rank_mode = args.operator_rank
    else:
        rank_mode = args.operator_rank
    return {
        'search_method': ('OperatorPcg' if active and not rank_profile else
                          'LegacyCompact' if args.profile == 'solve' else None),
        'sparse_backend': backend,
        'operator_rank_active': active,
        'operator_rank_mode': rank_mode,
        'resolved_rank_backend': resolved,
        'operator_rank_budget_seconds': args.operator_rank_seconds,
        'operator_rank_budget_entries': args.operator_rank_work_entries,
        'operator_rank_budget_workspace_bytes': args.operator_rank_workspace_mib * 1024**2,
        'operator_factor_representation': args.operator_factor_representation,
        'operator_factor_ownership': (args.operator_factor_ownership or 'reuse-accepted-copy-on-write')
            if args.profile == 'search' else None,
        'projected_reduction': args.projected_reduction,
        'preconditioner': args.preconditioner,
        'spqr_ordering': args.spqr_ordering.upper() if args.spqr_ordering else 'COLAMD',
        'schwarz_core_atoms': args.schwarz_core_atoms,
        'schwarz_overlap_hops': args.schwarz_overlap_hops,
        'schwarz_max_block_atoms': args.schwarz_max_block_atoms,
        'schwarz_storage_bytes': args.schwarz_storage_mib * 1024**2,
        'schwarz_scratch_bytes': args.schwarz_scratch_mib * 1024**2,
    }


def raw_elapsed(profile, raw, process_wall):
    if profile == 'prepare':
        return sum(raw.get(key, 0.) for key in ('construction_seconds', 'basis_seconds'))
    if profile == 'fixed':
        steps = raw.get('steps', [])
        return steps[0].get('fixed_step_wall_seconds', process_wall) if steps else process_wall
    if profile == 'solve':
        return sum(raw.get(key, 0.) for key in ('search_seconds', 'assessment_seconds'))
    if profile == 'search':
        return raw.get('search_seconds', process_wall)
    if profile == 'rank':
        return raw.get('rank_wall_seconds', process_wall)
    if profile == 'workflow':
        return raw.get('total_seconds', process_wall)
    if profile == 'postprocess':
        return sum(raw.get('phases', {}).values())
    return process_wall


def normalize_result(profile, raw):
    assessment = raw.get('returned_assessment') or {}
    rank_result = raw.get('rank_result') or {}
    convergence = assessment.get('runtime_convergence')
    qualified = convergence == 'passed' if convergence in ('passed', 'failed') else None
    if profile == 'search':
        qualified = None
    if profile == 'fixed' and isinstance(raw.get('valid'), bool):
        qualified = raw['valid']
    if profile == 'prepare' and isinstance(raw.get('raw_state_valid'), bool):
        qualified = raw['raw_state_valid']

    details = {}
    if profile == 'prepare':
        details = {key: raw.get(key) for key in ('generator', 'topology', 'memberships', 'raw_state_valid')}
    elif profile == 'fixed':
        details = {
            'valid': raw.get('valid'), 'reason': raw.get('reason'), 'mode': raw.get('mode'),
            'steps': [{key: step.get(key) for key in (
                'kind', 'valid', 'iterations', 'true_residual', 'predicted',
                'preparation_seconds', 'build_seconds', 'solve_seconds', 'fixed_step_wall_seconds')}
                      for step in raw.get('steps', [])],
        }
        for key in ('rank_diagnostics', 'preparation_work', 'partition'):
            if key in raw:
                details[key] = raw[key]
    elif profile == 'solve':
        search = raw.get('search') or {}
        primary = assessment.get('primary') or {}
        details = {
            'search_completed': raw.get('search_completed'), 'stop_reason': search.get('stop_reason'),
            'profile_evaluations': search.get('profile_evaluations'),
            'accepted_updates': search.get('accepted_updates'), 'endpoint_valid': primary.get('valid'),
            'assessment_execution': raw.get('assessment_execution'),
            'assessment_telemetry': raw.get('assessment_telemetry'),
            'derivative_reduction_micro_attribution': (raw.get('assessment_telemetry') or {}).get(
                'derivative_reduction_micro_attribution'),
            'assessment_stage': raw.get('assessment_stage'),
            'last_assessment_stage': raw.get('last_assessment_stage'),
            'active_assessment_stage': raw.get('active_assessment_stage'),
            'completed_assessment_stages': raw.get('completed_assessment_stages'),
            'completed_stage_seconds': raw.get('completed_stage_seconds'),
            'stage_dimensions': raw.get('stage_dimensions'),
        }
        if 'search_seconds' in raw:
            details['search_seconds'] = raw['search_seconds']
        if 'assessment_seconds' in raw:
            details['assessment_seconds'] = raw['assessment_seconds']
        if 'returned_state' in raw:
            details['state_available'] = raw['returned_state'] is not None
        for key in ('search_work', 'partition'):
            if key in raw:
                details[key] = raw[key]
    elif profile == 'search':
        search = raw.get('search') or {}
        details = {
            'measurement_scope': raw.get('measurement_scope'),
            'scope_description': raw.get('scope_description'),
            'assessment_execution': raw.get('assessment_execution'),
            'assessment_work': raw.get('assessment_work'),
            'search_completed': search.get('execution_complete'),
            'stop_reason': search.get('stop_reason'),
            'profile_evaluations': search.get('profile_evaluations'),
            'accepted_updates': search.get('accepted_updates'),
            'accepted_objective': search.get('accepted_objective'),
            'accepted_gradient_inf_norm': search.get('accepted_gradient_inf_norm'),
            'returned_search_state': search.get('returned_search_state'),
            'search_seconds': raw.get('search_seconds'),
            'free_columns': raw.get('free_columns'),
            'design_nonzeros': raw.get('design_nonzeros'),
            'search_work': raw.get('search_work'),
        }
        for key in ('failure_stage', 'active_search_stage', 'last_completed_search_stage',
                    'completed_search_stages', 'stage_seconds', 'stage_calls', 'stage_completed_calls',
                    'stage_dimensions', 'stage_nnz',
                    'spqr_factorization', 'factor_residency', 'operator_factor_representation',
                    'operator_factor_ownership'):
            if key in raw:
                details[key] = raw[key]
        if 'partition' in raw:
            details['partition'] = raw['partition']
    elif profile == 'rank':
        details = {key: rank_result.get(key) for key in (
            'status', 'reason', 'rank_lower', 'rank_upper', 'exact_rank', 'rank',
            'threshold', 'minimum_lower', 'maximum_upper', 'rank_backend', 'certificate', 'work_stage',
            'entries', 'seconds', 'workspace_bytes', 'estimated_total_entries',
            'estimated_remaining_entries', 'estimated_reconstruction_entries', 'design_nonzeros',
            'r_nonzeros', 'reflector_nonzeros', 'reflectors', 'threshold_lower', 'threshold_upper')}
        details['rows'] = raw.get('rank_rows', raw.get('rows'))
        details['columns'] = raw.get('free_columns')
        details['compact_extractions'] = raw.get('rank_compact_extractions')
        details['free_design_svds'] = (raw.get('work') or {}).get('free_design_svds')
        details['solver_policy'] = raw.get('solver_policy')
        details['local_witness'] = raw.get('local_witness')
        details['local_witness_seconds'] = raw.get('local_witness_seconds')
    elif profile in ('workflow', 'postprocess'):
        details = {
            'endpoint_count': len(raw.get('endpoints', [])),
            'target_count': len(raw.get('targets', [])),
            'objective': raw.get('objective'),
            'endpoint_statuses': [item.get('convergence') for item in raw.get('endpoints', [])],
        }
    return {
        'qualified': qualified,
        'scientific_status': ('search-only' if profile == 'search' else
                              convergence or rank_result.get('status') or 'not-assessed'),
        'details': details,
    }


def problem_result(raw):
    atoms = raw.get('atoms')
    endpoints = raw.get('endpoints')
    targets = raw.get('targets')
    if atoms is None and isinstance(targets, list):
        atoms = len(targets)
    components = raw.get('components')
    if components is None and isinstance(endpoints, list):
        components = len(endpoints)
    rows = raw.get('rows')
    return {'atoms': atoms, 'voxels': rows, 'components': components,
            'free_columns': raw.get('free_columns'), 'design_nonzeros': raw.get('design_nonzeros'),
            'parameters': 3 * atoms if isinstance(atoms, int) else None}


def numerics_result(raw):
    state = raw.get('returned_state') or raw.get('state_control') or {}
    assessment = raw.get('returned_assessment') or {}
    primary = assessment.get('primary') or {}
    rank = raw.get('rank_result')
    if isinstance(rank, dict):
        rank = rank.get('exact_rank')
        if rank is None:
            rank = raw.get('rank_result', {}).get('rank')
        if rank is None:
            rank = raw.get('rank_result', {}).get('rank_lower')
    if rank is None:
        rank = (assessment.get('design_spectrum') or {}).get('rank')
    if rank is None and raw.get('rank'):
        rank = raw['rank'][0].get('rank')
    gradient = raw.get('operator_gradient') or state.get('gradient')
    gradient_norm = (math.sqrt(sum(x * x for x in gradient))
                     if gradient and all(isinstance(x, (int, float)) and math.isfinite(x) for x in gradient)
                     else None)
    return {
        'rank': rank, 'objective': state.get('objective', raw.get('objective')),
        'gradient_norm': gradient_norm,
        'runtime_convergence': assessment.get('runtime_convergence'),
        'active_atoms': primary.get('active_atoms'),
    }


def execute_once(args, case, build, run_root, deadline):
    results = []
    raw = None
    rss = 0
    if args.profile == 'command' and (not args.model or not args.map_path or
                                      not args.model.is_file() or not args.map_path.is_file()):
        return {'status': 'unavailable', 'reason': 'command profile requires --model and --map',
                'elapsed_seconds': None, 'peak_rss_bytes': None, 'raw': None}

    if args.profile == 'fixed':
        frozen = run_root / 'frozen.json'
        command = sparse_command(build / 'bin/joint_sparse_benchmark', case, 'fixed', frozen,
                                 ('--fixed', 'freeze'), args.svd_mode)
        record, _ = run_process(command, run_root / 'setup', deadline, args.rss_limit, args.timeout)
        results.append({'role': 'state_setup', **record})
        if record['status'] != 'completed':
            return {'status': record['status'], 'elapsed_seconds': None, 'peak_rss_bytes': None,
                    'stages': results, 'raw': None}
        args.state_path = frozen

    try:
        if args.profile in ('workflow', 'postprocess'):
            driver_output = run_root / 'benchmark.sqlite'
        elif args.profile == 'command':
            driver_output = run_root / 'command-output'
            driver_output.mkdir()
        else:
            driver_output = run_root / 'driver.json'
        command = command_for_profile(args, case, driver_output, build)
    except ValueError as error:
        return {'status': 'unavailable', 'reason': str(error), 'elapsed_seconds': None,
                'peak_rss_bytes': None, 'stages': results, 'raw': None}
    if command is None:
        return {'status': 'unavailable', 'reason': 'command profile requires --model and --map',
                'elapsed_seconds': None, 'peak_rss_bytes': None, 'stages': results, 'raw': None}

    if args.profile == 'command':
        analysis, dump, database, export = command
        for role, stage_command in (('analysis', analysis), ('reload_export', dump)):
            directory = run_root / role
            record, _ = run_process(stage_command, directory, deadline, args.rss_limit, args.timeout)
            results.append({'role': role, **record})
            rss = max(rss, record.get('sampled_tree_peak_rss_bytes') or 0,
                      record.get('os_process_peak_rss_bytes') or 0)
            if record['status'] != 'completed':
                return {'status': record['status'], 'elapsed_seconds': None, 'peak_rss_bytes': None,
                        'stages': results, 'raw': None}
        exported_json = export / f'joint_result_{args.case}.json'
        exported_csv = export / f'joint_atoms_{args.case}.csv'
        persistence = database.is_file()
        raw = read(exported_json) if exported_json.is_file() else None
        export_ok = exported_json.is_file() and exported_csv.is_file()
        status = 'completed' if persistence and export_ok else 'process_error'
        return {'status': status, 'elapsed_seconds': sum(row.get('wall_seconds') or 0 for row in results),
                'peak_rss_bytes': rss, 'stages': results, 'persistence_result': persistence,
                'export_result': export_ok, 'export_bytes': exported_json.stat().st_size if export_ok else None,
                'raw': raw}

    directory = run_root / 'measured'
    record, _ = run_process(command, directory, deadline, args.rss_limit, args.timeout)
    results.append({'role': 'measurement', **record})
    rss = max(record.get('sampled_tree_peak_rss_bytes') or 0,
              record.get('os_process_peak_rss_bytes') or 0)
    if record['status'] != 'completed':
        stderr = (directory / 'stderr.txt').read_text(errors='replace') if (directory / 'stderr.txt').is_file() else ''
        reason = spqr_ordering_unavailable_reason(stderr)
        if reason:
            return {'status': 'unavailable', 'reason': reason, 'elapsed_seconds': record.get('wall_seconds'),
                    'process_wall_seconds': record.get('wall_seconds'), 'peak_rss_bytes': rss,
                    'stages': results, 'raw': None}
        raw = partial_driver_json(args.profile, directory,
                                  driver_output if args.profile not in ('workflow', 'postprocess') else None)
        raw_peak = (raw.get('peak_rss_bytes') or raw.get('process_peak_rss_bytes') or 0) if raw else 0
        return {'status': record['status'], 'reason': record.get('raw_status'),
                'elapsed_seconds': (raw_elapsed(args.profile, raw, record.get('wall_seconds') or 0)
                                    if raw else None),
                'process_wall_seconds': record.get('wall_seconds'),
                'peak_rss_bytes': max(rss, raw_peak), 'stages': results, 'raw': raw}
    try:
        raw = driver_json(directory, driver_output if args.profile not in ('workflow', 'postprocess') else None)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        return {'status': 'process_error', 'reason': str(error), 'elapsed_seconds': None,
                'peak_rss_bytes': None, 'stages': results, 'raw': None}
    return {'status': 'completed', 'elapsed_seconds': raw_elapsed(args.profile, raw, record.get('wall_seconds') or 0),
            'process_wall_seconds': record.get('wall_seconds'),
            'peak_rss_bytes': max(rss, raw.get('peak_rss_bytes') or raw.get('process_peak_rss_bytes') or 0),
            'stages': results, 'raw': raw}


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', choices=tuple(SCOPES), required=True)
    parser.add_argument('--case', required=True)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeat', type=int, default=1)
    parser.add_argument('--warmup', type=int, default=0)
    parser.add_argument('--timeout', type=float, default=600)
    parser.add_argument('--rss-limit', type=int, default=4 * 1024**3)
    parser.add_argument('--preconditioner', choices=('legacy', 'identity', 'diagonal', 'schwarz'), default='schwarz')
    parser.add_argument('--spqr-ordering', choices=('colamd', 'default', 'best', 'metis'))
    parser.add_argument('--operator-factor-representation', choices=('exported-fixed', 'native-qr'),
                        default='exported-fixed')
    parser.add_argument('--operator-factor-ownership', choices=('dedicated-fixed', 'dedicated-native',
                        'reuse-accepted-copy-on-write', 'reuse-accepted-handoff'),
                        default=None)
    parser.add_argument('--operator-rank', choices=('auto', 'dense', 'spqr-bounds'), default='auto')
    parser.add_argument('--operator-rank-seconds', type=float, default=120)
    parser.add_argument('--operator-rank-work-entries', type=int, default=100_000_000)
    parser.add_argument('--operator-rank-workspace-mib', type=int, default=256)
    parser.add_argument('--schwarz-core-atoms', type=int, default=128)
    parser.add_argument('--schwarz-overlap-hops', type=int, default=1)
    parser.add_argument('--schwarz-max-block-atoms', type=int, default=512)
    parser.add_argument('--schwarz-storage-mib', type=int, default=512)
    parser.add_argument('--schwarz-scratch-mib', type=int, default=256)
    parser.add_argument('--fixed-action', choices=('composed', 'normal'), default='normal')
    parser.add_argument('--svd-mode', choices=('legacy', 'values', 'auto'))
    parser.add_argument('--assessment-reduction', choices=('observation-tsqr', 'compact-stack-qr'),
                        default='observation-tsqr')
    parser.add_argument('--projected-reduction', choices=('observation-tiled-qr', 'structured-compact-qr',
                        'projected-tail-census'),
                        default='observation-tiled-qr')
    parser.add_argument('--rank-mode', choices=('prototype', 'oracle'), default='prototype')
    parser.add_argument('--cli', type=Path)
    parser.add_argument('--model', type=Path)
    parser.add_argument('--map', dest='map_path', type=Path)
    return parser


def validate_args(parser, args):
    if args.repeat < 1 or args.warmup < 0 or args.timeout <= 0 or args.rss_limit <= 0:
        parser.error('--repeat and limits must be positive; --warmup must be nonnegative')
    if (not math.isfinite(args.operator_rank_seconds) or args.operator_rank_seconds < 0 or
            args.operator_rank_work_entries < 0 or args.operator_rank_workspace_mib < 0):
        parser.error('operator rank budgets must be nonnegative and the time budget finite')
    if (args.schwarz_core_atoms <= 0 or args.schwarz_overlap_hops < 0 or
            args.schwarz_max_block_atoms < args.schwarz_core_atoms or
            args.schwarz_storage_mib <= 0 or args.schwarz_scratch_mib <= 0):
        parser.error('Schwarz core, max block, and memory limits must be positive; overlap must be nonnegative and max block must cover core')
    if args.profile == 'fixed' and args.preconditioner == 'legacy':
        parser.error('fixed profile requires identity, diagonal, or schwarz preconditioner')
    if args.profile == 'search' and args.preconditioner == 'legacy':
        parser.error('search profile requires identity, diagonal, or schwarz preconditioner')
    if args.spqr_ordering and args.profile != 'search':
        parser.error('--spqr-ordering is supported by the benchmark-only search profile')
    if args.operator_factor_representation != 'exported-fixed' and args.profile != 'search':
        parser.error('--operator-factor-representation is supported only by the benchmark-only search profile')
    if args.operator_factor_ownership is not None and args.profile != 'search':
        parser.error('--operator-factor-ownership is supported only by the benchmark-only search profile')
    if args.assessment_reduction != 'observation-tsqr' and args.profile != 'solve':
        parser.error('--assessment-reduction is supported only by the benchmark-only solve profile')
    if args.projected_reduction != 'observation-tiled-qr' and args.profile != 'solve':
        parser.error('--projected-reduction is supported only by the benchmark-only solve profile')


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)

    output = args.output.resolve()
    build = args.build_dir.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        metadata = build_metadata(build)
        if args.projected_reduction in ('structured-compact-qr', 'projected-tail-census') and metadata['backend'].upper() != 'SPQR':
            parser.error(f'{args.projected_reduction} requires the SPQR benchmark backend')
        commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, check=True,
                                capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        parser.error(str(error))
    metadata.update(commit=commit, profile=args.profile, case=args.case,
                    assessment_reduction=args.assessment_reduction,
                    projected_reduction=args.projected_reduction,
                    source_sha256=source_hash(ROOT), benchmark_sha256=sha(Path(__file__)))
    metadata['solver_policy'] = solver_policy_metadata(args, metadata['backend'])
    if args.profile in ('prepare', 'fixed', 'search', 'solve', 'rank'):
        driver = build / 'bin/joint_sparse_benchmark'
    elif args.profile in ('workflow', 'postprocess'):
        driver = build / 'bin/joint_postprocessing_benchmark'
    else:
        driver = args.cli or build / 'bin/RHBM-GEM'
    metadata['driver_sha256'] = sha(driver) if driver.is_file() else None

    if args.profile in ('workflow', 'postprocess') and args.case not in ('full', 'halo', 'multi'):
        parser.error('workflow and postprocess cases are full, halo, or multi')
    if args.profile == 'command' and not re.fullmatch(r'[A-Za-z0-9._-]+', args.case):
        parser.error('command case names may contain letters, digits, dot, underscore, and hyphen')
    if args.profile in ('prepare', 'fixed', 'search', 'solve', 'rank') and not synthetic_case(args.case) and ':' not in args.case:
        parser.error('Joint cases use chain-N, cube-N, or DATASET:CASE')

    work_parent = output.parent
    with tempfile.TemporaryDirectory(prefix='joint-benchmark-', dir=work_parent) as temporary:
        work = Path(temporary)
        try:
            case = case_input(args, work / 'inputs') if args.profile in ('prepare', 'fixed', 'search', 'solve', 'rank') else {
                'kind': 'synthetic', 'input_sha256': digest({
                    'driver': sha(build / 'bin/joint_postprocessing_benchmark')
                    if args.profile in ('workflow', 'postprocess') and (build / 'bin/joint_postprocessing_benchmark').is_file()
                    else sha(args.cli or build / 'bin/RHBM-GEM')
                    if args.profile == 'command' and (args.cli or build / 'bin/RHBM-GEM').is_file()
                    else None,
                    'case': args.case,
                })}
        except (KeyError, ValueError, OSError, RuntimeError) as error:
            parser.error(str(error))
        metadata['input_sha256'] = case['input_sha256']
        if 'fixture_sha256' in case:
            metadata['fixture_sha256'] = case['fixture_sha256']
        if (args.profile == 'command' and args.model and args.map_path and
                args.model.is_file() and args.map_path.is_file()):
            metadata['input_sha256'] = digest({'model': sha(args.model), 'map': sha(args.map_path)})
        records = []
        overall = 'completed'
        for kind, count in (('warmup', args.warmup), ('measurement', args.repeat)):
            for index in range(count):
                run_root = work / f'{kind}-{index + 1}'
                run_root.mkdir(parents=True)
                run_args = argparse.Namespace(**vars(args))
                result = execute_once(run_args, case, build, run_root,
                                      time.monotonic() + args.timeout)
                record = {'kind': kind, 'index': index + 1, **result}
                records.append(record)
                if result['status'] not in ('completed', 'unavailable'):
                    overall = result['status']
                if result['status'] == 'unavailable':
                    overall = 'unavailable'

        measurements = [record for record in records if record['kind'] == 'measurement']
        successful = [record for record in measurements if record['status'] == 'completed']
        all_measured = len(successful) == args.repeat
        elapsed = [record['elapsed_seconds'] for record in successful if record.get('elapsed_seconds') is not None]
        peaks = [record['peak_rss_bytes'] for record in successful if record.get('peak_rss_bytes') is not None]
        last_raw = next((record['raw'] for record in reversed(successful) if record.get('raw') is not None), None)
        run_summaries = []
        for record in records:
            run = {key: record.get(key) for key in (
                'kind', 'index', 'status', 'reason', 'elapsed_seconds',
                'process_wall_seconds', 'peak_rss_bytes', 'persistence_result',
                'export_result', 'export_bytes')}
            raw = record.get('raw') or {}
            run['problem'] = problem_result(raw)
            run['numerics'] = numerics_result(raw)
            run['result'] = normalize_result(args.profile, raw)
            run['driver_stage'] = raw.get('stage')
            run['measurement_scope'] = raw.get('measurement_scope')
            run['stages'] = [{key: stage.get(key) for key in (
                'role', 'status', 'exit_code', 'wall_seconds',
                'sampled_tree_peak_rss_bytes', 'os_process_peak_rss_bytes',
                'samples', 'maximum_sampling_gap_seconds')}
                             for stage in record.get('stages', [])]
            run_summaries.append(run)

        report = {
            'schema_version': 1, 'profile': args.profile, 'measurement_scope': SCOPES[args.profile],
            'case': args.case, 'metadata': metadata,
            'execution': {
                'status': overall,
                'returncode': 0 if all_measured and all(r['status'] == 'completed' for r in records) else None,
                'completed_repetitions': len(successful), 'failed_repetitions': args.repeat - len(successful),
                'runs': run_summaries,
            },
            'problem': problem_result(last_raw or {}),
            'numerics': numerics_result(last_raw or {}),
            'resources': {
                'elapsed_seconds': statistics.median(elapsed) if all_measured and elapsed else None,
                'peak_rss_bytes': max((r.get('peak_rss_bytes') or 0 for r in records), default=0) or None,
                'completed_measurement_peak_rss_bytes': max(peaks) if all_measured and peaks else None,
                'count': args.repeat, 'failed_count': args.repeat - len(successful),
                'samples': [r.get('elapsed_seconds') if r['status'] == 'completed' else None for r in measurements],
            },
            'result': normalize_result(args.profile, last_raw or {}),
        }
        if args.profile == 'command' and records:
            final = next((r for r in reversed(records) if r['kind'] == 'measurement'), {})
            report['result'].update(persistence=final.get('persistence_result'), export=final.get('export_result'),
                                    export_bytes=final.get('export_bytes'))
        write(output, report)
    print(json.dumps({'profile': args.profile, 'status': overall, 'output': output.name}))
    return 0 if overall in ('completed', 'unavailable') else 1


if __name__ == '__main__':
    raise SystemExit(main())
