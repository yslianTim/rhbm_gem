"""Run the current FixedNeighbor Joint benchmark profiles."""
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
from experiment_process import monitored
from experiment_provenance import source_hash
from joint_runtime_support import read


SCOPES = {
    'search': 'FixedNeighbor search-only with local LegacyCompact profile work',
    'solve': 'FixedNeighbor search, replay, endpoint certification, and RuntimeConvergence',
    'workflow': 'Joint workflow, postprocess, and persistence',
    'postprocess': 'Joint postprocess, uncertainty, and persistence',
    'command': 'CLI analysis, SQLite reload, and Joint export',
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
        'backend': 'EIGEN',
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


def synthetic_case(case):
    match = re.fullmatch(r'(chain|cube)-(\d+)', case)
    if not match:
        return None
    return match.group(1), int(match.group(2))


def case_input(case):
    synthetic = synthetic_case(case)
    if not synthetic:
        raise ValueError('Joint benchmark cases use chain-N or cube-N.')
    topology, atoms = synthetic
    return {
        'kind': 'synthetic',
        'topology': topology,
        'atoms': atoms,
        'input_sha256': digest({
            'generator': sha(ROOT / 'tests/support/JointSyntheticWorkload.cpp'),
            'topology': topology,
            'atoms': atoms,
        }),
    }


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


def fixed_neighbor_command(args, output, build):
    driver = build / 'bin/joint_fixed_neighbor_experiment'
    if not driver.is_file():
        raise ValueError('Build joint_fixed_neighbor_experiment with RHBM_GEM_BUILD_BENCHMARKS=ON')
    topology, atoms = synthetic_case(args.case)
    if args.profile == 'search':
        return [str(driver), '--scaling-only', str(output), topology, str(atoms)]
    return [str(driver), '--case', str(output), topology, str(atoms)]


def command_for_profile(args, output, build):
    if args.profile in ('search', 'solve'):
        return fixed_neighbor_command(args, output, build)
    if args.profile in ('workflow', 'postprocess'):
        driver = build / 'bin/joint_postprocessing_benchmark'
        if not driver.is_file():
            raise ValueError('Build joint_postprocessing_benchmark with RHBM_GEM_BUILD_BENCHMARKS=ON')
        return [str(driver), args.case, 'workflow' if args.profile == 'workflow' else 'post', str(output)]
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


def raw_elapsed(profile, raw, process_wall):
    fixed = raw.get('fixed_neighbor') or {}
    if profile == 'search':
        return fixed.get('search_seconds', process_wall)
    if profile == 'solve':
        return fixed.get('total_elapsed_seconds', process_wall)
    if profile == 'workflow':
        return raw.get('total_seconds', process_wall)
    if profile == 'postprocess':
        return sum(raw.get('phases', {}).values())
    return process_wall


def normalize_result(profile, raw):
    fixed = raw.get('fixed_neighbor') or {}
    if profile == 'search':
        return {
            'qualified': None,
            'scientific_status': 'search-only',
            'details': {
                'measurement_scope': 'fixed-neighbor-search-only',
                'search_completed': fixed.get('search_converged'),
                'stop_reason': fixed.get('search_reason'),
                'sweeps': fixed.get('sweeps'),
                'block_solves': fixed.get('block_solves'),
                'profile_evaluations': fixed.get('profile_evaluations'),
                'accepted_blocks': fixed.get('accepted_blocks'),
                'objective': fixed.get('objective'),
                'global_ac_kkt': fixed.get('final_global_ac_kkt'),
                'width_gradient_inf_norm': fixed.get('final_raw_width_gradient_inf_norm'),
                'search_seconds': fixed.get('search_seconds'),
                'local_profile_work': fixed.get('local_profile_work'),
            },
        }
    if profile == 'solve':
        runtime = fixed.get('runtime_convergence')
        return {
            'qualified': runtime == 'Passed' if runtime in ('Passed', 'Failed') else None,
            'scientific_status': runtime or 'not-assessed',
            'details': {
                'search_completed': fixed.get('search_converged'),
                'stop_reason': fixed.get('search_reason'),
                'endpoint_certified': fixed.get('endpoint_certified'),
                'runtime_convergence': runtime,
                'objective': fixed.get('objective'),
                'global_kkt': fixed.get('global_kkt'),
                'width_gradient_inf_norm': fixed.get('width_gradient_inf_norm'),
                'search_seconds': fixed.get('search_seconds'),
                'assessment_seconds': fixed.get('assessment_seconds'),
                'sweeps': fixed.get('sweeps'),
                'block_solves': sum(item.get('block_solves', 0) for item in fixed.get('sweep_telemetry', [])),
                'profile_evaluations': sum(item.get('profile_evaluations', 0)
                                           for item in fixed.get('sweep_telemetry', [])),
                'local_profile_work': fixed.get('local_profile_work'),
            },
        }
    if profile == 'workflow':
        return {
            'qualified': None,
            'scientific_status': 'workflow-complete',
            'details': {'endpoint_count': len(raw.get('endpoints', [])),
                        'target_count': len(raw.get('targets', [])),
                        'objective': raw.get('objective')},
        }
    if profile == 'postprocess':
        return {
            'qualified': None,
            'scientific_status': 'postprocess-complete',
            'details': {'endpoint_count': len(raw.get('endpoints', [])),
                        'target_count': len(raw.get('targets', [])),
                        'objective': raw.get('objective')},
        }
    return {'qualified': None, 'scientific_status': 'command-complete', 'details': {}}


def problem_result(raw):
    fixed = raw.get('fixed_neighbor') or {}
    atoms = raw.get('atoms', fixed.get('atoms'))
    rows = raw.get('rows', fixed.get('rows'))
    return {
        'atoms': atoms,
        'voxels': rows,
        'components': raw.get('components'),
        'free_columns': raw.get('free_columns'),
        'design_nonzeros': raw.get('design_nonzeros'),
        'parameters': raw.get('parameter_count', 3 * atoms if isinstance(atoms, int) else None),
    }


def numerics_result(raw):
    fixed = raw.get('fixed_neighbor') or {}
    if fixed:
        return {
            'rank': None,
            'objective': fixed.get('objective'),
            'gradient_norm': fixed.get('width_gradient_inf_norm',
                                       fixed.get('final_raw_width_gradient_inf_norm')),
            'runtime_convergence': fixed.get('runtime_convergence'),
            'active_atoms': None,
        }
    return {
        'rank': None,
        'objective': raw.get('objective'),
        'gradient_norm': None,
        'runtime_convergence': None,
        'active_atoms': None,
    }


def execute_once(args, build, run_root, deadline):
    results = []
    if args.profile == 'command' and (not args.model or not args.map_path or
                                      not args.model.is_file() or not args.map_path.is_file()):
        return {'status': 'unavailable', 'reason': 'command profile requires --model and --map',
                'elapsed_seconds': None, 'peak_rss_bytes': None, 'raw': None, 'stages': results}

    if args.profile in ('workflow', 'postprocess'):
        driver_output = run_root / 'benchmark.sqlite'
    elif args.profile == 'command':
        driver_output = run_root / 'command-output'
        driver_output.mkdir()
    else:
        driver_output = run_root / 'driver.json'
    try:
        command = command_for_profile(args, driver_output, build)
    except ValueError as error:
        return {'status': 'unavailable', 'reason': str(error), 'elapsed_seconds': None,
                'peak_rss_bytes': None, 'raw': None, 'stages': results}
    if command is None:
        return {'status': 'unavailable', 'reason': 'command profile requires --model and --map',
                'elapsed_seconds': None, 'peak_rss_bytes': None, 'raw': None, 'stages': results}

    if args.profile == 'command':
        analysis, dump, database, export = command
        rss = 0
        for role, stage_command in (('analysis', analysis), ('reload_export', dump)):
            record, _ = run_process(stage_command, run_root / role, deadline, args.rss_limit, args.timeout)
            results.append({'role': role, **record})
            rss = max(rss, record.get('sampled_tree_peak_rss_bytes') or 0,
                      record.get('os_process_peak_rss_bytes') or 0)
            if record['status'] != 'completed':
                return {'status': record['status'], 'elapsed_seconds': None,
                        'peak_rss_bytes': rss, 'stages': results, 'raw': None}
        exported_json = export / f'joint_result_{args.case}.json'
        exported_csv = export / f'joint_atoms_{args.case}.csv'
        persistence = database.is_file()
        raw = read(exported_json) if exported_json.is_file() else None
        export_ok = exported_json.is_file() and exported_csv.is_file()
        return {'status': 'completed' if persistence and export_ok else 'process_error',
                'elapsed_seconds': sum(row.get('wall_seconds') or 0 for row in results),
                'peak_rss_bytes': rss, 'stages': results, 'persistence_result': persistence,
                'export_result': export_ok,
                'export_bytes': exported_json.stat().st_size if export_ok else None, 'raw': raw}

    directory = run_root / 'measured'
    record, _ = run_process(command, directory, deadline, args.rss_limit, args.timeout)
    results.append({'role': 'measurement', **record})
    rss = max(record.get('sampled_tree_peak_rss_bytes') or 0,
              record.get('os_process_peak_rss_bytes') or 0)
    if record['status'] != 'completed':
        return {'status': record['status'], 'reason': record.get('raw_status'),
                'elapsed_seconds': None, 'process_wall_seconds': record.get('wall_seconds'),
                'peak_rss_bytes': rss, 'stages': results, 'raw': None}
    try:
        raw = driver_json(directory, driver_output if args.profile in ('search', 'solve') else None)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        return {'status': 'process_error', 'reason': str(error), 'elapsed_seconds': None,
                'peak_rss_bytes': rss, 'stages': results, 'raw': None}
    return {'status': 'completed', 'elapsed_seconds': raw_elapsed(args.profile, raw,
                                                                  record.get('wall_seconds') or 0),
            'process_wall_seconds': record.get('wall_seconds'),
            'peak_rss_bytes': max(rss, raw.get('peak_rss_bytes') or 0),
            'stages': results, 'raw': raw}


def solver_policy_metadata():
    return {
        'search_method': 'FixedNeighbor',
        'sparse_backend': 'EIGEN',
        'fixed_neighbor_core_atoms': 12,
        'fixed_neighbor_block_order': 'Forward',
        'fixed_neighbor_maximum_sweeps': 30,
        'fixed_neighbor_local_search': 'LegacyCompact',
        'fixed_neighbor_update_policy': 'OneAccepted',
    }


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
    parser.add_argument('--cli', type=Path)
    parser.add_argument('--model', type=Path)
    parser.add_argument('--map', dest='map_path', type=Path)
    return parser


def validate_args(parser, args):
    if args.repeat < 1 or args.warmup < 0 or args.timeout <= 0 or args.rss_limit <= 0:
        parser.error('--repeat and limits must be positive; --warmup must be nonnegative')
    if args.profile in ('workflow', 'postprocess') and args.case not in ('full', 'halo', 'multi'):
        parser.error('workflow and postprocess cases are full, halo, or multi')
    if args.profile == 'command' and not re.fullmatch(r'[A-Za-z0-9._-]+', args.case):
        parser.error('command case names may contain letters, digits, dot, underscore, and hyphen')
    if args.profile in ('search', 'solve') and not synthetic_case(args.case):
        parser.error('search and solve cases use chain-N or cube-N')


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)
    output = args.output.resolve()
    build = args.build_dir.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        metadata = build_metadata(build)
        commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, check=True,
                                capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        parser.error(str(error))
    metadata.update(commit=commit, profile=args.profile, case=args.case,
                    source_sha256=source_hash(ROOT), benchmark_sha256=sha(Path(__file__)))
    metadata['solver_policy'] = solver_policy_metadata()
    if args.profile in ('search', 'solve'):
        driver = build / 'bin/joint_fixed_neighbor_experiment'
    elif args.profile in ('workflow', 'postprocess'):
        driver = build / 'bin/joint_postprocessing_benchmark'
    else:
        driver = args.cli or build / 'bin/RHBM-GEM'
    metadata['driver_sha256'] = sha(driver) if driver.is_file() else None

    with tempfile.TemporaryDirectory(prefix='joint-benchmark-', dir=output.parent) as temporary:
        work = Path(temporary)
        case = case_input(args.case) if args.profile in ('search', 'solve') else {
            'kind': 'synthetic',
            'input_sha256': digest({'case': args.case, 'driver': metadata['driver_sha256']}),
        }
        metadata['input_sha256'] = case['input_sha256']
        records = []
        overall = 'completed'
        for kind, count in (('warmup', args.warmup), ('measurement', args.repeat)):
            for index in range(count):
                run_root = work / f'{kind}-{index + 1}'
                run_root.mkdir(parents=True)
                result = execute_once(args, build, run_root, time.monotonic() + args.timeout)
                record = {'kind': kind, 'index': index + 1, **result}
                records.append(record)
                if result['status'] not in ('completed', 'unavailable'):
                    overall = result['status']
                elif result['status'] == 'unavailable' and overall == 'completed':
                    overall = 'unavailable'

        measurements = [record for record in records if record['kind'] == 'measurement']
        successful = [record for record in measurements if record['status'] == 'completed']
        all_measured = len(successful) == args.repeat
        elapsed = [record['elapsed_seconds'] for record in successful if record.get('elapsed_seconds') is not None]
        peaks = [record['peak_rss_bytes'] for record in successful if record.get('peak_rss_bytes') is not None]
        last_raw = next((record['raw'] for record in reversed(successful) if record.get('raw') is not None), None)
        run_summaries = []
        for record in records:
            raw = record.get('raw') or {}
            run_summaries.append({
                key: record.get(key) for key in (
                    'kind', 'index', 'status', 'reason', 'elapsed_seconds',
                    'process_wall_seconds', 'peak_rss_bytes', 'persistence_result',
                    'export_result', 'export_bytes')
            } | {
                'problem': problem_result(raw),
                'numerics': numerics_result(raw),
                'result': normalize_result(args.profile, raw),
                'stages': [{key: stage.get(key) for key in (
                    'role', 'status', 'exit_code', 'wall_seconds',
                    'sampled_tree_peak_rss_bytes', 'os_process_peak_rss_bytes')}
                           for stage in record.get('stages', [])],
            })
        report = {
            'schema_version': 1,
            'profile': args.profile,
            'measurement_scope': SCOPES[args.profile],
            'case': args.case,
            'metadata': metadata,
            'execution': {
                'status': overall,
                'returncode': 0 if all_measured and all(r['status'] == 'completed' for r in records) else None,
                'completed_repetitions': len(successful),
                'failed_repetitions': args.repeat - len(successful),
                'runs': run_summaries,
            },
            'problem': problem_result(last_raw or {}),
            'numerics': numerics_result(last_raw or {}),
            'resources': {
                'elapsed_seconds': statistics.median(elapsed) if all_measured and elapsed else None,
                'peak_rss_bytes': max((r.get('peak_rss_bytes') or 0 for r in records), default=0) or None,
                'completed_measurement_peak_rss_bytes': max(peaks) if all_measured and peaks else None,
                'count': args.repeat,
                'failed_count': args.repeat - len(successful),
                'samples': [r.get('elapsed_seconds') if r['status'] == 'completed' else None
                            for r in measurements],
            },
            'result': normalize_result(args.profile, last_raw or {}),
        }
        if args.profile == 'command' and records:
            final = next((r for r in reversed(records) if r['kind'] == 'measurement'), {})
            report['result'].update(persistence=final.get('persistence_result'),
                                    export=final.get('export_result'),
                                    export_bytes=final.get('export_bytes'))
        write(output, report)
    print(json.dumps({'profile': args.profile, 'status': overall, 'output': output.name}))
    return 0 if overall in ('completed', 'unavailable') else 1


if __name__ == '__main__':
    raise SystemExit(main())
