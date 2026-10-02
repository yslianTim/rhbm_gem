"""Run a reproducible matrix of Joint Operator-PCG benchmark configurations."""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import subprocess
import sys
from pathlib import Path

from experiment_io import write


def build_configurations(args):
    configurations = []
    for topology in args.topologies:
        for atoms in args.atoms:
            for preconditioner in args.preconditioners:
                settings = [(None, None)] if preconditioner != 'schwarz' else [
                    (core, overlap) for core in args.cores for overlap in args.overlaps]
                for core, overlap in settings:
                    config = {
                        'topology': topology,
                        'atoms': atoms,
                        'preconditioner': preconditioner,
                        'core_atoms': core,
                        'overlap_hops': overlap,
                        'max_block_atoms': args.max_block_atoms if preconditioner == 'schwarz' else None,
                        'storage_mib': args.storage_mib if preconditioner == 'schwarz' else None,
                        'scratch_mib': args.scratch_mib if preconditioner == 'schwarz' else None,
                        'operator_rank': args.operator_rank,
                        'repeat': args.repeat,
                        'warmup': args.warmup,
                    }
                    configurations.append(config)
    keys = [configuration_key(config) for config in configurations]
    if len(keys) != len(set(keys)):
        raise ValueError('sweep dimensions contain duplicate configurations')
    return configurations


def configuration_key(config):
    fields = (config['topology'], config['atoms'], config['preconditioner'],
              config['core_atoms'], config['overlap_hops'], config['max_block_atoms'],
              config['operator_rank'], f"r{config.get('repeat', 1)}", f"w{config.get('warmup', 0)}")
    return '-'.join(str(value) for value in fields if value is not None)


def benchmark_command(script, config, args, output):
    command = [sys.executable, str(script), '--profile', 'solve',
               '--case', f"{config['topology']}-{config['atoms']}",
               '--build-dir', str(args.build_dir), '--output', str(output),
               '--repeat', str(args.repeat), '--warmup', str(args.warmup), '--timeout', str(args.timeout),
               '--rss-limit', str(args.rss_limit), '--preconditioner', config['preconditioner'],
               '--operator-rank', config['operator_rank']]
    if config['preconditioner'] == 'schwarz':
        command.extend(['--schwarz-core-atoms', str(config['core_atoms']),
                        '--schwarz-overlap-hops', str(config['overlap_hops']),
                        '--schwarz-max-block-atoms', str(config['max_block_atoms']),
                        '--schwarz-storage-mib', str(config['storage_mib']),
                        '--schwarz-scratch-mib', str(config['scratch_mib'])])
    return command


TIMING_FIELDS = (
    'partition_seconds', 'operator_prepare_seconds', 'metric_seconds', 'local_seconds',
    'factor_seconds', 'inverse_seconds', 'pcg_seconds', 'operator_normal_seconds',
    'operator_apply_seconds', 'operator_adjoint_seconds', 'rank_seconds',
    'operator_seconds', 'search_seconds', 'assessment_seconds', 'topology_setup_seconds',
    'linearization_setup_seconds', 'damping_setup_seconds', 'setup_seconds', 'iterative_seconds',
    'measured_seconds', 'wall_seconds', 'setup_fraction_of_search', 'pcg_fraction_of_search',
)
COUNTER_FIELDS = (
    'pcg_solves', 'pcg_iterations', 'operator_normals', 'operator_applications',
    'operator_adjoints', 'linearizations', 'damping_trials', 'accepted_updates',
    'profile_evaluations', 'local_builds', 'factor_builds', 'inverse_actions',
)
SUMMARY_FIELDS = (*TIMING_FIELDS, *COUNTER_FIELDS, 'peak_rss_bytes')


def distribution(values):
    values = [value for value in values if isinstance(value, (int, float))]
    if not values:
        return {'min': None, 'median': None, 'max': None}
    return {'min': min(values), 'median': statistics.median(values), 'max': max(values)}


def completed_measurement_runs(execution):
    return [item for item in execution.get('runs', [])
            if item.get('kind') == 'measurement' and item.get('status') == 'completed']


def summarize(config, report, individual_name, process_error=None):
    execution = report.get('execution', {}) if report else {}
    successful = completed_measurement_runs(execution)
    run = successful[-1] if successful else {}
    details = (run.get('result') or {}).get('details', {})
    work = details.get('search_work') or {}
    partition = details.get('partition') or {}
    policy = (report.get('metadata') or {}).get('solver_policy', {}) if report else {}
    numerics = report.get('numerics', {}) if report else {}
    status = execution.get('status', 'process_error')
    per_run_values = []
    per_run_iterations = []
    stop_reasons = []
    for item in successful:
        item_result = item.get('result') or {}
        item_details = item_result.get('details') or {}
        item_work = item_details.get('search_work') or {}
        search_seconds = item_details.get('search_seconds')
        assessment_seconds = item_details.get('assessment_seconds')
        pcg_seconds = item_work.get('pcg_seconds')
        setup_parts = [item_work.get(key) for key in (
            'partition_seconds', 'operator_prepare_seconds', 'metric_seconds',
            'local_seconds', 'factor_seconds')]
        topology_setup_seconds = item_work.get('partition_seconds')
        linearization_setup_seconds = (sum(setup_parts[1:4])
                                       if all(isinstance(value, (int, float)) for value in setup_parts[1:4])
                                       else None)
        damping_setup_seconds = item_work.get('factor_seconds')
        setup_seconds = (sum(setup_parts)
                         if all(isinstance(value, (int, float)) for value in setup_parts) else None)
        values = {key: item_work.get(key) for key in (*TIMING_FIELDS, *COUNTER_FIELDS)}
        operator_times = [item_work.get(key) for key in (
            'operator_prepare_seconds', 'operator_apply_seconds',
            'operator_adjoint_seconds', 'operator_normal_seconds')]
        values.update({
            'rank_seconds': item_work.get('operator_rank_seconds'),
            'search_seconds': search_seconds,
            'assessment_seconds': assessment_seconds,
            'topology_setup_seconds': topology_setup_seconds,
            'linearization_setup_seconds': linearization_setup_seconds,
            'damping_setup_seconds': damping_setup_seconds,
            'setup_seconds': setup_seconds,
            'iterative_seconds': pcg_seconds,
            'measured_seconds': (search_seconds + assessment_seconds
                                 if isinstance(search_seconds, (int, float)) and
                                 isinstance(assessment_seconds, (int, float)) else None),
            'wall_seconds': item.get('process_wall_seconds'),
            'peak_rss_bytes': item.get('peak_rss_bytes'),
            'setup_fraction_of_search': (setup_seconds / search_seconds
                                         if setup_seconds is not None and search_seconds else None),
            'pcg_fraction_of_search': (pcg_seconds / search_seconds
                                       if isinstance(pcg_seconds, (int, float)) and search_seconds else None),
            'operator_seconds': (sum(operator_times) if all(isinstance(value, (int, float))
                                                           for value in operator_times) else None),
        })
        values['accepted_updates'] = item_details.get('accepted_updates')
        values['profile_evaluations'] = item_details.get('profile_evaluations')
        per_run_values.append(values)
        stop_reasons.append(item_details.get('stop_reason'))
        counts = item_work.get('pcg_iteration_counts')
        if isinstance(counts, list):
            per_run_iterations.append(counts)

    measurements = {
        field: distribution([values.get(field) for values in per_run_values])
        for field in SUMMARY_FIELDS
    }
    repetition_medians = [statistics.median(counts) for counts in per_run_iterations if counts]
    pooled_iterations = [count for counts in per_run_iterations for count in counts]
    measurements['pcg_iterations_per_solve'] = distribution(repetition_medians)
    measurements['pcg_iterations_per_solve_pooled'] = distribution(pooled_iterations)
    pooled_summary = measurements['pcg_iterations_per_solve_pooled']
    pooled_mean = statistics.mean(pooled_iterations) if pooled_iterations else None
    row = {
        'configuration': config,
        'status': status,
        'reason': process_error or run.get('reason'),
        'individual_json': individual_name,
        'requested_repetitions': config.get('repeat', 1),
        'completed_repetitions': len(successful),
        'failed_repetitions': max(0, config.get('repeat', 1) - len(successful)),
        'measurements_complete': (status == 'completed' and len(successful) == config.get('repeat', 1)),
        'measurements': measurements,
        'pcg_iteration_counts': pooled_iterations,
        'pcg_iteration_counts_by_repetition': per_run_iterations,
        'pcg_iteration_repetitions': len(per_run_iterations),
        'pcg_iteration_valid_repetitions': len(repetition_medians),
        'pcg_iteration_telemetry_complete': (
            len(per_run_iterations) == len(successful) and len(repetition_medians) == len(successful)),
        'stop_reasons_by_repetition': stop_reasons,
        'pcg_iteration_budget_repetitions': sum(reason == 'pcg-iteration-budget' for reason in stop_reasons),
        'atoms': (report.get('problem') or {}).get('atoms') if report else config['atoms'],
        'rows': (report.get('problem') or {}).get('voxels') if report else None,
        'free_columns': work.get('operator_rank_columns'),
        'rank_backend': policy.get('resolved_rank_backend'),
        'sparse_backend': policy.get('sparse_backend'),
        'rank_status': work.get('operator_rank_status'),
        'rank_reason': work.get('operator_rank_reason'),
        'rank_seconds': measurements['rank_seconds']['median'],
        'rank_entries': work.get('operator_rank_entries'),
        'rank_workspace_bytes': work.get('operator_rank_workspace_bytes'),
        'compact_extractions': work.get('operator_rank_compact_extractions'),
        'free_design_svds': work.get('operator_rank_free_design_svds'),
        'blocks': partition.get('blocks'),
        'maximum_block_atoms': work.get('maximum_block_atoms'),
        'topology_bytes': work.get('topology_bytes'),
        'storage_bytes': work.get('storage_bytes'),
        'scratch_bytes_bound': work.get('scratch_bytes_bound'),
        'core_atom_counts': partition.get('core_atom_counts'),
        'block_atom_counts': partition.get('block_atom_counts'),
        'atom_membership': partition.get('atom_membership'),
        'atom_membership_histogram': partition.get('atom_membership_histogram'),
        'pcg_solves': measurements['pcg_solves']['median'],
        'pcg_iterations': measurements['pcg_iterations']['median'],
        'pcg_iterations_min': pooled_summary['min'],
        'pcg_iterations_median': pooled_summary['median'],
        'pcg_iterations_max': pooled_summary['max'],
        'pcg_iterations_mean': pooled_mean,
        'pcg_iterations_per_solve': measurements['pcg_iterations_per_solve']['median'],
        'pcg_iterations_per_solve_pooled_median': measurements['pcg_iterations_per_solve_pooled']['median'],
        'damping_trials': measurements['damping_trials']['median'],
        'local_builds': measurements['local_builds']['median'],
        'factor_builds': measurements['factor_builds']['median'],
        'inverse_actions': measurements['inverse_actions']['median'],
        'iterations_per_solve': measurements['pcg_iterations_per_solve']['median'],
        'operator_normals': measurements['operator_normals']['median'],
        'operator_applications': measurements['operator_applications']['median'],
        'operator_adjoints': measurements['operator_adjoints']['median'],
        'partition_seconds': measurements['partition_seconds']['median'],
        'local_seconds': measurements['local_seconds']['median'],
        'factor_seconds': measurements['factor_seconds']['median'],
        'inverse_seconds': measurements['inverse_seconds']['median'],
        'operator_seconds': measurements['operator_seconds']['median'],
        'operator_prepare_seconds': measurements['operator_prepare_seconds']['median'],
        'pcg_seconds': measurements['pcg_seconds']['median'],
        'search_seconds': measurements['search_seconds']['median'],
        'assessment_seconds': measurements['assessment_seconds']['median'],
        'topology_setup_seconds': measurements['topology_setup_seconds']['median'],
        'linearization_setup_seconds': measurements['linearization_setup_seconds']['median'],
        'damping_setup_seconds': measurements['damping_setup_seconds']['median'],
        'setup_seconds': measurements['setup_seconds']['median'],
        'iterative_seconds': measurements['iterative_seconds']['median'],
        'measured_seconds': measurements['measured_seconds']['median'],
        'setup_fraction_of_search': measurements['setup_fraction_of_search']['median'],
        'pcg_fraction_of_search': measurements['pcg_fraction_of_search']['median'],
        'wall_seconds': measurements['wall_seconds']['median'],
        'peak_rss_bytes': measurements['peak_rss_bytes']['max'],
        'objective': numerics.get('objective'),
        'runtime_convergence': numerics.get('runtime_convergence'),
        'state_available': details.get('state_available'),
        'stop_reason': details.get('stop_reason'),
    }
    return row


CSV_BASE_FIELDS = (
    'configuration', 'status', 'reason', 'individual_json', 'atoms', 'rows', 'free_columns',
    'sparse_backend', 'rank_backend', 'rank_status', 'rank_reason', 'rank_seconds', 'rank_entries',
    'rank_workspace_bytes', 'compact_extractions', 'free_design_svds', 'blocks',
    'maximum_block_atoms', 'topology_bytes', 'storage_bytes', 'scratch_bytes_bound',
    'requested_repetitions', 'completed_repetitions', 'failed_repetitions', 'measurements_complete',
    'pcg_iteration_repetitions', 'pcg_iteration_valid_repetitions',
    'pcg_iteration_telemetry_complete', 'pcg_iteration_budget_repetitions',
    'pcg_solves', 'pcg_iterations', 'pcg_iterations_per_solve',
    'pcg_iterations_per_solve_pooled_median', 'operator_normals', 'operator_applications',
    'operator_adjoints', 'damping_trials', 'local_builds', 'factor_builds',
    'inverse_actions', 'pcg_iteration_counts', 'pcg_iterations_min',
    'pcg_iterations_median', 'pcg_iterations_max', 'pcg_iterations_mean',
    'iterations_per_solve',
    'partition_seconds', 'local_seconds', 'factor_seconds', 'inverse_seconds',
    'operator_seconds', 'operator_prepare_seconds', 'pcg_seconds', 'search_seconds',
    'assessment_seconds', 'topology_setup_seconds', 'linearization_setup_seconds',
    'damping_setup_seconds', 'setup_seconds', 'iterative_seconds', 'measured_seconds',
    'setup_fraction_of_search', 'pcg_fraction_of_search', 'wall_seconds', 'peak_rss_bytes', 'objective',
    'runtime_convergence', 'state_available', 'stop_reason',
)
CSV_MEASUREMENT_FIELDS = tuple(
    f'{field}_{stat}' for field in SUMMARY_FIELDS for stat in ('min', 'median', 'max')) + tuple(
    f'{field}_{stat}' for field in ('pcg_iterations_per_solve', 'pcg_iterations_per_solve_pooled')
    for stat in ('min', 'median', 'max'))
CSV_FIELDS = (*CSV_BASE_FIELDS, *CSV_MEASUREMENT_FIELDS)


def write_csv(path, rows):
    with path.open('w', newline='') as output:
        writer = csv.DictWriter(output, fieldnames=CSV_FIELDS, extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            flattened = dict(row)
            flattened['configuration'] = json.dumps(row['configuration'], sort_keys=True)
            for field, summary in row.get('measurements', {}).items():
                for stat, value in summary.items():
                    flattened[f'{field}_{stat}'] = value
            flattened['pcg_iteration_counts'] = json.dumps(row.get('pcg_iteration_counts', []))
            writer.writerow(flattened)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='Aggregate JSON output path')
    parser.add_argument('--csv', type=Path, help='Optional CSV summary path')
    parser.add_argument('--topologies', nargs='+', choices=('chain', 'cube'), default=['chain'])
    parser.add_argument('--atoms', nargs='+', type=int, default=[64, 128])
    parser.add_argument('--cores', nargs='+', type=int, default=[32, 64])
    parser.add_argument('--overlaps', nargs='+', type=int, default=[0, 1])
    parser.add_argument('--preconditioners', nargs='+', choices=('identity', 'diagonal', 'schwarz'),
                        default=['schwarz'])
    parser.add_argument('--operator-rank', choices=('auto', 'dense', 'spqr-bounds'), default='auto')
    parser.add_argument('--max-block-atoms', type=int, default=512)
    parser.add_argument('--storage-mib', type=int, default=512)
    parser.add_argument('--scratch-mib', type=int, default=256)
    parser.add_argument('--repeat', type=int, default=1)
    parser.add_argument('--warmup', type=int, default=0)
    parser.add_argument('--timeout', type=float, default=600)
    parser.add_argument('--rss-limit', type=int, default=4 * 1024**3)
    return parser


def validate_args(parser, args):
    if args.repeat < 1 or args.warmup < 0 or args.timeout <= 0 or args.rss_limit <= 0:
        parser.error('--repeat and limits must be positive; --warmup must be nonnegative')
    if (not args.atoms or any(value <= 0 for value in args.atoms) or
            not args.cores or any(value <= 0 for value in args.cores) or
            not args.overlaps or any(value < 0 for value in args.overlaps)):
        parser.error('atoms and cores must be positive; overlaps must be nonnegative')
    if (args.max_block_atoms <= 0 or args.storage_mib <= 0 or args.scratch_mib <= 0):
        parser.error('Schwarz block and memory limits must be positive')
    for core in args.cores:
        if core > args.max_block_atoms:
            parser.error('max block atoms must cover every requested core')
    if (len(args.topologies) != len(set(args.topologies)) or
            len(args.atoms) != len(set(args.atoms)) or
            len(args.cores) != len(set(args.cores)) or
            len(args.overlaps) != len(set(args.overlaps)) or
            len(args.preconditioners) != len(set(args.preconditioners))):
        parser.error('sweep dimensions must not contain duplicate values')


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)
    try:
        configurations = build_configurations(args)
    except ValueError as error:
        parser.error(str(error))

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    run_directory = output.parent / f'{output.stem}_runs'
    run_directory.mkdir(exist_ok=True)
    script = Path(__file__).with_name('joint_benchmark.py').resolve()
    rows = []
    for config in configurations:
        name = configuration_key(config) + '.json'
        individual = run_directory / name
        command = benchmark_command(script, config, args, individual)
        error = None
        try:
            process = subprocess.run(command, check=False, capture_output=True, text=True,
                                     timeout=max(120, args.timeout * (args.repeat + 1)))
            if process.returncode:
                error = f'joint_benchmark exited {process.returncode}'
        except subprocess.TimeoutExpired:
            error = 'joint_benchmark process timed out'
        try:
            report = json.loads(individual.read_text()) if individual.is_file() else None
        except (OSError, json.JSONDecodeError) as failure:
            report = None
            error = error or f'cannot read individual result: {failure}'
        row = summarize(config, report, str(Path(run_directory.name) / name), error)
        if error and row['status'] == 'completed':
            row['status'] = 'timeout' if 'timed out' in error else 'process_error'
        elif error and row['status'] == 'process_error' and 'timed out' in error:
            row['status'] = 'timeout'
        rows.append(row)

    aggregate = {
        'schema_version': 1,
        'tool': 'joint_schwarz_sweep',
        'measurement_profile': 'solve',
        'configuration_count': len(configurations),
        'configurations': rows,
    }
    write(output, aggregate)
    if args.csv:
        csv_path = args.csv.resolve()
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        write_csv(csv_path, rows)
    print(json.dumps({'status': 'completed', 'configurations': len(rows), 'output': output.name}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
