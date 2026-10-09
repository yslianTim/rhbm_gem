"""Matched LegacyCompact, OperatorPcg, and FixedNeighbor route frontier campaign."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from experiment_io import ROOT, read, sha, write
from experiment_provenance import source_hash

ROUTES = {
    'LegacyCompact': 'legacy',
    'OperatorPcg': 'schwarz',
    'FixedNeighbor': 'fixed-neighbor',
}
FIXED_NEIGHBOR_POLICY = {
    'core_atoms': 12,
    'block_order': 'forward',
    'local_search': 'LegacyCompact',
    'maximum_sweeps': 30,
}
SCOPES = {'search-only': 'search', 'full-endpoint': 'solve'}
DEFAULT_CASES = tuple(
    f'{topology}-{atoms}'
    for atoms in (512, 768, 1024)
    for topology in ('chain', 'cube')
)


def _measurement(report):
    runs = (report.get('execution') or {}).get('runs') or []
    return next((run for run in reversed(runs) if run.get('kind') == 'measurement'), {})


def _complete(report, scope):
    run = _measurement(report)
    if run.get('status') != 'completed':
        return False
    details = (run.get('result') or {}).get('details') or {}
    if details.get('search_completed') is not True:
        return False
    if scope == 'full-endpoint':
        if details.get('endpoint_certified') is not None:
            return details.get('endpoint_certified') is True and details.get('runtime_convergence') == 'Passed'
        return (report.get('result') or {}).get('qualified') is True
    return True


def result_summary(report, scope, timeout_seconds):
    run = _measurement(report)
    details = (run.get('result') or {}).get('details') or {}
    work = details.get('search_work') or {}
    numerics = report.get('numerics') or {}
    resources = report.get('resources') or {}
    status = run.get('status') or (report.get('execution') or {}).get('status')
    completed = _complete(report, scope)
    elapsed = resources.get('elapsed_seconds')
    observed = run.get('process_wall_seconds')
    search_seconds = details.get('search_seconds')
    if search_seconds is None:
        search_seconds = work.get('search_seconds')
    objective = (details.get('objective') if details.get('objective') is not None else
                 details.get('accepted_objective', numerics.get('objective')))
    failure_stage = None
    if not completed:
        failure_stage = run.get('failure_stage') or run.get('driver_stage')
        if failure_stage == 'complete':
            failure_stage = 'assessment' if details.get('search_completed') is True else 'search'
    failure_reason = run.get('reason') if status != 'completed' else (
        details.get('stop_reason') if not completed else None)
    assessment_seconds = details.get('assessment_seconds')
    if assessment_seconds is None and scope == 'full-endpoint' and not completed and observed is not None:
        assessment_seconds = max(0., observed - (search_seconds or 0.))
    return {
        'status': status,
        'completed': completed,
        'timeout_seconds': timeout_seconds,
        'search_seconds': search_seconds,
        'assessment_seconds': assessment_seconds,
        'total_seconds': elapsed if elapsed is not None else observed,
        'observed_process_seconds': observed,
        'peak_rss_mb': ((run.get('peak_rss_bytes') or 0) / 1024**2
                        if run.get('peak_rss_bytes') is not None else None),
        'search_completed': details.get('search_completed'),
        'failure_stage': failure_stage,
        'failure_reason': failure_reason,
        'objective': objective,
        'final_gradient_norm': (details.get('accepted_gradient_inf_norm')
                                if details.get('accepted_gradient_inf_norm') is not None
                                else numerics.get('gradient_norm')),
        'final_global_ac_kkt': details.get('global_ac_kkt'),
        'final_width_gradient_inf_norm': details.get('width_gradient_inf_norm'),
        'endpoint_certified': (details.get('endpoint_certified')
                               if details.get('endpoint_certified') is not None
                               else (report.get('result') or {}).get('qualified')),
        'runtime_convergence': details.get('runtime_convergence', numerics.get('runtime_convergence')),
        'endpoint_inner': details.get('endpoint_inner'),
        'endpoint_gradient': details.get('endpoint_gradient'),
        'endpoint_local': details.get('endpoint_local'),
        'endpoint_identified': details.get('endpoint_identified'),
        'endpoint_trust': details.get('endpoint_trust'),
        'factorization_count': (work.get('factor_builds') if work.get('factor_builds') is not None
                                else work.get('fixed_factor', {}).get('calls')),
        'factor_seconds': (details.get('factor_seconds') if details.get('factor_seconds') is not None
                           else work.get('factor_seconds')),
        'profile_evaluations': (details.get('profile_evaluations')
                                if details.get('profile_evaluations') is not None
                                else work.get('profile_evaluations')),
        'pcg_solves': work.get('pcg_solves'),
        'pcg_iterations': work.get('pcg_iterations'),
        'block_solves': details.get('block_solves'),
        'sweeps': details.get('sweeps'),
    }


def record_from_report(report, topology, atoms, route, scope, timeout_seconds, report_path):
    problem = report.get('problem') or {}
    return {
        'topology': topology,
        'atoms': atoms,
        'rows': problem.get('voxels'),
        'parameter_count': problem.get('parameters', 3 * atoms),
        'route': route,
        'measurement_scope': scope,
        'input_sha256': (report.get('metadata') or {}).get('input_sha256'),
        'report_path': str(report_path.relative_to(ROOT)) if report_path.is_relative_to(ROOT) else str(report_path),
        'result': result_summary(report, scope, timeout_seconds),
    }


def compare_group(records):
    complete = [record for record in records if (record.get('diagnostic_result') or {}).get('completed')]
    time_key = 'search_seconds' if records and records[0]['measurement_scope'] == 'search-only' else 'total_seconds'
    timed = [record for record in complete if (record['diagnostic_result'].get(time_key) is not None)]
    fastest = min(timed, key=lambda record: record['diagnostic_result'][time_key])['route'] if timed else None
    rss_rows = []
    for record in records:
        diagnostic = record.get('diagnostic_result') or {}
        formal = record.get('formal_result') or {}
        if diagnostic.get('peak_rss_mb') is not None:
            rss_rows.append((record, diagnostic['peak_rss_mb']))
        elif formal.get('peak_rss_mb') is not None:
            rss_rows.append((record, formal['peak_rss_mb']))
    least_rss = min(rss_rows, key=lambda item: item[1])[0]['route'] if rss_rows else None
    formal_passes = [record['route'] for record in records
                     if (record.get('formal_result') or {}).get('completed')]
    input_hashes = {record.get('input_sha256') for record in records if record.get('input_sha256')}
    row_counts = {record.get('rows') for record in records if record.get('rows') is not None}
    parameter_counts = {record.get('parameter_count') for record in records
                        if record.get('parameter_count') is not None}
    return {
        'topology': records[0]['topology'] if records else None,
        'atoms': records[0]['atoms'] if records else None,
        'measurement_scope': records[0]['measurement_scope'] if records else None,
        'fastest_route': fastest,
        'fastest_metric': time_key if fastest else None,
        'least_rss_route': least_rss,
        'formal_envelope_pass_routes': formal_passes,
        'matched_input': len(input_hashes) <= 1 and len(row_counts) <= 1 and len(parameter_counts) <= 1,
        'rows': next(iter(row_counts)) if len(row_counts) == 1 else None,
        'parameter_count': next(iter(parameter_counts)) if len(parameter_counts) == 1 else None,
        'routes': [{
            'route': record['route'],
            'formal_envelope': record.get('formal_result'),
            'diagnostic': record.get('diagnostic_result'),
        } for record in records],
    }


def classify_position(groups):
    by_topology = {}
    resource_advantage = False
    wins = {}
    for group in groups:
        winner = group.get('fastest_route')
        if winner:
            wins[winner] = wins.get(winner, 0) + 1
            by_topology.setdefault(group['topology'], []).append(winner)
        formal = {route['route']: route.get('formal_envelope') or {}
                  for route in group.get('routes', [])}
        fixed_pass = formal.get('FixedNeighbor', {}).get('completed') is True
        memory_failed = any(formal.get(route, {}).get('status') == 'rss_limit'
                          for route in ('LegacyCompact', 'OperatorPcg'))
        resource_advantage |= fixed_pass and memory_failed
    topology_winners = {}
    for topology, items in by_topology.items():
        counts = {route: items.count(route) for route in items}
        topology_winners[topology] = max(counts, key=counts.get)
    if resource_advantage:
        return 'FixedNeighbor-memory-advantaged'
    if len(set(topology_winners.values())) > 1:
        return 'mixed-by-topology'
    if wins.get('FixedNeighbor', 0) > len([g for g in groups if g.get('fastest_route')]) / 2:
        return 'FixedNeighbor-time-dominant'
    if wins.get('OperatorPcg', 0) >= wins.get('LegacyCompact', 0) and wins.get('OperatorPcg', 0):
        return 'OperatorPcg-dominant'
    if wins:
        return 'mixed-by-topology' if len(set(topology_winners.values())) > 1 else 'no-measured-benefit'
    return 'no-measured-benefit'


def analyze(records):
    grouped = {}
    for record in records:
        key = (record['measurement_scope'], record['topology'], record['atoms'])
        grouped.setdefault(key, []).append(record)
    groups = []
    for key in sorted(grouped, key=lambda item: (item[0], item[2], item[1])):
        groups.append(compare_group(sorted(grouped[key], key=lambda row: row['route'])))
    first_failures = {}
    for route in ROUTES:
        for topology in ('chain', 'cube'):
            rows = sorted((row for row in records if row['route'] == route and
                           row['topology'] == topology and row['measurement_scope'] == 'search-only'),
                          key=lambda item: item['atoms'])
            first = next((row for row in rows if not (row.get('formal_result') or {}).get('completed')), None)
            first_failures[f'{topology}:{route}'] = first['atoms'] if first else None
    return {
        'route_position': classify_position(groups),
        'groups': groups,
        'first_formal_search_failure_atoms': first_failures,
    }


def write_outputs(output_dir, manifest, records):
    output_dir.mkdir(parents=True, exist_ok=True)
    write(output_dir / 'campaign-manifest.json', manifest)
    analysis = analyze(records)
    write(output_dir / 'analysis.json', analysis)
    columns = [
        'topology', 'atoms', 'rows', 'parameter_count', 'route', 'measurement_scope',
        'formal_status', 'formal_pass', 'formal_search_seconds', 'formal_total_seconds',
        'formal_peak_rss_mb', 'formal_failure_stage', 'formal_failure_reason',
        'diagnostic_status', 'diagnostic_completed', 'diagnostic_search_seconds',
        'diagnostic_assessment_seconds', 'diagnostic_total_seconds', 'diagnostic_peak_rss_mb',
        'diagnostic_reason',
        'objective', 'factorization_count', 'factor_seconds', 'profile_evaluations',
        'pcg_solves', 'pcg_iterations', 'block_solves', 'sweeps',
    ]
    with (output_dir / 'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator='\n')
        writer.writeheader()
        for record in records:
            formal = record.get('formal_result') or {}
            diagnostic = record.get('diagnostic_result') or {}
            writer.writerow({
                'topology': record['topology'], 'atoms': record['atoms'], 'rows': record.get('rows'),
                'parameter_count': record.get('parameter_count'), 'route': record['route'],
                'measurement_scope': record['measurement_scope'],
                'formal_status': formal.get('status'), 'formal_pass': formal.get('completed'),
                'formal_search_seconds': formal.get('search_seconds'),
                'formal_total_seconds': formal.get('total_seconds'),
                'formal_peak_rss_mb': formal.get('peak_rss_mb'),
                'formal_failure_stage': formal.get('failure_stage'),
                'formal_failure_reason': formal.get('failure_reason'),
                'diagnostic_status': diagnostic.get('status'),
                'diagnostic_completed': diagnostic.get('completed'),
                'diagnostic_search_seconds': diagnostic.get('search_seconds'),
                'diagnostic_assessment_seconds': diagnostic.get('assessment_seconds'),
                'diagnostic_total_seconds': diagnostic.get('total_seconds'),
                'diagnostic_peak_rss_mb': diagnostic.get('peak_rss_mb'),
                'diagnostic_reason': diagnostic.get('reason'),
                'objective': diagnostic.get('objective'),
                'factorization_count': diagnostic.get('factorization_count'),
                'factor_seconds': diagnostic.get('factor_seconds'),
                'profile_evaluations': diagnostic.get('profile_evaluations'),
                'pcg_solves': diagnostic.get('pcg_solves'),
                'pcg_iterations': diagnostic.get('pcg_iterations'),
                'block_solves': diagnostic.get('block_solves'),
                'sweeps': diagnostic.get('sweeps'),
            })
    readme = (
        '# Matched Joint route frontier\n\n'
        'This campaign measures LegacyCompact, OperatorPcg, and FixedNeighbor in separate processes '
        'with the same frozen synthetic workload, initial widths, SPQR backend, one Eigen thread, '
        'resource monitor, and 4 GiB RSS limit. Search-only is the primary comparison. Full endpoint '
        f'runs are limited to 512 and 768 atoms. Formal runs use a {manifest["formal_envelope"]["wall_seconds"]}-second cap; '
        f'diagnostic reruns use a {manifest["diagnostic_envelope"]["wall_seconds"]}-second wall cap with the same RSS limit; '
        f'retry policy: {manifest.get("diagnostic_retry_policy", "retry every formal non-pass")}.\n\n'
        f'P1 route position: **{analysis["route_position"]}**.\n'
    )
    (output_dir / 'README.md').write_text(readme)
    return analysis


def _case(case):
    topology, size = case.split('-', 1)
    atoms = int(size)
    if topology not in ('chain', 'cube') or atoms not in (512, 768, 1024, 2048):
        raise ValueError(f'Unsupported frontier case: {case}')
    return topology, atoms


def _run_case(args, profile, case, route, timeout, report_path):
    command = [
        sys.executable, str(ROOT / 'tests/integration/joint_benchmark.py'),
        '--profile', profile, '--case', case,
        '--build-dir', str(args.build_dir), '--output', str(report_path),
        '--timeout', str(timeout), '--rss-limit', str(args.rss_limit),
        '--preconditioner', ROUTES[route], '--repeat', '1', '--warmup', '0',
        '--operator-rank', 'auto', '--operator-rank-seconds', '120',
        '--operator-rank-work-entries', '100000000', '--operator-rank-workspace-mib', '256',
        '--schwarz-core-atoms', '128', '--schwarz-overlap-hops', '1',
        '--schwarz-max-block-atoms', '512', '--schwarz-storage-mib', '512',
        '--schwarz-scratch-mib', '256',
    ]
    if route == 'FixedNeighbor':
        command.extend(['--fixed-core-atoms', '12'])
    print(f'Running {case} {route} {profile} with {timeout:g}s cap', flush=True)
    completed = subprocess.run(command, cwd=ROOT, check=False)
    if not report_path.is_file():
        raise RuntimeError(f'Benchmark runner did not write {report_path} (exit {completed.returncode})')
    report = read(report_path)
    return record_from_report(report, *_case(case), route,
                              'search-only' if profile == 'search' else 'full-endpoint',
                              timeout, report_path)


def run_campaign(args):
    output = args.output_dir
    individual = output / 'individual-results'
    cases = list(args.cases)
    if args.include_2048:
        missing_lower = set(DEFAULT_CASES) - set(cases)
        if missing_lower:
            raise ValueError('2048 frontier requires all 512, 768, and 1024 lower-size cases first')
        cases.extend(('chain-2048', 'cube-2048'))
    scopes = [('search-only', 'search', cases)]
    full_cases = [case for case in cases if _case(case)[1] in (512, 768)]
    if not args.search_only:
        scopes.append(('full-endpoint', 'solve', full_cases))
    manifest = {
        'schema_version': 1,
        'phase': 'P1 matched route performance frontier',
        'source_revision': subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                           check=True, capture_output=True, text=True).stdout.strip(),
        'source_sha256': source_hash(ROOT),
        'benchmark_source_hashes': {
            'route_runner': sha(Path(__file__)),
            'shared_runner': sha(ROOT / 'tests/integration/joint_benchmark.py'),
            'global_driver': sha(ROOT / 'tests/experiments/joint_sparse_benchmark.cpp'),
            'fixed_neighbor_driver': sha(ROOT / 'tests/experiments/joint_fixed_neighbor.cpp'),
        },
        'branch': subprocess.run(['git', 'branch', '--show-current'], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip(),
        'backend': 'SPQR',
        'eigen_threads': 1,
        'topologies': ['chain', 'cube'],
        'sizes': sorted({_case(case)[1] for case in cases}),
        'routes': list(ROUTES),
        'measurement_scopes': [name for name, _, _ in scopes],
        'formal_envelope': {'wall_seconds': args.formal_timeout, 'rss_bytes': args.rss_limit},
        'diagnostic_envelope': {'wall_seconds': args.diagnostic_timeout, 'rss_bytes': args.rss_limit},
        'diagnostic_retry_policy': 'deferred' if args.skip_diagnostics else 'retry every formal non-pass',
        'route_policy': {
            'LegacyCompact': 'global LegacyCompact search',
            'OperatorPcg': 'global OperatorPcg with Schwarz, 128 core atoms, one overlap hop',
            'FixedNeighbor': 'serial forward Gauss-Seidel, 12 atom core, at most one trusted accepted local update per block visit, local LegacyCompact',
        },
        'fixed_neighbor_policy': FIXED_NEIGHBOR_POLICY,
        'formal_results': [],
    }
    cache = (args.build_dir / 'CMakeCache.txt').read_text()
    if 'RHBM_GEM_JOINT_SPARSE_BACKEND:STRING=SPQR' not in cache:
        raise ValueError('Matched route frontier requires the SPQR benchmark build')
    records = []
    for scope_name, profile, scope_cases in scopes:
        for case in scope_cases:
            topology, atoms = _case(case)
            for route in ROUTES:
                formal_path = individual / scope_name / f'{case}-{route}-formal.json'
                formal_path.parent.mkdir(parents=True, exist_ok=True)
                if formal_path.is_file():
                    formal_report = read(formal_path)
                    formal_record = record_from_report(formal_report, topology, atoms, route,
                                                       scope_name, args.formal_timeout, formal_path)
                else:
                    formal_record = _run_case(args, profile, case, route, args.formal_timeout, formal_path)
                record = {
                    'topology': topology, 'atoms': atoms, 'rows': formal_record.get('rows'),
                    'parameter_count': formal_record.get('parameter_count'),
                    'route': route, 'measurement_scope': scope_name,
                    'formal_result': formal_record['result'],
                    'diagnostic_result': None,
                }
                manifest['formal_results'].append({
                    'topology': topology, 'atoms': atoms, 'route': route, 'scope': scope_name,
                    'status': formal_record['result']['status'],
                    'completed': formal_record['result']['completed'],
                    'path': (str(formal_path.relative_to(ROOT)) if formal_path.is_relative_to(ROOT)
                             else str(formal_path)),
                })
                if formal_record['result']['completed']:
                    record['diagnostic_result'] = dict(formal_record['result'], source='formal-envelope')
                else:
                    diagnostic_path = individual / scope_name / f'{case}-{route}-diagnostic.json'
                    if diagnostic_path.is_file():
                        diagnostic_report = read(diagnostic_path)
                        diagnostic_record = record_from_report(
                            diagnostic_report, topology, atoms, route, scope_name,
                            args.diagnostic_timeout, diagnostic_path)
                        record['diagnostic_result'] = dict(diagnostic_record['result'],
                                                           source='diagnostic-extended-envelope')
                        record['diagnostic_report_path'] = (
                            str(diagnostic_path.relative_to(ROOT)) if diagnostic_path.is_relative_to(ROOT)
                            else str(diagnostic_path))
                    elif args.skip_diagnostics:
                        record['diagnostic_result'] = {
                            'status': 'not-run', 'completed': False,
                            'reason': 'extended diagnostic deferred until formal frontier is recorded',
                            'source': 'deferred-by-policy',
                        }
                    else:
                        diagnostic_record = _run_case(args, profile, case, route,
                                                      args.diagnostic_timeout, diagnostic_path)
                        record['diagnostic_result'] = dict(diagnostic_record['result'],
                                                           source='diagnostic-extended-envelope')
                        record['diagnostic_report_path'] = (
                            str(diagnostic_path.relative_to(ROOT)) if diagnostic_path.is_relative_to(ROOT)
                            else str(diagnostic_path))
                records.append(record)
                write_outputs(output, manifest, records)
    manifest['analysis'] = analyze(records)
    write_outputs(output, manifest, records)
    return manifest['analysis']


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--formal-timeout', type=float, default=600.)
    parser.add_argument('--diagnostic-timeout', type=float, default=7200.)
    parser.add_argument('--rss-limit', type=int, default=4 * 1024**3)
    parser.add_argument('--include-2048', action='store_true')
    parser.add_argument('--search-only', action='store_true')
    parser.add_argument('--skip-diagnostics', action='store_true',
                        help='record formal results first and defer extended retries')
    parser.add_argument('--cases', nargs='+', default=DEFAULT_CASES)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.build_dir = args.build_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    if args.formal_timeout <= 0 or args.diagnostic_timeout < args.formal_timeout or args.rss_limit <= 0:
        parser.error('resource limits must be positive and diagnostic timeout must cover formal timeout')
    for case in args.cases:
        try:
            _case(case)
        except ValueError as error:
            parser.error(str(error))
        if _case(case)[1] == 2048 and not args.include_2048:
            parser.error('2048 frontier requires --include-2048')
    if args.include_2048 and set(DEFAULT_CASES) - set(args.cases):
        parser.error('--include-2048 requires all six 512/768/1024 lower-size cases')
    try:
        analysis = run_campaign(args)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.error(str(error))
    print(json.dumps({'route_position': analysis['route_position'],
                      'output_dir': str(args.output_dir)}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
