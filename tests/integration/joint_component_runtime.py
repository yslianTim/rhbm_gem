#!/usr/bin/env python3
"""Self-contained runtime regressions and fresh Map/Model entry points."""
import argparse
import json
from statistics import median
from pathlib import Path
import subprocess
import tempfile
from joint_runtime_support import read, write, require, scientific, differences, unpack

import joint_fixture_records as records
from joint_offline_support import replay_passed
import numpy as np

CATALOG = Path(__file__).resolve().parents[1]/'fixtures/joint_component/catalog.json'


def execute(executable, *args):
    subprocess.run([str(executable), *map(str, args)], check=True)


def runtime_expected(expected):
    out = {k: v for k, v in expected.items() if k not in ('qualification_checks', 'qualification_failure')}
    checks = expected.get('qualification_checks')
    if checks is not None:
        out['runtime_checks'] = {k: v for k, v in checks.items() if k != 'derivative'}
    required = ('inner', 'b_gradient', 'local_correction', 'identified')
    values = [checks.get(k) for k in required] if checks is not None else []
    out['runtime_convergence'] = ('unavailable' if not expected['usable_state'] or checks is None
                                  else 'failed' if any(v is False for v in values)
                                  else 'passed' if all(v is True for v in values) else 'unavailable')
    return out


def backend_differences(expected, actual, scale):
    """Fixed backend contract; historical search traces remain in differences()."""
    delta = []
    for key in ('usable_state', 'runtime_convergence', 'runtime_checks'):
        if (key in expected) != (key in actual) or expected.get(key) != actual.get(key):
            delta.append('/'+key)
    if not isinstance(actual.get('search_success'), bool) or not actual.get('stop_reason'):
        delta.append('/search-termination')
    for key, maximum in (('accepted_updates', 100), ('profile_evaluations', 200)):
        value = actual.get(key)
        if type(value) is not int or not 0 <= value <= maximum:
            delta.append('/'+key)
    for trial in actual.get('trials', []):
        if trial['accepted'] and (not trial['valid'] or trial['trust_passed'] is not True):
            delta.append('/untrusted-accepted-trial')
    left, right = expected.get('last_trusted_state'), actual.get('last_trusted_state')
    if (left is None) != (right is None):
        return delta+['/last_trusted_state']
    if left is None:
        return delta
    same_face = left.get('active_atoms') == right.get('active_atoms')
    endpoint_flags = ('valid', 'feasible', 'kkt_passed') + (('free_rank',) if same_face else ())
    for key in endpoint_flags:
        if key not in left or key not in right or left[key] != right[key]:
            delta.append('/last_trusted_state/'+key)
    a, b = left['objective']/scale**2, right['objective']/scale**2
    converged = expected['runtime_convergence'] == 'passed'
    if not np.isfinite(b) or (abs(a-b) if converged else b-a) > 1e-12:
        delta.append('/normalized-objective')
    if converged:
        for key in ('beta', 'eta', 'b'):
            x, y = np.asarray(left[key]), np.asarray(right[key])
            if x.shape != y.shape or not np.all(np.isfinite(y)) or np.any(np.abs(x-y) > 1e-10*(1+np.maximum(np.abs(x), np.abs(y)))):
                delta.append('/last_trusted_state/'+key)
    return delta


def regression(args):
    catalog = read(args.catalog)
    results = []
    for dataset, entry in catalog['datasets'].items():
        selected = entry['extended_cases'] if args.extended else entry['default_cases']
        if args.dataset and dataset != args.dataset:
            continue
        if args.case:
            selected = [args.case]
        elif args.all_starts:
            selected = [start+'-'+precision for precision in catalog['precisions'] for start in catalog['starts']]
        if not selected:
            continue
        source = unpack(args.catalog, dataset, args.work_dir/'fixtures')
        cases = read(source/'cases.json')
        data = records.load(source)
        for case in selected:
            require(case in cases, 'Unknown fixture case: '+case)
            output = args.work_dir/'results'/dataset/(case+'.json')
            execute(args.executable, 'fixture', source, case, output)
            actual = read(output)
            records.verify_census(data, actual['census'])
            state = actual['record']['last_trusted_state']
            if state is not None:
                y = data['y64' if case.endswith('double') else 'y32']
                require(replay_passed(data, y, state, max(1., np.linalg.norm(y))), 'Independent endpoint replay failed')
            expected = runtime_expected(cases[case]['expected'])
            scale = max(1., np.linalg.norm(data['y64' if case.endswith('double') else 'y32']))
            delta = (differences(expected, actual['record']) if args.strict_history else
                     backend_differences(expected, actual['record'], scale))
            if not args.strict_history:
                for parity in actual['backend_parity'].values():
                    require(parity['status'] != 'failed', 'Same-state backend parity failed: '+str(parity))
                    if expected['runtime_convergence'] == 'passed':
                        require(parity['status'] == 'passed', 'Missing backend parity for a converged endpoint')
            require(actual['api_contract_passed'] and not delta, dataset+'/'+case+': '+str(delta[:20]))
            old_state = expected.get('last_trusted_state')
            results.append(dict(dataset=dataset, case=case, passed=True,
                                endpoint_same_active_face=(old_state['active_atoms'] == state['active_atoms']) if old_state and state else None))
    require(results, 'No regression cases selected')
    result = dict(passed=True, cases=results, oracle='pre-extraction frozen component states', comparison='strict-history' if args.strict_history else 'backend-numerical')
    write(args.work_dir/'summary.json', result)
    return result


def summarize(root):
    files = sorted(Path(root).rglob('runtime.json'))
    require(files, 'No runtime results')
    for path in files:
        value = read(path)
        require(value['regular_certificate'] == 'not-run', 'Runtime claimed an offline certificate')
        require(value['prediction_available'] and value['prediction_reassembly_difference'] == 0 and value['monolithic_parity_passed'],
                'Incomplete or inconsistent fresh-input result: '+str(path))
    return dict(passed=True, cases=len(files))


def compare(left, right):
    a = {p.relative_to(left): scientific(read(p)) for p in left.rglob('*.json')}
    b = {p.relative_to(right): scientific(read(p)) for p in right.rglob('*.json')}
    require(a and b, 'Missing comparison records')
    delta = differences(a, b)
    require(not delta, 'Runtime comparison differs: '+str(delta[:20]))
    return dict(passed=True, records=len(a))


def physical_smoke(executable):
    with tempfile.TemporaryDirectory(prefix='joint-runtime-') as temporary:
        root = Path(temporary)
        execute(executable, 'physical', root/'physical')
        summarize(root/'physical')
        execute(executable, 'run', root/'physical/model.cif', root/'physical/map.mrc', root/'loaded')
        summarize(root/'loaded')
    return dict(passed=True, cases=3)


def benchmark(args):
    samples = []
    for _ in range(3):
        result = subprocess.run([str(args.executable), str(args.input), args.case],
                                check=True, capture_output=True, text=True)
        sample = json.loads(result.stdout)
        require(sample['available'], 'Benchmark has no assembled state')
        samples.append(sample)
    keys = ('construction_seconds', 'search_seconds', 'assessment_seconds', 'assembly_seconds',
            'total_seconds', 'process_peak_rss_bytes')
    result = dict(samples=samples, median={key: median(row[key] for row in samples) for key in keys})
    write(args.output, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    reg = commands.add_parser('regression')
    reg.add_argument('--catalog', type=Path, default=CATALOG)
    reg.add_argument('--executable', type=Path, required=True)
    reg.add_argument('--work-dir', type=Path, required=True)
    reg.add_argument('--extended', action='store_true')
    reg.add_argument('--strict-history', action='store_true', help='Also require historical search decisions and counts')
    reg.add_argument('--all-starts', action='store_true')
    reg.add_argument('--dataset', choices=tuple(read(CATALOG)['datasets']))
    reg.add_argument('--case')
    smoke = commands.add_parser('physical-smoke')
    smoke.add_argument('--executable', type=Path, required=True)
    for command in ('run', 'physical'):
        sub = commands.add_parser(command)
        sub.add_argument('--executable', type=Path, required=True)
        sub.add_argument('--output', type=Path, required=True)
        if command == 'run':
            sub.add_argument('--model', type=Path, required=True)
            sub.add_argument('--map', type=Path, required=True)
    sub = commands.add_parser('benchmark')
    for name in ('executable', 'input', 'output'): sub.add_argument('--'+name, type=Path, required=True)
    sub.add_argument('--case', default='first-stage-float32')
    sub = commands.add_parser('summarize'); sub.add_argument('--run', type=Path, required=True)
    sub = commands.add_parser('compare')
    for name in ('left', 'right', 'output'): sub.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'regression' and args.case and not args.dataset:
        parser.error('--case requires --dataset')
    if args.command == 'regression': result = regression(args)
    elif args.command == 'benchmark': result = benchmark(args)
    elif args.command == 'physical-smoke': result = physical_smoke(args.executable)
    elif args.command == 'summarize': result = summarize(args.run)
    elif args.command == 'compare':
        result = compare(args.left, args.right); write(args.output, result)
    else:
        inputs = (args.model, args.map, args.output) if args.command == 'run' else (args.output,)
        execute(args.executable, args.command, *inputs); result = summarize(args.output)
    print(result)


if __name__ == '__main__': main()
