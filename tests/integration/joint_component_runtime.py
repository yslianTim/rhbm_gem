#!/usr/bin/env python3
"""Self-contained runtime regressions and fresh Map/Model entry points."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
from joint_runtime_support import read, write, require, scientific, differences, unpack

import joint_fixture_records as records
from joint_fixture_records import CATALOG
from joint_numerical_reference import replay, replay_passed
import numpy as np


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


def runtime_differences(expected, actual, scale):
    """Compare optimization evidence while allowing the certification contract to change."""
    delta = []
    checks=actual.get('runtime_checks')
    required_checks=('ac_feasibility','conditional_ac_optimality','independent_ac_optimality',
                     'independent_kkt_consistency','independent_replay','width_gradient_consistency',
                     'width_stationarity','design_identifiability','local_correction','numerical_identifiability')
    if not isinstance(checks,dict) or any(type(checks.get(name)) is not bool for name in required_checks):
        if actual.get('endpoint_verification_status')!='not-run' or actual.get('runtime_convergence')!='unavailable':
            delta.append('/runtime_checks/reference-free-contract')
    elif actual.get('runtime_convergence')=='passed' and not all(checks[name] for name in required_checks):
        delta.append('/runtime_convergence-with-failed-check')
    statuses=('passed','failed','unavailable','not-run')
    endpoint=actual.get('endpoint_verification_status')
    qualification=actual.get('numerical_qualification_status')
    usable=actual.get('usable_state')
    search_success=actual.get('search_success')
    if endpoint not in statuses:
        delta.append('/endpoint-verification-status')
    if qualification not in statuses:
        delta.append('/numerical-qualification-status')
    if type(usable) is not bool:
        delta.append('/trusted-state-availability')
    if type(search_success) is not bool:
        delta.append('/search-success')
    if qualification=='passed':
        if (usable is not True or search_success is not True or endpoint!='passed' or
                actual.get('runtime_convergence')!='passed' or not isinstance(checks,dict) or
                any(checks.get(name) is not True for name in required_checks)):
            delta.append('/qualification-without-complete-trusted-state')
    elif usable is True or search_success is True or actual.get('runtime_convergence')=='passed':
        delta.append('/unqualified-scientific-state')
    if actual.get('reference_solve_calls',0)!=0:
        delta.append('/reference_solve_calls')
    candidate=actual.get('candidate_state')
    if type(actual.get('candidate_available')) is not bool or actual['candidate_available'] != isinstance(candidate,dict):
        delta.append('/candidate-availability')
    if isinstance(candidate,dict):
        for key in ('beta','eta','b','objective'):
            value=candidate.get(key)
            if not np.all(np.isfinite(np.asarray(value,dtype=float))):
                delta.append('/candidate_state/'+key)
    if not isinstance(actual.get('search_success'), bool) or not actual.get('stop_reason'):
        delta.append('/search-termination')
    for key, maximum in (('accepted_updates', 100),):
        value = actual.get(key)
        if type(value) is not int or not 0 <= value <= maximum:
            delta.append('/'+key)
    profile_evaluations=actual.get('profile_evaluations')
    if type(profile_evaluations) is not int or profile_evaluations < 0:
        delta.append('/profile_evaluations')
    for trial in actual.get('trials', []):
        if trial['accepted'] and (not trial['valid'] or trial['trust_passed'] is not True):
            delta.append('/untrusted-accepted-trial')
    left, right = expected.get('last_trusted_state'), candidate
    if left is None:
        return delta
    if right is None:
        return delta+['/candidate_state']
    a, b = left['objective']/scale**2, right['objective']/scale**2
    objective_error=abs(a-b) if expected['runtime_convergence']=='passed' else b-a
    if not np.isfinite(b) or objective_error > 1e-12:
        delta.append('/normalized-objective')
    for key in ('beta','eta','b'):
        x, y = np.asarray(left[key]), np.asarray(right[key])
        if x.shape != y.shape or not np.all(np.isfinite(y)):
            delta.append('/candidate_state/'+key)
    return delta


def candidate_comparison(data, y, expected, actual, scale):
    old=expected.get('last_trusted_state')
    new=actual.get('candidate_state')
    if new is None:
        return None
    replay_new=replay(data,y,new,scale)
    result={
        'candidate_available':actual.get('candidate_available'),
        'primary_ac_kkt':new.get('projected_kkt'),
        'primary_width_gradient_inf':float(np.max(np.abs(new['b_gradient']))) if isinstance(new.get('b_gradient'),list) and new['b_gradient'] else None,
        'independent_python_candidate_ac_kkt':replay_new['projected_kkt'],
        'independent_python_candidate_width_gradient_inf':float(np.max(np.abs(replay_new['b_gradient']))) if replay_new['b_gradient'].size else 0.,
        'candidate_residual_norm':float(np.linalg.norm(replay_new['residual'])),
        'production_independent_replay':actual.get('independent_replay'),
        'search_accepted_updates':actual.get('accepted_updates'),
        'legacy_accepted_updates':expected.get('accepted_updates'),
        'search_profile_evaluations':actual.get('profile_evaluations'),
        'legacy_profile_evaluations':expected.get('profile_evaluations'),
        'search_stop_reason':actual.get('stop_reason'),
        'legacy_stop_reason':expected.get('stop_reason'),
        'legacy_runtime_convergence':expected.get('runtime_convergence'),
        'reference_free_runtime_convergence':actual.get('runtime_convergence'),
        'parameter_comparison_basis':'historical_endpoint_may_be_a_different_sweep',
        'endpoint_verification_status':actual.get('endpoint_verification_status'),
        'numerical_qualification_status':actual.get('numerical_qualification_status'),
        'rank_assessment_status':actual.get('rank_assessment_status'),
        'verification_reason':actual.get('verification_reason'),
        'early_stop':actual.get('accepted_updates',0)<expected.get('accepted_updates',0)}
    if old is None:
        result['legacy_endpoint_available']=False
        return result
    replay_old=replay(data,y,old,scale)
    result.update(legacy_endpoint_available=True,
        candidate_objective_difference=0.5*(np.dot(replay_new['residual'],replay_new['residual'])-
                                             np.dot(replay_old['residual'],replay_old['residual']))/(scale*scale),
        ac_kkt_difference=replay_new['projected_kkt']-replay_old['projected_kkt'],
        width_gradient_difference_inf=float(np.max(np.abs(replay_new['b_gradient']-replay_old['b_gradient']))),
        residual_difference_norm=float(np.linalg.norm(replay_new['residual']-replay_old['residual'])),
        candidate_parameter_max_abs={key:float(np.max(np.abs(np.asarray(new[key])-np.asarray(old[key]))))
                                     for key in ('beta','eta','b')})
    return result


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
            y = data['y64' if case.endswith('double') else 'y32']
            state = actual['record']['candidate_state']
            if state is not None and state.get('valid') and all(key in state for key in ('projected_kkt','b_gradient','rss')):
                require(replay_passed(data, y, state, max(1., np.linalg.norm(y))), 'Independent candidate replay failed')
            expected = runtime_expected(cases[case]['expected'])
            scale = max(1., np.linalg.norm(data['y64' if case.endswith('double') else 'y32']))
            delta = (differences(expected, actual['record']) if args.strict_history else
                     runtime_differences(expected, actual['record'], scale))
            if not args.strict_history:
                for parity in actual['derivative_parity'].values():
                    require(parity['status'] != 'failed', 'Same-state derivative parity failed: '+str(parity))
                    if expected['runtime_convergence'] == 'passed':
                        require(parity['status'] == 'passed', 'Missing derivative parity for a converged endpoint')
            require(actual['api_contract_passed'] and not delta, dataset+'/'+case+': '+str(delta[:20]))
            old_state = expected.get('last_trusted_state')
            comparison=candidate_comparison(data,y,expected,actual['record'],scale)
            results.append(dict(dataset=dataset, case=case, passed=True,
                                legacy_runtime_convergence=expected['runtime_convergence'],
                                reference_free_runtime_convergence=actual['record']['runtime_convergence'],
                                convergence_status_changed=expected['runtime_convergence']!=actual['record']['runtime_convergence'],
                                endpoint_same_active_face=(old_state['active_atoms'] == state['active_atoms']) if old_state and state and 'active_atoms' in state else None,
                                numerical_comparison=comparison))
    require(results, 'No regression cases selected')
    result = dict(passed=True, cases=results, oracle='pre-extraction frozen component states',
                  comparison='strict-history' if args.strict_history else 'reference-free-runtime-numerical',
                  certification_contract='reference-free-v1')
    write(args.work_dir/'summary.json', result)
    return result


def summarize(root):
    files = sorted(Path(root).rglob('runtime.json'))
    require(files, 'No runtime results')
    for path in files:
        value = read(path)
        require(value['regular_certificate'] == 'not-run', 'Runtime claimed an offline certificate')
        require(value['prediction_available'] and value['prediction_reassembly_passed'] and value['monolithic_parity_passed'],
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
    sub = commands.add_parser('summarize'); sub.add_argument('--run', type=Path, required=True)
    sub = commands.add_parser('compare')
    for name in ('left', 'right', 'output'): sub.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'regression' and args.case and not args.dataset:
        parser.error('--case requires --dataset')
    if args.command == 'regression': result = regression(args)
    elif args.command == 'physical-smoke': result = physical_smoke(args.executable)
    elif args.command == 'summarize': result = summarize(args.run)
    elif args.command == 'compare':
        result = compare(args.left, args.right); write(args.output, result)
    else:
        inputs = (args.model, args.map, args.output) if args.command == 'run' else (args.output,)
        execute(args.executable, args.command, *inputs); result = summarize(args.output)
    print(result)


if __name__ == '__main__': main()
