#!/usr/bin/env python3
"""Run selected offline certificates without a historical experiment directory."""
import argparse
from pathlib import Path
import subprocess
from joint_runtime_support import read, write, require, unpack, differences
from joint_component_runtime import CATALOG
import joint_fixture_records as records
from joint_offline_support import certificate, validate_audit, replay_passed
import numpy as np


def certify_local(data, case, directory):
    fit, audit, scope, context = (read(directory/name) for name in ('fit.json', 'audit.json', 'scope.json', 'context.json'))
    validate_audit(audit, fit)
    require(scope['parent_snapshot_sha256'] == data['hash'], 'Changed parent snapshot')
    if fit['usable_state']:
        state = fit['last_trusted_state']
        require(scope['state'] == {key: state[key] for key in ('eta', 'beta')}, 'Replaced actual component state')
        require(all(fit['primary'][key] == state[key] for key in ('eta', 'beta')), 'Assessment changed coefficients')
    if scope['weak_direction_available']:
        directions = context['audit']['directions']
        require(all(abs(np.linalg.norm(d)-1) <= 1e-12 for d in directions), 'Non-unit local direction')
        require(directions[2] == fit['width_spectrum']['weak_directions'][0], 'Wrong local weak direction')
    y = data['y64' if case.endswith('double') else 'y32']
    for trial in fit.get('trials', []):
        if trial['accepted']:
            require(trial.get('trust', {}).get('passed') is True, 'Untrusted accepted trial')
            require(replay_passed(data, y, trial, fit['residual_scale']), 'Independent accepted-state replay failed')
    result = certificate(fit, audit, replay_passed(data, y, fit.get('primary'), fit['residual_scale']))
    result['scope'] = 'component-local'
    result['checks'].update(local_weak_direction=scope['weak_direction_available'],
                            actual_state_profile=fit.get('assembled_profile_agrees', False), usable_state=fit['usable_state'])
    result['failures'] = [k for k, v in result['checks'].items() if v is False]
    result['unavailable'] = [k for k, v in result['checks'].items() if v is None]
    result['regular_qualified'] = all(v is True for v in result['checks'].values())
    write(directory/'certificate.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--catalog', type=Path, default=CATALOG)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--two-step-only', action='store_true')
    parser.add_argument('--dataset', choices=tuple(read(CATALOG)['datasets']))
    parser.add_argument('--case', default='first-stage-double')
    args = parser.parse_args()
    cases = [(args.dataset, args.case)] if args.dataset else [
        ('baseline', 'first-stage-double'), ('near-0.02', 'narrower-double'),
        ('weak-1e-4', 'first-stage-double'), ('active-a', 'first-stage-double')]
    if args.two_step_only and not args.dataset:
        cases = [(dataset, case) for dataset, entry in read(args.catalog)['datasets'].items() for case in entry['default_cases']]
    summary = []
    import tempfile
    for dataset, case in cases:
        source = unpack(args.catalog, dataset, args.work_dir/'fixtures')
        data = records.load(source)
        # Every invocation owns fresh audit output and fresh precision caches.
        args.work_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=dataset+'-', dir=args.work_dir) as temp:
            output = Path(temp)/'audit'
            if args.two_step_only:
                subprocess.run([str(args.executable), 'two-step-fixture', str(source), case, str(output)], check=True)
                expected = read(source/'cases.json')[case]['expected']
                expected = {k: v for k, v in expected.items() if k in ('qualification_checks', 'qualification_failure')}
                actual = read(output/'0.json')
                require(not differences(expected, actual), 'Historical two-step evidence changed: '+dataset+'/'+case)
                summary.append(dict(dataset=dataset, case=case, passed=True))
                continue
            subprocess.run([str(args.executable), 'fixture', str(source), case, str(output)], check=True)
            cert = certify_local(data, case, output/'0')
            expected = read(source/'cases.json')[case]['offline_expected']
            require({k: cert[k] for k in expected} == expected, 'Local certificate changed: '+dataset+'/'+case)
            write(args.work_dir/dataset/(case+'-certificate.json'), cert)
            summary.append(dict(dataset=dataset, case=case, regular=cert['regular_qualified']))
    write(args.work_dir/'summary.json', dict(passed=True, scopes=summary))


if __name__ == '__main__': main()
