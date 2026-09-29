"""Report-only compact receipt verification shared with reference campaigns."""
import itertools
import statistics
from pathlib import Path
import experiment_io as v
from joint_validation_checks import sparse_parity as parity, scaled, svd_parity, audit_parity


def summarize_compact_receipt(r, root, *, cases, modes, backends, command_cases, performance_ratio_limits):
    out = dict(fixed={}, audits={}, replay={}, commands={}, baseline_controls={})
    numerical = bool(r.get('matrix_files'))
    for path,digest in r.get('matrix_files', {}).items():
        if v.sha(root/path) != digest:
            raise RuntimeError('Captured matrix hash mismatch')
    for backend in backends:
        out['fixed'][backend] = {}
        for case in cases:
            group = r['fixed'].get(backend, {}).get(case, {})
            records = out['fixed'][backend][case] = {}
            for mode in modes:
                runs = group.get(mode, [])
                complete = len(runs) == 3 and all(x['process']['status'] == 'completed' for x in runs)
                numerical &= complete
                if not complete:
                    records[mode] = dict(complete=False)
                    continue
                states = [x['result'] for x in runs]
                records[mode] = dict(complete=True, reference_seconds=statistics.median(s['reference_seconds'] for s in states),
                    derivative_seconds=statistics.median(s['derivative_seconds'] for s in states),
                    combined_seconds=statistics.median(s['reference_seconds']+s['derivative_seconds'] for s in states),
                    reference_svd_seconds=statistics.median(s['work']['reference_svd_seconds'] for s in states),
                    free_design_svd_seconds=statistics.median(s['work']['free_design_svd_seconds'] for s in states),
                    peak_rss_bytes=max(max(x['process'].get('os_process_peak_rss_bytes') or 0,
                                          x['process'].get('sampled_tree_peak_rss_bytes', 0)) for x in runs))
            pairs = []
            for mode in (mode for mode in modes if mode != 'legacy'):
                for a,b in itertools.product(group.get('legacy', []), group.get(mode, [])):
                    if 'result' not in a or 'result' not in b:
                        numerical = False
                        continue
                    check = parity(a['result'], b['result'])
                    passed = (all(x.get('passed', False) for x in check.values()) and
                              a['result']['initial_b'] == b['result']['initial_b'] and
                              a['result']['trust']['passed'] and b['result']['trust']['passed'] and
                              a['result']['derivative_valid'] and b['result']['derivative_valid'])
                    pairs.append(dict(mode=mode, passed=passed, checks=check)); numerical &= passed
            records['comparisons'] = pairs
            audits = r['audits'].get(backend, {}).get(case, {})
            checks = {}
            for mode in (mode for mode in modes if mode != 'legacy'):
                a,b = audits.get('legacy', {}), audits.get(mode, {})
                if a.get('process', {}).get('status') == b.get('process', {}).get('status') == 'completed':
                    if any(v.sha(root/row['output']) != row['output_sha256'] for row in (a,b)):
                        raise RuntimeError('Audit output hash mismatch')
                    checks[mode] = audit_parity(v.read(root/a['output']), v.read(root/b['output']))
                else:
                    checks[mode] = dict(passed=False, unavailable=True)
                numerical &= checks[mode]['passed']
            out['audits'][backend+'/'+case] = checks
            control = r['baseline_controls'].get(backend, {}).get(case, {})
            if control.get('process', {}).get('status') == 'completed' and group.get('legacy'):
                check = parity(control['result'], group['legacy'][0]['result'])
                passed = all(x['passed'] for x in check.values())
                out['baseline_controls'][backend+'/'+case] = dict(passed=passed, checks=check)
                numerical &= passed
            else:
                numerical = False
    for name, modes in r.get('replay', {}).items():
        checks = {}
        for mode in (mode for mode in modes if mode != 'legacy'):
            a,b = modes.get('legacy', {}), modes.get(mode, {})
            checks[mode] = svd_parity(a['result'], b['result']) if 'result' in a and 'result' in b else dict(passed=False)
            numerical &= checks[mode]['passed']
        out['replay'][name] = checks
    numerical &= all(any(k.startswith(case+'/'+role+'-') for k in r.get('replay', {}))
                     for case in cases for role in ('reference','free-design'))
    performance = True
    for case in cases:
        group = out['fixed']['spqr'][case]; a,b = group['legacy'],group['auto']
        if not a['complete'] or not b['complete']:
            performance = False
            continue
        keys = ('reference_svd_seconds', 'free_design_svd_seconds', 'combined_seconds') if case == 'single-512' else ('combined_seconds',)
        ratios = {k: b[k]/a[k] for k in keys}
        limit = performance_ratio_limits[case]
        group['performance_gate'] = dict(ratios=ratios, maximum_ratio=limit, passed=all(x <= limit for x in ratios.values()))
        performance &= group['performance_gate']['passed']
    for case in command_cases:
        groups = r.get('commands', {}).get(case, {})
        out['commands'][case] = {}
        for role in ('baseline','candidate'):
            runs = groups.get(role, [])
            complete = len(runs) == 3 and all(x['status'] == 'completed' for x in runs)
            out['commands'][case][role] = dict(complete=complete,
                converged=complete and all(x.get('runtime_convergence') == 'passed' for x in runs),
                median_seconds=statistics.median(x['total_command_seconds'] for x in runs) if complete else None,
                peak_rss_bytes=max((max(x.get('os_process_peak_rss_bytes') or 0, x.get('sampled_tree_peak_rss_bytes',0)) for x in runs), default=0))
    numerical &= out['commands']['single-128']['baseline']['converged'] and out['commands']['single-128']['candidate']['converged']
    command_pairs=[]
    for a,b in itertools.product(r.get('command_endpoints', {}).get('baseline', []), r.get('command_endpoints', {}).get('candidate', [])):
        x,y=a['assembled_state'],b['assembled_state']
        errors={k:scaled(x[k],y[k]) for k in ('ac','b')}
        objective=abs(x['objective']/a['observation_scale']**2-y['objective']/b['observation_scale']**2)
        passed=a['atom_ids']==b['atom_ids'] and all(e is not None and e<=1e-10 for e in errors.values()) and objective<=1e-12
        command_pairs.append(dict(errors=errors,normalized_objective_difference=objective,passed=passed)); numerical &= passed
    numerical &= len(command_pairs)==9
    out['commands']['single-128']['endpoint_comparisons']=command_pairs
    out.update(numerical_passed=bool(numerical),performance_passed=bool(performance),
               compact_gate_passed=bool(numerical and performance),
               complete_512_passed=out['commands']['single-512']['candidate']['converged'])
    return out
