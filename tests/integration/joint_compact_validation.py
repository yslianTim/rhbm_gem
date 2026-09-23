"""Bounded three-way compact SVD acceptance; report-only rechecks saved evidence."""
from __future__ import annotations
import argparse
import itertools
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time
from types import SimpleNamespace

import numpy as np
import joint_validation as v
from joint_runtime_support import unpack
from joint_sparse_validation import parity

MODES = ('legacy', 'values', 'auto')
CASES = ('single-128', 'heterogeneous-168', 'single-512')
BASELINE_SOURCE = '9718531067e72cd3ecf180d12ae8033b8e39360b7fe1d49d24ae847de4240fa0'


def scaled(a, b):
    x, y = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if x.shape != y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        return None
    return float(np.max(np.abs(x-y)/(1+np.maximum(np.abs(x), np.abs(y))), initial=0))


def svd_parity(a, b):
    x, y = np.asarray(a['singular_values'], dtype=float), np.asarray(b['singular_values'], dtype=float)
    same = bool(a['valid'] == b['valid'] and a['rank'] == b['rank'] and x.shape == y.shape
                and np.isfinite(x).all() and np.isfinite(y).all())
    maximum = float(x[0]) if x.size else 0.
    error = float(np.max(np.abs(x-y), initial=0)/max(maximum, np.finfo(float).tiny)) if same else None
    solution = scaled(a['solution'], b['solution'])
    threshold = abs(a['threshold']-b['threshold']) <= 1e-10*max(abs(a['threshold']), np.finfo(float).tiny)
    # Retain the actual weak-end values; a global norm alone does not verify rank.
    near = [dict(index=k, baseline=float(s), candidate=float(y[k]), threshold=a['threshold'])
            for k,s in enumerate(x) if same and (s <= 2*a['threshold'] or k >= len(x)-3)]
    return dict(same_rank_status=same, spectrum_error=error, solution_error=solution,
                threshold_agrees=threshold, weak_values=near,
                passed=same and a['valid'] and error <= 1e-10 and threshold and solution is not None and solution <= 1e-10)


def audit_parity(a, b):
    endpoints = parity(a, b)
    checks = [svd_parity(x, y) for x,y in zip(a.get('svd_records', []), b.get('svd_records', []))]
    same_calls = (len(a.get('svd_records', [])) == len(b.get('svd_records', [])) > 0
                  and all(all(x.get(k)==y.get(k) for k in ('role','rows','columns','relative_threshold','absolute_override'))
                          for x,y in zip(a.get('svd_records', []), b.get('svd_records', []))))
    left, right = a.get('derivative_audit', {}), b.get('derivative_audit', {})
    arrays = {k: scaled(left.get(k, []), right.get(k, []))
              for k in ('coefficients', 'correction', 'projected', 'jacobian', 'response')}
    passed = (all(c['passed'] for c in endpoints.values()) and a['trust']['passed'] and b['trust']['passed']
              and a['initial_b'] == b['initial_b'] and same_calls and all(c['passed'] for c in checks)
              and left.get('valid') is True and right.get('valid') is True
              and a['derivative_reason'] == b['derivative_reason']
              and all(k in left and k in right and len(left[k]) > 0 for k in arrays)
              and all(e is not None and e <= 1e-10 for e in arrays.values()))
    return dict(endpoints=endpoints, spectra=checks, derivative_errors=arrays, passed=passed)


def summary(r, root=Path('.')):
    out = dict(fixed={}, audits={}, replay={}, commands={}, baseline_controls={})
    numerical = bool(r.get('matrix_files'))
    for path,digest in r.get('matrix_files', {}).items():
        if v.sha(root/path) != digest:
            raise RuntimeError('Captured matrix hash mismatch')
    for backend in ('eigen', 'spqr'):
        out['fixed'][backend] = {}
        for case in CASES:
            group = r['fixed'].get(backend, {}).get(case, {})
            records = out['fixed'][backend][case] = {}
            for mode in MODES:
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
            for mode in ('values', 'auto'):
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
            for mode in ('values', 'auto'):
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
        for mode in ('values', 'auto'):
            a,b = modes.get('legacy', {}), modes.get(mode, {})
            checks[mode] = svd_parity(a['result'], b['result']) if 'result' in a and 'result' in b else dict(passed=False)
            numerical &= checks[mode]['passed']
        out['replay'][name] = checks
    numerical &= all(any(k.startswith(case+'/'+role+'-') for k in r.get('replay', {}))
                     for case in CASES for role in ('reference','free-design'))
    performance = True
    for case in CASES:
        group = out['fixed']['spqr'][case]; a,b = group['legacy'],group['auto']
        if not a['complete'] or not b['complete']:
            performance = False
            continue
        keys = ('reference_svd_seconds', 'free_design_svd_seconds', 'combined_seconds') if case == 'single-512' else ('combined_seconds',)
        ratios = {k: b[k]/a[k] for k in keys}
        limit = .7 if case == 'single-512' else 1.1
        group['performance_gate'] = dict(ratios=ratios, maximum_ratio=limit, passed=all(x <= limit for x in ratios.values()))
        performance &= group['performance_gate']['passed']
    for case in ('single-128','single-512'):
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


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('spqr','eigen','baseline-spqr','baseline-eigen','work-dir'):
        parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--report-only',action='store_true')
    args=parser.parse_args()
    for key,value in vars(args).items():
        if isinstance(value,Path):setattr(args,key,value.resolve())
    receipt=args.work_dir/'receipt.json'
    if args.report_only:
        v.write(args.work_dir/'summary.json',summary(v.read(receipt),args.work_dir));return
    if receipt.exists():raise SystemExit('Use a fresh work directory; preserve all previous attempts.')
    if not (args.spqr/'bin/joint_validation').is_file():
        raise SystemExit('Build the joint_validation target in the candidate SPQR build before measuring.')
    v.process_tree_rss(os.getpid())
    started=time.monotonic();deadline=started+4800
    r=dict(platform=platform.platform(),limits=dict(campaign_seconds=4800,process_or_pipeline_seconds=600,rss_bytes=4*1024**3),
           fixed={},audits={},replay={},baseline_controls={},commands={},command_endpoints={},inputs={},builds={},finished=False)
    r['sources']={name:v.sha(v.ROOT/name) for name in ('tests/experiments/joint_sparse_benchmark.cpp','tests/experiments/joint_validation.cpp',
        'tests/integration/joint_compact_validation.py','tests/integration/joint_sparse_validation.py',
        'tests/integration/joint_validation.py','tests/integration/joint_runtime_support.py')}
    roots=dict(eigen=args.eigen,spqr=args.spqr,baseline_eigen=args.baseline_eigen,baseline_spqr=args.baseline_spqr)
    for name,root in roots.items():
        binaries=[root/'bin/RHBM-GEM',root/'bin/joint_sparse_benchmark',root/'bin/joint_validation',*list((root/'src').glob('librhbm_gem.*'))]
        r['builds'][name]=dict(binaries={str(p):v.sha(p) for p in binaries if p.is_file()},
            configuration=(root/'generated/Release/SimulationConfiguration-CXX.txt').read_text(),
            source_fingerprint=(root/'generated/Release/SimulationBuildInfo.hpp').read_text())
        if name.startswith('baseline_') and BASELINE_SOURCE not in r['builds'][name]['source_fingerprint']:
            raise SystemExit('Baseline must be the frozen ce58c897 production sources.')
        if platform.system()=='Darwin':
            r['builds'][name]['linked_libraries']=subprocess.check_output(['otool','-L',str(root/'src/librhbm_gem.dylib')],text=True)
    def save():
        r['wall_seconds']=time.monotonic()-started;v.write(receipt,r)
    def run(command,directory):
        output=directory/'state.json'; process=v.monitored([*command,output],directory,deadline)
        return dict(process=process,result=v.read(output) if output.exists() else {})
    inputs={}
    save()
    for count in (128,512):
        case=f'single-{count}';directory=args.work_dir/'inputs'/case
        generated=v.monitored([args.spqr/'bin/joint_validation','single',count,directory],directory/'generation',deadline)
        if generated['status']!='completed':raise RuntimeError('Input generation failed')
        inputs[case]=(directory/'input.cif',directory/'input.map',directory/'widths.json')
        prepared=v.monitored([args.spqr/'bin/joint_sparse_benchmark','prepare',*inputs[case]],directory/'prepare',deadline)
        if prepared['status']!='completed':raise RuntimeError('Initialization failed')
        r['inputs'][case]={p.name:v.sha(p) for p in inputs[case]}
    dataset=unpack(v.ROOT/'tests/fixtures/joint_component/catalog.json','heterogeneous-168',args.work_dir/'inputs')
    r['inputs']['heterogeneous-168']={p.name:v.sha(p) for p in dataset.iterdir() if p.is_file()}
    def command(root,case):
        return [root/'bin/joint_sparse_benchmark',*(['fixture',dataset,'first-stage-float32'] if case=='heterogeneous-168' else ['initial',*inputs[case]])]
    # All formal timing precedes matrix capture and replay diagnostics.
    for backend in ('spqr','eigen'):
        r['fixed'][backend]={};r['baseline_controls'][backend]={}
        for case in CASES:
            groups=r['fixed'][backend][case]={}
            for mode in MODES:
                runs=groups[mode]=[]
                for repetition in (1,2,3):
                    directory=args.work_dir/'fixed'/backend/case/mode/str(repetition);output=directory/'state.json'
                    process=v.monitored([*command(roots[backend],case),output,'--svd-mode',mode],directory,deadline)
                    runs.append(dict(process=process,**({'result':v.read(output)} if output.exists() else {})));save()
                    print(f'fixed {backend} {case} {mode} {repetition}: {process["status"]}',flush=True)
                    if process['status']!='completed':break
            r['baseline_controls'][backend][case]=run(command(roots['baseline_'+backend],case),args.work_dir/'baseline-controls'/backend/case);save()
    for case in ('single-128','single-512'):
        groups=r['commands'][case]={}
        for role,root in [('baseline',args.baseline_spqr),('candidate',args.spqr)]:
            runs=groups[role]=[]
            cmd_args=SimpleNamespace(work_dir=args.work_dir,stage_dir='commands/'+role,cli=root/'bin/RHBM-GEM')
            for repetition in (1,2,3):
                record=v.resource_run(cmd_args,case,*inputs[case][:2],repetition,deadline);runs.append(record)
                if case=='single-128' and record['status']=='completed':
                    path=args.work_dir/cmd_args.stage_dir/case/f'run-{repetition}'/'joint_result_validation.json'
                    outcome=v.read(path)
                    r['command_endpoints'].setdefault(role,[]).append({k:outcome[k] for k in ('atom_ids','observation_scale','assembled_state')})
                save();print(f'command {case} {role} {repetition}: {record["status"]}',flush=True)
                if record['status']!='completed':break
    for backend in ('spqr','eigen'):
        r['audits'][backend]={}
        for case in CASES:
            records=r['audits'][backend][case]={}
            for mode in MODES:
                directory=args.work_dir/'audits'/backend/case/mode;output=directory/'state.json'
                options=['--audit']
                if backend=='spqr' and mode=='auto':options=['--capture',directory/'matrices']
                process=v.monitored([*command(roots[backend],case),output,'--svd-mode',mode,*options],directory,deadline)
                records[mode]=dict(process=process,output=str(output.relative_to(args.work_dir)),
                    output_sha256=v.sha(output) if output.exists() else None);save()
    for case in CASES:
        directory=args.work_dir/'audits/spqr'/case/'auto/matrices'
        for path in sorted(directory.glob('*.json')):
            key=case+'/'+path.stem;modes=r['replay'][key]={}
            for mode in MODES:
                modes[mode]=run([args.spqr/'bin/joint_sparse_benchmark','replay',path,mode],args.work_dir/'replay'/case/path.stem/mode);save()
    r['matrix_files']={str(p.relative_to(args.work_dir)):v.sha(p) for p in (args.work_dir/'audits/spqr').glob('*/auto/matrices/*') if p.is_file()}
    for build in r['builds'].values():
        if any(v.sha(path)!=digest for path,digest in build['binaries'].items()):
            raise RuntimeError('A measured binary changed during the campaign.')
    r['finished']=True;save()
    v.write(args.work_dir/'summary.json',summary(r,args.work_dir))


if __name__=='__main__':main()
