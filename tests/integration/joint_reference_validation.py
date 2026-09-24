"""Incremental endpoint-reference acceptance; historical compact evidence stays immutable."""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import itertools
import os
from pathlib import Path
import shutil
import sqlite3
import statistics
import subprocess
import tarfile
import time
from types import SimpleNamespace

import joint_compact_validation as compact
import joint_validation as v

BASELINE = 'f1ac45c36945e8e89e588af5c305158a7723b268a21c5c7ab0b22a3ae676a755'
CANDIDATE = '62d178d17643a4a8195aceac17a5e4ced5af59cfb79a7e288f87bdd5d82c2de1'
EXPORTS = ('result.sqlite', 'joint_result_validation.json', 'joint_atoms_validation.csv')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source_proof(ref):
    archive = tarfile.open(fileobj=io.BytesIO(subprocess.check_output(['git', 'archive', ref], cwd=v.ROOT)))
    names = sorted(m.name for m in archive if m.isfile() and (
        m.name == 'CMakeLists.txt' or
        m.name.startswith('src/') and m.name.endswith(('.cpp', '.hpp', '/CMakeLists.txt')) or
        m.name.startswith('include/') and m.name.endswith('.hpp') or
        m.name.startswith('cmake/') and m.name.endswith('.cmake')))
    files = {n: hashlib.sha256(archive.extractfile(n).read()).hexdigest() for n in names}
    digest = hashlib.sha256(''.join(n+'\n'+h+'\n' for n,h in files.items()).encode()).hexdigest()
    return dict(commit=subprocess.check_output(['git','rev-parse',ref],cwd=v.ROOT,text=True).strip(),
                production_source_sha256=digest, files=files)


def build_info(root, expected):
    header = (root/'generated/Release/SimulationBuildInfo.hpp').read_text()
    require(f'SOURCE_SHA256 "{expected}"' in header, 'Build source fingerprint mismatch: '+str(root))
    files = [root/'bin/RHBM-GEM', root/'bin/joint_sparse_benchmark', *sorted((root/'src').glob('librhbm_gem.*'))]
    return dict(source_fingerprint=header, binaries={str(p):v.sha(p) for p in files},
                configuration=(root/'generated/Release/SimulationConfiguration-CXX.txt').read_text(),
                linked_libraries=subprocess.check_output(['otool','-L',str(root/'src/librhbm_gem.dylib')],text=True))


def endpoint_comparison(a, b):
    x,y = a.get('assembled_state'),b.get('assembled_state')
    if not x or not y:
        return dict(passed=False, reason='missing-state')
    errors = {k:compact.scaled(x[k],y[k]) for k in ('ac','b')}
    # JointState.objective is already divided by observation_scale squared.
    objective = abs(x['objective']-y['objective'])
    ranks = lambda r: [[(e['name'],e['scope'],e['rank']) for e in c['ranks']] for c in r['components']]
    face = lambda r: [k for k,a in enumerate(r['assembled_state']['ac'][::2]) if a > 0]
    evidence = lambda r: [[(e['name'],e['scope'],e['status']) for e in c['evidence']] for c in r['components']]
    same = a['atom_ids']==b['atom_ids'] and ranks(a)==ranks(b) and face(a)==face(b) and evidence(a)==evidence(b)
    return dict(parameter_errors=errors, normalized_objective_difference=objective,
                same_identity_rank_face_evidence=same,
                passed=same and all(e is not None and e<=1e-10 for e in errors.values()) and objective<=1e-12)


def exports(root, expected_widths=None):
    require(all((root/n).is_file() for n in EXPORTS), 'Missing completed-command export: '+str(root))
    outcome = v.read(root/EXPORTS[1])
    require(outcome['schema_version'] in (3,4,5) and outcome['assembled_state'] is not None, 'Invalid exported state')
    with (root/EXPORTS[2]).open() as f:
        require(len(list(csv.reader(f)))==len(outcome['atom_ids'])+1, 'CSV atom count mismatch')
    with sqlite3.connect((root/EXPORTS[0]).resolve().as_uri()+'?mode=ro',uri=True) as db:
        require(db.execute('PRAGMA quick_check').fetchone()==('ok',), 'Invalid SQLite export')
    if expected_widths is not None:
        require(outcome['initialization']['b']==expected_widths, 'Production initialization changed')
    return outcome


def historical(root):
    r = v.read(root/'receipt.json')
    require(r['finished'], 'Historical campaign did not finish')
    for name,b in r['builds'].items():
        expected = compact.BASELINE_SOURCE if name.startswith('baseline_') else BASELINE
        require(f'SOURCE_SHA256 "{expected}"' in b['source_fingerprint'], 'Historical source mismatch')
    for case,files in r['inputs'].items():
        for name,digest in files.items():
            require(v.sha(root/'inputs'/case/name)==digest, 'Historical input hash mismatch')
    for backend,cases in r['fixed'].items():
        for case,modes in cases.items():
            for mode,runs in modes.items():
                for i,run in enumerate(runs,1):
                    path=root/'fixed'/backend/case/mode/str(i)
                    require(v.read(path/'state.json')==run['result'], 'Historical fixed sample mismatch')
                    require(v.read(path/'process.json')==run['process'], 'Historical process mismatch')
    for backend,cases in r['baseline_controls'].items():
        for case,run in cases.items():
            require(v.read(root/'baseline-controls'/backend/case/'state.json')==run['result'], 'Historical control mismatch')
    for key,modes in r['replay'].items():
        for mode,run in modes.items():
            require(v.read(root/'replay'/key/mode/'state.json')==run['result'], 'Historical replay mismatch')
    out = compact.summary(r,root)  # Also verifies captured matrices and audit hashes.
    endpoints={}
    for case,roles in r['commands'].items():
        endpoints[case]={}
        for role,runs in roles.items():
            states=endpoints[case][role]=[]
            for run in runs:
                path=root/'commands'/role/case/f'run-{run["repetition"]}'
                require(v.read(path/'receipt.json')==run, 'Historical command receipt mismatch')
                if run['status']=='completed':
                    state=exports(path,v.read(root/'inputs'/case/'widths.json')['b'])
                    expected=compact.BASELINE_SOURCE if role=='baseline' else BASELINE
                    require(state['metadata']['software']['source_sha256']==expected, 'Historical export source mismatch')
                    states.append(state)
    comparisons=[endpoint_comparison(a,b) for a,b in itertools.product(
        endpoints['single-128']['baseline'],endpoints['single-128']['candidate'])]
    require(len(comparisons)==9 and all(c['passed'] for c in comparisons), 'Historical normalized endpoint comparison failed')
    out['corrected_128_endpoint_comparisons']=comparisons
    require(out['compact_gate_passed'], 'Historical compact gate incomplete or failed')
    return out


def command_group(rows, root, candidate=False):
    complete = len(rows)==3 and sorted(x.get('repetition',0) for x in rows)==[1,2,3] and len({x.get('directory') for x in rows})==3 and all(x['status']=='completed' and x['total_command_seconds']<=600 and
        max(x.get('os_process_peak_rss_bytes') or 0,x.get('sampled_tree_peak_rss_bytes',0))<=4*1024**3 for x in rows)
    outcomes=[]
    for row in rows:
        if row['status']!='completed':
            continue
        directory=root/row['directory']
        for name,digest in row.get('export_hashes',{}).items():
            require(v.sha(directory/name)==digest, 'Export hash mismatch')
        require(set(row.get('export_hashes',{}))==set(EXPORTS), 'Missing export fingerprint')
        outcomes.append(exports(directory,v.read(root/'inputs'/row['case']/'widths.json')['b']))
    converged=complete and all(x['runtime_convergence']=='passed' for x in outcomes)
    zero_references=bool(outcomes) and all(x['costs']['search_reference_seconds']==0 and
        all(c['reference_evaluations']==0 for c in x['components']) for x in outcomes)
    return dict(complete=complete,converged=converged,search_reference_zero=zero_references,
        passed=converged and (not candidate or zero_references),
        median_seconds=statistics.median(x['total_command_seconds'] for x in rows) if complete else None,
        peak_rss_bytes=max((max(x.get('os_process_peak_rss_bytes') or 0,x.get('sampled_tree_peak_rss_bytes',0)) for x in rows),default=0),
        median_costs={k:statistics.median(x['costs'][k] for x in outcomes) for k in outcomes[0]['costs']} if complete else {}),outcomes


def summarize(r, root):
    # Manifest covers historical evidence and every new measured/audited artifact.
    for path,digest in r.get('files',{}).items():
        require(v.sha(root/path)==digest, 'Evidence hash mismatch: '+path)
    history=historical(root/'historical')
    require(history==r['historical'], 'Historical aggregation changed')
    out=dict(historical=history,commands={},audits={},finished=r['finished'])
    numerical=True
    for case in ('single-128','single-512'):
        groups=r['commands'].get(case,{})
        a,aa=command_group(groups.get('baseline',[]),root)
        b,bb=command_group(groups.get('candidate',[]),root,True)
        comparisons=[endpoint_comparison(x,y) for x,y in itertools.product(aa,bb)]
        row=dict(baseline=a,candidate=b,endpoint_comparisons=comparisons,
                 speedup=a['median_seconds']/b['median_seconds'] if a['complete'] and b['complete'] else None,
                 rss_ratio=b['peak_rss_bytes']/a['peak_rss_bytes'] if a['complete'] and b['complete'] else None)
        out['commands'][case]=row
        if case=='single-128':numerical &= a['passed'] and b['passed'] and len(comparisons)==9 and all(x['passed'] for x in comparisons)
    for backend in ('spqr','eigen'):
        for case in compact.CASES:
            modes=r['audits'].get(backend,{}).get(case,{})
            checks={}
            for mode in compact.MODES:
                entry=modes.get(mode,{})
                if entry.get('process',{}).get('status')!='completed':
                    checks[mode]=dict(passed=False,reason='incomplete-audit');numerical=False;continue
                current=v.read(root/entry['output'])
                previous=v.read(root/'historical/audits'/backend/case/mode/'state.json')
                historic=compact.audit_parity(previous,current)
                legacy=modes.get('legacy',{})
                local=compact.audit_parity(v.read(root/legacy['output']),current) if legacy.get('process',{}).get('status')=='completed' else dict(passed=False)
                checks[mode]=dict(historical=historic,latest_legacy=local,passed=historic['passed'] and local['passed'])
                numerical &= checks[mode]['passed']
            out['audits'][backend+'/'+case]=checks
    intact=r['finished'] and bool(r.get('files'))
    out.update(latest_numerical_passed=bool(intact and numerical),
        historical_compact_passed=bool(intact and r['historical']['compact_gate_passed']),
        latest_512_passed=bool(intact and numerical and out['commands']['single-512']['candidate']['passed']))
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('prior-work-dir','baseline-spqr','spqr','eigen','work-dir'):
        parser.add_argument('--'+key,type=Path,required=key=='work-dir')
    parser.add_argument('--report-only',action='store_true')
    args=parser.parse_args();root=args.work_dir.resolve()
    if args.report_only:
        v.write(root/'summary.json',summarize(v.read(root/'receipt.json'),root));return
    require(all(getattr(args,k) for k in ('prior_work_dir','baseline_spqr','spqr','eigen')), 'Supply all build and prior-work directories')
    require(not root.exists(), 'Use a fresh work directory')
    old=args.prior_work_dir.resolve();history=historical(old)
    proofs={ref:source_proof(ref) for ref in ('ce58c897','e32919f3','c9e0f8c9')}
    for ref,expected in zip(proofs,(compact.BASELINE_SOURCE,BASELINE,CANDIDATE)):
        require(proofs[ref]['production_source_sha256']==expected, 'Commit source proof mismatch')
    prior=v.read(old/'receipt.json')
    for name,digest in prior['sources'].items():
        content=subprocess.check_output(['git','show','e32919f3:'+name],cwd=v.ROOT)
        require(hashlib.sha256(content).hexdigest()==digest, 'Historical harness source mismatch: '+name)
    builds={k:getattr(args,k).resolve() for k in ('baseline_spqr','spqr','eigen')}
    info={k:build_info(p,BASELINE if k=='baseline_spqr' else CANDIDATE) for k,p in builds.items()}
    config=lambda text:dict(line.split('=',1) for line in text.splitlines() if '=' in line)
    settings=('version','compiler','system','configuration','flags','release_flags','boost','eigen','joint_sparse_backend','spqr','cholmod',
              'suitesparse_config','spqr_transitive_link','spqr_blas_version','openmp','root','dependency_provider','compile_features','compile_options','compile_definitions','link_options')
    left,right=(config(info[k]['configuration']) for k in ('baseline_spqr','spqr'))
    require(all(left.get(k)==right.get(k) for k in settings),'Baseline/candidate build settings differ')
    shutil.copytree(old,root/'historical');shutil.copytree(old/'inputs',root/'inputs')
    r=dict(finished=False,historical=history,source_proofs=proofs,builds=info,commands={},audits={},
           limits=dict(campaign_seconds=4800,pipeline_seconds=600,rss_bytes=4*1024**3),
           sources={n:v.sha(v.ROOT/n) for n in ('tests/integration/joint_reference_validation.py','tests/integration/joint_compact_validation.py','tests/integration/joint_validation.py')})
    v.process_tree_rss(os.getpid())
    started=time.monotonic();deadline=started+4800
    def save():
        r['wall_seconds']=time.monotonic()-started;v.write(root/'receipt.json',r)
    save()
    # Frozen serialized initialization is verified before timing complete commands.
    for case in ('single-128','single-512'):
        inputs=root/'inputs'/case;expected=v.read(inputs/'widths.json')
        for role in ('baseline_spqr','spqr'):
            path=root/'initialization'/case/role
            p=v.monitored([builds[role]/'bin/joint_sparse_benchmark','prepare',inputs/'input.cif',inputs/'input.map',path/'widths.json'],path,deadline)
            require(p['status']=='completed' and v.read(path/'widths.json')==expected,'Initialization preflight failed')
    for case in ('single-128','single-512'):
        groups=r['commands'][case]={'baseline':[],'candidate':[]};stopped=set();inputs=root/'inputs'/case
        for repetition in (1,2,3):
            for role,key in (('baseline','baseline_spqr'),('candidate','spqr')):
                if role in stopped:continue
                a=SimpleNamespace(work_dir=root,stage_dir='commands/'+role,cli=builds[key]/'bin/RHBM-GEM')
                row=v.resource_run(a,case,inputs/'input.cif',inputs/'input.map',repetition,deadline)
                row.update(case=case,directory=f'commands/{role}/{case}/run-{repetition}')
                if row['status']=='completed':
                    path=root/row['directory'];outcome=exports(path,v.read(inputs/'widths.json')['b'])
                    require(outcome['metadata']['software']['source_sha256']==(BASELINE if role=='baseline' else CANDIDATE), 'Command export source mismatch')
                    row['export_hashes']={n:v.sha(path/n) for n in EXPORTS}
                    row['components']=[{k:c[k] for k in ('id','stop_reason','search_completed','profile_evaluations','reference_evaluations','accepted_updates','runtime_convergence')} for c in outcome['components']]
                else:stopped.add(role)
                groups[role].append(row);save();print('command',case,role,repetition,row['status'],flush=True)
    for backend in ('spqr','eigen'):
        cases=r['audits'][backend]={}
        for case in compact.CASES:
            modes=cases[case]={};inputs=root/'inputs'/case
            command=[builds[backend]/'bin/joint_sparse_benchmark',*(['fixture',inputs,'first-stage-float32'] if case=='heterogeneous-168' else ['initial',inputs/'input.cif',inputs/'input.map',inputs/'widths.json'])]
            for mode in compact.MODES:
                directory=root/'audits'/backend/case/mode;output=directory/'state.json'
                p=v.monitored([*command,output,'--svd-mode',mode,'--audit'],directory,deadline)
                modes[mode]=dict(process=p,output=str(output.relative_to(root)));save()
                print('audit',backend,case,mode,p['status'],flush=True)
    for b in info.values():
        require(all(v.sha(p)==h for p,h in b['binaries'].items()),'A measured binary changed')
    r['files']={str(p.relative_to(root)):v.sha(p) for p in sorted(root.rglob('*')) if p.is_file() and p!=root/'receipt.json'}
    r['finished']=True;save();v.write(root/'summary.json',summarize(r,root))


if __name__=='__main__':main()
