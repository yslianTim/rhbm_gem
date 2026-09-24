"""Four bounded PR2/PR3 campaigns. Resource or numerical failures never promote search."""
from __future__ import annotations
import argparse
import hashlib
import io
import statistics
import subprocess
import tarfile
import time
from pathlib import Path
import joint_validation as v
from joint_operator_validation import fingerprint
from joint_runtime_support import unpack

BASELINE='c6c869cd4f4939104dfde96984ae17f736a8ce4a'
KINDS=('schwarz','diagonal','identity')


def require_current_build(build):
    cache=(build/'CMakeCache.txt').read_text()
    source=Path(next(line.split('=',1)[1] for line in cache.splitlines() if line.startswith('CMAKE_HOME_DIRECTORY:')))
    libraries=list((build/'src').glob('librhbm_gem.*'))
    outputs=[build/'bin/joint_sparse_benchmark',*libraries]
    if not libraries or any(not p.is_file() for p in outputs): raise ValueError(f'Build is incomplete: {build}')
    sources=[p for folder in ('src','include','cmake') for p in (source/folder).rglob('*') if p.is_file() and p.suffix in ('.cpp','.hpp','.h','.cmake','.txt')]
    sources += [source/'CMakeLists.txt']
    core_time=max(p.stat().st_mtime_ns for p in sources)
    driver_time=(source/'tests/experiments/joint_sparse_benchmark.cpp').stat().st_mtime_ns
    if core_time>min(p.stat().st_mtime_ns for p in libraries) or max(core_time,driver_time)>(build/'bin/joint_sparse_benchmark').stat().st_mtime_ns:
        raise ValueError(f'Sources are newer than measured binaries; rebuild before starting a campaign: {build}')


def frozen_hash():
    raw=subprocess.check_output(['git','archive',BASELINE,'src','include','cmake','CMakeLists.txt'],cwd=v.ROOT)
    with tarfile.open(fileobj=io.BytesIO(raw)) as archive:
        return v.digest({m.name:hashlib.sha256(archive.extractfile(m).read()).hexdigest() for m in archive
                         if m.isfile() and Path(m.name).suffix in ('.cpp','.hpp','.h','.cmake','.txt')})


def scientific_parity(left,right):
    """Trajectories may differ; returned evidence and qualified endpoints may not."""
    a,b=left.get('returned_assessment'),right.get('returned_assessment')
    if (a is None)!=(b is None): return dict(passed=False,reason='state-availability')
    if a is None:
        a,b=left.get('search',{}).get('initial',{}),right.get('search',{}).get('initial',{})
        same=a.get('valid') is False and b.get('valid') is False and a.get('reason')==b.get('reason')
        return dict(passed=same,reason='equivalent-explicit-unavailability' if same else 'uncertified-states')
    failures=[]
    for key in ('runtime_convergence','runtime_checks','runtime_failure'):
        if a.get(key)!=b.get(key): failures.append(key)
    for key in ('design_spectrum','width_spectrum','profile_jacobian_spectrum'):
        if a.get(key,{}).get('rank')!=b.get(key,{}).get('rank'): failures.append(key+'/rank')
    sa,sb=left['returned_state'],right['returned_state']
    if sa.get('active_atoms')!=sb.get('active_atoms'): failures.append('active-face')
    if left.get('search_completed') and not right.get('search_completed'): failures.append('search-incomplete')
    if a.get('runtime_convergence')=='passed':
        for key in ('beta','b'):
            x,y=sa.get(key,[]),sb.get(key,[])
            if len(x)!=len(y) or any(abs(u-w)>1e-10*(1+max(abs(u),abs(w))) for u,w in zip(x,y)): failures.append(key)
        oa=sa['objective']/left['search']['residual_scale']**2
        ob=sb['objective']/right['search']['residual_scale']**2
        if abs(oa-ob)>1e-12: failures.append('normalized-objective')
    return dict(passed=not failures,reason='same-evidence-and-qualified-endpoint',differences=failures)


def compare(root):
    reports={stage:v.read(root/(stage+'.json')) for stage in ('baseline','eigen','spqr','large-local')}
    if any(b['source_sha256']!=frozen_hash() for b in reports['baseline']['builds'].values()): raise ValueError('Wrong baseline source')
    out=dict(comparisons={},performance={},large={},promote=False)
    numerical=True; complete=all(r.get('finished') is True for r in reports.values())
    required={'single-128','single-512','heterogeneous-168','chain-128','cube-128','chain-512','cube-512'}
    catalog=v.read(v.ROOT/'tests/fixtures/joint_component/catalog.json')
    required.update(name+'/'+case for name,item in catalog['datasets'].items() if name!='heterogeneous-168' for case in item['default_cases'])
    for backend in ('eigen','spqr'):
        base=reports['baseline']; candidate=reports[backend]
        present={key.split('/',2)[2] for key in base['runs'] if key.startswith(backend+'/legacy/')}
        complete &= required<=present
        if base['inputs']!=candidate['inputs']: raise ValueError('Different input fingerprints')
        for key,runs in base['runs'].items():
            if not key.startswith(backend+'/legacy/'): continue
            case=key.split('/',2)[2]
            for kind in KINDS:
                actual=candidate['runs'].get(f'{backend}/{kind}/{case}',[]); checks=[]
                for x,y in zip(runs,actual):
                    if x['process']['status']!='completed' or y['process']['status']!='completed' or x.get('result',{}).get('stage')!='complete' or y.get('result',{}).get('stage')!='complete':
                        checks.append(dict(passed=None,reason='incomplete-measurement'))
                    else: checks.append(scientific_parity(x['result'],y['result']))
                done=len(checks)==3 and all(c['passed'] is not None for c in checks)
                passed=all(c['passed'] for c in checks) if done else None
                out['comparisons'][f'{backend}/{kind}/{case}']=dict(complete=done,passed=passed,checks=checks)
                numerical &= passed is True; complete &= done
        def median(kind,case):
            report=base if kind=='legacy' else candidate
            runs=report['runs'].get(f'{backend}/{kind}/{case}',[])
            if len(runs)!=3 or any(r['process']['status']!='completed' or r.get('result',{}).get('stage')!='complete' or not r['result'].get('search_completed') for r in runs): return None
            return statistics.median(r['result']['search_seconds'] for r in runs)
        legacy,new=median('legacy','single-512'),median('schwarz','single-512')
        wins=[]
        for case in ('heterogeneous-168','chain-128','cube-128','chain-512','cube-512'):
            a,b=median('schwarz',case),median('diagonal',case)
            if a is not None and b is not None and a<b and out['comparisons'][f'{backend}/schwarz/{case}']['passed'] is True and out['comparisons'][f'{backend}/diagonal/{case}']['passed'] is True: wins.append(case)
        out['performance'][backend]=dict(legacy_512=legacy,schwarz_512=new,
            within_ten_percent=None if legacy is None or new is None else new<=1.1*legacy,schwarz_wins=wins)
    for key,runs in reports['large-local']['runs'].items():
        row=runs[0]; x=row.get('result',{}); counts=x.get('work',{})
        passed=row['process']['status']=='completed' and x.get('local_available') and x.get('local_passed') and all(counts.get(k)==0 for k in ('numeric','reference_solves','free_design_svds','derivative_preparations'))
        out['large'][key]=dict(passed=bool(passed),process_status=row['process']['status'],reason=x.get('local_reason'))
    e,s=out['performance']['eigen'],out['performance']['spqr']
    expected_large={f'{backend}/{topology}-{size}' for backend in ('eigen','spqr') for topology in ('chain','cube') for size in (128,512,2000,5000,10000)}
    complete &= expected_large==set(out['large'])
    out['comparator_sha256']=v.sha(Path(__file__))
    out['numerical_passed']=numerical; out['comparisons_complete']=complete
    out['performance_passed']=e['within_ten_percent'] is True and s['within_ten_percent'] is True and bool(set(e['schwarz_wins'])&set(s['schwarz_wins']))
    # Test/regression evidence must be reviewed separately before changing the default.
    out['benchmark_promotion_eligible']=complete and numerical and out['performance_passed'] and all(x['passed'] for x in out['large'].values())
    out['promote']=False
    out['promotion_note']='Legacy remains default until both backend regressions and all benchmark gates are verified.'
    v.write(root/'comparison.json',out); return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',required=True,choices=('baseline','eigen','spqr','large-local','compare'))
    p.add_argument('--work-dir',type=Path,required=True); p.add_argument('--input-dir',type=Path,required=True)
    p.add_argument('--eigen',type=Path); p.add_argument('--spqr',type=Path)
    p.add_argument('--model',type=Path,default=Path('/Users/yslian/data/6z6u.cif'))
    p.add_argument('--map',type=Path,default=Path('/Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map'))
    args=p.parse_args(); root=args.work_dir.resolve(); inputs=args.input_dir.resolve()
    if args.stage=='compare':
        result=compare(root); raise SystemExit(0 if result['benchmark_promotion_eligible'] else 2)
    if args.eigen is None or args.spqr is None: p.error('Both build directories are required')
    receipt=root/(args.stage+'.json')
    if receipt.exists(): raise SystemExit('Keep prior receipts; refusing to restart an existing campaign.')
    builds={'eigen':args.eigen.resolve(),'spqr':args.spqr.resolve()}
    selected=('eigen','spqr') if args.stage in ('baseline','large-local') else (args.stage,)
    for backend in selected: require_current_build(builds[backend])
    root.mkdir(parents=True,exist_ok=True); deadline=time.monotonic()+4800
    report=dict(baseline_commit=BASELINE,limits=dict(process_seconds=600,rss_bytes=4*1024**3,campaign_seconds=4800),
                builds={b:fingerprint(builds[b]) for b in selected},runs={},controls={},inputs={},
                driver_sha256=v.sha(v.ROOT/'tests/experiments/joint_sparse_benchmark.cpp'),
                generator_sha256=v.sha(v.ROOT/'tests/support/JointOperatorWorkload.cpp'),
                runner_sha256=v.sha(Path(__file__)),numeric_environment={k:v.ENV[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')})
    if args.stage=='baseline' and any(b['source_sha256']!=frozen_hash() for b in report['builds'].values()): raise ValueError('Baseline production sources differ')
    def save(): v.write(receipt,report)
    def run(command,directory,output):
        process=v.monitored(command,directory,deadline); row=dict(process=process)
        if output.is_file(): row['result']=v.read(output)
        print(f'{directory.relative_to(root)}: {process["status"]}',flush=True); return row
    save()
    if args.stage=='large-local':
        for backend in selected:
            for topology in ('chain','cube'):
                for size in (128,512,2000,5000,10000):
                    key=f'{backend}/{topology}-{size}'; d=root/args.stage/key; d.mkdir(parents=True,exist_ok=True); output=d/'state.json'
                    report['runs'][key]=[run([builds[backend]/'bin/joint_sparse_benchmark','synthetic',topology,size,'local',output,'--resources'],d,output)]; save()
        report['finished']=True; save(); return
    catalog=v.read(v.ROOT/'tests/fixtures/joint_component/catalog.json'); cases={}
    for name,item in catalog['datasets'].items():
        wanted=item['default_cases'] if name!='heterogeneous-168' else ['first-stage-float32']
        if not wanted: continue
        folder=unpack(v.ROOT/'tests/fixtures/joint_component/catalog.json',name,inputs)
        report['inputs'][name]={path.name:v.sha(path) for path in folder.iterdir() if path.is_file()}
        for case in wanted: cases[name if name=='heterogeneous-168' else name+'/'+case]=['fixture',folder,case]
    for size in (128,512):
        folder=inputs/f'single-{size}'
        if not (folder/'widths.json').exists():
            if args.stage!='baseline': raise ValueError('Generate frozen inputs during baseline first')
            generated=run([builds['eigen']/'bin/joint_validation','single',size,folder],root/'generation'/str(size),folder/'unused.json')
            report['controls'][f'generate-{size}']=generated
            if generated['process']['status']=='completed':
                report['controls'][f'prepare-{size}']=run([builds['eigen']/'bin/joint_sparse_benchmark','prepare',folder/'input.cif',folder/'input.map',folder/'widths.json'],root/'prepare'/str(size),folder/'widths.json')
        paths=[folder/'input.cif',folder/'input.map',folder/'widths.json']
        report['inputs'][f'single-{size}']={path.name:v.sha(path) for path in paths if path.is_file()}
        cases[f'single-{size}']=['initial',*paths]
    for topology in ('chain','cube'):
        for size in (128,512): cases[f'{topology}-{size}']=['synthetic',topology,size,'fixed']
    save()
    # Put the declared performance gate first; campaign exhaustion is explicit.
    ordered=['single-512','heterogeneous-168','single-128','chain-128','cube-128','chain-512','cube-512']
    ordered += sorted(set(cases)-set(ordered))
    for case in ordered:
        for backend in selected:
            for kind in (('legacy',) if args.stage=='baseline' else KINDS):
                key=f'{backend}/{kind}/{case}'; report['runs'][key]=[]
                for repeat in range(1,4):
                    d=root/args.stage/key/str(repeat); d.mkdir(parents=True,exist_ok=True); output=d/'state.json'
                    row=run([builds[backend]/'bin/joint_sparse_benchmark',*cases[case],output,'--search',kind,'--resources'],d,output)
                    report['runs'][key].append(row); save()
                    if row['process']['status']!='completed': break
    for backend in selected:
        d=root/args.stage/backend/'6z6u'; d.mkdir(parents=True,exist_ok=True)
        if args.model.is_file() and args.map.is_file():
            report['controls'][backend+'/6z6u']=run([builds[backend]/'bin/joint_component_runtime','run',args.model,args.map,d/'output'],d,d/'output/completion.json')
            report['controls'][backend+'/6z6u']['scope']='unchanged public workflow; production search remains legacy until promotion'
        else: report['controls'][backend+'/6z6u']=dict(process=dict(status='not-run-missing-input'))
        save()
    report['finished']=True; save()


if __name__=='__main__': main()
