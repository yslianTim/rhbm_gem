"""One 30-minute fixed-state/fixed-step campaign; never resumes or promotes search."""
from __future__ import annotations
import argparse
import math
import statistics
import time
from pathlib import Path
import joint_validation as v
from joint_operator_validation import fingerprint
from joint_search_validation import finite, require_current_build
from joint_runtime_support import unpack

BASELINE='ec24f3056427408cb1978cff004ed477908fb0d0'
MODES=('A','B','C')
CASES=('chain-8','single-128','heterogeneous-168','single-512')


def vector_close(a,b,relative=1e-8,absolute=1e-12):
    if not isinstance(a,list) or not isinstance(b,list) or not a or len(a)!=len(b): return False
    if not all(finite(x) for x in a+b): return False
    return math.hypot(*(x-y for x,y in zip(a,b)))<=absolute+relative*math.hypot(*a)


def parity(a,b):
    failures=[]
    for key in ('valid','reason'):
        if key not in a or a.get(key)!=b.get(key): failures.append(key)
    sa,sb=a.get('state_control',{}),b.get('state_control',{})
    for key in ('eta','beta','free_columns','scale','rank_rows'):
        if key not in sa or sa.get(key)!=sb.get(key): failures.append('state/'+key)
    for key in ('residual','gradient'):
        if not vector_close(sa.get(key),sb.get(key),relative=0,absolute=1e-12): failures.append('state/'+key)
    oa,ob=sa.get('objective'),sb.get('objective')
    if not finite(oa) or not finite(ob) or abs(oa-ob)>1e-12: failures.append('state/objective')
    ra,rb=a.get('rank',[]),b.get('rank',[])
    if len(ra)!=1 or len(rb)!=1: failures.append('rank/missing')
    else:
        for key in ('valid','rank','rows','columns','relative_threshold','absolute_override'):
            if key not in ra[0] or ra[0].get(key)!=rb[0].get(key): failures.append('rank/'+key)
        x,y=ra[0].get('singular_values'),rb[0].get('singular_values')
        if not isinstance(x,list) or not isinstance(y,list) or not x or len(x)!=len(y) or not all(finite(z) for z in x+y) or max(abs(u-w) for u,w in zip(x,y))>1e-10*max(x): failures.append('rank/spectrum')
        x,y=ra[0].get('threshold'),rb[0].get('threshold')
        if not finite(x) or not finite(y) or abs(x-y)>1e-10*max(abs(x),abs(y)): failures.append('rank/threshold')
    if a.get('valid') is not True or b.get('valid') is not True: failures.append('unavailable-operator')
    for key in ('apply','adjoint','normal'):
        if not vector_close(a.get(key),b.get(key)): failures.append(key)
    x,y=a.get('operator_gradient',[]),b.get('operator_gradient',[])
    if not x or len(x)!=len(y) or not all(finite(z) for z in x+y) or any(abs(u-w)>1e-13+2e-9*abs(u) for u,w in zip(x,y)): failures.append('operator-gradient')
    for result in (a,b):
        if result.get('normal_q_actions')!=(2 if result.get('mode')=='normal' else 6): failures.append('normal-q-actions')
        if any(result.get('work',{}).get(k)!=0 for k in ('reference_solves','derivative_preparations')): failures.append('forbidden-work')
    aa,bb=a.get('steps',[]),b.get('steps',[])
    if len(aa)!=3 or len(bb)!=3: failures.append('steps/missing')
    for index,kind in enumerate(('identity','diagonal','schwarz')):
        if len(aa)<=index or len(bb)<=index: continue
        x,y=aa[index],bb[index]
        if x.get('kind')!=kind or y.get('kind')!=kind or x.get('valid') is not True or y.get('valid') is not True: failures.append(kind+'/unavailable'); continue
        if not vector_close(x.get('step'),y.get('step'),1e-10,1e-10): failures.append(kind+'/step')
        if not finite(x.get('predicted')) or not finite(y.get('predicted')) or abs(x['predicted']-y['predicted'])>1e-12: failures.append(kind+'/predicted')
        if any(not finite(r.get('true_residual')) or r['true_residual']>1e-10 for r in (x,y)): failures.append(kind+'/true-residual')
    return dict(passed=not failures,differences=failures)


def improvement(a,b,qualified):
    if len(a)!=3 or len(b)!=3 or not all(finite(x) and x>=0 for x in a+b): return dict(complete=False,observed=False)
    ma,mb=statistics.median(a),statistics.median(b)
    return dict(complete=True,baseline_median=ma,candidate_median=mb,baseline_range=[min(a),max(a)],candidate_range=[min(b),max(b)],
                ratio=mb/ma if ma else None,observed=qualified and mb<ma and all(y<x for x,y in zip(a,b)))


def total_bounds(result,index):
    """Old receipts omitted gradient setup; summed adjoint time bounds that action.

    Retain their subtotals unchanged. New receipts explicitly time gradient setup.
    These are sums of tracked numerical phases, not full process wall time.
    """
    steps=result.get('steps',[])
    if len(steps)<=index: return None,None
    total=steps[index].get('total_seconds')
    if not finite(total) or total<0: return None,None
    if result.get('total_includes_gradient') is True: return total,total
    adjoint=result.get('search_work',{}).get('operator_adjoint_seconds')
    if not finite(adjoint) or adjoint<0: return total,None
    return total,total+adjoint


def compare(root):
    report=v.read(root/'campaign.json'); out=dict(comparisons={},performance={},promote=False)
    for case in CASES:
        for backend in ('eigen','spqr'):
            for left,right in (('A','B'),('B','C'),('A','C')):
                x,y=(report['runs'].get(f'{case}/{backend}/{mode}',[]) for mode in (left,right))
                checks=[]
                for a,b in zip(x,y):
                    if any(r['process']['status']!='completed' or r.get('result',{}).get('stage')!='complete' for r in (a,b)):
                        checks.append(dict(passed=None,reason='incomplete-measurement'))
                    elif a.get('state_sha256')!=report.get('states',{}).get(case,{}).get('sha256') or a.get('state_sha256')!=b.get('state_sha256'):
                        checks.append(dict(passed=False,reason='state-fingerprint'))
                    else: checks.append(parity(a['result'],b['result']))
                complete=len(x)==len(y)==3 and all(c['passed'] is not None for c in checks)
                passed=complete and all(c['passed'] is True for c in checks)
                key=f'{case}/{backend}/{left}-{right}'
                out['comparisons'][key]=dict(complete=complete,passed=passed,checks=checks)
                if not complete: continue
                fields={'compact':lambda r:r.get('preparation_work',{}).get('operator_compact_seconds'),
                        'preparation':lambda r:r.get('preparation_work',{}).get('operator_preparation_seconds'),
                        'normal':lambda r:r.get('normal_action_seconds')}
                for index,kind in enumerate(('identity','diagonal','schwarz')):
                    fields[kind+'-total']=lambda r,i=index:r['steps'][i].get('total_seconds') if len(r.get('steps',[]))>i else None
                out['performance'][key]={name:improvement([get(r['result']) for r in x],[get(r['result']) for r in y],passed) for name,get in fields.items()}
                for index,kind in enumerate(('identity','diagonal','schwarz')):
                    baseline_bounds=[total_bounds(r['result'],index) for r in x]
                    candidate_bounds=[total_bounds(r['result'],index) for r in y]
                    out['performance'][key][kind+'-total']['scope']='recorded subtotal; old receipts exclude gradient setup'
                    out['performance'][key][kind+'-gradient-inclusive-bounds']=dict(
                        baseline_bounds=baseline_bounds,candidate_bounds=candidate_bounds,
                        conservative_improvement=improvement([t[0] for t in baseline_bounds],[t[1] for t in candidate_bounds],passed))
    out['comparator_sha256']=v.sha(Path(__file__))
    v.write(root/'comparison.json',out); return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work-dir',type=Path,required=True); p.add_argument('--input-dir',type=Path)
    for name in ('baseline-eigen','baseline-spqr','eigen','spqr'): p.add_argument('--'+name,type=Path)
    p.add_argument('--compare',action='store_true')
    args=p.parse_args(); root=args.work_dir.resolve()
    if args.compare: compare(root); return
    builds={key:getattr(args,key.replace('-','_')) for key in ('baseline-eigen','baseline-spqr','eigen','spqr')}
    if args.input_dir is None or any(x is None for x in builds.values()): p.error('All four builds and input-dir are required')
    if (root/'campaign.json').exists(): raise SystemExit('Refusing to replace or resume a prior campaign')
    builds={k:b.resolve() for k,b in builds.items()}
    for build in builds.values(): require_current_build(build)
    v.process_tree_rss(0)  # Fail before starting the budget if monitoring is unavailable.
    root.mkdir(parents=True,exist_ok=True)
    report=dict(baseline_commit=BASELINE,limits=dict(campaign_seconds=1800,process_seconds=300,rss_bytes=4*1024**3),
                builds={k:fingerprint(b) for k,b in builds.items()},runner_sha256=v.sha(Path(__file__)),
                driver_sha256=v.sha(v.ROOT/'tests/experiments/joint_sparse_benchmark.cpp'),
                diagnostic_sha256=v.sha(v.ROOT/'tests/support/JointFixedDiagnostic.hpp'),
                generator_sha256=v.sha(v.ROOT/'tests/support/JointOperatorWorkload.cpp'),
                numeric_environment={k:v.ENV[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
                inputs={},states={},preparation={},runs={},finished=False)
    expected=v.read(v.ROOT/'docs/developer/figures/joint-fixed-actions/baseline.json')
    if expected['commit']!=BASELINE or any(report['builds']['baseline-'+backend]['source_sha256']!=expected['source_sha256'] for backend in ('eigen','spqr')):
        raise ValueError('Baseline source differs from frozen instrumentation-only source')
    for key in builds:
        source=Path(report['builds'][key]['source_root'])
        if v.sha(source/'tests/support/JointOperatorWorkload.cpp')!=report['generator_sha256'] or v.sha(source/'tests/support/JointFixedDiagnostic.hpp')!=report['diagnostic_sha256'] or v.sha(source/'tests/experiments/joint_sparse_benchmark.cpp')!=report['driver_sha256']:
            raise ValueError('Diagnostic wrappers differ')
    start=time.monotonic(); deadline=start+1800
    def save(): report['elapsed_seconds']=time.monotonic()-start; v.write(root/'campaign.json',report)
    def run(command,directory,output):
        try: row=dict(process=v.monitored(command,directory,deadline,seconds=300))
        except KeyboardInterrupt:
            report['interruption']=dict(status='user-cancelled',directory=str(directory),command=list(map(str,command)))
            save(); raise
        if output.is_file(): row['result']=v.read(output)
        print(f'{directory.relative_to(root)}: {row["process"]["status"]}',flush=True)
        return row
    save()
    for case in CASES:
        if time.monotonic()>=deadline:
            report['preparation'][case]=dict(process=dict(status='not-run-budget')); save(); continue
        if case=='chain-8': command=['synthetic','chain','8','fixed']; files=[]
        elif case=='heterogeneous-168':
            folder=unpack(v.ROOT/'tests/fixtures/joint_component/catalog.json',case,args.input_dir)
            command=['fixture',folder,'first-stage-float32']; files=sorted(folder.iterdir())
        else:
            size=case.split('-')[1]; folder=args.input_dir/case
            files=[folder/'input.cif',folder/'input.map',folder/'widths.json']
            if not all(f.is_file() for f in files):
                generation=run([builds['baseline-eigen']/'bin/joint_validation','single',size,folder],root/'preparation'/case/'generate',root/'unused.json')
                report['preparation'][case+'/generate']=generation; save()
                if generation['process']['status']!='completed': continue
                prepared=run([builds['baseline-eigen']/'bin/joint_sparse_benchmark','prepare',*files],root/'preparation'/case/'initialize',root/'unused.json')
                report['preparation'][case+'/initialize']=prepared; save()
                if prepared['process']['status']!='completed': continue
            command=['initial',*files]
        report['inputs'][case]={f.name:v.sha(f) for f in files if f.is_file()}
        state=root/'states'/f'{case}.json'; state.parent.mkdir(exist_ok=True)
        frozen=run([builds['baseline-eigen']/'bin/joint_sparse_benchmark',*command,state,'--fixed','freeze'],root/'preparation'/case/'freeze',state)
        report['preparation'][case]=frozen; save()
        if frozen['process']['status']!='completed': continue
        report['states'][case]=dict(sha256=v.sha(state),path=str(state))
        stopped=set()
        for repeat in range(3):
            for backend in ('eigen','spqr'):
                for mode in MODES[repeat:]+MODES[:repeat]:
                    key=f'{case}/{backend}/{mode}'; runs=report['runs'].setdefault(key,[])
                    if key in stopped: continue
                    d=root/'runs'/key/str(repeat+1); output=d/'result.json'
                    build=builds[('baseline-' if mode=='A' else '')+backend]
                    row=run([build/'bin/joint_sparse_benchmark',*command,output,'--fixed','normal' if mode=='C' else 'composed','--state',state,'--resources'],d,output)
                    row['state_sha256']=v.sha(state); runs.append(row)
                    if row['process']['status']!='completed': stopped.add(key)
                    save()
    report['finished']=True; save(); compare(root)


if __name__=='__main__': main()
