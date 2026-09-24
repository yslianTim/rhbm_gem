"""One serial 60-minute campaign: fixed (5), search (45), rank (10). Never resumes/promotes."""
from __future__ import annotations
import argparse
import statistics
import sys
import time
from pathlib import Path
import joint_validation as v
from joint_operator_validation import fingerprint
from joint_search_validation import finite, require_current_build, scientific_parity
from joint_fixed_validation import parity as fixed_parity
from joint_runtime_support import unpack

BACKENDS=('eigen','spqr')
STAGE_SECONDS={'fixed':300,'search':2700,'rank':600}
RESOURCE_STOPS={'time-limit','rss-limit','rss-limit-observed-after-exit','not-run-budget'}


def search_schedule():
    for case in ('single-128','heterogeneous-168'):
        for backend in BACKENDS:
            for kind in ('legacy','identity','diagonal','schwarz'): yield backend,case,kind,1
    for repeat in range(1,4):
        for index,backend in enumerate(BACKENDS if repeat%2 else BACKENDS[::-1]):
            for kind in (('legacy','schwarz') if (repeat+index)%2 else ('schwarz','legacy')):
                yield backend,'single-512',kind,repeat


def completed(row):
    return row.get('process',{}).get('status')=='completed' and row.get('result',{}).get('stage')=='complete'


def qualified(row):
    r=row.get('result',{})
    return completed(row) and r.get('search_completed') is True and (r.get('returned_assessment') or {}).get('runtime_convergence')=='passed' and scientific_parity(r,r)['passed']


def status(checks):
    failed=any(x is False for x in checks); incomplete=not checks or any(x is None for x in checks)
    return dict(status='failed' if failed else 'incomplete' if incomplete else 'passed',
                failed=failed,incomplete=incomplete,exit_code=1 if failed else 3 if incomplete else 0)


def fixed_check(row,kinds):
    if not completed(row): return None
    r=row['result']; steps=r.get('steps',[])
    if r.get('total_includes_gradient') is not True or r.get('valid') is not True or [s.get('kind') for s in steps]!=list(kinds): return False
    for step in steps:
        fields=['fixed_step_wall_seconds','total_seconds','preparation_seconds','gradient_seconds','metric_seconds','partition_seconds','build_seconds','solve_seconds']
        if step.get('valid') is not True or any(not finite(step.get(k)) or step[k]<0 for k in fields): return False
        if not finite(step.get('true_residual')) or step['true_residual']>1e-10 or not finite(step.get('predicted')): return False
        if step['fixed_step_wall_seconds']+1e-9<step['total_seconds']: return False
    return r.get('work',{}).get('reference_solves')==0 and r.get('work',{}).get('derivative_preparations')==0


def rank_check(proto,oracles):
    if not completed(proto): return dict(passed=None,reason='incomplete-prototype')
    r=proto['result']; evidence=r.get('rank_result',{}); decision=evidence.get('status')
    if r.get('rank_compact_extractions')!=0 or any(r.get('work',{}).get(k)!=0 for k in ('free_design_svds','reference_solves','derivative_preparations')):
        return dict(passed=False,reason='forbidden-dense-work')
    if decision=='unavailable': return dict(passed=None,reason=evidence.get('reason'))
    if decision not in ('full-rank','deficient'): return dict(passed=False,reason='invalid-rank-decision')
    fields=('minimum_lower','maximum_lower','maximum_upper','threshold_lower','threshold_upper')
    if any(not finite(evidence.get(k)) for k in fields): return dict(passed=False,reason='nonfinite-rank-bounds')
    if any(type(evidence.get(k)) is not int for k in ('rank_lower','rank_upper')) or not 0<=evidence['rank_lower']<=evidence['rank_upper']<=r.get('free_columns',-1): return dict(passed=False,reason='invalid-rank-bounds')
    if decision=='full-rank' and (evidence['rank_lower']!=r['free_columns'] or evidence['minimum_lower']<=evidence['threshold_upper'] or any(not finite(evidence.get(k)) for k in ('reconstruction_error','orthogonal_minimum'))): return dict(passed=False,reason='missing-full-rank-certificate')
    if any(not completed(o) for o in oracles): return dict(passed=None,reason='incomplete-oracle')
    for oracle in oracles:
        data=oracle['result']; actual=data.get('rank_result',{}); values=actual.get('singular_values',[])
        if data.get('input_sha256')!=r.get('input_sha256') or data.get('free_columns')!=r.get('free_columns'): return dict(passed=False,reason='rank-input-fingerprint')
        if actual.get('valid') is not True or not values or not all(finite(x) for x in values): return dict(passed=None,reason='unavailable-oracle')
        rank=actual.get('rank')
        if type(rank) is not int or not evidence['rank_lower']<=rank<=evidence['rank_upper']: return dict(passed=False,reason='rank-oracle-disagreement')
        if decision=='full-rank' and rank!=r['free_columns']: return dict(passed=False,reason='false-full-rank')
        if decision=='deficient' and rank>=r['free_columns']: return dict(passed=False,reason='false-deficiency')
        if not evidence['threshold_lower']<=actual['threshold']<=evidence['threshold_upper']: return dict(passed=False,reason='threshold-oracle-disagreement')
    return dict(passed=True,reason='bounded-oracle-match' if oracles else 'certified-rank-without-dense-oracle')


def compare(root):
    report=v.read(root/'campaign.json'); runs=report['runs']; checks=[]
    source_hashes=[report.get('builds',{}).get(b,{}).get('source_sha256') for b in BACKENDS]
    if not all(isinstance(h,str) and h for h in source_hashes) or len(set(source_hashes))!=1: raise ValueError('Missing or different source fingerprints')
    out=dict(fixed={},search={},performance={},rank={},promote=False,rank_search_backend='dense')
    for case,kinds in (('chain-8',('identity','diagonal','schwarz')),('single-128',('schwarz',))):
        rows=[runs.get(f'fixed/{b}/{case}',{}) for b in BACKENDS]
        values=[fixed_check(r,kinds) for r in rows]
        match=fixed_parity(rows[0]['result'],rows[1]['result'],kinds)['passed'] if all(x is True for x in values) else None
        expected=report.get('states',{}).get(case)
        identity=all(r.get('state_sha256')==expected for r in rows) if expected else None
        out['fixed'][case]=dict(checks=values,cross_backend_parity=match,state_fingerprint=identity); checks.extend([*values,identity])
        if match is not None: checks.append(match)
    fixed_end=len(checks)
    for backend,case,kind,repeat in search_schedule():
        if kind=='legacy': continue
        left=runs.get(f'search/{backend}/{case}/legacy/{repeat}',{})
        right=runs.get(f'search/{backend}/{case}/{kind}/{repeat}',{})
        key=f'{backend}/{case}/{kind}/{repeat}'
        if not completed(left) or not completed(right): result=dict(passed=None,reason='incomplete-measurement')
        elif left.get('input_sha256')!=right.get('input_sha256'): result=dict(passed=False,reason='input-fingerprint')
        else: result=scientific_parity(left['result'],right['result'])
        result['qualified_endpoints']=qualified(left) and qualified(right)
        out['search'][key]=result
        checks.append(False if result['passed'] is False else True if result['passed'] is True and result['qualified_endpoints'] else None)
    for backend in BACKENDS:
        medians={}
        for kind in ('legacy','schwarz'):
            rows=[runs.get(f'search/{backend}/single-512/{kind}/{i}',{}) for i in range(1,4)]
            values=[r.get('result',{}).get('search_seconds') for r in rows]
            medians[kind]=statistics.median(values) if all(qualified(r) for r in rows) and all(finite(x) and x>0 for x in values) else None
        paired=all(out['search'][f'{backend}/single-512/schwarz/{i}']['passed'] is True for i in range(1,4))
        passed=medians['schwarz']<=1.1*medians['legacy'] if paired and all(x is not None for x in medians.values()) else None
        out['performance'][backend]=dict(medians=medians,within_ten_percent=passed); checks.append(passed)
    search_end=len(checks)
    for size in (128,512,2000):
        for topology in ('chain','cube'):
            case=f'{topology}-{size}'
            oracles=[runs.get(f'rank/{b}/{case}/oracle',{}) for b in BACKENDS] if size<=512 else []
            result=rank_check(runs.get(f'rank/spqr/{case}/prototype',{}),oracles)
            out['rank'][case]=result; checks.append(result['passed'])
    out['sections']={name:status(values) for name,values in (('fixed',checks[:fixed_end]),('search',checks[fixed_end:search_end]),('rank',checks[search_end:]))}
    out.update(status(checks)); out['comparator_sha256']=v.sha(Path(__file__))
    v.write(root/'comparison.json',out); return out


def validate_builds(builds):
    for b in builds.values(): require_current_build(b)
    records={b:fingerprint(path) for b,path in builds.items()}
    if len({r['source_sha256'] for r in records.values()})!=1: raise ValueError('Different source fingerprints')
    keys=('CMAKE_BUILD_TYPE','CMAKE_CXX_COMPILER','CMAKE_CXX_FLAGS','CMAKE_CXX_FLAGS_RELEASE','BUILD_TESTING','RHBM_GEM_ENABLE_UMAP','RHBM_GEM_ROOT_MODE','RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS')
    configs={b:{line.split(':',1)[0]:line.split('=',1)[1] for line in r['cache'].splitlines() if ':' in line and '=' in line and not line.startswith('//')} for b,r in records.items()}
    if any(configs['eigen'].get(k)!=configs['spqr'].get(k) for k in keys): raise ValueError('Different build settings')
    for b,c in configs.items():
        if c.get('RHBM_GEM_JOINT_SPARSE_BACKEND')!=b.upper() or c.get('CMAKE_BUILD_TYPE')!='Release': raise ValueError('Wrong backend/build type')
    return records


def campaign(root,inputs,builds):
    if root.exists(): raise ValueError('Refusing to replace or resume an existing campaign directory')
    records=validate_builds(builds); v.process_tree_rss(0)
    root.mkdir(parents=True)
    report=dict(schema_version=1,base_commit='5f61bb6e',builds=records,inputs={},states={},runs={},stages={},finished=False,
                limits=dict(total_seconds=3600,stages=STAGE_SECONDS,search_process_seconds=600,other_process_seconds=180,rss_bytes=4*1024**3),
                numeric_environment={k:v.ENV[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
                wrappers={str(p.relative_to(v.ROOT)):v.sha(p) for p in [Path(__file__),v.ROOT/'tests/integration/joint_validation.py',v.ROOT/'tests/integration/joint_search_validation.py',v.ROOT/'tests/integration/joint_fixed_validation.py',v.ROOT/'tests/integration/joint_operator_validation.py',v.ROOT/'tests/integration/joint_runtime_support.py',v.ROOT/'tests/experiments/joint_sparse_benchmark.cpp',*sorted((v.ROOT/'tests/support').glob('Joint*Diagnostic.hpp')),v.ROOT/'tests/support/JointOperatorWorkload.cpp']})
    for case in ('chain-8','single-128'):
        for b in BACKENDS: report['runs'][f'fixed/{b}/{case}']=dict(process=dict(status='not-run'))
    for b,c,k,r in search_schedule(): report['runs'][f'search/{b}/{c}/{k}/{r}']=dict(process=dict(status='not-run'))
    for size in (128,512,2000):
        for topology in ('chain','cube'):
            for b,kind in [('spqr','prototype')]+([(b,'oracle') for b in BACKENDS] if size<=512 else []):
                report['runs'][f'rank/{b}/{topology}-{size}/{kind}']=dict(process=dict(status='not-run'))
    start=time.monotonic(); total_deadline=start+3600; current_key=None
    def save():
        report['elapsed_seconds']=time.monotonic()-start; v.write(root/'campaign.json',report)
    def run(key,command,deadline,seconds=180,output=None):
        nonlocal current_key
        current_key=key; directory=root/'runs'/key; output=output or directory/'result.json'
        row=dict(process=dict(status='running')); report['runs'][key]=row; save()
        row['process']=v.monitored(command,directory,deadline,seconds=seconds)
        if output.is_file(): row['result']=v.read(output)
        print(f'{key}: {row["process"]["status"]}',flush=True); save(); current_key=None
        return row
    prepared={}
    def prepare(case,deadline):
        if case in prepared: return prepared[case]
        if time.monotonic()>=deadline: return None
        if case=='chain-8': command=['synthetic','chain',8,'fixed']; files=[]
        elif case=='heterogeneous-168':
            folder=unpack(v.ROOT/'tests/fixtures/joint_component/catalog.json',case,inputs)
            command=['fixture',folder,'first-stage-float32']; files=sorted(folder.iterdir())
        else:
            folder=inputs/case; files=[folder/'input.cif',folder/'input.map',folder/'widths.json']
            if not all(p.is_file() for p in files):
                generated=run(f'prepare/{case}/generate',[builds['eigen']/'bin/joint_validation','single',case.split('-')[1],folder],deadline)
                if generated['process']['status']!='completed': return None
                initialized=run(f'prepare/{case}/initialize',[builds['eigen']/'bin/joint_sparse_benchmark','prepare',*files],deadline)
                if initialized['process']['status']!='completed': return None
            command=['initial',*files]
        report['inputs'][case]={p.name:v.sha(p) for p in files}; prepared[case]=command; save(); return command
    try:
        save()
        for stage,budget in STAGE_SECONDS.items():
            began=time.monotonic(); deadline=min(total_deadline,began+budget)
            report['stages'][stage]=dict(start_seconds=began-start,budget_seconds=budget); save()
            if stage=='fixed':
                for case in ('chain-8','single-128'):
                    command=prepare(case,deadline)
                    if command is None: continue
                    state=root/'states'/f'{case}.json'; state.parent.mkdir(exist_ok=True)
                    row=run(f'freeze/{case}',[builds['eigen']/'bin/joint_sparse_benchmark',*command,state,'--fixed','freeze'],deadline,output=state)
                    if row['process']['status']!='completed': continue
                    report['states'][case]=v.sha(state)
                    for b in BACKENDS:
                        key=f'fixed/{b}/{case}'; output=root/'runs'/key/'result.json'
                        args=[builds[b]/'bin/joint_sparse_benchmark',*command,output,'--fixed','normal','--state',state,'--resources']
                        if case=='single-128': args+=['--fixed-preconditioner','schwarz']
                        row=run(key,args,deadline,output=output); row['state_sha256']=v.sha(state); save()
            elif stage=='search':
                stopped=set()
                for b,case,kind,repeat in search_schedule():
                    key=f'search/{b}/{case}/{kind}/{repeat}'
                    if (b,case,kind) in stopped:
                        report['runs'][key]['process']['status']='not-run-group-stopped'; continue
                    command=prepare(case,deadline)
                    if command is None: continue
                    output=root/'runs'/key/'result.json'
                    row=run(key,[builds[b]/'bin/joint_sparse_benchmark',*command,output,'--search',kind,'--resources'],deadline,600,output)
                    row['input_sha256']=v.digest(report['inputs'][case]); save()
                    if row['process']['status']!='completed': stopped.add((b,case,kind))
            else:
                stopped={}
                for size in (128,512,2000):
                    for topology in ('chain','cube'):
                        for b,kind in [('spqr','prototype')]+([(b,'oracle') for b in BACKENDS] if size<=512 else []):
                            key=f'rank/{b}/{topology}-{size}/{kind}'
                            if topology in stopped and size>stopped[topology]:
                                report['runs'][key]['process']['status']='not-run-topology-stopped'; continue
                            output=root/'runs'/key/'result.json'
                            row=run(key,[builds[b]/'bin/joint_sparse_benchmark','synthetic',topology,size,'rank' if kind=='prototype' else 'rank-oracle',output,'--resources'],deadline,output=output)
                            reason=row.get('result',{}).get('rank_result',{}).get('reason','')
                            if row['process']['status'] in RESOURCE_STOPS or reason in ('rank-time-budget','rank-work-budget','rank-memory-budget'): stopped.setdefault(topology,size)
            report['stages'][stage]['elapsed_seconds']=time.monotonic()-began
            for key,row in report['runs'].items():
                if key.startswith(stage+'/') and row['process']['status']=='not-run': row['process']['status']='not-run-budget-or-preparation'
            save()
        report['finished']=True
    except KeyboardInterrupt:
        report['interruption']='user-cancelled'
        if current_key: report['runs'][current_key]['process']['status']='user-cancelled'
        for row in report['runs'].values():
            if row['process']['status']=='not-run': row['process']['status']='not-run-user-cancelled'
        raise
    finally: save()
    # A mutation while measuring invalidates the campaign instead of comparing mixed binaries.
    if any(v.sha(v.ROOT/path)!=digest for path,digest in report['wrappers'].items()): raise ValueError('Wrappers changed during campaign')
    if records!=validate_builds(builds): raise ValueError('Source or binaries changed during campaign')
    return compare(root)['exit_code']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True); parser.add_argument('--input-dir',type=Path)
    parser.add_argument('--eigen',type=Path); parser.add_argument('--spqr',type=Path); parser.add_argument('--compare',action='store_true')
    args=parser.parse_args(); root=args.work_dir.resolve()
    if args.compare: return compare(root)['exit_code']
    if any(x is None for x in (args.input_dir,args.eigen,args.spqr)): parser.error('input-dir and both builds required')
    return campaign(root,args.input_dir.resolve(),{b:getattr(args,b).resolve() for b in BACKENDS})


if __name__=='__main__':
    try: raise SystemExit(main())
    except KeyboardInterrupt: raise SystemExit(3)
    except (OSError,ValueError,KeyError,TypeError,IndexError) as error:
        print(f'Invalid campaign input: {error}',file=sys.stderr); raise SystemExit(2)
