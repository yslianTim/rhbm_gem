"""Bounded offline validation. Production algorithms and convergence policies are unchanged."""
from __future__ import annotations
import argparse
import collections
import hashlib
import json
import math
import os
import platform
import re
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SEEDS = list(range(20260921, 20260941))
ENV = dict(os.environ, **{k: '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')})
LIMITS = dict(process_seconds=600, rss_bytes=4*1024**3, sample_seconds=.1, stage_seconds={'a':1200,'b':1200,'c':4800}, total_seconds=7200)


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def digest(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',',':'), allow_nan=False).encode()).hexdigest()


def process_tree_rss(pid):
    result = subprocess.run(['ps','-axo','pid=,ppid=,rss='], capture_output=True, text=True, check=True)
    rows = [tuple(map(int, line.split())) for line in result.stdout.splitlines()]
    descendants = {pid}
    while True:
        found = {p for p, parent, _ in rows if parent in descendants}
        if found <= descendants:
            break
        descendants |= found
    return sum(rss*1024 for p, _, rss in rows if p in descendants)


def monitored(command, directory, deadline, *, rss_limit=LIMITS['rss_bytes'], seconds=600):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    available = deadline-start
    if available <= 0:
        return dict(status='not-run-budget', command=list(map(str,command)), wall_seconds=0)
    expires = min(deadline, start+seconds)
    peak = 0; samples = 0; gap = 0.; last = start; next_update = start+30
    record = dict(command=list(map(str,command)), status='running', sample_interval_seconds=.1)
    # The parent owns the entire process group, including time and all CLI children.
    timed = ['/usr/bin/time', '-l' if sys.platform=='darwin' else '-v', *map(str,command)]
    with (directory/'stdout.txt').open('w') as stdout, (directory/'stderr.txt').open('w') as stderr:
        process = subprocess.Popen(timed, stdout=stdout, stderr=stderr, env=ENV, start_new_session=True)
        try:
            while process.poll() is None:
                now = time.monotonic(); gap=max(gap,now-last); last=now
                peak=max(peak,process_tree_rss(process.pid)); samples+=1
                reason = 'rss-limit' if peak > rss_limit else ('time-limit' if now >= expires else None)
                if reason:
                    record['status']=reason
                    os.killpg(process.pid, signal.SIGTERM)
                    try: process.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL); process.wait()
                    break
                if now >= next_update:
                    print(f"running {directory.name}: {now-start:.0f}s, sampled peak {peak/1024**2:.0f} MiB",flush=True)
                    next_update=now+30
                time.sleep(max(0,min(.1-(time.monotonic()-now),expires-time.monotonic())))
        except BaseException:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL); process.wait()
            raise
    if record['status']=='running': record['status']='completed' if process.returncode==0 else 'process-failure'
    text=(directory/'stderr.txt').read_text()
    match = re.search(r'(\d+)\s+maximum resident set size',text) if sys.platform=='darwin' else re.search(r'Maximum resident set size \(kbytes\): (\d+)',text)
    record.update(wall_seconds=time.monotonic()-start, exit_code=process.returncode,
                  sampled_tree_peak_rss_bytes=peak, samples=samples, maximum_sampling_gap_seconds=gap,
                  os_process_peak_rss_bytes=int(match[1])*(1 if sys.platform=='darwin' else 1024) if match else None)
    if record['status']=='completed' and (record['os_process_peak_rss_bytes'] or 0)>rss_limit:
        record['status']='rss-limit-observed-after-exit'
    write(directory/'process.json',record)
    return record


def run_a(args, deadline):
    root=args.work_dir/args.stage_dir; root.mkdir(exist_ok=True)
    records=[]; hashes=[]
    for start in range(3):
        output=root/f'start-{start}.json'
        run=monitored([args.executable,'weak',start,output],root/f'process-{start}',deadline)
        run['start']=start
        if output.exists():
            data=read(output); h=digest(data['snapshot']); hashes.append(h)
            data['snapshot_sha256']=h; write(output,data)
            run.update(audit_complete=data.get('audit_complete',False),case_complete=data.get('audit_complete',False) or data['outcome']['assembled_state'] is None,snapshot_sha256=h,output_sha256=sha(output))
        records.append(run)
    if len(set(hashes))>1: raise AssertionError('Weak snapshots changed across starts')
    write(root/'receipt.json',dict(runs=records, complete=all(x.get('case_complete') for x in records)))


AXIS=np.arange(25)*.5-6
Z,Y,X=np.meshgrid(AXIS,AXIS,AXIS,indexing='ij')
TARGET=(X*X+Y*Y+Z*Z<=6.25).ravel()
TRUTH=np.array([[2.,.5,.2],[2.3,.55,.15]])


def clean_map(distance, shift):
    values=np.zeros(X.shape)
    for x,a,b,c in [(0,2,.5,.2),(distance+shift,2.3,.55,.15)]:
        square=(X-x)**2+Y**2+Z**2
        radius=np.sqrt(square); keep=square<=6.25
        gaussian=np.exp(-square/(2*b*b))/(2*math.pi*b*b)**1.5
        charge=np.zeros(X.shape); center=radius<1e-5
        charge[center]=math.sqrt(2/math.pi)/b
        nonzero=~center
        charge[nonzero]=np.fromiter((math.erf(float(t)/(math.sqrt(2)*b))/float(t) for t in radius[nonzero]),float)
        values[keep]+=(a*gaussian+c*charge)[keep]
    return values.ravel()


def noise_field(seed, kind):
    rng=np.random.Generator(np.random.PCG64(seed))
    if kind=='iid': return rng.normal(size=(25,25,25)).ravel()
    field=rng.normal(size=(31,31,31))
    kernel=np.exp(-.5*np.arange(-3,4,dtype=float)**2)
    kernel/=np.sqrt(np.sum(kernel*kernel)) # unit marginal variance, no sample normalization
    for axis in range(3): field=np.apply_along_axis(lambda row:np.convolve(row,kernel,mode='valid'),axis,field)
    return field.ravel()


def conditions():
    for distance in (1.2,.6):
        yield distance,0.,'none',0.,None,['fixed','production']
        for kind in ('iid','correlated'):
            for sigma in (.01,.05,.1):
                for seed in SEEDS:
                    yield distance,0.,kind,sigma,seed,['fixed','production'] if sigma==.05 else ['fixed']
        for shift in (.05,.15):
            modes=['fixed','production'] if shift==.15 else ['fixed']
            yield distance,shift,'none',0.,None,modes
            for seed in SEEDS: yield distance,shift,'iid',.05,seed,modes


def metrics(outcome, prediction, values, clean, rows):
    result=dict(state_available=outcome['assembled_state'] is not None, runtime_convergence=outcome['runtime_convergence'],
                initialization_valid=outcome['initialization']['valid'], stop_reasons=[c['stop_reason'] for c in outcome['components']],
                failure_reasons=[e['name']+':'+e['status'] for c in outcome['components'] for e in c['evidence'] if e['status'] in ('failed','unavailable')],
                parameters=None, errors=None, prediction_rmse=None, residual_neighbor_correlation=None)
    if result['state_available']:
        state=outcome['assembled_state']; ac=np.array(state['ac']).reshape(2,2)
        parameters=np.column_stack([ac[:,0],state['b'],ac[:,1]])
        prediction=np.array(prediction); residual=prediction-values[rows]
        result.update(parameters=parameters.tolist(),errors=(parameters-TRUTH).tolist(),
                      prediction_rmse=float(np.sqrt(np.mean((prediction-clean[rows])**2))), residual_rmse=float(np.sqrt(np.mean(residual**2))))
        lookup={int(row):i for i,row in enumerate(rows)}; left=[]; right=[]
        for row,i in lookup.items():
            for stride,coordinate in ((1,row%25),(25,(row//25)%25),(625,row//625)):
                if coordinate<24 and row+stride in lookup: left.append(residual[i]); right.append(residual[lookup[row+stride]])
        if len(left)>2 and np.std(left)>0 and np.std(right)>0:
            result['residual_neighbor_correlation']=float(np.corrcoef(left,right)[0,1])
    return result


def run_b(args, deadline):
    root=args.work_dir/args.stage_dir;root.mkdir(exist_ok=True)
    records=[]; maps={}
    for index,(distance,shift,kind,sigma,seed,modes) in enumerate(conditions()):
        case=root/f'{index:03d}';case.mkdir(exist_ok=True)
        record=dict(index=index,distance=distance,shift=shift,noise=kind,sigma_fraction=sigma,seed=seed,modes=modes)
        if time.monotonic()>=deadline:
            record['status']='not-run-budget';records.append(record);continue
        key=(distance,shift)
        if key not in maps: maps[key]=clean_map(*key)
        clean=maps[key]; reference=maps.setdefault((distance,0.),clean_map(distance,0.))
        scale=float(np.sqrt(np.mean(reference[TARGET]**2)))
        values=clean.copy() if seed is None else clean+sigma*scale*noise_field(seed,kind)
        spec=dict(distance=distance,shift=shift,values=values.tolist(),clean=clean.tolist(),modes=modes)
        write(case/'input.json',spec)
        record.update(input_sha256=sha(case/'input.json'),observation_sha256=hashlib.sha256(values.astype('<f8').tobytes()).hexdigest(),noise_sigma=sigma*scale)
        run=monitored([args.executable,'stat',case/'input.json',case/'output.json'],case,deadline)
        record.update(status=run['status'],process=run)
        if run['status']=='completed':
            result=read(case/'output.json'); rows=np.array(result['row_ids'],dtype=int)
            assert np.array_equal(rows,np.flatnonzero(TARGET)), 'Changed target domain'
            record['generator_max_error']=result['generator_max_error'];record['census']=result['census']
            record['results']={mode:metrics(result[mode]['outcome'],result[mode]['prediction'],values,clean,rows) for mode in modes}
            record['software']=result['fixed']['outcome']['metadata']['software']
            record['output_sha256']=sha(case/'output.json')
        records.append(record)
        write(root/'receipt.json',dict(seeds=SEEDS,rng='numpy.PCG64',numpy_version=np.__version__,records=records))
        if index%40==0: print(f'B {index+1}/326 input cases',flush=True)
    write(root/'receipt.json',dict(seeds=SEEDS,rng='numpy.PCG64',numpy_version=np.__version__,records=records))


def resource_run(args, case, model, map_path, repetition, deadline):
    root=args.work_dir/args.stage_dir/case/f'run-{repetition}';root.mkdir(parents=True,exist_ok=True)
    database=root/'result.sqlite';started=time.monotonic(); pipeline_deadline=min(deadline,started+600)
    analysis=monitored([args.cli,'potential_analysis','--estimator','joint-components','-a',model,'-m',map_path,
        '-d',database,'-k','validation','--exclude-hydrogen','true','--asymmetry','false','--only-backbone','false',
        '--map-normalization','false','-j','1','-v','0'],root/'analysis',pipeline_deadline)
    result=dict(repetition=repetition,analysis=analysis,status=analysis['status'])
    if analysis['status']=='completed':
        export=monitored([args.cli,'result_dump','--printer','joint','-d',database,'-k','validation','-o',root,'-v','0'],root/'export',pipeline_deadline)
        result['export']=export;result['status']=export['status']
        if export['status']=='completed':
            outcome=read(root/'joint_result_validation.json')
            assert outcome['schema_version'] in (3,4) and outcome['metadata']['map_normalization']['divisor']==1
            csv=(root/'joint_atoms_validation.csv').read_text().splitlines()
            assert len(csv)==len(outcome['atom_ids'])+1
            result.update(runtime_convergence=outcome['runtime_convergence'], state_available=outcome['assembled_state'] is not None,
                          initialization_valid=outcome['initialization']['valid'],
                          invalid_initialization_reasons=dict(collections.Counter(a['reason'] for a in outcome['initialization']['atoms'] if a['reason']!='valid-width')),
                          stop_reasons=dict(collections.Counter(c['stop_reason'] for c in outcome['components'])),
                          costs=outcome['costs'],software=outcome['metadata']['software'],
                          component_statuses=dict(collections.Counter(c['runtime_convergence'] for c in outcome['components'])),
                          failure_reasons=dict(collections.Counter(e['name']+':'+e['status'] for c in outcome['components'] for e in c['evidence'] if e['status'] in ('failed','unavailable'))),
                          input_hashes={k:outcome['metadata'][k] for k in ('model_sha256','map_sha256')},
                          output_bytes=sum(p.stat().st_size for p in root.iterdir() if p.is_file()))
    # External measurement ends when export exits, before parsing its results.
    result['total_command_seconds']=analysis['wall_seconds']+result.get('export',{}).get('wall_seconds',0)
    if 'costs' in result:
        estimator=sum(result['costs'][key] for key in ('construction_seconds','initialization_seconds','search_seconds','assessment_seconds','assembly_seconds'))
        result['outside_estimator_seconds']=result['total_command_seconds']-estimator
    stages=[analysis]+([result['export']] if 'export' in result else [])
    result['sampled_tree_peak_rss_bytes']=max(s.get('sampled_tree_peak_rss_bytes',0) for s in stages)
    result['os_process_peak_rss_bytes']=max((s['os_process_peak_rss_bytes'] for s in stages if s.get('os_process_peak_rss_bytes') is not None),default=None)
    write(root/'receipt.json',result)
    return result


def run_c(args, deadline):
    root=args.work_dir/args.stage_dir;root.mkdir(exist_ok=True)
    cases=[('catalogue',1000),('catalogue',10000),('catalogue',37406),('single',128),('single',512),('user',2192),('single',2167),('multi',2167)]
    report={};blocked=set();inputs={}
    for kind,count in cases:
        name=f'{kind}-{count}';directory=root/name;directory.mkdir(exist_ok=True)
        row=dict(kind=kind,count=count,runs=[]);report[name]=row
        if kind in blocked or time.monotonic()>=deadline:
            row['status']='not-run-larger-family-limit' if kind in blocked else 'not-run-budget';continue
        if kind=='user': model,map_path=args.model,args.map
        else:
            generation=monitored([args.executable,kind,count,directory],directory/'generation',deadline)
            row['generation']=generation
            if generation['status']!='completed':row['status']='generation-'+generation['status'];continue
            model,map_path=directory/'input.cif',directory/'input.map'
        inputs[name]=(model,map_path)
        row['input_hashes']=dict(model=sha(model),map=sha(map_path))
        inspection=monitored([args.executable,'inspect',model,map_path,directory/'census.json'],directory/'inspect',deadline)
        row['inspection']=inspection
        if inspection['status']!='completed': row['status']='inspection-'+inspection['status'];continue
        row['census']=read(directory/'census.json')
        write(root/'receipt.json',report)
        if kind=='user':
            assert row['census']==dict(targets=1559,contributors=2192,rows=262801,memberships=866216,component_sizes=[2167,25],catalogue=37406,model_atoms=70190,map_size=[258]*3)
        elif kind=='catalogue':
            controls=[r['census'] for r in report.values() if r.get('kind')=='catalogue' and 'census' in r]
            assert all({k:v for k,v in c.items() if k not in ('catalogue','model_atoms')}=={k:v for k,v in controls[0].items() if k not in ('catalogue','model_atoms')} for c in controls)
        elif kind=='single': assert row['census']['component_sizes']==[count]
        else: assert max(row['census']['component_sizes'])<=32 and sum(row['census']['component_sizes'])==2167
        run=resource_run(args,name,model,map_path,1,deadline);row['runs'].append(run);row['status']=run['status']
        if run['status'].startswith('rss-limit') or run['status']=='time-limit':blocked.add(kind)
        write(root/'receipt.json',report)
        print(f'C {name}: {run["status"]}, {run["total_command_seconds"]:.1f}s',flush=True)
    for repetition in (2,3):
        for name,row in report.items():
            if row['status']!='completed' or time.monotonic()>=deadline:continue
            run=resource_run(args,name,*inputs[name],repetition,deadline);row['runs'].append(run)
            if run['status']!='completed':row['status']=run['status']
            write(root/'receipt.json',report)
            print(f'C {name} repeat {repetition}: {run["status"]}',flush=True)
    write(root/'receipt.json',report)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resume-stage','--resume-c',action='store_true',help='Preserve earlier evidence and charge only the remaining stage budget after a diagnosed tooling/product fix.')
    parser.add_argument('--stage-dir','--resource-dir',dest='stage_dir')
    parser.add_argument('--stage',choices=('a','b','c'),required=True)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--cli',type=Path,default=ROOT/'build/joint-v1-default/bin/RHBM-GEM')
    parser.add_argument('--model',type=Path,default=Path('/Users/yslian/data/6z6u.cif'))
    parser.add_argument('--map',type=Path,default=Path('/Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map'))
    args=parser.parse_args()
    args.stage_dir=args.stage_dir or args.stage
    for key in ('work_dir','executable','cli','model','map'):setattr(args,key,getattr(args,key).resolve())
    args.work_dir.mkdir(parents=True,exist_ok=True)
    receipt=args.work_dir/'campaign.json'
    campaign=read(receipt) if receipt.exists() else dict(limits=LIMITS,platform=platform.platform(),stages={},base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    previous=campaign['stages'].get(args.stage)
    if previous and not (args.resume_stage and not (args.work_dir/args.stage_dir).exists()):
        raise SystemExit('Stage already attempted; recovery requires a fresh stage directory and retains its consumed budget.')
    process_tree_rss(os.getpid()) # fail before starting experiments if monitoring permission is unavailable
    used=sum(r['wall_seconds'] for r in campaign['stages'].values())
    prior_seconds=previous['wall_seconds'] if previous else 0
    budget=min(LIMITS['stage_seconds'][args.stage]-prior_seconds,LIMITS['total_seconds']-used)
    start=time.monotonic()
    campaign['stages'][args.stage]=dict(status='running',wall_seconds=prior_seconds+budget,started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),budget_seconds=budget)
    if previous:campaign['stages'][args.stage]['previous_attempt']=previous
    campaign['stages'][args.stage]['stage_directory']=args.stage_dir
    campaign['stages'][args.stage]['resource_directory']=args.stage_dir if args.stage=='c' else None
    write(receipt,campaign) # conservatively charge the full budget if the orchestrator is interrupted
    try:
        globals()['run_'+args.stage](args,start+budget)
        campaign['stages'][args.stage]['status']='finished'
    except BaseException:
        campaign['stages'][args.stage]['status']='interrupted-or-failed'
        raise
    finally:
        campaign['stages'][args.stage]['wall_seconds']=prior_seconds+time.monotonic()-start
        sources=[ROOT/'tests/experiments/joint_validation.cpp',Path(__file__),ROOT/'tests/support/JointPrecisionAudit.cpp',ROOT/'tests/support/JointPrecisionAudit.hpp']
        campaign['stages'][args.stage]['source_hashes']={str(p.relative_to(ROOT)):sha(p) for p in sources}
        campaign['stages'][args.stage]['executable_sha256']=sha(args.executable)
        campaign['stages'][args.stage]['cli_sha256']=sha(args.cli)
        write(receipt,campaign)
    print(json.dumps(campaign['stages'][args.stage],indent=2),flush=True)

if __name__=='__main__':main()
