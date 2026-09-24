"""PR0/PR1 bounded baseline, operator parity and preparation-only scaling.

Run baseline before candidate; each stage retains its own immutable receipt.
Use compare after both stages. Large cases never enter a numerical solve.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import platform
import statistics
import subprocess
import tarfile
import time
from pathlib import Path

import joint_validation as v
from joint_runtime_support import unpack, differences, scientific
from joint_sparse_validation import parity

BASELINE = '15e71e3fc230205e6efd7f60c8a915eb1989733b'


def source_hash(root):
    files = sorted(p for folder in ('src', 'include', 'cmake') for p in (root/folder).rglob('*')
                   if p.is_file() and p.suffix in ('.cpp', '.hpp', '.h', '.cmake', '.txt'))
    files += [root/'CMakeLists.txt']
    return v.digest({str(p.relative_to(root)): v.sha(p) for p in files})


def baseline_source_hash():
    data=subprocess.check_output(['git','archive',BASELINE,'src','include','cmake','CMakeLists.txt'],cwd=v.ROOT)
    with tarfile.open(fileobj=io.BytesIO(data)) as archive:
        return v.digest({m.name:hashlib.sha256(archive.extractfile(m).read()).hexdigest()
                         for m in archive if m.isfile() and Path(m.name).suffix in ('.cpp','.hpp','.h','.cmake','.txt')})


def fingerprint(root):
    cache = (root/'CMakeCache.txt').read_text()
    source = Path(next(line.split('=', 1)[1] for line in cache.splitlines()
                       if line.startswith('CMAKE_HOME_DIRECTORY:')))
    binaries = [root/'bin/joint_sparse_benchmark', root/'bin/joint_component_runtime',
                root/'bin/joint_validation', root/'bin/RHBM-GEM', *sorted((root/'src').glob('librhbm_gem.*'))]
    return dict(source_root=str(source), source_sha256=source_hash(source), cache=cache,
                binaries={str(p.relative_to(root)): v.sha(p) for p in binaries if p.is_file()},
                configuration=(root/'generated/Release/SimulationConfiguration-CXX.txt').read_text(),
                linked_libraries=subprocess.check_output(
                    ['otool', '-L', str(root/'src/librhbm_gem.dylib')] if platform.system() == 'Darwin'
                    else ['ldd', str(root/'src/librhbm_gem.so')], text=True))


def compare(root):
    baseline, candidate = v.read(root/'baseline.json'), v.read(root/'candidate.json')
    expected=baseline_source_hash()
    if any(b['source_sha256']!=expected for b in baseline['builds'].values()):
        raise ValueError('Baseline source does not match the frozen commit')
    def inputs(report):
        return {case:sorted((Path(name).name,digest) for name,digest in files.items())
                for case,files in report.get('inputs',{}).items()}
    if inputs(baseline)!=inputs(candidate):
        raise ValueError('Baseline and candidate input fingerprints differ')
    out = {}
    for key, runs in baseline['fixed'].items():
        actual = candidate['fixed'].get(key, [])
        if runs and actual and runs[0]['process']['status']==actual[0]['process']['status']=='not-run-driver-null-widths':
            out[key]=dict(passed=None, complete=False, reason='Covered by real public-runtime control'); continue
        checks = []
        for x, y in zip(runs, actual):
            a, b = x.get('result'), y.get('result')
            if not a or not b or x['process']['status'] != 'completed' or y['process']['status'] != 'completed':
                checks.append(dict(passed=None, reason='incomplete-measurement')); continue
            checks_for_roles = parity(a, b)
            for role, check in checks_for_roles.items():
                if not check.get('available'):
                    check['expected_unavailable'] = (a.get(role, {}).get('valid') is False and
                        b.get(role, {}).get('valid') is False and not differences(a.get(role), b.get(role)))
            passed = all(c.get('passed') or c.get('expected_unavailable') for c in checks_for_roles.values())
            op = b.get('operator')
            checks.append(dict(passed=passed and (op is None or op.get('passed', False)), roles=checks_for_roles,
                               operator=op, reason='same-state'))
        complete=len(checks)==3 and all(c['passed'] is not None for c in checks)
        out[key] = dict(checks=checks, complete=complete,
                        passed=all(c['passed'] for c in checks) if complete else None)
    for backend in ('eigen','spqr'):
        for control in ('physical','real-runtime'):
            left=root/'baseline'/backend/control/'output'; right=root/'candidate'/backend/control/'output'
            files=sorted(left.rglob('runtime.json')) if left.exists() else []
            delta=[] if {p.relative_to(left) for p in files}=={p.relative_to(right) for p in right.rglob('runtime.json')} else ['output-files-differ']
            for path in files:
                other=right/path.relative_to(left)
                delta += differences(scientific(v.read(path)),scientific(v.read(other))) if other.is_file() else ['missing-output']
            key=backend+'/'+control
            processes=[r['controls'].get(key,{}).get('process',{}).get('status') for r in (baseline,candidate)]
            complete=bool(files) and processes==['completed','completed']
            out[key]=dict(passed=not delta if complete else None, complete=complete,
                          differences=delta, compared_files=len(files), process_statuses=processes)
    v.write(root/'comparison.json', out)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('baseline', 'candidate', 'large', 'compare'), required=True)
    parser.add_argument('--eigen', type=Path)
    parser.add_argument('--spqr', type=Path)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--model', type=Path, default=Path('/Users/yslian/data/6z6u.cif'))
    parser.add_argument('--map', type=Path, default=Path('/Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map'))
    args = parser.parse_args(); root=args.work_dir.resolve()
    if args.stage == 'compare':
        result=compare(root)
        if any(r['passed'] is False for r in result.values()): raise SystemExit(1)
        if any(r['passed'] is None and r.get('reason')!='Covered by real public-runtime control' for r in result.values()): raise SystemExit(2)
        return
    if args.eigen is None or args.spqr is None: parser.error('--eigen and --spqr are required')
    receipt=root/(args.stage+'.json')
    if receipt.exists(): raise SystemExit('Keep prior receipts; select a new campaign directory.')
    root.mkdir(parents=True, exist_ok=True)
    builds={'eigen': args.eigen.resolve(), 'spqr': args.spqr.resolve()}
    deadline=time.monotonic()+4800
    report=dict(baseline_commit=BASELINE, platform=platform.platform(),
                limits=dict(process_seconds=600, rss_bytes=4*1024**3, campaign_seconds=4800),
                numeric_environment={k:v.ENV[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
                builds={name:fingerprint(path) for name,path in builds.items()}, fixed={}, inputs={}, controls={})
    if args.stage=='baseline' and any(b['source_sha256']!=baseline_source_hash() for b in report['builds'].values()):
        raise ValueError('Baseline source does not match the frozen commit')
    def save(): v.write(receipt, report)
    def run(command, directory, output):
        process=v.monitored(command, directory, deadline)
        row=dict(process=process)
        if output.is_file(): row['result']=v.read(output)
        print(f'{directory.relative_to(root)}: {process["status"]}', flush=True)
        return row
    save()
    if args.stage == 'large':
        report['generator_sha256']=v.sha(v.ROOT/'tests/support/JointOperatorWorkload.cpp')
        for backend, build in builds.items():
            for topology in ('chain','cube'):
                for size in (128,512,2000,5000,10000):
                    key=f'{backend}/{topology}-{size}'; directory=root/'large'/key; output=directory/'state.json'
                    directory.mkdir(parents=True,exist_ok=True)
                    report['fixed'][key]=[run([build/'bin/joint_sparse_benchmark','synthetic',topology,size,'prepare',output,'--resources'],directory,output)]
                    result=report['fixed'][key][0].get('result')
                    if result:
                        report['fixed'][key][0]['manifest_sha256']=v.digest({k:x for k,x in result.items() if k not in ('resources','work','peak_rss_bytes') and not k.endswith('seconds')})
                    save()
        report['finished']=True; save(); return
    fixture=unpack(v.ROOT/'tests/fixtures/joint_component/catalog.json','heterogeneous-168',root/'inputs')
    report['inputs']['heterogeneous-168']={p.name:v.sha(p) for p in fixture.iterdir() if p.is_file()}
    inputs={}
    for size in (128,512):
        directory=root/'inputs'/f'single-{size}'
        if not (directory/'widths.json').exists():
            generated=run([builds['eigen']/'bin/joint_validation','single',size,directory],directory/'generation',directory/'unused.json')
            report['controls'][f'generate-{size}']=generated; save()
            if generated['process']['status']!='completed': continue
            prepared=run([builds['eigen']/'bin/joint_sparse_benchmark','prepare',directory/'input.cif',directory/'input.map',directory/'widths.json'],directory/'prepare',directory/'widths.json')
            report['controls'][f'prepare-{size}']=prepared; save()
            if prepared['process']['status']!='completed': continue
        inputs[f'single-{size}']=(directory/'input.cif',directory/'input.map',directory/'widths.json')
    if args.model.is_file() and args.map.is_file():
        directory=root/'inputs'/'6z6u'; directory.mkdir(parents=True,exist_ok=True)
        if not (directory/'widths.json').exists():
            report['controls']['prepare-6z6u']=run([builds['eigen']/'bin/joint_sparse_benchmark','prepare',args.model,args.map,directory/'widths.json'],directory/'prepare',directory/'widths.json'); save()
        if (directory/'widths.json').exists(): inputs['6z6u']=(args.model,args.map,directory/'widths.json')
    else: report['controls']['6z6u']=dict(status='not-run-missing-input')
    # The old fixed-state driver cannot decode null widths for analytically
    # eliminated halo contributors. Exercise the unchanged public workflow for
    # this case; never replace missing widths with invented initial values.
    real_runtime = '6z6u' in inputs and any(x is None for x in v.read(inputs['6z6u'][2])['b'])
    for name, paths in inputs.items(): report['inputs'][name]={str(p):v.sha(p) for p in paths}
    save()
    for backend, build in builds.items():
        for name in ('single-128','heterogeneous-168','single-512','6z6u'):
            key=backend+'/'+name; runs=report['fixed'][key]=[]
            if name=='6z6u' and real_runtime:
                runs.append(dict(process=dict(status='not-run-driver-null-widths'),
                                 reason='Covered by real public-runtime control')); save(); continue
            for repetition in (1,2,3):
                directory=root/args.stage/key/str(repetition); directory.mkdir(parents=True,exist_ok=True); output=directory/'state.json'
                if name=='heterogeneous-168': command=[build/'bin/joint_sparse_benchmark','fixture',fixture,'first-stage-float32',output]
                elif name in inputs: command=[build/'bin/joint_sparse_benchmark','initial',*inputs[name],output]
                else: runs.append(dict(process=dict(status='not-run-missing-input'))); break
                if args.stage=='candidate': command+=['--resources','--operator']
                row=run(command,directory,output); runs.append(row); save()
                if row['process']['status']!='completed': break
        control=root/args.stage/backend/'physical'; control.mkdir(parents=True,exist_ok=True)
        report['controls'][backend+'/physical']=run([build/'bin/joint_component_runtime','physical',control/'output'],control,control/'output/completion.json'); save()
        if real_runtime:
            directory=root/args.stage/backend/'real-runtime'; directory.mkdir(parents=True,exist_ok=True)
            report['controls'][backend+'/real-runtime']=run([build/'bin/joint_component_runtime','run',args.model,args.map,directory/'output'],directory,directory/'output/completion.json')
            report['controls'][backend+'/real-runtime']['outputs']={str(p.relative_to(directory)):v.sha(p) for p in (directory/'output').rglob('*.json')}
            save()
        if args.stage=='candidate':
            directory=root/args.stage/backend/'workflow'; directory.mkdir(parents=True,exist_ok=True)
            report['controls'][backend+'/workflow']=run([build/'bin/joint_sparse_benchmark','synthetic','cube',128,'workflow',directory/'state.json','--resources'],directory,directory/'state.json'); save()
    report['medians']={key:{field:statistics.median(r['result'][field] for r in runs)
                                for field in ('primary_seconds','reference_seconds','derivative_seconds','peak_rss_bytes')}
                       for key,runs in report['fixed'].items() if len(runs)==3 and all(r['process']['status']=='completed' for r in runs)}
    report['finished']=True; save()


if __name__=='__main__': main()
