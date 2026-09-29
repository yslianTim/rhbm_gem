"""Bounded three-way compact SVD acceptance; report-only rechecks saved evidence."""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import platform
import subprocess
import time
from types import SimpleNamespace

import joint_validation_support as v
from joint_runtime_support import unpack
from joint_compact_support import summarize_compact_receipt

MODES = ('legacy', 'values', 'auto')
CASES = ('single-128', 'heterogeneous-168', 'single-512')
BASELINE_SOURCE = '9718531067e72cd3ecf180d12ae8033b8e39360b7fe1d49d24ae847de4240fa0'





def summary(r, root=Path('.')):
    return summarize_compact_receipt(r, root, cases=CASES, modes=MODES, backends=('eigen', 'spqr'),
        command_cases=('single-128', 'single-512'),
        performance_ratio_limits={'single-128': 1.1, 'heterogeneous-168': 1.1, 'single-512': .7})


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
