"""Shared infrastructure and command measurement for Joint validations."""
import collections
import time
from pathlib import Path

from experiment_io import ROOT, digest, read, sha, write
from experiment_process import ENV, RSS_LIMIT_BYTES, monitored, process_tree_rss


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
            assert outcome['schema_version'] in (3,4,5) and outcome['metadata']['map_normalization']['divisor']==1
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
