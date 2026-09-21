"""Summarize saved bounded-validation records without running or reclassifying fits."""
from __future__ import annotations
import argparse
import collections
import math
from pathlib import Path
import numpy as np
from joint_validation import read, write, sha, digest


def wilson(failures, total):
    if not total:return None
    z=1.959963984540054;p=failures/total;den=1+z*z/total
    center=(p+z*z/(2*total))/den;half=z*math.sqrt(p*(1-p)/total+z*z/(4*total*total))/den
    return [max(0,center-half),min(1,center+half)]


def parameter_stats(results):
    if not results:return dict(n=0,bias=None,rmse=None,variance=None,bias_mcse=None,mse_mcse=None)
    errors=np.array([r['errors'] for r in results]);n=len(errors)
    variance=errors.var(axis=0,ddof=1) if n>1 else None
    squared=errors**2;rmse=np.sqrt(squared.mean(axis=0))
    mse_mcse=squared.std(axis=0,ddof=1)/math.sqrt(n) if n>1 else None
    return dict(n=n,bias=errors.mean(axis=0).tolist(),rmse=rmse.tolist(),variance=variance.tolist() if variance is not None else None,
                bias_mcse=np.sqrt(variance/n).tolist() if variance is not None else None,
                mse_mcse=mse_mcse.tolist() if mse_mcse is not None else None)


def statistical_summary(records):
    groups=collections.defaultdict(list);noiseless={}
    for record in records:
        for mode,result in record.get('results',{}).items():
            if record['noise']=='none':noiseless[record['distance'],record['shift'],mode]=result
        for mode in record['modes']:
            key=(record['distance'],record['shift'],record['noise'],record['sigma_fraction'],mode)
            groups[key].append((record,record.get('results',{}).get(mode)))
    summary=[]
    for key,attempts in sorted(groups.items()):
        results=[v for _,v in attempts if v is not None];available=[r for r in results if r['state_available']]
        converged=[r for r in available if r['runtime_convergence']=='passed']
        failed=[r for r in available if r['runtime_convergence']!='passed']
        failure=sum(r['runtime_convergence']!='passed' for r in available)
        unavailable=sum(not r['state_available'] for r in results)
        row=dict(distance=key[0],shift=key[1],noise=key[2],sigma_fraction=key[3],mode=key[4],attempts=len(attempts),
                 recorded=len(results),process_incomplete=len(attempts)-len(results),unavailable=unavailable,unavailable_interval95=wilson(unavailable,len(results)),
                 unavailable_rate=unavailable/len(results) if results else None,
                 nonconverged_available=failure,nonconverged_interval95=wilson(failure,len(available)),
                 nonconverged_available_rate=failure/len(available) if available else None,
                 convergence_counts=dict(collections.Counter(r['runtime_convergence'] for r in results)),
                 all_available=parameter_stats(available),converged_only=parameter_stats(converged),nonconverged_only=parameter_stats(failed),
                 failure_reasons=dict(collections.Counter(reason for r in results for reason in r['failure_reasons'])))
        for metric in ('prediction_rmse','residual_rmse','residual_neighbor_correlation'):
            values=[r[metric] for r in available if r.get(metric) is not None]
            row[metric]=dict(n=len(values),mean=float(np.mean(values)) if values else None,
                             mcse=float(np.std(values,ddof=1)/math.sqrt(len(values))) if len(values)>1 else None)
        base=noiseless.get((key[0],key[1],key[4]))
        if key[1]>0 and base and base['parameters'] is not None:
            differences=[np.array(r['parameters'])-np.array(base['parameters']) for r in available]
            row['difference_from_noiseless_reference']=dict(n=len(differences),mean=np.mean(differences,axis=0).tolist() if differences else None,
                interpretation='paired procedural reference, not a proven global optimum or physical truth')
        summary.append(row)
    return summary


def compact_a(work, out):
    root=work/read(work/'campaign.json')['stages']['a'].get('stage_directory','a')
    data=[]
    for path in sorted(root.glob('start-*.json')):
        record=read(path);snapshot=record.pop('snapshot');snapshot_hash=digest(snapshot)
        assert record['snapshot_sha256']==snapshot_hash
        write(out/'weak-snapshot.json',snapshot)
        record['raw_record_sha256']=sha(path)
        state=record['outcome']['assembled_state']
        record['actual_state_identity_verified']=None
        if 'actual_state_assessment' in record:
            assessed=record['actual_state_assessment']['primary']
            assert state is not None and assessed['beta']==state['ac'] and assessed['eta']==state['log_b']
            record['actual_state_identity_verified']=True
        record['parameter_errors']=None
        if state is not None:
            ac=np.array(state['ac']).reshape(2,2)
            parameters=np.column_stack([ac[:,0],state['b'],ac[:,1]])
            truth=np.column_stack([snapshot['truth_a'],snapshot['truth_b'],snapshot['truth_c']])
            record['parameter_errors']=dict(roles=['target','halo'],order=['A','B','C'],values=(parameters-truth).tolist())
        # Accepted/rejected trial counts and the final diagnostic state are retained; full traces remain local.
        restart=record.get('diagnostic_restart')
        if restart:
            record['diagnostic_restart']={k:restart.get(k) for k in ('runtime_convergence','stop_reason','lm_status','profile_evaluations',
                'accepted_updates','last_trusted_state','primary','local_correction_inf','runtime_checks','settings','usable_state')}
        record['diagnostics_applicable']=record['outcome']['assembled_state'] is not None
        record['case_complete']=record.get('audit_complete',False) or (record['outcome']['assembled_state'] is None)
        data.append(record)
    assert len({r['snapshot_sha256'] for r in data})<=1, 'Weak inputs differ between starts'
    receipt=read(root/'receipt.json')
    original_receipt=read(work/'a/receipt.json') if root!=work/'a' else receipt
    write(out/'weak-halo.json',dict(runs=data,receipt=receipt,original_receipt=original_receipt,
        collection_complete=len(data)==3 and all(r['case_complete'] for r in data) and all(r['status']=='completed' for r in receipt['runs']),
        completion_note='No endpoint audit applies to a missing trusted state; this collection flag does not change numerical statuses.'))
    return data


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--work-dir',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();work=args.work_dir;out=args.output;out.mkdir(parents=True,exist_ok=True)
    a=compact_a(work,out)
    b_root=work/read(work/'campaign.json')['stages']['b'].get('stage_directory','b')
    b=read(b_root/'receipt.json');summary=statistical_summary(b['records'])
    compact=[]
    for record in b['records']:
        row={k:v for k,v in record.items() if k!='process'}
        if 'process'in record:row['measurement']={k:v for k,v in record['process'].items() if k!='command'}
        compact.append(row)
    write(out/'noise-runs.json',dict(seeds=b['seeds'],rng=b['rng'],numpy_version=b['numpy_version'],records=compact))
    write(out/'noise-summary.json',dict(parameter_order=['A','B','C'],roles=['target','halo'],coefficient_scale='input-map; normalization disabled',
        intervals='95% Wilson binomial intervals; MCSE for bias, mean squared error and residual metrics; n=1 variance/MCSE unavailable',groups=summary))
    resource=work/read(work/'campaign.json')['stages']['c'].get('resource_directory','c')/'receipt.json'
    if resource.exists():
        measurements=read(resource)
        for name,case in measurements.items():
            for run in case['runs']:
                if run['status']!='completed':
                    directory=resource.parent/name/f"run-{run['repetition']}"
                    run['partial_output_bytes']=sum(p.stat().st_size for p in directory.iterdir() if p.is_file() and p.name!='receipt.json')
                    run['runtime_convergence']=None # No exported numerical result; this is not a numerical success.
        write(out/'resources.json',measurements)
    if resource!=work/'c/receipt.json':write(out/'resource-preflight.json',read(work/'c/receipt.json'))
    profiles=[]
    markers=('joint_component::','Eigen::SparseQR<','Eigen::JacobiSVD<','RunLocalAlphaTraining','TrainAlphaR')
    for path in sorted((work/'profiles').glob('*/receipt.json')):
        profile=read(path);profile['selected_stack_frames']={}
        for name,expected in profile['stack_sha256'].items():
            stack=path.parent/name
            assert sha(stack)==expected, 'Diagnostic stack changed after capture'
            profile['selected_stack_frames'][name]=[line.strip() for line in stack.read_text().splitlines()
                if any(marker in line for marker in markers)][:20]
        profiles.append(profile)
    if profiles:write(out/'diagnostic-profiles.json',profiles)
    write(out/'campaign.json',read(work/'campaign.json'))
    # A small complete table; all uncertainties and both analysis populations remain in JSON.
    lines=['| Distance Å | Shift Å | Noise / σ | Init | Available / attempts | Converged | Target RMSE A / B / C | Halo RMSE A / B / C |',
           '| --- | --- | --- | --- | --- | --- | --- | --- |']
    for row in summary:
        stats=row['all_available'];rmse=stats['rmse'];fmt=lambda values:' / '.join(f'{v:.3g}' for v in values) if values is not None else 'unavailable'
        lines.append(f"| {row['distance']} | {row['shift']} | {row['noise']} / {row['sigma_fraction']} | {row['mode']} | {stats['n']} / {row['attempts']} | {row['converged_only']['n']} | {fmt(rmse[0] if rmse else None)} | {fmt(rmse[1] if rmse else None)} |")
    (out/'noise-table.md').write_text('\n'.join(lines)+'\n')
    files=sorted(p for p in out.iterdir() if p.is_file() and p.name!='manifest.json')
    write(out/'manifest.json',dict(schema_version=1,files={p.name:sha(p) for p in files},
        raw_output_directory=str(work.resolve()),report_sources={str(Path(__file__).name):sha(__file__)}))
    print(f'{len(a)} weak starts; {len(compact)} statistical inputs; {len(summary)} statistical groups')

if __name__=='__main__':main()
