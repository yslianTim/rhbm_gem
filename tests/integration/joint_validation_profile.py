"""A bounded diagnostic sampling run, excluded from end-to-end timing samples (macOS)."""
import argparse
import os
from pathlib import Path
import signal
import subprocess
import time
from joint_validation import ENV, LIMITS, process_tree_rss, read, sha, write


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--case',required=True)
    parser.add_argument('--cli',type=Path,required=True)
    args=parser.parse_args();work=args.work_dir.resolve();campaign=read(work/'campaign.json')
    stage=campaign['stages']['c']
    if stage['status']!='finished':raise SystemExit('Finish the serial measurement campaign before profiling.')
    root=work/'profiles'/args.case
    if root.exists():raise SystemExit('Profile already attempted; keep the existing evidence.')
    remaining=min(4800-stage['wall_seconds'],7200-sum(s['wall_seconds']for s in campaign['stages'].values()))
    if remaining<32:raise SystemExit('Insufficient remaining campaign budget for this diagnostic.')
    source=work/stage['resource_directory']/args.case
    case=read(work/stage['resource_directory']/'receipt.json')[args.case]
    if case['kind']=='user':
        model=Path('/Users/yslian/data/6z6u.cif')
        map_path=Path('/Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map')
    else:model,map_path=source/'input.cif',source/'input.map'
    if case['input_hashes']!=dict(model=sha(model),map=sha(map_path)):
        raise SystemExit('Profile inputs no longer match the measured case.')
    if sha(args.cli)!=stage['cli_sha256']:
        raise SystemExit('Profile executable no longer matches the measured command.')
    root.mkdir(parents=True)
    command=[str(args.cli.resolve()),'potential_analysis','--estimator','joint-components','-a',str(model),'-m',str(map_path),
             '-d',str(root/'diagnostic.sqlite'),'-k','diagnostic','--exclude-hydrogen','true','--asymmetry','false',
             '--only-backbone','false','--map-normalization','false','-j','1','-v','0']
    start=time.monotonic();record=dict(case=args.case,command=command,status='running',wall_seconds=32,
        scope='diagnostic sampling only; not an uninstrumented performance repetition',cli_sha256=sha(args.cli),
        model_sha256=sha(model),map_sha256=sha(map_path),source_sha256=sha(__file__),
        runner_sha256=sha(Path(__file__).with_name('joint_validation.py')),sample_at_seconds=[5,15,25],sample_duration_seconds=1)
    stage['wall_seconds']+=32
    campaign.setdefault('diagnostic_profiles',[]).append(record);write(work/'campaign.json',campaign)
    samplers=[];peak=0
    with (root/'stdout.txt').open('w') as stdout,(root/'stderr.txt').open('w') as stderr:
        process=subprocess.Popen(command,stdout=stdout,stderr=stderr,env=ENV,start_new_session=True)
        try:
            times=[5,15,25]
            while process.poll() is None:
                now=time.monotonic();elapsed=now-start;peak=max(peak,process_tree_rss(process.pid))
                if peak>LIMITS['rss_bytes'] or elapsed>=32:
                    record['status']='rss-limit' if peak>LIMITS['rss_bytes'] else 'intentional-diagnostic-stop'
                    break
                if times and elapsed>=times[0]:
                    second=times.pop(0)
                    sampler=subprocess.Popen(['/usr/bin/sample',str(process.pid),'1','-file',str(root/f'stack-{second}.txt')],
                        stdout=stderr,stderr=stderr,start_new_session=True)
                    samplers.append(sampler)
                time.sleep(max(0,.1-(time.monotonic()-now)))
            if record['status']=='running':record['status']='command-completed' if process.returncode==0 else 'process-failure'
        finally:
            if process.poll() is None:
                os.killpg(process.pid,signal.SIGTERM)
                try:process.wait(timeout=1)
                except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
            for sampler in samplers:
                if sampler.poll() is None:os.killpg(sampler.pid,signal.SIGKILL)
                sampler.wait()
            elapsed=time.monotonic()-start;stage['wall_seconds']+=elapsed-32
            record.update(wall_seconds=elapsed,exit_code=process.returncode,sampled_tree_peak_rss_bytes=peak,
                          sample_exit_codes=[s.returncode for s in samplers],
                          stack_sha256={p.name:sha(p) for p in root.glob('stack-*.txt')})
            write(root/'receipt.json',record);write(work/'campaign.json',campaign)
    print(record)

if __name__=='__main__':main()
