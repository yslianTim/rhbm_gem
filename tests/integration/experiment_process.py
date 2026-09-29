"""Process-only watchdog and resource measurements for experiment runners."""
from __future__ import annotations
import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path

from experiment_io import write

ENV = dict(os.environ, **{k: '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')})
RSS_LIMIT_BYTES = 4 * 1024**3


def process_tree_rss(pid):
    result = subprocess.run(['ps', '-axo', 'pid=,ppid=,rss='], capture_output=True, text=True, check=True)
    rows = [tuple(map(int, line.split())) for line in result.stdout.splitlines()]
    descendants = {pid}
    while True:
        found = {p for p, parent, _ in rows if parent in descendants}
        if found <= descendants:
            break
        descendants |= found
    return sum(rss * 1024 for p, _, rss in rows if p in descendants)


def monitored(command, directory, deadline, *, rss_limit=RSS_LIMIT_BYTES, seconds=600):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    available = deadline - start
    if available <= 0:
        return dict(status='not-run-budget', command=list(map(str, command)), wall_seconds=0)
    expires = min(deadline, start + seconds)
    peak = 0; samples = 0; gap = 0.; last = start; next_update = start + 30
    record = dict(command=list(map(str, command)), status='running', sample_interval_seconds=.1)
    # The parent owns the entire process group, including time and all CLI children.
    timed = ['/usr/bin/time', '-l' if sys.platform == 'darwin' else '-v', *map(str, command)]
    with (directory / 'stdout.txt').open('w') as stdout, (directory / 'stderr.txt').open('w') as stderr:
        process = subprocess.Popen(timed, stdout=stdout, stderr=stderr, env=ENV, start_new_session=True)
        try:
            while process.poll() is None:
                now = time.monotonic(); gap = max(gap, now - last); last = now
                peak = max(peak, process_tree_rss(process.pid)); samples += 1
                reason = 'rss-limit' if peak > rss_limit else ('time-limit' if now >= expires else None)
                if reason:
                    record['status'] = reason
                    os.killpg(process.pid, signal.SIGTERM)
                    try: process.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL); process.wait()
                    break
                if now >= next_update:
                    print(f"running {directory.name}: {now-start:.0f}s, sampled peak {peak/1024**2:.0f} MiB", flush=True)
                    next_update = now + 30
                time.sleep(max(0, min(.1 - (time.monotonic() - now), expires - time.monotonic())))
        except BaseException:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL); process.wait()
            raise
    if record['status'] == 'running': record['status'] = 'completed' if process.returncode == 0 else 'process-failure'
    text = (directory / 'stderr.txt').read_text()
    match = re.search(r'(\d+)\s+maximum resident set size', text) if sys.platform == 'darwin' else re.search(r'Maximum resident set size \(kbytes\): (\d+)', text)
    record.update(wall_seconds=time.monotonic() - start, exit_code=process.returncode,
                  sampled_tree_peak_rss_bytes=peak, samples=samples, maximum_sampling_gap_seconds=gap,
                  os_process_peak_rss_bytes=int(match[1]) * (1 if sys.platform == 'darwin' else 1024) if match else None)
    if record['status'] == 'completed' and (record['os_process_peak_rss_bytes'] or 0) > rss_limit:
        record['status'] = 'rss-limit-observed-after-exit'
    write(directory / 'process.json', record)
    return record
