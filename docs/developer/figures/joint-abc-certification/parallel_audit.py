"""Repeat the frozen endpoint audit in independent single-threaded processes.

Only the expensive active-a dataset is split by case. Its exact endpoint cache
keys are checked for cross-case matches before splitting. All other datasets
retain the original audit order and reuse policy. Search is never rerun.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import shutil
import sys
import time

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO/'tests/integration'))
import joint_abc_certification as cert


def main():
    root, executable = (Path(p).resolve() for p in sys.argv[1:])
    before = cert.fixed.experiment.provenance(executable, REPO)
    cert.require(before == cert.read(root/'provenance.json'), 'Repeat audit must use the frozen search sources/build.')
    cert.require(cert.read(root/'search-execution.json')['returncode'] == 0, 'Search is incomplete.')
    parts = root.parent/(root.name+'-audit-parts')
    cert.require(not parts.exists(), 'Use a fresh partition directory.'); parts.mkdir()
    active = root/'datasets/active-a'; keys = {}
    for path in sorted(active.glob('*/fits/*.json')):
        fit = cert.read(path)
        key = json.dumps([path.stem.endswith('double'), fit.get('primary'), fit['joint_qualified'], fit.get('width_spectrum')], sort_keys=True)
        cert.require(key not in keys or keys[key] == path.stem, 'Cross-case endpoint reuse prevents this partition.')
        keys[key] = path.stem
    jobs = []
    for dataset in sorted((root/'datasets').iterdir()):
        labels = cert.coverage.CASES if dataset.name == 'active-a' else [None]
        for label in labels:
            part = parts/(dataset.name+('-'+label if label else '')); part.mkdir()
            for name in ('inputs.json', 'input-hashes.json'): shutil.copyfile(root/name, part/name)
            target = part/'datasets'/dataset.name; target.mkdir(parents=True)
            for name in ('voxels.csv', 'contributors.csv', 'snapshot.json', 'scoring-truth.json'):
                shutil.copyfile(dataset/name, target/name)
            for variant in cert.VARIANTS:
                fits = target/variant/'fits'; fits.mkdir(parents=True)
                (target/variant/'audits').mkdir()
                for path in sorted((dataset/variant/'fits').glob('*.json')):
                    if label is None or path.stem == label: shutil.copyfile(path, fits/path.name)
            jobs.append(part)
    start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=8) as workers:
        pending = {workers.submit(cert.audit, part, executable): part for part in jobs}
        for future in as_completed(pending):
            future.result(); print('Completed '+pending[future].name, flush=True)
    cert.require(before == cert.fixed.experiment.provenance(executable, REPO), 'Audit source/build changed.')
    copied = {}; duplicate_generating_seconds = 0.
    for part in jobs:
        for path in sorted((part/'datasets').rglob('*.json')):
            relative = path.relative_to(part)
            if 'audits' not in relative.parts and not path.name.startswith('truth-boundary-'): continue
            value = cert.read(path)
            if relative in copied:
                cert.require(cert.coverage.joint.scientific(value) == cert.coverage.joint.scientific(copied[relative]), 'Partitioned generating-model audits disagree.')
                duplicate_generating_seconds += value['seconds']
                continue
            target = root/relative; target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target); copied[relative] = value
    execution = [cert.read(part/'audit-execution.json') for part in jobs]
    cert.write(root/'audit-execution.json', {'returncode': 0, 'seconds': time.perf_counter()-start,
        'sum_partition_seconds': sum(r['seconds'] for r in execution),
        'duplicate_generating_scan_seconds': duplicate_generating_seconds,
        'process_peak_rss_bytes': max(r['process_peak_rss_bytes'] for r in execution),
        'memory_scope': 'Maximum child process high-water RSS; concurrent aggregate RSS is not measured.',
        'policy': '8 independent processes; single numerical thread per process; active-a split by case after verifying no cross-case cache identity.',
        'partitions': [p.name for p in jobs]})
    cert.paths_and_manifest(root)
    cert.write(root/'execution-status.json', {'complete': True, 'inputs_and_sources_stable': True})
    cert.summarize(root)


if __name__ == '__main__': main()
