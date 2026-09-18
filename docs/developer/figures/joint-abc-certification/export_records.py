"""Export reviewed certification results, including immutable rows, without duplicates.

Usage: python3 export_records.py RUN_DIRECTORY OUTPUT_DIRECTORY
The archive contains each observation snapshot once. Per-branch residual tables
are reproducible from those rows and the retained coefficients; their original
hashes remain in the raw inventory.
"""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import sys
import tarfile


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def table(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    root, output = map(Path, sys.argv[1:]); output.mkdir(parents=True, exist_ok=True)
    result = read(root/'results.json')
    evidence, counts, costs, controls, starts, representatives = [], [], [], [], [], []
    for group in result['results']:
        dataset, variant = group['dataset'], group['variant']
        path = root/'datasets'/dataset/variant
        counts.append({'dataset': dataset, 'variant': variant,
                       'legacy_qualified': sum(c['joint_qualified'] for c in group['cases']),
                       'regular_qualified': sum(c['certificate']['regular_qualified'] for c in group['cases']),
                       'recovered_double': sum(c['oracle_recovered'] is True for c in group['cases']),
                       'accurate_float32': sum(c['float32_accuracy_passed'] is True for c in group['cases']),
                       'multistart_consistent': group['multistart_consistent']})
        representatives.append({'dataset': dataset, 'variant': variant,
                                **{'legacy_'+k: v for k, v in group['representatives'].items()},
                                **{'regular_'+k: v for k, v in group['regular_representatives'].items()}})
        controls.extend({'dataset': dataset, 'variant': variant, **r} for r in group['fixed_b_comparisons'])
        starts.extend({'dataset': dataset, 'variant': variant, **r} for r in group['pairwise_comparisons'])
        for case in group['cases']:
            label = case['case']; fit = read(path/'fits'/f'{label}.json')
            audit = read(path/'audits'/f'{label}.json')
            for row in case['certificate']['evidence']:
                evidence.append({'dataset': dataset, 'variant': variant, 'case': label, **row})
            for phase, value in fit.get('resources', {}).items():
                costs.append({'dataset': dataset, 'variant': variant, 'case': label, 'phase': phase,
                              'seconds': value['seconds'], 'process_peak_rss_bytes': value['process_peak_rss_bytes'],
                              'included_in': 'fit total'})
            costs.append({'dataset': dataset, 'variant': variant, 'case': label, 'phase': 'search reference',
                          'seconds': fit.get('search_reference_seconds', 0),
                          'process_peak_rss_bytes': fit.get('resources', {}).get('search', {}).get('process_peak_rss_bytes'),
                          'included_in': 'search; RSS is search-process high-water mark'})
            for phase, value in [('new endpoint audit', audit), ('high precision', audit.get('precision', {}))]:
                costs.append({'dataset': dataset, 'variant': variant, 'case': label, 'phase': phase,
                              'seconds': value.get('seconds', 0), 'process_peak_rss_bytes': value.get('process_peak_rss_bytes'),
                              'included_in': 'new endpoint audit' if phase == 'high precision' else 'audit process'})
            boundary = path/'audits'/f'{label}-boundary.json'
            if boundary.exists():
                value = read(boundary)
                costs.append({'dataset': dataset, 'variant': variant, 'case': label, 'phase': 'boundary scan',
                              'seconds': value['seconds'], 'process_peak_rss_bytes': value['process_peak_rss_bytes'],
                              'included_in': 'audit process'})
    for dataset in sorted((root/'datasets').iterdir()):
        value = read(dataset/'initialization.json')['resources']
        costs.append({'dataset': dataset.name, 'variant': 'shared', 'case': 'all', 'phase': 'initialization',
                      'seconds': value['seconds'], 'process_peak_rss_bytes': value['process_peak_rss_bytes'], 'included_in': 'preparation'})
        for p in sorted(dataset.glob('truth-boundary-*.json')):
            value = read(p)
            costs.append({'dataset': dataset.name, 'variant': 'generating-model', 'case': p.stem, 'phase': 'boundary scan',
                          'seconds': value['seconds'], 'process_peak_rss_bytes': value['process_peak_rss_bytes'], 'included_in': 'audit process'})
    for name, rows in [('failure-matrix-evidence', evidence), ('dataset-counts', counts), ('phase-costs', costs),
                       ('fixed-b-comparison', controls), ('multistart', starts), ('representatives', representatives)]:
        table(output/f'{name}.csv', rows)
    for name in ['results.json', 'failure-matrix.csv', 'estimates.csv', 'input-hashes.json', 'provenance.json',
                 'audit-provenance.json', 'audit-revision.json', 'reproducibility.json', 'scientific-validation.json',
                 'execution-status.json', 'search-execution.json', 'audit-execution.json', 'trajectory-audit.json',
                 'boundary-summary.csv', 'domain-equivalence.json', 'search-reproducibility.json']:
        if (root/name).exists(): shutil.copyfile(root/name, output/name)
    inventory = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in sorted(root.rglob('*')) if p.is_file() and p.name not in ('artifact-index.json', 'reproducibility.json')}
    write(output/'raw-artifact-index.json', inventory)
    selected = []
    for p in sorted(root.rglob('*')):
        relative = p.relative_to(root)
        if not p.is_file() or p.suffix not in ('.json', '.csv'): continue
        if 'residuals' in relative.parts: continue
        if len(relative.parts) > 3 and relative.parts[0] == 'datasets':
            if relative.parts[2] in ('legacy', 'guarded', 'guarded-log') and relative.name in (
                    'voxels.csv', 'dataset.json', 'snapshot.json', 'initialization.json', 'scoring-truth.json',
                    'preparation-resources.json', 'forward-status.json'):
                continue
            if relative.parts[2] in ('fits', 'controls', 'weak-directions'): continue
        selected.append((p, str(relative)))
    with (output/'scientific-records.tar.gz').open('wb') as stream, gzip.GzipFile(fileobj=stream, mode='wb', mtime=0, filename='') as compressed:
        with tarfile.open(fileobj=compressed, mode='w|') as archive:
            for path, name in selected:
                data = path.read_bytes(); info = tarfile.TarInfo(name)
                info.size = len(data); info.mode = 0o644; archive.addfile(info, io.BytesIO(data))
    write(output/'bundle-description.json', {
        'archive': 'scientific-records.tar.gz', 'files': len(selected),
        'contents': 'All branch fits/traces/audits, controls, independent model evidence and one immutable voxel/contributor snapshot per dataset.',
        'omissions': 'Duplicate observation/initialization copies and per-branch residual tables; reconstruct predictions from retained rows and fitted coefficients. Raw original hashes are retained.',
        'memory_semantics': 'Process high-water RSS, not isolated per-phase allocation. Nested reference and precision times must not be added twice.'})


if __name__ == '__main__': main()
