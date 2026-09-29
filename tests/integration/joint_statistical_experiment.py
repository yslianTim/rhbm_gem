"""Run the fixed-seed Joint noise and position-mismatch research experiment."""
from __future__ import annotations

import argparse
import collections
import hashlib
import math
import os
import time
from pathlib import Path

import numpy as np
from experiment_io import read, sha, write
from experiment_process import monitored, process_tree_rss

SEEDS = list(range(20260921, 20260941))
PROCESS_SECONDS = 600
RSS_BYTES = 4 * 1024**3
STAGE_SECONDS = 1200

AXIS = np.arange(25) * .5 - 6
Z, Y, X = np.meshgrid(AXIS, AXIS, AXIS, indexing='ij')
TARGET = (X * X + Y * Y + Z * Z <= 6.25).ravel()
TRUTH = np.array([[2., .5, .2], [2.3, .55, .15]])


def clean_map(distance, shift):
    values = np.zeros(X.shape)
    for x, a, b, c in [(0, 2, .5, .2), (distance + shift, 2.3, .55, .15)]:
        square = (X - x) ** 2 + Y ** 2 + Z ** 2
        radius = np.sqrt(square)
        keep = square <= 6.25
        gaussian = np.exp(-square / (2 * b * b)) / (2 * math.pi * b * b) ** 1.5
        charge = np.zeros(X.shape)
        center = radius < 1e-5
        charge[center] = math.sqrt(2 / math.pi) / b
        nonzero = ~center
        charge[nonzero] = np.fromiter(
            (math.erf(float(t) / (math.sqrt(2) * b)) / float(t) for t in radius[nonzero]), float)
        values[keep] += (a * gaussian + c * charge)[keep]
    return values.ravel()


def noise_field(seed, kind):
    rng = np.random.Generator(np.random.PCG64(seed))
    if kind == 'iid':
        return rng.normal(size=(25, 25, 25)).ravel()
    field = rng.normal(size=(31, 31, 31))
    kernel = np.exp(-.5 * np.arange(-3, 4, dtype=float) ** 2)
    kernel /= np.sqrt(np.sum(kernel * kernel))  # unit marginal variance, no sample normalization
    for axis in range(3):
        field = np.apply_along_axis(lambda row: np.convolve(row, kernel, mode='valid'), axis, field)
    return field.ravel()


def conditions():
    for distance in (1.2, .6):
        yield distance, 0., 'none', 0., None, ['fixed', 'production']
        for kind in ('iid', 'correlated'):
            for sigma in (.01, .05, .1):
                for seed in SEEDS:
                    yield distance, 0., kind, sigma, seed, ['fixed', 'production'] if sigma == .05 else ['fixed']
        for shift in (.05, .15):
            modes = ['fixed', 'production'] if shift == .15 else ['fixed']
            yield distance, shift, 'none', 0., None, modes
            for seed in SEEDS:
                yield distance, shift, 'iid', .05, seed, modes


def metrics(outcome, prediction, values, clean, rows):
    result = dict(state_available=outcome['assembled_state'] is not None,
                  runtime_convergence=outcome['runtime_convergence'],
                  initialization_valid=outcome['initialization']['valid'],
                  stop_reasons=[c['stop_reason'] for c in outcome['components']],
                  failure_reasons=[e['name'] + ':' + e['status'] for c in outcome['components']
                                   for e in c['evidence'] if e['status'] in ('failed', 'unavailable')],
                  parameters=None, errors=None, prediction_rmse=None, residual_rmse=None,
                  residual_neighbor_correlation=None)
    if result['state_available']:
        state = outcome['assembled_state']
        ac = np.array(state['ac']).reshape(2, 2)
        parameters = np.column_stack([ac[:, 0], state['b'], ac[:, 1]])
        prediction = np.array(prediction)
        residual = prediction - values[rows]
        result.update(parameters=parameters.tolist(), errors=(parameters - TRUTH).tolist(),
                      prediction_rmse=float(np.sqrt(np.mean((prediction - clean[rows]) ** 2))),
                      residual_rmse=float(np.sqrt(np.mean(residual ** 2))))
        lookup = {int(row): i for i, row in enumerate(rows)}
        left = []
        right = []
        for row, i in lookup.items():
            for stride, coordinate in ((1, row % 25), (25, (row // 25) % 25), (625, row // 625)):
                if coordinate < 24 and row + stride in lookup:
                    left.append(residual[i])
                    right.append(residual[lookup[row + stride]])
        if len(left) > 2 and np.std(left) > 0 and np.std(right) > 0:
            result['residual_neighbor_correlation'] = float(np.corrcoef(left, right)[0, 1])
    return result


def wilson(failures, total):
    if not total:
        return None
    z = 1.959963984540054
    p = failures / total
    den = 1 + z * z / total
    center = (p + z * z / (2 * total)) / den
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / den
    return [max(0, center - half), min(1, center + half)]


def parameter_stats(results):
    if not results:
        return dict(n=0, bias=None, rmse=None, variance=None, bias_mcse=None, mse_mcse=None)
    errors = np.array([r['errors'] for r in results])
    n = len(errors)
    variance = errors.var(axis=0, ddof=1) if n > 1 else None
    squared = errors ** 2
    rmse = np.sqrt(squared.mean(axis=0))
    mse_mcse = squared.std(axis=0, ddof=1) / math.sqrt(n) if n > 1 else None
    return dict(n=n, bias=errors.mean(axis=0).tolist(), rmse=rmse.tolist(),
                variance=variance.tolist() if variance is not None else None,
                bias_mcse=np.sqrt(variance / n).tolist() if variance is not None else None,
                mse_mcse=mse_mcse.tolist() if mse_mcse is not None else None)


def statistical_summary(records):
    groups = collections.defaultdict(list)
    noiseless = {}
    for record in records:
        for mode, result in record.get('results', {}).items():
            if record['noise'] == 'none':
                noiseless[record['distance'], record['shift'], mode] = result
        for mode in record['modes']:
            key = (record['distance'], record['shift'], record['noise'], record['sigma_fraction'], mode)
            groups[key].append((record, record.get('results', {}).get(mode)))
    summary = []
    for key, attempts in sorted(groups.items()):
        results = [value for _, value in attempts if value is not None]
        available = [result for result in results if result['state_available']]
        converged = [result for result in available if result['runtime_convergence'] == 'passed']
        failed = [result for result in available if result['runtime_convergence'] != 'passed']
        unavailable = sum(not result['state_available'] for result in results)
        execution_counts = collections.Counter(record.get('status', 'not-run') for record, _ in attempts)
        row = dict(distance=key[0], shift=key[1], noise=key[2], sigma_fraction=key[3], mode=key[4],
                   attempted_cases=len(attempts),
                   completed_cases=execution_counts['completed'],
                   qualified_cases=len(converged), failed_cases=len(failed), unavailable_cases=unavailable,
                   execution_status_counts=dict(execution_counts),
                   numerical_qualification_counts=dict(collections.Counter(
                       result['runtime_convergence'] if result['state_available'] else 'unavailable'
                       for result in results)),
                   attempts=len(attempts), recorded=len(results),
                   process_incomplete=len(attempts) - execution_counts['completed'],
                   unavailable=unavailable, unavailable_interval95=wilson(unavailable, len(results)),
                   unavailable_rate=unavailable / len(results) if results else None,
                   nonconverged_available=len(failed),
                   nonconverged_interval95=wilson(len(failed), len(available)),
                   nonconverged_available_rate=len(failed) / len(available) if available else None,
                   convergence_counts=dict(collections.Counter(r['runtime_convergence'] for r in results)),
                   all_available=parameter_stats(available), converged_only=parameter_stats(converged),
                   nonconverged_only=parameter_stats(failed),
                   failure_reasons=dict(collections.Counter(reason for r in results for reason in r['failure_reasons'])))
        for metric in ('prediction_rmse', 'residual_rmse', 'residual_neighbor_correlation'):
            values = [r[metric] for r in available if r.get(metric) is not None]
            row[metric] = dict(n=len(values), mean=float(np.mean(values)) if values else None,
                               mcse=float(np.std(values, ddof=1) / math.sqrt(len(values))) if len(values) > 1 else None)
        base = noiseless.get((key[0], key[1], key[4]))
        if key[1] > 0 and base and base['parameters'] is not None:
            differences = [np.array(r['parameters']) - np.array(base['parameters']) for r in available]
            row['difference_from_noiseless_reference'] = dict(
                n=len(differences), mean=np.mean(differences, axis=0).tolist() if differences else None,
                interpretation='paired procedural reference, not a proven global optimum or physical truth')
        summary.append(row)
    return summary


def run_experiment(args, cases, deadline, receipt):
    root = args.work_dir
    root.mkdir(parents=True, exist_ok=True)
    records = []
    maps = {}
    for index, (distance, shift, kind, sigma, seed, modes) in enumerate(cases):
        case = root / f'{index:03d}'
        case.mkdir(exist_ok=True)
        record = dict(index=index, distance=distance, shift=shift, noise=kind,
                      sigma_fraction=sigma, seed=seed, modes=modes)
        if time.monotonic() >= deadline:
            record['status'] = 'not-run-budget'
            records.append(record)
            continue
        key = (distance, shift)
        if key not in maps:
            maps[key] = clean_map(*key)
        clean = maps[key]
        reference = maps.setdefault((distance, 0.), clean_map(distance, 0.))
        scale = float(np.sqrt(np.mean(reference[TARGET] ** 2)))
        values = clean.copy() if seed is None else clean + sigma * scale * noise_field(seed, kind)
        spec = dict(distance=distance, shift=shift, values=values.tolist(), clean=clean.tolist(), modes=modes)
        write(case / 'input.json', spec)
        record.update(input_sha256=sha(case / 'input.json'),
                      observation_sha256=hashlib.sha256(values.astype('<f8').tobytes()).hexdigest(),
                      noise_sigma=sigma * scale)
        run = monitored([args.executable, 'stat', case / 'input.json', case / 'output.json'], case, deadline,
                        rss_limit=RSS_BYTES, seconds=PROCESS_SECONDS)
        record.update(status=run['status'], process=run)
        if run['status'] == 'completed':
            result = read(case / 'output.json')
            rows = np.array(result['row_ids'], dtype=int)
            assert np.array_equal(rows, np.flatnonzero(TARGET)), 'Changed target domain'
            record['generator_max_error'] = result['generator_max_error']
            record['census'] = result['census']
            record['results'] = {mode: metrics(result[mode]['outcome'], result[mode]['prediction'],
                                               values, clean, rows) for mode in modes}
            record['software'] = result['fixed']['outcome']['metadata']['software']
            record['output_sha256'] = sha(case / 'output.json')
        records.append(record)
        write(receipt, dict(seeds=SEEDS, rng='numpy.PCG64', numpy_version=np.__version__,
                            scope='smoke' if args.smoke else 'full', records=records))
        if index % 40 == 0:
            print(f'B {index + 1}/{len(cases)} input cases', flush=True)
    result = dict(seeds=SEEDS, rng='numpy.PCG64', numpy_version=np.__version__,
                  scope='smoke' if args.smoke else 'full', records=records)
    write(receipt, result)
    return result


def write_report(receipt, output):
    output.mkdir(parents=True, exist_ok=True)
    compact = []
    for record in receipt['records']:
        row = {key: value for key, value in record.items() if key != 'process'}
        if 'process' in record:
            row['process'] = {key: value for key, value in record['process'].items()
                              if key not in ('command', 'os_process_peak_rss_bytes', 'samples',
                                             'maximum_sampling_gap_seconds')}
        compact.append(row)
    write(output / 'noise-runs.json', dict(seeds=receipt['seeds'], rng=receipt['rng'],
         numpy_version=receipt['numpy_version'], scope=receipt['scope'], records=compact))
    summary = statistical_summary(receipt['records'])
    write(output / 'noise-summary.json', dict(
        parameter_order=['A', 'B', 'C'], roles=['target', 'halo'],
        coefficient_scale='input-map; normalization disabled',
        intervals='95% Wilson binomial intervals; MCSE for bias, mean squared error and residual metrics; n=1 variance/MCSE unavailable',
        scope=receipt['scope'], groups=summary))
    lines = ['| Distance Å | Shift Å | Noise / σ | Init | Available / attempts | Converged | Target RMSE A / B / C | Halo RMSE A / B / C |',
             '| --- | --- | --- | --- | --- | --- | --- | --- |']
    for row in summary:
        stats = row['all_available']
        rmse = stats['rmse']
        fmt = lambda values: ' / '.join(f'{value:.3g}' for value in values) if values is not None else 'unavailable'
        lines.append(f"| {row['distance']} | {row['shift']} | {row['noise']} / {row['sigma_fraction']} | {row['mode']} | {stats['n']} / {row['attempted_cases']} | {row['qualified_cases']} | {fmt(rmse[0] if rmse else None)} | {fmt(rmse[1] if rmse else None)} |")
    (output / 'noise-table.md').write_text('\n'.join(lines) + '\n')
    files = sorted(path for path in output.iterdir() if path.is_file() and path.name != 'manifest.json')
    write(output / 'manifest.json', dict(schema_version=1, files={path.name: sha(path) for path in files},
                                         report_source_sha256=sha(Path(__file__))))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--smoke', action='store_true', help='Run the first two fixed cases only; not a research result.')
    args = parser.parse_args()
    args.work_dir = args.work_dir.resolve()
    args.executable = args.executable.resolve()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    if args.smoke:
        cases = list(conditions())[:2]
    else:
        cases = list(conditions())
        if len(cases) != 326 or sum(len(case[-1]) for case in cases) != 450:
            raise AssertionError('The fixed statistical experiment matrix changed.')
    process_tree_rss(os.getpid())
    started = time.monotonic()
    receipt_path = args.work_dir / 'receipt.json'
    receipt = run_experiment(args, cases, started + STAGE_SECONDS, receipt_path)
    summary = write_report(receipt, args.work_dir / 'report')
    if args.smoke:
        if len(receipt['records']) != 2 or any(row['status'] != 'completed' for row in receipt['records']):
            raise SystemExit('The research smoke must complete both fixed input cases.')
        if receipt['records'][1]['observation_sha256'] != \
                'c5147ebda1fe72d1012b49d91838743ce4e8ca00abf4b3f8d3dcfa233fba02e1':
            raise AssertionError('The fixed-seed smoke input changed.')
        if sum(row['attempted_cases'] for row in summary) != 3 or \
                sum(row['completed_cases'] for row in summary) != 3:
            raise AssertionError('The smoke summary denominator is incomplete.')
    print(f"{receipt['scope']}: {len(receipt['records'])} attempted inputs; {len(summary)} statistical groups")


if __name__ == '__main__':
    main()
