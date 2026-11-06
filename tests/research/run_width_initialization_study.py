#!/usr/bin/env python3
"""One FirstStage fit per dataset; sequential orthogonal Width seed/bound study."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import statistics as stats
import subprocess
from run_first_stage_sampling_study import shared_equal, summarize

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / 'docs/research/first-stage-sampling'


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def errors(values):
    return {'count': len(values), 'bias': stats.mean(values),
            'rmse': math.sqrt(stats.mean(x*x for x in values)),
            'mae': stats.mean(abs(x) for x in values)} if values else {'count': 0}


def quantiles(values):
    values = sorted(x for x in values if x is not None and math.isfinite(x) and x > 0)
    def q(p):
        t = (len(values)-1)*p
        return values[math.floor(t)] + (values[math.ceil(t)]-values[math.floor(t)])*(t-math.floor(t))
    return {str(p): q(p) for p in (.5, .9, .95, .99)} if values else {}


def truth_summary(atoms):
    out = {}
    for source in ('ols', 'mdpde'):
        out[source] = {label: errors([a[source][k]-a['truth'][k] for a in atoms
                      if a[source][k] is not None]) for label, k in (('width', 1), ('amplitude', 0))}
    paired = [{'phase': a['phase'], 'id': a['id'], 'alpha': a['alpha'],
               'ols_error': abs(a['ols'][1]-a['truth'][1]), 'mdpde_error': abs(a['mdpde'][1]-a['truth'][1]),
               'delta': abs(a['ols'][1]-a['truth'][1])-abs(a['mdpde'][1]-a['truth'][1])}
              for a in atoms if a['ols'][1] is not None and a['mdpde'][1] is not None]
    ds = [p['delta'] for p in paired]
    out['paired'] = {'mean': stats.mean(ds), 'median': stats.median(ds),
                     'improvement_fraction': stats.mean(d < 0 for d in ds),
                     'unchanged_fraction': stats.mean(d == 0 for d in ds),
                     'tie_tolerance_A': 1e-12,
                     'improvement_fraction_beyond_tolerance': stats.mean(d < -1e-12 for d in ds),
                     'unchanged_fraction_within_tolerance': stats.mean(abs(d) <= 1e-12 for d in ds),
                     'worst_OLS_regression': max(ds), 'worst_MDPDE_regression': -min(ds),
                     'outliers': sorted(paired, key=lambda x: x['delta'], reverse=True)[:5],
                     'all_deltas': ds}
    return out


def joint_summary(run, atoms):
    if 'joint' not in run:
        return run
    out = {k: v for k, v in run.items() if k not in ('joint', 'initial_tau', 'bounded_tau', 'first_sweeps')}
    out['initial_tau_sha256'] = digest(run['initial_tau'])
    out['bounded_tau_sha256'] = digest(run['bounded_tau'])
    fit = run['joint']
    checks = fit['checks'] + [x for c in fit['components'] for x in c['checks']]
    out['joint'] = {k: fit[k] for k in ('objective', 'runtime_qualification', 'target_qualification', 'tau', 'ac', 'eta_gradient') if k in fit}
    out['joint']['checks'] = fit['checks']
    out['joint']['failed_checks'] = [c for c in checks if c['status'] != 0]
    for key, names in [('KKT', ('conditional-ac-optimality', 'width-stationarity')),
                       ('Replay', ('independent-replay',)), ('Endpoint', ('endpoint-verification',))]:
        out['joint'][key] = all(any(c['name'] == n for c in checks) for n in names) and all(
            c['status'] == 0 for c in checks if c['name'] in names)
    out['joint']['sweeps'] = [c['sweeps'] for c in fit['components']]
    out['joint']['sweeps_max'] = max(out['joint']['sweeps'], default=0)
    out['joint']['sweeps_sum'] = sum(out['joint']['sweeps'])
    out['joint']['stop_reasons'] = sorted({c['stop'] for c in fit['components']})
    assert not out['joint']['failed_checks'], (run['id'], out['joint']['failed_checks'])
    assert all(out['joint'][k] for k in ('KKT', 'Replay', 'Endpoint'))
    if atoms and 'truth' in atoms[0] and 'tau' in fit:
        out['initial_width_truth'] = errors([t-a['truth'][1] for t, a in zip(run['bounded_tau'], atoms, strict=True)])
        out['final_width_truth'] = errors([t-a['truth'][1] for t, a in zip(fit['tau'], atoms, strict=True)])
        out['final_amplitude_truth'] = errors([fit['ac'][2*i]-a['truth'][0] for i, a in enumerate(atoms)])
    return out


def validate_factorial(data):
    runs = data['runs']
    for a, b in (('MM', 'OM'), ('MO', 'OO')):
        if 'bound' in runs[a] and 'bound' in runs[b]:
            assert runs[a]['bound'] == runs[b]['bound']
    for a, b in (('MM', 'MO'), ('OM', 'OO')):
        if 'initial_tau' in runs[a] and 'initial_tau' in runs[b]:
            assert runs[a]['initial_tau'] == runs[b]['initial_tau']
    if 'MM-common' in runs:
        for a, b in (('MM-common', 'MO-common'), ('OM-common', 'OO-common')):
            if 'joint' in runs[a] and 'joint' in runs[b]:
                assert runs[a]['bounded_tau'] == runs[b]['bounded_tau']
                assert runs[a]['joint'] == runs[b]['joint'], 'Repeated identical run changed endpoint/evidence'


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--build', type=Path, default=Path('build/release'))
    parser.add_argument('--real-data', type=Path, default=Path('/Users/yslian/data'))
    parser.add_argument('--scope', choices=('all', 'synthetic', 'real'), default='all')
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    build = (ROOT / args.build).resolve()
    directory = build / 'width-initialization-study'
    directory.mkdir(exist_ok=True)
    output_path = DOCS / 'ols-vs-mdpde-initialization-results.json'
    output = json.loads(output_path.read_text()) if output_path.exists() else {'synthetic': {}, 'real': {}}
    output['source_commit'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    output['design'] = {'cells': {'MM': ['MDPDE', 'MDPDE'], 'OM': ['OLS', 'MDPDE'],
                                 'MO': ['MDPDE', 'OLS'], 'OO': ['OLS', 'OLS']},
                        'qualification': 'shared MDPDE FirstStage gate plus source-specific finite positive raw Width; no independent OLS qualification',
                        'phases': 'Three deterministic sensitivity cases; no independent-trial confidence intervals',
                        'runtime': 'SolverConfigured through return, including search/qualification/assembly; excludes sampling/fitting/adapter/calibration/projection/downstream',
                        'rss': 'Cumulative process high-water RSS; same process runs all cells serially; cannot infer cell-specific peak allocation'}
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_THREAD_LIMIT='4', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    original = json.loads((DOCS / 'results.json').read_text())
    previous_neighbor = json.loads((DOCS / 'neighbor-selection-results.json').read_text())
    datasets = ([] if args.scope == 'real' else ['fine', 'coarse', 'noisy', 'crowded']) + ([] if args.scope == 'synthetic' else ['6z6u', '8rqb'])
    for scenario in datasets:
        real = scenario in ('6z6u', '8rqb')
        paths = ()
        section = {'datasets': {}, 'atoms': []}
        if real:
            paths = (args.real_data / f'{scenario}.cif', args.real_data /
                     ('emd_11103_unsharpened.map' if scenario == '6z6u' else 'emd_19436_unsharpened.map'))
            if not all(p.is_file() for p in paths):
                output['real'][scenario] = {'status': 'source-unavailable'}
                continue
            section['sources'] = {p.name: {'bytes': p.stat().st_size, 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()} for p in paths}
            assert section['sources'] == original['real_data'][f'{scenario}.cif'], 'Real source SHA-256 mismatch'
        for phase in range(1 if real else 3):
            name = scenario if real else f'{scenario}-r{phase}'
            path = directory / f'{name}.json'
            if not args.summarize_only:
                print(name, flush=True)
                with path.with_suffix('.log').open('w') as log:
                    subprocess.run([str(build / 'bin/first_stage_sampling_study'), str(path),
                                    'real' if real else f'{scenario}:{phase}', 'width', 'production', *map(str, paths)],
                                   cwd=ROOT, env=env, stdout=log, stderr=log, check=True, timeout=900)
            data = json.loads(path.read_text())
            validate_factorial(data)
            parity = {}
            # Committed evidence also supports fresh checkouts without old build artifacts.
            prior_atoms = previous_neighbor['real' if real else 'synthetic'][scenario]['F0']['atoms']
            expected_atoms = [{k: v for k, v in a.items() if k != '_phase'}
                              for a in prior_atoms if a['_phase'] == phase]
            assert shared_equal(expected_atoms, data['atoms']), f'Committed FirstStage parity {name}'
            mm = data['runs']['MM']
            summary_input = dict(data, policy='F0', calibration=data['calibrations']['MDPDE'],
                                 joint=mm['joint'], joint_seconds=mm['solver_seconds'])
            expected_summary = previous_neighbor['runs'][f'{name}-F0-production']
            assert summarize(summary_input)['joint'] == expected_summary['joint'], f'Committed JCE summary parity {name}'
            parity['committed-FirstStage-and-JCE-summary'] = True
            for reference_path in (build / f'neighbor-selection-study/{name}-F0-production.json',
                                   build / f'first-stage-study/{name}-F0-production.json',
                                   DOCS / ('baseline.json' if name == 'coarse-r0' else 'baseline-6z6u.json' if name == '6z6u' else 'no-baseline')):
                if not reference_path.is_file():
                    continue
                reference = json.loads(reference_path.read_text())
                assert shared_equal(reference['atoms'], data['atoms']), f'FirstStage parity {name}'
                assert shared_equal(reference['joint'], data['runs']['MM']['joint']), f'MM parity {name}'
                assert shared_equal(reference['calibration'], data['calibrations']['MDPDE']), f'Calibration parity {name}'
                assert reference['qscore'] == data['qscore']
                parity[str(reference_path.relative_to(ROOT))] = True
            assert parity, f'No exact baseline reference for {name}'
            atom_keys = ('id', 'truth', 'ols', 'mdpde', 'alpha', 'qualified', 'qualification', 'native_status',
                         'original_tau', 'used_tau', 'seed_source', 'reason', 'retained', 'fitting_candidates', 'target')
            atoms = [dict({k: a[k] for k in atom_keys if k in a}, phase=phase) for a in data['atoms']]
            for i, atom in enumerate(atoms):
                atom['sources'] = {s: data['selections'][s]['atoms'][i] for s in ('OLS', 'MDPDE')}
            section['atoms'].extend(atoms)
            keys = ('sampling_seconds', 'alpha_training_and_fitting_seconds', 'adapter_seconds', 'calibration_seconds',
                    'alpha_training_seconds', 'first_stage_fitting_seconds', 'first_stage_bookkeeping_seconds',
                    'initialization_seconds', 'initialization_peak_rss_bytes', 'calibrations', 'donor_cohorts',
                    'matched_cohort_calibrations', 'paired_endpoint_audits', 'common_bound', 'observation_scale', 'initial_objective', 'qscore', 'grid_spacing', 'total_seconds', 'peak_rss_bytes')
            summary = {k: data[k] for k in keys if k in data}
            summary['exact_baseline_parity'] = parity
            summary['factorial_validated'] = True
            summary['observation_sha256'] = digest(data['observations'])
            summary['raw_sampling_sha256'] = digest([a['neighbor_selection']['observations'] for a in data['atoms']])
            summary['runs'] = {k: joint_summary(v, atoms) for k, v in data['runs'].items()}
            section['datasets'][name] = summary
        for source in ('ols', 'mdpde'):
            section[f'{source}_width_quantiles'] = quantiles([a[source][1] for a in section['atoms']])
        section['validity_pairs'] = {label: sum((a['sources']['OLS']['width_valid'], a['sources']['MDPDE']['width_valid']) == pair for a in section['atoms'])
                                     for label, pair in {'both': (True, True), 'OLS_only': (True, False), 'MDPDE_only': (False, True), 'neither': (False, False)}.items()}
        if not real:
            section['truth'] = truth_summary(section['atoms'])
            section['by_alpha'] = {str(alpha): truth_summary([a for a in section['atoms'] if a['alpha'] == alpha])
                                   for alpha in sorted({a['alpha'] for a in section['atoms']})}
            for source in ('ols', 'mdpde'):
                assert math.isclose(section['truth'][source]['width']['rmse'], previous_neighbor['synthetic'][scenario]['F0'][f'{source}_width_rmse'], rel_tol=1e-14)
        else:
            section['largest_estimator_discrepancies'] = sorted([
                {'id': a['id'], 'ols': a['ols'][1], 'mdpde': a['mdpde'][1], 'difference': a['ols'][1]-a['mdpde'][1]}
                for a in section['atoms'] if a['ols'][1] is not None and a['mdpde'][1] is not None], key=lambda a: abs(a['difference']), reverse=True)[:10]
        output['real' if real else 'synthetic'][scenario] = section
    output_path.write_text(json.dumps(output, separators=(',', ':'), allow_nan=False) + '\n')
    print('Saved', output_path, flush=True)


if __name__ == '__main__':
    main()
