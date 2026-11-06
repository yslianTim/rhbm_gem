#!/usr/bin/env python3
"""Sequential fixed-Fibonacci selection study; independent of the P0–P6 ledger."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import statistics as stats
import subprocess
from run_first_stage_sampling_study import shared_equal, summarize

POLICIES = ('F0', 'FO', 'FQ')
ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / 'docs/research/first-stage-sampling'


def average(values):
    return stats.mean(values) if values else None


def contamination(rows):
    if not rows:
        return {'count': 0}
    # Each row carries its atom's epsilon: 1e-12 times its largest raw target signal.
    ratios = [abs(r[12]) / (abs(r[11]) + r[13]) for r in rows]
    return {'count': len(rows), 'neighbor_abs_mean': average([abs(r[12]) for r in rows]),
            'target_abs_mean': average([abs(r[11]) for r in rows]),
            'rho_mean': average(ratios), 'rho_median': stats.median(ratios),
            'rho_max': max(ratios), 'small_denominator_fraction': average([abs(r[11]) <= r[13] for r in rows])}


def selection_summary(atoms):
    rows = []
    universe = []
    per_atom = []
    for atom in atoms:
        observations = atom['neighbor_selection']['observations']
        epsilon = 1e-12 * max(abs(r[11]) for r in observations) if 'truth' in atom else None
        rows.extend([*r, epsilon] if epsilon is not None else r for r in observations)
        universe.extend(atom['neighbor_selection']['shells'])
        if epsilon is not None:
            selected = [[*r, epsilon] for r in observations if r[0] <= 1. and r[{'F0': 5, 'FO': 6, 'FQ': 7}[atom['_policy']]]]
            per_atom.append({'id': atom['id'], 'phase': atom['_phase'], 'width_error': atom['mdpde'][1] - atom['truth'][1],
                             'fitting_contamination': contamination(selected)})
    def counts(rr):
        return {'raw': len(rr), 'f0': sum(r[5] for r in rr), 'fo': sum(r[6] for r in rr),
                'fq': sum(r[7] for r in rr), 'ownership_rejected': sum(r[9] for r in rr),
                'cone_rejected': sum(r[10] for r in rr), 'qscore_rejected': sum(not r[7] for r in rr),
                'fo_new': sum(r[6] and not r[5] for r in rr), 'fq_new': sum(r[7] and not r[5] for r in rr),
                'fo_new_rejections': sum(r[5] and not r[6] for r in rr),
                'fq_new_rejections': sum(r[5] and not r[7] for r in rr),
                'fo_fq_differences': sum(r[6] != r[7] for r in rr),
                'fq_rejected_vs_fo': sum(r[6] and not r[7] for r in rr),
                'shared_qscore_vs_full_differences': sum(r[7] != r[8] for r in rr),
                'shared_inclusion_violations': sum((not r[8] and r[6]) or (not r[6] and r[5]) for r in rr)}
    out = {domain: counts([r for r in rows if domain == 'raw' or r[0] <= 1.]) for domain in ('raw', 'fitting')}
    assert out['raw']['shared_inclusion_violations'] == 0
    assert out['fitting']['fq_new_rejections'] == 0
    assert out['fitting']['fq_rejected_vs_fo'] == 0
    assert out['fitting']['shared_qscore_vs_full_differences'] == 0
    out['per_atom_counts'] = [{'id': a['id'], 'phase': a['_phase'],
        'raw': counts(a['neighbor_selection']['observations']),
        'fitting': counts([r for r in a['neighbor_selection']['observations'] if r[0]<=1.])}
        for a in atoms]
    out['neighbor_universe'] = {domain: {
        'shells': len(ss), 'different_sets': sum(set(s['fixed_ids']) != set(s['qscore_ids']) for s in ss),
        'fixed_only_candidates': sum(len(set(s['fixed_ids']) - set(s['qscore_ids'])) for s in ss),
        'qscore_only_candidates': sum(len(set(s['qscore_ids']) - set(s['fixed_ids'])) for s in ss)}
        for domain in ('raw', 'fitting')
        for ss in [[s for s in universe if domain == 'raw' or s['radius'] <= 1.]]}
    groups = {'fo_retained': lambda r: r[6], 'fq_retained': lambda r: r[7], 'f0_retained': lambda r: r[5], 'fo_new': lambda r: r[6] and not r[5],
              'fq_new': lambda r: r[7] and not r[5], 'fo_rejected': lambda r: not r[6],
              'fq_rejected': lambda r: not r[7], 'fq_new_vs_fo': lambda r: r[7] and not r[6]}
    if rows and len(rows[0]) == 14:
        out['contamination'] = {domain: {name: contamination([r for r in rows if predicate(r) and (domain == 'raw' or r[0] <= 1.)])
            for name, predicate in groups.items()} for domain in ('raw', 'fitting')}
        out['per_atom'] = per_atom
        x = [a['fitting_contamination']['rho_mean'] for a in per_atom]
        y = [abs(a['width_error']) for a in per_atom]
        out['contamination_absolute_width_error_correlation'] = stats.correlation(x, y) if len(x)>1 and len(set(x))>1 and len(set(y))>1 else None
    out['shells'] = []
    for radius in sorted({r[0] for r in rows}):
        rr = [r for r in rows if r[0] == radius]
        shell = {'radius': radius, **counts(rr)}
        if len(rr[0]) == 14:
            shell['contamination'] = {name: contamination([r for r in rr if predicate(r)]) for name, predicate in groups.items()}
        out['shells'].append(shell)
    return out


def paired(reference, candidate):
    differences = []
    for a, b in zip(reference, candidate, strict=True):
        assert (a['_phase'], a['id'], a['truth']) == (b['_phase'], b['id'], b['truth'])
        for x, y in zip(a['neighbor_selection']['observations'], b['neighbor_selection']['observations'], strict=True):
            assert x == y, 'Fixed geometry/response/mask diagnostics changed between policies'
        error_a = abs(a['mdpde'][1] - a['truth'][1])
        error_b = abs(b['mdpde'][1] - b['truth'][1])
        differences.append({'phase': a['_phase'], 'id': a['id'], 'delta': error_b - error_a,
                            'baseline_error': error_a, 'candidate_error': error_b})
    delta = [d['delta'] for d in differences]
    return {'mean': average(delta), 'median': stats.median(delta),
            'improvement_fraction': average([d < 0. for d in delta]),
            'unchanged_fraction': average([d == 0. for d in delta]), 'worst_regression': max(delta),
            'outliers': sorted(differences, key=lambda d: d['delta'], reverse=True)[:8]}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--build', type=Path, default=Path('build/release'))
    parser.add_argument('--real-data', type=Path, default=Path('/Users/yslian/data'))
    parser.add_argument('--real-only', action='store_true')
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    build = (ROOT / args.build).resolve()
    directory = build / 'neighbor-selection-study'
    directory.mkdir(exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_THREAD_LIMIT='4', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    binary = build / 'bin/first_stage_sampling_study'
    old = json.loads((DOCS / 'results.json').read_text())
    output_path = DOCS / 'neighbor-selection-results.json'
    output = json.loads(output_path.read_text()) if output_path.exists() else {'runs': {}, 'synthetic': {}, 'real': {}}
    output['epsilon_definition'] = '1e-12 * max absolute analytic target contribution among the same atom raw samples; small denominator means abs(target)<=epsilon'
    output['phases'] = 'Three deterministic sensitivity phases; no random-sample confidence intervals.'
    output['initial_head'] = '0e9661a81f3480629549361101b1e3609b97686b'
    output['source_commit'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()

    def run(name, scenario, policy, mode='production', paths=()):
        path = directory / f'{name}-{policy}-{mode}.json'
        if not args.summarize_only:
            print(path.name, flush=True)
            with path.with_suffix('.log').open('w') as log:
                subprocess.run([str(binary), str(path), scenario, policy, str(mode), *map(str, paths)],
                               cwd=ROOT, env=env, stdout=log, stderr=log, check=True, timeout=600)
        data = json.loads(path.read_text())
        if policy == 'F0' and mode == 'production':
            baseline = build / 'first-stage-study' / path.name
            if baseline.is_file():
                previous = json.loads(baseline.read_text())
                assert all(shared_equal(previous[k], data[k]) for k in ('atoms', 'calibration', 'joint', 'qscore')), f'F0 parity: {name}'
        if policy == 'F0' and mode == 'production' and baseline.is_file():
            output.setdefault('baseline_exact_parity', {})[name] = True
        summary = summarize(data)
        assert not summary['joint']['failed_checks'], (name, policy, summary['joint']['failed_checks'])
        assert all(summary['joint'][key] for key in ('kkt_passed', 'replay_passed', 'endpoint_passed'))
        summary['initialization_peak_rss_bytes'] = data['initialization_peak_rss_bytes']
        summary['jce_peak_rss_bytes'] = data['jce_peak_rss_bytes']
        def quantiles(values):
            values = sorted(v for v in values if v is not None)
            def quantile(p):
                position = p * (len(values)-1)
                lo = math.floor(position)
                hi = math.ceil(position)
                return values[lo] + (values[hi]-values[lo])*(position-lo)
            return {str(p): quantile(p) for p in (.5, .9, .95, .99)}
        summary['initial_width_quantiles'] = quantiles(data.get('initial_tau', [a['mdpde'][1] for a in data['atoms']]))
        summary['final_width_quantiles'] = quantiles(data['joint']['tau'])
        summary['bounded_projection_count'] = data.get('bounded_projection_count', data['calibration']['initial_violations'])
        output['runs'][path.stem] = summary
        return data

    scenarios = ('6z6u', '8rqb') if args.real_only else ('fine', 'coarse', 'noisy', 'crowded')
    for scenario in scenarios:
        paths = ()
        if args.real_only:
            map_name = 'emd_11103_unsharpened.map' if scenario == '6z6u' else 'emd_19436_unsharpened.map'
            paths = (args.real_data / f'{scenario}.cif', args.real_data / map_name)
            if not all(p.is_file() for p in paths):
                output['real'][scenario] = {'status': 'unavailable'}
                continue
            provenance = {p.name: {'bytes': p.stat().st_size, 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()} for p in paths}
            assert provenance == old['real_data'][f'{scenario}.cif'], 'Real-map source mismatch'
        all_atoms = {p: [] for p in POLICIES}
        for phase in range(1 if args.real_only else 3):
            name = scenario if args.real_only else f'{scenario}-r{phase}'
            for policy in POLICIES:
                data = run(name, 'real' if args.real_only else f'{scenario}:{phase}', policy, paths=paths)
                for atom in data['atoms']:
                    atom['_phase'] = phase
                    atom['_policy'] = policy
                all_atoms[policy].extend(data['atoms'])
            if phase == 0 and (args.real_only or scenario == 'crowded'):
                bound = max(output['runs'][f'{name}-{p}-production']['calibration']['bound'] for p in POLICIES)
                for policy in POLICIES:
                    run(name, 'real' if args.real_only else f'{scenario}:{phase}', policy, bound, paths)
        section = {}
        for policy, atoms in all_atoms.items():
            runs = [output['runs'][f'{scenario if args.real_only else scenario+"-r"+str(phase)}-{policy}-production'] for phase in range(1 if args.real_only else 3)]
            cell = {'selection': selection_summary(atoms), 'initialization_seconds_median': stats.median(r['initialization_seconds'] for r in runs),
                    'raw_mean_per_phase': average([r['raw'] for r in runs]),
                    'selected_fitting_mean_per_phase': average([r['fitting_candidates'] for r in runs]),
                    'positive_fitting_mean_per_phase': average([r['retained'] for r in runs]),
                    'working_width_variance_median': stats.median(a['working_information']['amplitude_width_covariance'][3] for a in atoms if a['working_information']),
                    'donors_per_phase': [r['calibration']['donors'] for r in runs],
                    'fallback_per_phase': [r['calibration']['fallback'] for r in runs],
                    'atoms': [{k: a[k] for k in ('id', '_phase', 'ols', 'mdpde', 'alpha', 'rank', 'retained', 'distinct_fitting_radii', 'qualified', 'reason', 'seed_source')} for a in atoms]}
            if not args.real_only:
                for estimator in ('ols', 'mdpde'):
                    for coordinate, index in (('amplitude', 0), ('width', 1)):
                        errors = [a[estimator][index] - a['truth'][index] for a in atoms]
                        cell[f'{estimator}_{coordinate}_bias'] = average(errors)
                        cell[f'{estimator}_{coordinate}_rmse'] = math.sqrt(average([e*e for e in errors]))
                if policy != 'F0':
                    cell['paired_vs_f0'] = paired(all_atoms['F0'], atoms)
            else:
                for a, b in zip(all_atoms['F0'], atoms, strict=True):
                    assert a['id'] == b['id']
                    assert a['neighbor_selection'] == b['neighbor_selection'], 'Real fixed geometry/response parity'
            section[policy] = cell
        # Store policy-independent masks/counts/contamination groups once per dataset.
        common_selection = section['F0']['selection'].copy()
        policy_fields = ('per_atom', 'contamination_absolute_width_error_correlation')
        for key in policy_fields:
            common_selection.pop(key, None)
        for policy in POLICIES:
            diagnostic = section[policy]['selection']
            section[policy]['selection'] = {k: diagnostic[k] for k in policy_fields if k in diagnostic}
        section['selection_diagnostics'] = common_selection
        section['Q0_external_reference'] = [old['runs'][f'{scenario if args.real_only else scenario+"-r"+str(phase)}-Q0-production'] for phase in range(1 if args.real_only else 3)]
        if args.real_only:
            section['provenance'] = provenance
            audit_path = directory / f'audit-{scenario}.json'
            if not args.summarize_only:
                subprocess.run([str(binary), str(audit_path), 'audit', scenario, str(directory), *map(str, paths)], env=env, cwd=ROOT, check=True, timeout=600, stdout=subprocess.DEVNULL)
            section['endpoint_audit'] = json.loads(audit_path.read_text())
        output['real' if args.real_only else 'synthetic'][scenario] = section
    with output_path.open('w') as stream:
        stream.write('{\n')
        for index, (key, value) in enumerate(output.items()):
            if index:
                stream.write(',\n')
            stream.write('  ' + json.dumps(key) + ': ')
            if key in ('runs', 'synthetic', 'real'):
                stream.write('{\n' + ',\n'.join('    ' + json.dumps(k) + ': ' + json.dumps(v, separators=(',', ':'), allow_nan=False) for k, v in value.items()) + '\n  }')
            else:
                stream.write(json.dumps(value, allow_nan=False))
        stream.write('\n}\n')
    print('Saved', output_path, flush=True)


if __name__ == '__main__':
    main()
