"""Summarize the preregistered acceptance gates from independently audited results."""
import csv
import json
from pathlib import Path
import sys


def read(path):
    return json.loads(path.read_text())


def main():
    root = Path(sys.argv[1]); result = read(root/'results.json')
    gates = {'complete_216_branches': result['case_count'] == 216 and result['matrix_complete'],
             'legacy_frozen_parity': not result['legacy_parity_failures'],
             'legacy_32_qualified': result['counts']['legacy']['legacy_qualified'] == 32,
             'original_qualified_branches_preserved': not result['qualified_baseline_regressions'],
             'guarded_accepted_states_replay': not result['untrusted_accepted_states'],
             'exact_legacy_domain_equivalence': read(root/'domain-equivalence.json')['passed']}
    qualified_replays, main_cases, budgets, negative_controls, original_regular = [], [], [], [], []
    narrower, boundary, precision = [], [], []
    for group in result['results']:
        path = root/'datasets'/group['dataset']/group['variant']
        if group['dataset'] == 'heterogeneous-168':
            main_cases.append(group['recovery_passed'] and group['multistart_consistent'] and
                              all(c['certificate']['regular_qualified'] for c in group['cases']))
        if group['dataset'] in ('zero-signal', 'duplicate'):
            negative_controls.append(group['nonidentifiability_evidence_verified'])
        for case in group['cases']:
            label = case['case']; fit = read(path/'fits'/f'{label}.json')
            audit = read(path/'audits'/f'{label}.json')
            old = read(Path(__file__).resolve().parents[1]/'joint-abc-coverage/datasets'/group['dataset']/'fits'/f'{label}.json')
            if group['variant'] != 'legacy' and old['joint_qualified']:
                original_regular.append(case['certificate']['regular_qualified'])
            if case['certificate']['regular_qualified']:
                raw = read(path/'audits'/f'{label}-raw-replay.json')
                qualified_replays.append(all(v is not None and v['passed'] for v in raw.values()))
            budgets.append(fit.get('profile_evaluations', 0) <= 200 and fit.get('accepted_updates', 0) <= 100)
            if group['dataset'] == 'near-0.02' and label.startswith('narrower'):
                narrower.append({'variant': group['variant'], 'case': label,
                                 'recovered': case['oracle_recovered'] if label.endswith('double') else case['float32_accuracy_passed'],
                                 'regular_qualified': case['certificate']['regular_qualified'],
                                 'profile_evaluations': fit['profile_evaluations'], 'accepted_updates': fit['accepted_updates'],
                                 'stop_reason': fit['stop_reason'],
                                 'untrusted_accepted_evaluations': [t['evaluation'] for t in fit['trials'] if t['accepted'] and not t['trust']['passed']]})
            if 'precision' in audit:
                p = audit['precision']
                precision.append({'dataset': group['dataset'], 'variant': group['variant'], 'case': label,
                                  'agreement_passed': p['agreement_passed'], 'fixed_face_feasible': p.get('fixed_face_feasible'),
                                  'derivative_status': audit['derivative_status'], 'local_correction_inf': p.get('local_correction_inf'),
                                  'local_correction_passed': p.get('local_correction_passed')})
    gates.update({'all_regular_endpoints_independently_replayed': bool(qualified_replays) and all(qualified_replays),
                  'original_32_guarded_branches_have_regular_certificates': len(original_regular) == 64 and all(original_regular),
                  'main_24_branches_unchanged_thresholds': len(main_cases) == 3 and all(main_cases),
                  'all_searches_within_budget': all(budgets),
                  'negative_controls_have_independent_evidence': len(negative_controls) == 6 and all(negative_controls),
                  'guarded_narrower_never_accept_untrusted': all(not r['untrusted_accepted_evaluations'] for r in narrower if r['variant'] != 'legacy')})
    for p in sorted((root/'datasets/active-a').rglob('*boundary*.json')):
        value = read(p)
        for atom in (1, 5, 9):
            rows = [r for r in value['precision100']['rows'] if r['atom'] == atom]
            last = rows[-1]
            boundary.append({'source': str(p.relative_to(root)), 'atom': atom+1,
                             'precision_agreement': value['agreement_passed'], 'steps': len(rows),
                             'valid_constrained_profiles': sum(r['profile_valid'] for r in rows),
                             'smallest_step_prediction_slope': last['prediction_slope'],
                             'first_order_prediction_norm': last['first_order_prediction_norm'],
                             'near_zero_rows': last['near_zero_rows'],
                             'near_zero_first_order_max': last['near_zero_first_order_max']})
    gates['boundary_evidence_complete'] = len(boundary) == 78 and all(r['steps'] == 21 and r['precision_agreement'] for r in boundary)
    reproducibility = root/'reproducibility.json'
    gates['exact_scientific_repeat'] = read(reproducibility)['passed'] if reproducibility.exists() else None
    value = {'passed': all(v is True for v in gates.values()), 'gates': gates,
             'narrower_outcomes': narrower, 'precision_endpoints': precision,
             'interpretation': 'Trustworthy unresolved derivatives/nonregular or unidentifiable controls are retained; 72/72 qualification is not a gate.'}
    (root/'scientific-validation.json').write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    with (root/'boundary-summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(boundary[0])); writer.writeheader(); writer.writerows(boundary)
    print(json.dumps({'passed': value['passed'], 'gates': gates}, indent=2))
    if any(v is False for v in gates.values()): raise SystemExit('Scientific acceptance gate failed.')


if __name__ == '__main__': main()
