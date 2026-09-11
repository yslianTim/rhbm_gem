#!/usr/bin/env python3
"""Parse observation-only phase audit schema 1; never compare across domains."""
import argparse
import json
import math
import re
from pathlib import Path

PATTERN = re.compile(r'Second-stage phase audit(?P<counters> counters)?: schema=(\d+), payload=(\{.*\})')


def parse(text):
    events, counters = [], []
    for line in text.splitlines():
        match = PATTERN.search(line)
        if not match:
            continue
        if match[2] != '1':
            raise ValueError(f'Unsupported phase audit schema: {match[2]}')
        record = json.loads(match[3])
        (counters if match['counters'] else events).append(record)
    by_id = {e['candidate_id']: e for e in events}
    if len(by_id) != len(events):
        raise ValueError('Duplicate candidate ID')
    for event in events:
        parent_id = event['parent_id']
        if not parent_id:
            continue
        if parent_id not in by_id:
            raise ValueError(f'Missing parent: {parent_id}')
        parent = by_id[parent_id]
        if (event['domain_id'], event['attempt']) != (parent['domain_id'], parent['attempt']):
            raise ValueError('Cross-domain or cross-attempt comparison')
        if event['objective'] is not None and parent['objective'] is not None:
            expected = event['objective']['total'] - parent['objective']['total']
            if abs(expected - event['delta_parent']) > 1e-12 * max(1, abs(expected)):
                raise ValueError('Parent objective delta mismatch')
    return events, counters


def compare_sample(sample, reference, reproduced):
    result = dict(status='unavailable', objective_delta=None, p99_delta=None,
                  residual_max=None, residual_max_delta=None,
                  objective_and_max_improve=None, objective_and_all_improve=None)
    if reference is None:
        return result
    if sample.get('objective') is not None and reference.get('objective') is not None:
        result['objective_delta'] = sample['objective']['total'] - reference['objective']['total']
    def qualified(operator):
        return (operator and operator.get('status') == 'available' and
                operator.get('complete') and operator.get('solver_qualified') and
                operator.get('p99') is not None and len(operator['p99']) == 3 and
                all(x is not None and math.isfinite(x) for x in operator['p99']))
    op, ref_op = sample.get('operator'), reference.get('operator')
    if not reproduced or not qualified(op) or not qualified(ref_op):
        return result
    delta = result['objective_delta']
    if delta is None or not math.isfinite(delta):
        return result
    changes = [a-b for a, b in zip(op['p99'], ref_op['p99'])]
    max_delta = max(op['p99']) - max(ref_op['p99'])
    result.update(status='available', p99_delta=changes, residual_max=max(op['p99']),
                  residual_max_delta=max_delta, objective_and_max_improve=delta < 0 and max_delta < 0,
                  objective_and_all_improve=delta < 0 and all(x < 0 for x in changes))
    return result


def summarize(events):
    by_id = {e['candidate_id']: e for e in events}
    baselines = {(e['attempt'], e['domain_id']): e for e in events if e['stage'] == 'baseline'}
    attempts = {}
    for event in events:
        attempt = attempts.setdefault(event['attempt'], {
            'attempt': event['attempt'], 'domain_id': event['domain_id'],
            'operator_reproduced': event['operator_reproduced'], 'stages': {},
            'opposite_changes': [], 'oversized_steps': [], 'no_sampled_descent': [],
            'gate_blocked_corrections': [], 'local_improvement_assembly_worsens': False, 'correction_samples': []})
        attempt['stages'].setdefault(event['stage'], []).append(event['candidate_id'])
        samples = event['direction_samples']
        if samples and all(s['delta_parent'] is not None for s in samples):
            if samples[0]['delta_parent'] >= 0 and any(s['delta_parent'] < 0 for s in samples[1:]):
                attempt['oversized_steps'].append(event['candidate_id'])
            if all(s['delta_parent'] >= 0 for s in samples):
                attempt['no_sampled_descent'].append(event['candidate_id'])
        parent = by_id.get(event['parent_id'])
        if event['stage'] == 'boundary-correction' and 5 <= event['attempt'] <= 8:
            baseline = baselines.get((event['attempt'], event['domain_id']))
            for sample in samples:
                attempt['correction_samples'].append(dict(
                    candidate_id=event['candidate_id'], attempt=event['attempt'], domain_id=event['domain_id'],
                    alpha=sample['alpha'], objective=sample.get('objective'), operator=sample.get('operator'),
                    member_gates=sample.get('member_gates'),
                    parent=compare_sample(sample, parent, event['operator_reproduced']),
                    baseline=compare_sample(sample, baseline, event['operator_reproduced'])))
        op, parent_op = event['operator'], parent['operator'] if parent else None
        if (event['operator_reproduced'] and op and parent_op and
                op['status'] == parent_op['status'] == 'available' and event['delta_parent'] is not None):
            differences = [a-b for a, b in zip(op['p99'], parent_op['p99'])]
            opposite = [i for i, delta in enumerate(differences) if delta * event['delta_parent'] < 0]
            if opposite:
                attempt['opposite_changes'].append({'candidate_id': event['candidate_id'],
                    'objective_delta': event['delta_parent'], 'p99_delta': differences, 'coordinates': opposite})
        if ('correction' in event['stage'] and event['disposition'] == 'rejected' and
                event['delta_parent'] is not None and event['delta_parent'] < 0):
            attempt['gate_blocked_corrections'].append({'candidate_id': event['candidate_id'],
                'reason': event['reason'], 'delta_parent': event['delta_parent'],
                'delta_baseline': event['delta_baseline']})
    for attempt in attempts.values():
        local = [by_id[i] for i in attempt['stages'].get('local-polish', [])]
        assembly = [by_id[i] for i in attempt['stages'].get('assembly-after-polish', [])]
        attempt['local_improvement_assembly_worsens'] = (
            any(e['disposition'] == 'accepted' and e['delta_parent'] is not None and e['delta_parent'] < 0 for e in local)
            and any(e['delta_parent'] is not None and e['delta_parent'] > 0 for e in assembly))
    return list(attempts.values())


def report(events, counters, summaries):
    by_id = {e['candidate_id']: e for e in events}
    def label(identifier):
        e = by_id[identifier]
        key = e['key']
        return e['stage'] + (f" (first={key[0]}, n={len(key)})" if key else '') + ' #' + identifier.rsplit('/', 1)[-1]

    lines = ['# Second-stage phase audit', '',
        'Objective deltas are within one attempt and frozen domain. Residual coordinates: log peak height, log width, offset-to-peak ratio.',
        'Facts: objectives and isolated T(x) − x. Direction samples are observations, not a proof of descent or stationarity.',
        'Final retention means member parameters equal the final selection; it does not mean a candidate was committed independently.', '',
        f'Candidates: {len(events)}; diagnostic time: {sum(c["elapsed_ms"] for c in counters)/1000:.3f} s; '
        f'failures: {sum(c["failures"] for c in counters)}.', '']
    probes = [p for summary in summaries for p in summary['correction_samples']]
    if probes:
        probes.sort(key=lambda p: (not (p['attempt'] in (7, 8) and p['alpha'] == .5),
                                   p['attempt'], p['candidate_id'], -p['alpha']))
        fmt = lambda n: 'unavailable' if n is None else f'{n:.9g}'
        flag = lambda value: 'unavailable' if value is None else 'yes' if value else 'no'
        lines += ['## Boundary correction direction samples', '',
                  'R∞ = max of the three p99 coordinates. Positive gate margin is slack; strict improvement requires margin > 0.',
                  'Missing optional probe evidence in older logs means not measured (未量測). Member gates are not complete production acceptance.', '',
                  '| Attempt / α | Objective | P99 (H, W, C) | Members | Global strict | Blocked or unavailable members |',
                  '|---|---:|---|---|---|---|']
        for probe in probes:
            gates = probe['member_gates']
            blockers = []
            if gates:
                for member in gates['members']:
                    for name in ('previous_gate', 'best_gate'):
                        gate = member[name]
                        if gate['status'] not in ('pass', 'not-applicable'):
                            blockers.append(f"{member['key']} {name}={gate['status']} margin={fmt(gate['margin'])}")
            op = probe['operator']
            residual = ', '.join(fmt(x) for x in op['p99']) if op and op['status'] == 'available' else '未量測/unavailable'
            lines.append(f"| {probe['attempt']} / {probe['alpha']} | {fmt(probe['objective']['total'] if probe['objective'] else None)} | "
                         f"{residual} | {gates['status'] if gates else '未量測'} | "
                         f"{gates['global_strict_improvement']['status'] if gates else '未量測'} | " + '; '.join(blockers) + ' |')
        lines += ['', '| Attempt / α | Reference | Δ objective | Δ p99 (H, W, C) | Δ R∞ | J and R∞ improve | J and all improve |',
                  '|---|---|---:|---|---:|---|---|']
        for probe in probes:
            for name in ('parent', 'baseline'):
                comparison = probe[name]
                changes = ', '.join(fmt(x) for x in comparison['p99_delta']) if comparison['p99_delta'] is not None else 'unavailable'
                lines.append(f"| {probe['attempt']} / {probe['alpha']} | {name} | {fmt(comparison['objective_delta'])} | {changes} | "
                             f"{fmt(comparison['residual_max_delta'])} | {flag(comparison['objective_and_max_improve'])} | "
                             f"{flag(comparison['objective_and_all_improve'])} |")
        lines += ['']
    for summary in sorted(summaries, key=lambda s: (not 5 <= s['attempt'] <= 8, s['attempt'])):
        lines += [f'## Attempt {summary["attempt"]} (domain {summary["domain_id"]})', '',
            f'Operator reproduced: {summary["operator_reproduced"]}.', '']
        if not summary['operator_reproduced']:
            lines += ['**Evidence insufficient for objective/fixed-point consistency: baseline operator did not reproduce or was unavailable.**', '']
        lines += ['| Stage | Total | Δ baseline | Δ parent | P99 (H, W, R) | Disposition |',
                  '|---|---:|---:|---:|---|---|']
        for stage in ('baseline', 'post-joint-offset', 'unrestricted-proposal', 'production-proposal', 'assembly-before-polish',
                      'assembly-after-polish', 'boundary-endpoint', 'boundary-correction', 'rescue-correction',
                      'boundary-final', 'final-selection'):
            for candidate_id in summary['stages'].get(stage, []):
                e = by_id[candidate_id]
                fmt = lambda n: 'unavailable' if n is None else f'{n:.10g}'
                residual = ', '.join(fmt(n) for n in e['operator']['p99']) if e['operator'] and e['operator']['status'] == 'available' else 'unavailable'
                lines.append(f'| {stage} | {fmt(e["objective"]["total"] if e["objective"] else None)} | '
                    f'{fmt(e["delta_baseline"])} | {fmt(e["delta_parent"])} | {residual} | {e["disposition"]} {e["reason"]} |')
        lines += ['', f'Local improvement with assembly deterioration: {summary["local_improvement_assembly_worsens"]}.', '',
                  'Sampled oversized steps: ' + ', '.join(map(label, summary['oversized_steps'])) + '.', '',
                  'No sampled descent (not a mathematical proof): ' + ', '.join(map(label, summary['no_sampled_descent'])) + '.', '',
                  'Objective/residual changes in opposite directions:', '']
        lines += ['- ' + label(case['candidate_id']) + ': Δ objective=' + f"{case['objective_delta']:.6g}" + ', Δ P99=' + ', '.join(f'{x:.6g}' for x in case['p99_delta']) for case in summary['opposite_changes']] or ['None with qualified, reproduced evidence.']
        lines += ['', 'Rejected corrections improving the global objective relative to their parent:', '']
        lines += ['- ' + label(case['candidate_id']) + ': ' + case['reason'] + f"; Δ parent={case['delta_parent']:.6g}, Δ baseline={case['delta_baseline']:.6g}" for case in summary['gate_blocked_corrections']] or ['None.']
        lines += ['']
    lines += ['## Interpretation limits', '',
        'A smaller sampled step improving the objective supports investigating step size. No sampled descent supports inspecting the update objective, but does not establish a non-descent direction.',
        'Local improvements followed by a worse assembly support investigating coupling between clusters. A rejected globally improving correction warrants inspecting the recorded member gates; it does not justify relaxing them.',
        'When the baseline operator does not reproduce, resolve that discrepancy before interpreting objective/residual consistency.']
    return '\n'.join(lines) + '\n'


def compatibility(text, events):
    by_id = {e['candidate_id']: e for e in events}
    pattern = re.compile(r'Second-stage (?P<kind>solver|compatibility) audit: schema=1, payload=(\{.*\})')
    solvers, observations = [], []
    for match in pattern.finditer(text):
        record = json.loads(match[2])
        if match['kind'] == 'solver':
            solvers.append(record)
            continue
        event = by_id.get(record['candidate_id'])
        if not event:
            raise ValueError('Compatibility record has no phase event')
        record['stage'] = event['stage']
        record['key'] = event['key']
        record['operator_reproduced'] = event['operator_reproduced']
        if record['kind'] == 'direction':
            parent = by_id.get(event['parent_id'])
            record['objective'] = event['objective']
            record['parent_objective'] = parent['objective'] if parent else None
            record['delta_parent'] = event['delta_parent']
            record['operator'] = event['operator']
            record['parent_operator'] = parent['operator'] if parent else None
            classification = record['derivative']['classification']
            delta = event['delta_parent']
            record['sampled_oversized_step'] = classification == 'descent' and delta is not None and delta > 0
            record['audit_ascent_direction'] = classification == 'ascent'
        observations.append(record)
    gradients = []
    for attempt in (5, 8):
        rows = [r for r in observations if r['kind'] == 'gradient' and r['attempt'] == attempt]
        if not rows:
            continue
        stable = [r for r in rows if r['derivative']['classification'] in ('ascent', 'descent')]
        stable.sort(key=lambda r: abs(r['derivative']['samples'][-1]['slope']), reverse=True)
        event = by_id[rows[0]['candidate_id']]
        gradients.append(dict(attempt=attempt, components=len(rows), stable_components=len(stable),
                              uncertain_or_unavailable=len(rows)-len(stable),
                              largest_components=[dict(atom_index=r['atom_index'], coordinate=r['coordinate'],
                                  slope=r['derivative']['samples'][-1]['slope']) for r in stable[:10]],
                              operator=event['operator'], operator_reproduced=event['operator_reproduced']))
    return solvers, dict(observations=observations, gradients=gradients)


def compatibility_report(solvers, data):
    lines = ['# Objective / fixed-point compatibility diagnostics', '',
             'Coordinates: log height, log width, physical offset divided by the fixed parent peak height.',
             'Central differences use h=1e-3, 3e-4, 1e-4. No invalid perturbation is clipped.',
             'Stability requires common sign, <=10% spread, and signal above the roundoff estimate at all three h.',
             'A solved frozen surrogate does not imply nonlinear stationarity or an audit minimum.', '']
    if not data['observations']:
        return '\n'.join(lines + ['Not measured (未量測).']) + '\n'
    lines += ['## Directions', '', '| Attempt | Stage / key first | Δ objective | Small-step derivative | Slope (smallest h) | Oversized step |',
              '|---|---|---:|---|---:|---|']
    for row in data['observations']:
        if row['kind'] != 'direction':
            continue
        derivative = row['derivative']
        samples = derivative['samples']
        slope = samples[-1]['slope'] if samples else None
        lines.append(f"| {row['attempt']} | {row['stage']} / {row['key'][:1]} | {row['delta_parent']} | "
                     f"{derivative['classification']} | {slope} | {row['sampled_oversized_step']} |")
    lines += ['', '## Full gradients', '']
    for gradient in data['gradients']:
        lines += [f"### Attempt {gradient['attempt']}", '',
                  f"Components: {gradient['components']}; stable: {gradient['stable_components']}; "
                  f"uncertain/unavailable: {gradient['uncertain_or_unavailable']}.", '',
                  '| Atom index | Coordinate | Derivative |', '|---:|---:|---:|']
        for row in gradient['largest_components']:
            lines.append(f"| {row['atom_index']} | {row['coordinate']} | {row['slope']:.9g} |")
        lines += ['', 'Operator top atoms: ' + json.dumps(gradient['operator']['top_atoms']) + '.', '']
    lines += ['## Solver observations', '', f'Records: {len(solvers)}. Raw equations and support observations are in solver_audit.json.',
              'Sources named production are actual iteration solves; final-selection sources are isolated T(x) evaluations.',
              'Used normal residuals divide by max(1, ||right-hand side||₂). MDPDE refreshed residuals only diagnose the returned parameters; they never update them.',
              'Variance updates and shape-support changes must be checked separately from beta solver status.']
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('log', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    log_text = args.log.read_text()
    events, counters = parse(log_text)
    if not events:
        parser.error('No phase audit events (requires trace build, Debug logging, non-quiet mode)')
    summaries = summarize(events)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, data in [('candidates.json', events), ('attempts.json', summaries), ('counters.json', counters),
                       ('direction_samples.json', [p for s in summaries for p in s['correction_samples']])]:
        (args.output_dir / name).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    (args.output_dir / 'report.md').write_text(report(events, counters, summaries))
    solvers, diagnostics = compatibility(log_text, events)
    for name, data in [('solver_audit.json', solvers), ('compatibility.json', diagnostics)]:
        (args.output_dir / name).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    (args.output_dir / 'compatibility.md').write_text(compatibility_report(solvers, diagnostics))


if __name__ == '__main__':
    main()
