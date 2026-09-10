#!/usr/bin/env python3
"""Parse observation-only phase audit schema 1; never compare across domains."""
import argparse
import json
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


def summarize(events):
    by_id = {e['candidate_id']: e for e in events}
    attempts = {}
    for event in events:
        attempt = attempts.setdefault(event['attempt'], {
            'attempt': event['attempt'], 'domain_id': event['domain_id'],
            'operator_reproduced': event['operator_reproduced'], 'stages': {},
            'opposite_changes': [], 'oversized_steps': [], 'no_sampled_descent': [],
            'gate_blocked_corrections': [], 'local_improvement_assembly_worsens': False})
        attempt['stages'].setdefault(event['stage'], []).append(event['candidate_id'])
        samples = event['direction_samples']
        if samples and all(s['delta_parent'] is not None for s in samples):
            if samples[0]['delta_parent'] >= 0 and any(s['delta_parent'] < 0 for s in samples[1:]):
                attempt['oversized_steps'].append(event['candidate_id'])
            if all(s['delta_parent'] >= 0 for s in samples):
                attempt['no_sampled_descent'].append(event['candidate_id'])
        parent = by_id.get(event['parent_id'])
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
    for summary in sorted(summaries, key=lambda s: (not 5 <= s['attempt'] <= 8, s['attempt'])):
        lines += [f'## Attempt {summary["attempt"]} (domain {summary["domain_id"]})', '',
            f'Operator reproduced: {summary["operator_reproduced"]}.', '']
        if not summary['operator_reproduced']:
            lines += ['**Evidence insufficient for objective/fixed-point consistency: baseline operator did not reproduce or was unavailable.**', '']
        lines += ['| Stage | Total | Δ baseline | Δ parent | P99 (H, W, R) | Disposition |',
                  '|---|---:|---:|---:|---|---|']
        for stage in ('baseline', 'unrestricted-proposal', 'production-proposal', 'assembly-before-polish',
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('log', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    events, counters = parse(args.log.read_text())
    if not events:
        parser.error('No phase audit events (requires trace build, Debug logging, non-quiet mode)')
    summaries = summarize(events)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, data in [('candidates.json', events), ('attempts.json', summaries), ('counters.json', counters)]:
        (args.output_dir / name).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    (args.output_dir / 'report.md').write_text(report(events, counters, summaries))


if __name__ == '__main__':
    main()
