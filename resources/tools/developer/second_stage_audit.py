#!/usr/bin/env python3
"""Summarize passive second-stage decision records (schema 1)."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re

MARKER = re.compile(r'Second-stage audit: schema=(\d+), payload=(.*)$')


def reject_constant(value):
    raise ValueError(f'Non-JSON numeric constant: {value}')


def validate_batch(record):
    anomalies = record['anomalies']
    details = anomalies['details']
    if len(details) > 5 or anomalies['shown'] != len(details):
        raise ValueError('Audit detail limit/count mismatch')
    if anomalies['total'] != anomalies['shown'] + anomalies['omitted']:
        raise ValueError('Audit total/shown/omitted mismatch')
    if anomalies['total'] != sum(value['total'] for value in anomalies['categories'].values()):
        raise ValueError('Audit category totals do not match')
    if any(type(value) is not int or value < 0 for value in
           [anomalies['total'], anomalies['shown'], anomalies['omitted'], *(count for category in anomalies['categories'].values() for count in category.values())]):
        raise ValueError('Audit counts must be nonnegative integers')
    for category in anomalies['categories'].values():
        if category['total'] != category['shown'] + category['omitted']:
            raise ValueError('Audit category shown/omitted mismatch')
    for stage in record['stages'].values():
        if stage['total'] != stage['accepted'] + stage['rejected'] + stage['skipped']:
            raise ValueError('Stage totals do not match')


def parse(text):
    records = []
    for line_number, line in enumerate(text.splitlines(), 1):
        match = MARKER.search(line)
        if not match:
            continue
        if match[1] != '1':
            raise ValueError(f'Line {line_number}: unsupported audit schema {match[1]}')
        try:
            record = json.loads(match[2], parse_constant=reject_constant)
            if record['kind'] not in ('start', 'iteration', 'terminal'):
                raise ValueError('Unknown record kind')
            if record['kind'] != 'start':
                validate_batch(record)
            if record['kind'] == 'iteration':
                for audit in record['selection_audit'].values():
                    if type(audit['executed']) is not bool:
                        raise ValueError('Selection audit execution must be explicit')
                    if audit['result'] not in ('skipped', 'passed', 'rejected', 'unavailable', 'empty_after_salvage'):
                        raise ValueError('Unknown selection audit result')
                    if not audit['executed'] and audit['result'] in ('passed', 'rejected', 'empty_after_salvage'):
                        raise ValueError('An unexecuted gate cannot pass or reject')
                if record['convergence']['reference'] != 'iteration_previous':
                    raise ValueError('Outer operator reference must be iteration_previous')
            if record['kind'] == 'terminal':
                polish = record['final_polish']
                if polish['applied'] and not (polish['attempted'] and polish['objective_accepted'] and polish['operator_certified'] is True):
                    raise ValueError('Applied polish requires both production acceptance and certification')
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f'Line {line_number}: {error}') from error
        records.append(record)
    if not records:
        raise ValueError('No current audit records; enable SECOND_STAGE_AUDIT, Debug logging, and non-quiet execution')
    if records[0]['kind'] != 'start' or sum(r['kind'] == 'start' for r in records) != 1:
        raise ValueError('Expected one run starting with a start record')
    terminals = [r for r in records if r['kind'] == 'terminal']
    if len(terminals) > 1 or (terminals and records[-1]['kind'] != 'terminal'):
        raise ValueError('Terminal record must appear once, at the end')
    iterations = [r for r in records if r['kind'] == 'iteration']
    attempts = [r['attempt'] for r in iterations]
    if attempts != sorted(set(attempts)):
        raise ValueError('Iteration attempts must be unique and ordered')
    return {'schema': 1, 'complete': bool(terminals), 'start': records[0],
            'iterations': iterations, 'terminal': terminals[0] if terminals else None}


def report(audit):
    lines = ['# Second-stage decision audit', '',
             'Only evidence calculated by production is recorded. Objectives belong to the stated scope, reference, and frozen environment; this report does not compare different environments.', '',
             'Run status: ' + ('complete' if audit['complete'] else 'incomplete (no terminal record; interruption or observation failure is possible)'), '',
             '| Attempt | Accepted / rejected clusters | Candidate objective | Ordinary audit | After rescue | Certificate | Anomalies shown / total |',
             '|---|---:|---:|---|---|---|---:|']
    for row in audit['iterations']:
        selection = row['selection_audit']
        def status(item):
            return ('executed' if item['executed'] else 'skipped') + ': ' + item['result']
        anomaly = row['anomalies']
        score = row.get('score', {}).get('candidate', {})
        objective = score.get('value')
        total = objective.get('total') if isinstance(objective, dict) else None
        score_text = format(total, '.6g') if total is not None else 'unavailable: ' + str(score.get('reason', 'not_evaluated'))
        lines.append(f"| {row['attempt']} | {row['accepted_clusters']} / {row['rejected_clusters']} | {score_text} | "
                     f"{status(selection['ordinary'])} | {status(selection['after_rescue'])} | "
                     f"{row['convergence']['status']} | {anomaly['shown']} / {anomaly['total']} |")
    terminal = audit['terminal']
    if terminal:
        polish = terminal['final_polish']
        lines += ['', f"Stop reason: {terminal['stop_reason']}. Final state: {terminal['final_state_source']}.",
                  'Final polish: ' + ', '.join(f"{name}={polish[name]}" for name in
                     ('attempted', 'objective_accepted', 'operator_certified', 'applied')) + '.']
    lines += ['', '## Convergence and score references', '',
              '| Attempt | Score reference / objective revision | Accepted-active p99 | Nominal-operator p99 | Complete / qualified | Blockers |',
              '|---|---|---|---|---|---|']
    for row in audit['iterations']:
        certificate = row['convergence']
        def residual(name):
            values = certificate.get(name)
            return 'not_evaluated' if values is None else ', '.join('null' if x is None else format(x, '.4g') for x in values)
        blockers = certificate.get('blockers')
        blocker_text = 'not_evaluated' if blockers is None else (', '.join(k for k, v in blockers.items() if v) or 'none')
        lines.append(f"| {row['attempt']} | {row.get('score', {}).get('reference', 'unavailable')} / "
                     f"{row.get('objective_revision', 'unavailable')} | {residual('accepted_active_p99')} | "
                     f"{residual('operator_nominal_p99')} | {certificate.get('complete', 'not_evaluated')} / "
                     f"{certificate.get('qualified', 'not_evaluated')} | {blocker_text} |")
    for row in [*audit['iterations'], *([terminal] if terminal else [])]:
        details = row['anomalies']['details']
        if details:
            lines += ['', f"## {row['kind']} {row['attempt']}: bounded details", '']
            for event in details:
                lines.append(f"- {event['stage']}, first atom {event['first_atom']}, trial {event['trial']}: "
                             f"{event['outcome']} ({event['reason']}); scope={event['scope']}, reference={event['reference']}.")
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('log', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    try:
        audit = parse(args.log.read_text())
    except ValueError as error:
        parser.error(str(error))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'audit.json').write_text(json.dumps(audit, indent=2, allow_nan=False) + '\n')
    (args.output_dir / 'report.md').write_text(report(audit))


if __name__ == '__main__':
    main()
