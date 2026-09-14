#!/usr/bin/env python3
"""Compare test-only numerical probe logs: BASELINE OFF ON (no external datasets)."""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import shlex

SOURCE = Path(__file__).resolve().parents[2] / 'resources/tools/developer/second_stage_audit.py'
spec = importlib.util.spec_from_file_location('audit', SOURCE)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def read(path):
    result = {}
    for line in path.read_text().splitlines():
        if not line.startswith('NUMERICAL '):
            continue
        parts = shlex.split(line)
        key = tuple(parts[1:4])
        work = tuple(map(int, parts[5:9]))
        payload = parts[9:]
        if 'SUMMARY' in payload:
            pos = payload.index('SUMMARY') + 1
            # The old baseline contains additional informational timing records.
            # Compare the six stable basic final-summary fields only.
            payload[pos] = re.findall(r' - (accepted_iterations|best_iteration|stop_reason|best_audit_objective|final_uses_polish|final_state_source) = (\S+)', payload[pos])
        result[key] = work, payload
    if len(result) != 12:
        raise ValueError(f'{path}: expected 12 probe cases, got {len(result)}')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline', 'off', 'on'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    baseline, off, on = map(read, (args.baseline, args.off, args.on))
    if baseline.keys() != off.keys() or off.keys() != on.keys():
        raise ValueError('Probe input matrices differ')
    for key in baseline:
        if not baseline[key][1] == off[key][1] == on[key][1]:
            raise ValueError(f'Numerical decisions, summary or persisted state differ: {key}')
        if off[key][0] != on[key][0]:
            raise ValueError(f'ON added solver/operator/objective/snapshot work: {key}')
    if any(audit.MARKER.search(line) for line in args.off.read_text().splitlines()):
        raise ValueError('OFF emitted extra audit records')
    runs, current = [], []
    for line in args.on.read_text().splitlines():
        match = audit.MARKER.search(line)
        if not match:
            continue
        row = json.loads(match[2])
        if row['kind'] == 'start' and current:
            runs.append(audit.parse('\n'.join(current))); current = []
        current.append(line)
    if current:
        runs.append(audit.parse('\n'.join(current)))
    if len(runs) != 6 or not all(run['complete'] for run in runs):
        raise ValueError('Expected six complete non-quiet ON audit runs')
    # Debug scheduling and fixed key order make these independent of worker completion.
    for left, right in ((0,1), (2,3), (4,5)):
        if runs[left]['iterations'] != runs[right]['iterations']:
            raise ValueError('Parallel audit ordering or totals differ')
    print(json.dumps({'cases': 12, 'baseline_off_on_decisions_and_persistence': 'exact',
                      'on_off_work_counts': 'equal', 'complete_audit_runs': len(runs),
                      'parallel_iteration_records': 'equal'}, indent=2))


if __name__ == '__main__':
    main()
