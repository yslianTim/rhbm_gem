"""Collect the fixed weak-halo identifiability diagnosis across three starts."""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
from experiment_io import digest, read, sha, write


def collect_start(executable, work_dir, start):
    raw = work_dir / f'start-{start}.json'
    process = subprocess.run([str(executable), 'weak', str(start), str(raw)],
                             capture_output=True, text=True)
    run = dict(start=start, execution=dict(status='completed' if process.returncode == 0 else 'process_error',
                                           exit_code=process.returncode))
    if process.returncode != 0:
        run['stderr'] = process.stderr
        return run, None

    record = read(raw)
    snapshot = record.pop('snapshot')
    snapshot_hash = digest(snapshot)
    record['snapshot_sha256'] = snapshot_hash
    record['raw_record_sha256'] = sha(raw)
    state = record['outcome']['assembled_state']
    record['actual_state_identity_verified'] = None
    if 'actual_state_assessment' in record:
        assessed = record['actual_state_assessment']['primary']
        assert state is not None and assessed['beta'] == state['ac'] and assessed['eta'] == state['log_b']
        record['actual_state_identity_verified'] = True
    record['parameter_errors'] = None
    if state is not None:
        ac = np.array(state['ac']).reshape(2, 2)
        parameters = np.column_stack([ac[:, 0], state['b'], ac[:, 1]])
        truth = np.column_stack([snapshot['truth_a'], snapshot['truth_b'], snapshot['truth_c']])
        record['parameter_errors'] = dict(roles=['target', 'halo'], order=['A', 'B', 'C'],
                                          values=(parameters - truth).tolist())
    restart = record.get('diagnostic_restart')
    if restart:
        record['diagnostic_restart'] = {key: restart.get(key) for key in (
            'runtime_convergence', 'stop_reason', 'lm_status', 'profile_evaluations',
            'accepted_updates', 'last_trusted_state', 'primary', 'local_correction_inf',
            'runtime_checks', 'settings', 'usable_state')}
    record['execution'] = run['execution']
    record['diagnostics_applicable'] = state is not None
    record['diagnostic_complete'] = record.get('audit_complete', False) or state is None
    return record, snapshot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    work_dir = args.work_dir.resolve()
    executable = args.executable.resolve()
    output = args.output.resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    if any((work_dir / f'start-{start}.json').exists() for start in range(3)):
        raise SystemExit('Use a fresh work directory to preserve existing diagnostic output.')

    records = []
    snapshots = []
    for start in range(3):
        record, snapshot = collect_start(executable, work_dir, start)
        records.append(record)
        if snapshot is not None:
            snapshots.append(snapshot)
    hashes = {digest(snapshot) for snapshot in snapshots}
    if len(hashes) > 1:
        raise AssertionError('Weak diagnostic inputs changed between starts.')
    if snapshots:
        write(output.with_name('weak-snapshot.json'), snapshots[0])
    executions_complete = len(records) == 3 and all(r['execution']['status'] == 'completed' for r in records)
    cases_complete = executions_complete and all(r.get('diagnostic_complete', False) for r in records)
    report = dict(schema_version=1, diagnostic='weak-halo-identifiability',
                  metadata=dict(executable_sha256=sha(executable),
                                source_sha256=sha(Path(__file__)),
                                native_source_sha256=sha(Path(__file__).parents[1] / 'experiments' /
                                                         'joint_offline_diagnostic.cpp')),
                  execution=dict(status='completed' if executions_complete else 'process_error',
                                 completed_starts=sum(r['execution']['status'] == 'completed' for r in records),
                                 attempted_starts=len(records)),
                  snapshot_sha256=next(iter(hashes), None), runs=records,
                  collection_complete=cases_complete,
                  completion_note='Collection completion is separate from runtime convergence and identifiability.')
    write(output, report)
    print(json.dumps(report['execution'], indent=2))
    if not executions_complete:
        raise SystemExit('One or more offline diagnostic processes did not complete.')


if __name__ == '__main__':
    main()
