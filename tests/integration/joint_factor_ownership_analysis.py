"""Summarize search-only accepted-factor ownership campaigns."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def _close(left, right):
    return (isinstance(left, (int, float)) and isinstance(right, (int, float)) and
            math.isclose(left, right, rel_tol=1e-10, abs_tol=1e-12))


def _run(report):
    details = (report.get('result') or {}).get('details') or {}
    work = details.get('search_work') or {}
    residency = work.get('factor_residency') or details.get('factor_residency') or {}
    policy = (report.get('metadata') or {}).get('solver_policy') or {}
    ownership = (details.get('operator_factor_ownership') or policy.get('operator_factor_ownership'))
    if ownership is None:
        representation = details.get('operator_factor_representation') or policy.get('operator_factor_representation')
        ownership = {'native-qr': 'dedicated-native'}.get(representation, 'dedicated-fixed')
    runs = report.get('execution', {}).get('runs') or []
    stages = [stage for run in runs for stage in run.get('stages', [])]
    rss = [stage.get('sampled_tree_peak_rss_bytes') for stage in stages
           if stage.get('sampled_tree_peak_rss_bytes') is not None]
    factorization = details.get('spqr_factorization') or {}
    symbolic = (factorization.get('symbolic') or {}).get('calls')
    numeric = (factorization.get('numeric') or {}).get('calls')
    return {
        'case': report.get('case'),
        'ownership': ownership,
        'status': report.get('execution', {}).get('status'),
        'search_seconds': details.get('search_seconds'),
        'process_wall_seconds': max((run.get('process_wall_seconds', 0) for run in runs), default=None),
        'peak_rss_bytes': max(rss, default=None),
        'maximum_concurrent_factor_count': residency.get('maximum_concurrent_factor_count'),
        'maximum_concurrent_owned_bytes_estimate': residency.get('maximum_concurrent_owned_bytes_estimate'),
        'maximum_concurrent_factor_stage': residency.get('maximum_concurrent_factor_stage'),
        'symbolic_calls': symbolic,
        'symbolic_reuses': (factorization.get('symbolic') or {}).get('reuses'),
        'symbolic_seconds': (factorization.get('symbolic') or {}).get('seconds'),
        'derived_symbolic_reuses': numeric - symbolic if isinstance(numeric, int) and isinstance(symbolic, int) else None,
        'numeric_calls': numeric,
        'numeric_seconds': (factorization.get('numeric') or {}).get('seconds'),
        'operator_fixed_factor_seconds': work.get('operator_fixed_factor_seconds'),
        'operator_prepare_seconds': work.get('operator_prepare_seconds'),
        'operator_rank_seconds': work.get('operator_rank_seconds'),
        'operator_normal_seconds': work.get('operator_normal_seconds'),
        'preconditioner_factor_seconds': work.get('factor_seconds'),
        'reuse_attempts': work.get('accepted_factor_reuse_attempts'),
        'reuse_accepted': work.get('accepted_factor_reuse_accepted'),
        'reuse_fallbacks': work.get('accepted_factor_reuse_fallbacks'),
        'stop_reason': details.get('stop_reason'),
        'profile_evaluations': details.get('profile_evaluations'),
        'accepted_updates': details.get('accepted_updates'),
        'pcg_iteration_counts': work.get('pcg_iteration_counts'),
        'rank_status': work.get('operator_rank_status'),
        'rank_certificate': work.get('operator_rank_certificate'),
        'rank_reason': work.get('operator_rank_reason'),
        'returned_search_state': details.get('returned_search_state'),
        'accepted_objective': details.get('accepted_objective'),
        'accepted_gradient_inf_norm': details.get('accepted_gradient_inf_norm'),
    }


def trajectory_parity(reference, candidate):
    exact = ('status', 'stop_reason', 'profile_evaluations', 'accepted_updates',
             'pcg_iteration_counts', 'rank_status', 'rank_certificate', 'rank_reason')
    checks = {key: reference.get(key) == candidate.get(key) for key in exact}
    for key in ('accepted_objective', 'accepted_gradient_inf_norm'):
        a, b = reference.get(key), candidate.get(key)
        checks[key] = ((a is None and b is None) or _close(a, b))
    left, right = reference.get('returned_search_state'), candidate.get('returned_search_state')
    checks['returned_search_state'] = (isinstance(left, list) and isinstance(right, list) and
        len(left) == len(right) and all(_close(a, b) for a, b in zip(left, right)))
    return {'passed': all(checks.values()), 'checks': checks}


def analyze(paths):
    rows = [_run(json.loads(Path(path).read_text())) for path in paths]
    comparisons = []
    for case in sorted({row['case'] for row in rows}):
        group = [row for row in rows if row['case'] == case]
        baseline = next((row for row in group if row['ownership'] == 'dedicated-fixed'), None)
        for row in group:
            if row is baseline:
                continue
            comparisons.append({'case': case, 'ownership': row['ownership'],
                                'trajectory_parity': (trajectory_parity(baseline, row)
                                                      if baseline else None)})
    return {'schema_version': 1, 'rows': rows, 'comparisons': comparisons}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    result = analyze(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'rows': len(result['rows']), 'comparisons': len(result['comparisons']),
                      'output': str(args.output)}))


if __name__ == '__main__':
    main()
