"""Classify fixed-policy PCG iteration scaling from Joint sweep results."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from experiment_io import write


GROUP_FIELDS = (
    'topology', 'preconditioner', 'core_atoms', 'overlap_hops',
    'max_block_atoms', 'storage_mib', 'scratch_mib', 'operator_rank',
    'operator_rank_seconds', 'operator_rank_work_entries',
    'operator_rank_workspace_mib', 'sparse_backend', 'measurement_scope',
)
SUPPORT_FIELDS = (
    'pcg_solves', 'linearizations', 'damping_trials', 'operator_normals',
    'operator_applications', 'operator_adjoints', 'setup_seconds', 'pcg_seconds',
    'search_seconds', 'wall_seconds', 'peak_rss_bytes',
)
MINIMUM_POINTS = 4
GROWTH_RATIO_THRESHOLD = 1.5
GROWTH_SLOPE_THRESHOLD = .25


def configuration_key(row):
    config = row.get('configuration') or {}
    return {
        'topology': config.get('topology', row.get('topology')),
        'preconditioner': config.get('preconditioner', row.get('preconditioner')),
        'core_atoms': config.get('core_atoms', row.get('core_atoms')),
        'overlap_hops': config.get('overlap_hops', row.get('overlap_hops')),
        'max_block_atoms': config.get('max_block_atoms', row.get('max_block_atoms')),
        'operator_rank': config.get('operator_rank', row.get('operator_rank')),
        'sparse_backend': row.get('sparse_backend'),
        'storage_mib': config.get('storage_mib', row.get('storage_mib')),
        'scratch_mib': config.get('scratch_mib', row.get('scratch_mib')),
        'operator_rank_seconds': config.get('operator_rank_seconds', row.get('operator_rank_seconds')),
        'operator_rank_work_entries': config.get('operator_rank_work_entries', row.get('operator_rank_work_entries')),
        'operator_rank_workspace_mib': config.get('operator_rank_workspace_mib', row.get('operator_rank_workspace_mib')),
        'measurement_scope': config.get('measurement_scope', row.get('measurement_scope')),
    }


def iterations_per_solve(row):
    if row.get('pcg_iteration_telemetry_complete') is not True:
        return None
    measurements = row.get('measurements') or {}
    per_solve = measurements.get('pcg_iterations_per_solve') or {}
    if per_solve.get('median') is not None:
        return per_solve['median']
    return row.get('pcg_iterations_per_solve', row.get('iterations_per_solve'))


def measurement_point(row):
    config = row.get('configuration') or {}
    atoms = row.get('atoms') if row.get('atoms') is not None else config.get('atoms')
    counts = row.get('stop_reasons_by_repetition') or []
    budget_repetitions = row.get('pcg_iteration_budget_repetitions')
    if budget_repetitions is None:
        budget_repetitions = sum(reason == 'pcg-iteration-budget' for reason in counts)
        if not counts and row.get('stop_reason') == 'pcg-iteration-budget':
            budget_repetitions = 1
    point = {'atoms': atoms, 'iterations_per_solve': iterations_per_solve(row),
             'block_count': row.get('blocks', row.get('block_count')),
             'pcg_iteration_budget_repetitions': budget_repetitions,
             'stop_reasons_by_repetition': counts}
    point.update({field: row.get(field) for field in SUPPORT_FIELDS})
    point['valid_iteration_point'] = (
        isinstance(atoms, (int, float)) and atoms > 0 and
        isinstance(point['iterations_per_solve'], (int, float)) and point['iterations_per_solve'] > 0)
    return point


def loglog_slope(points, value_field):
    values = [(math.log(point['atoms']), math.log(point[value_field])) for point in points
              if isinstance(point.get('atoms'), (int, float)) and point['atoms'] > 0 and
              isinstance(point.get(value_field), (int, float)) and point[value_field] > 0]
    if len(values) < 2:
        return None
    mean_x = sum(x for x, _ in values) / len(values)
    mean_y = sum(y for _, y in values) / len(values)
    denominator = sum((x - mean_x) ** 2 for x, _ in values)
    if denominator == 0:
        return None
    return sum((x - mean_x) * (y - mean_y) for x, y in values) / denominator


def endpoint_growth_ratio(points, value_field):
    values = [point[value_field] for point in points
              if isinstance(point.get('atoms'), (int, float)) and point['atoms'] > 0 and
              isinstance(point.get(value_field), (int, float)) and point[value_field] > 0]
    return values[-1] / values[0] if len(values) >= 2 else None


def adjacent_scaling(points):
    out = []
    for left, right in zip(points, points[1:]):
        p0, p1 = left['atoms'], right['atoms']
        i0, i1 = left['iterations_per_solve'], right['iterations_per_solve']
        if p1 <= p0:
            continue
        ratio = i1 / i0
        out.append({
            'from_atoms': p0,
            'to_atoms': p1,
            'iteration_growth_ratio': ratio,
            'iteration_loglog_slope': math.log(ratio) / math.log(p1 / p0),
            'approximately_doubling': abs(p1 / p0 / 2 - 1) <= .15,
        })
    return out


def classify_group(key, rows):
    points = []
    excluded = []
    for row in rows:
        if row.get('_evidence_excluded'):
            excluded.append({
                'atoms': row.get('atoms'),
                'reason': 'evidence-ineligible',
                'evidence_exclusion_reasons': row.get('evidence_exclusion_reasons') or [],
            })
            continue
        config = row.get('configuration') or {}
        atoms = row.get('atoms') if row.get('atoms') is not None else config.get('atoms')
        complete = (
            row.get('status') == 'completed' and row.get('measurements_complete') is True and
            row.get('requested_repetitions') == row.get('completed_repetitions'))
        if not complete:
            excluded.append({'atoms': atoms, 'reason': 'incomplete-measurements'})
            continue
        point = measurement_point(row)
        points.append(point)
    points.sort(key=lambda point: (
        not isinstance(point['atoms'], (int, float)),
        point['atoms'] if isinstance(point['atoms'], (int, float)) else str(point['atoms'])))
    valid_points = [point for point in points if point['valid_iteration_point']]
    atom_counts = [point['atoms'] for point in points if isinstance(point['atoms'], (int, float))]
    duplicate_atoms = len(atom_counts) != len(set(atom_counts))
    missing_policy = [field for field in GROUP_FIELDS if key.get(field) is None]
    not_comparable = (key.get('preconditioner') != 'schwarz' or bool(missing_policy) or duplicate_atoms)

    iteration_slope = loglog_slope(valid_points, 'iterations_per_solve')
    iteration_growth = endpoint_growth_ratio(valid_points, 'iterations_per_solve')
    supporting = {
        field: {
            'loglog_slope': loglog_slope(points, field),
            'growth_ratio': endpoint_growth_ratio(points, field),
        }
        for field in ('pcg_solves', 'linearizations', 'damping_trials')
    }
    adjacent = adjacent_scaling(valid_points)
    larger_points = [point for point in points if atom_counts and
                     isinstance(point['atoms'], (int, float)) and point['atoms'] > min(atom_counts)]
    budget_limited = any(point['pcg_iteration_budget_repetitions'] for point in larger_points)

    if not_comparable:
        gate = 'not-comparable'
    elif budget_limited:
        gate = 'pcg-budget-limited'
    elif len(valid_points) < MINIMUM_POINTS:
        gate = 'insufficient-evidence'
    elif iteration_slope >= GROWTH_SLOPE_THRESHOLD and iteration_growth >= GROWTH_RATIO_THRESHOLD:
        gate = 'growth-observed'
    else:
        gate = 'stable'

    diagnostics = []
    nonlinear_growth = any((summary['growth_ratio'] or 0) >= GROWTH_RATIO_THRESHOLD
                           for summary in supporting.values())
    if gate == 'stable' and nonlinear_growth:
        diagnostics.append('nonlinear-work-growth')
    cost_growth = any((endpoint_growth_ratio(points, field) or 0) >= GROWTH_RATIO_THRESHOLD
                      for field in ('wall_seconds', 'peak_rss_bytes'))
    if gate == 'stable' and cost_growth:
        diagnostics.append('cost-growth-with-stable-krylov')
    if duplicate_atoms:
        diagnostics.append('duplicate-global-size-points')
    if missing_policy:
        diagnostics.append('missing-comparison-policy:' + ','.join(missing_policy))
    if budget_limited:
        diagnostics.append('coarse-correction-investigation-warranted')

    if gate in ('growth-observed', 'pcg-budget-limited'):
        interpretation = 'coarse correction investigation is warranted'
    elif gate == 'stable':
        interpretation = 'the available points show no clear PCG iteration growth'
    elif gate == 'insufficient-evidence':
        interpretation = 'more complete global sizes are needed for a scaling classification'
    else:
        interpretation = 'these rows do not share a comparable one-level Schwarz policy'

    multi_block_rows = [point for point in points
                        if isinstance(point.get('block_count'), (int, float)) and point['block_count'] >= 2]
    multi_valid_points = [point for point in multi_block_rows if point['valid_iteration_point']]
    multi_atoms = [point['atoms'] for point in multi_block_rows
                   if isinstance(point.get('atoms'), (int, float))]
    multi_larger = [point for point in multi_block_rows
                    if multi_atoms and isinstance(point.get('atoms'), (int, float)) and
                    point['atoms'] > min(multi_atoms)]
    multi_budget_limited = any(point['pcg_iteration_budget_repetitions'] for point in multi_larger)
    multi_slope = loglog_slope(multi_valid_points, 'iterations_per_solve')
    multi_growth = endpoint_growth_ratio(multi_valid_points, 'iterations_per_solve')
    if not_comparable:
        multi_gate = 'not-comparable'
    elif multi_budget_limited:
        multi_gate = 'pcg-budget-limited'
    elif len(multi_valid_points) < MINIMUM_POINTS:
        multi_gate = 'insufficient-evidence'
    elif multi_slope >= GROWTH_SLOPE_THRESHOLD and multi_growth >= GROWTH_RATIO_THRESHOLD:
        multi_gate = 'growth-observed'
    else:
        multi_gate = 'stable'

    return {
        **key,
        'points': points,
        'excluded_rows': excluded,
        'valid_point_count': len(valid_points),
        'iteration_loglog_slope': iteration_slope,
        'iteration_growth_ratio': iteration_growth,
        'adjacent_scaling': adjacent,
        'doubling_ratios': [item for item in adjacent if item['approximately_doubling']],
        'supporting_scaling': supporting,
        'coarse_gate': gate,
        'multi_block_valid_point_count': len(multi_valid_points),
        'multi_block_point_sizes': [point['atoms'] for point in multi_block_rows],
        'multi_block_counts': [point['block_count'] for point in multi_block_rows],
        'multi_block_iteration_loglog_slope': multi_slope,
        'multi_block_iteration_growth_ratio': multi_growth,
        'multi_block_gate': multi_gate,
        'diagnostics': diagnostics,
        'interpretation': interpretation,
    }


def analyze(document):
    if not isinstance(document, dict):
        raise ValueError('input must be a JSON object')
    if document.get('tool') != 'joint_schwarz_sweep':
        raise ValueError('input must be a joint_schwarz_sweep aggregate')
    rows = document.get('configurations')
    if not isinstance(rows, list):
        raise ValueError('input is missing the configurations list')
    grouped = {}
    keys_by_group = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        row = dict(row)
        if row.get('measurement_scope') is None:
            row['measurement_scope'] = ('search-only' if document.get('measurement_profile') == 'search' else
                                        'joint_search_and_returned_state_assessment')
        key = configuration_key(row)
        identity = tuple(key[field] for field in GROUP_FIELDS)
        if row.get('evidence_eligible') is False:
            config = row.get('configuration') or {}
            atoms = row.get('atoms') if row.get('atoms') is not None else config.get('atoms')
            excluded_row = dict(row)
            excluded_row['_evidence_excluded'] = True
            grouped.setdefault(identity, [])
            keys_by_group[identity] = key
            excluded_row['atoms'] = atoms
            grouped[identity].append(excluded_row)
            continue
        grouped.setdefault(identity, []).append(row)
        keys_by_group[identity] = key
    groups = [classify_group(keys_by_group[identity], grouped[identity])
              for identity in sorted(grouped, key=lambda item: tuple(str(value) for value in item))]
    return {'schema_version': 1, 'tool': 'joint_scaling_analysis', 'groups': groups}


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True, help='joint_schwarz_sweep aggregate JSON')
    parser.add_argument('--output', type=Path, required=True, help='Analysis JSON output path')
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        document = json.loads(args.input.read_text())
        result = analyze(document)
    except (OSError, json.JSONDecodeError, ValueError) as error:
        parser.error(str(error))
    write(args.output.resolve(), result)
    print(json.dumps({'status': 'completed', 'groups': len(result['groups']),
                      'output': args.output.name}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
