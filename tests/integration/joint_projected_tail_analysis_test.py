#!/usr/bin/env python3
"""Validate the sparse projected-tail census campaign and recorded limits."""
from __future__ import annotations

import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / 'docs/developer/figures/joint-projected-tail-r1'
CASES = ('chain-128', 'cube-128', 'chain-256', 'cube-256', 'chain-512', 'cube-512')
RSS_LIMIT_BYTES = 4 * 1024**3
TIME_LIMIT_SECONDS = 600


def load_case(case):
    report = json.loads((CAMPAIGN / f'{case}-census.json').read_text())
    run = report['execution']['runs'][0]
    details = run['result']['details']
    projected = details['assessment_telemetry']['derivative_reduction_micro_attribution']['projected_reduction']
    return report, run, details, projected


def close(left, right, tolerance=1e-12):
    return abs(left - right) <= tolerance * (1. + abs(left))


def analyze_campaign():
    rows = []
    for case in CASES:
        report, run, details, projected = load_case(case)
        atoms = int(case.rsplit('-', 1)[1])
        n = projected['observations']
        p = projected['free_design_columns']
        m = projected['width_columns']
        q_nnz = projected['q_transformed_nonzeros']
        tail_nnz = projected['tail_nonzeros']
        if (run['status'] != 'completed' or details.get('assessment_execution') != 'completed' or
                report['metadata'].get('assessment_reduction') != 'compact-stack-qr' or
                report['metadata'].get('backend', '').upper() != 'SPQR'):
            raise ValueError(f'{case}: census did not complete with the configured SPQR compact assessment')
        if (run.get('process_wall_seconds', TIME_LIMIT_SECONDS + 1) > TIME_LIMIT_SECONDS or
                run.get('peak_rss_bytes', RSS_LIMIT_BYTES + 1) >= RSS_LIMIT_BYTES):
            raise ValueError(f'{case}: census exceeded the production resource envelope')
        if (projected.get('projected_reduction_kind') != 'projected-tail-census' or
                projected.get('projected_reduction_attempts') != 1 or
                projected.get('projected_reduction_accepted') != 0 or
                projected.get('projected_reduction_fallbacks') != 1 or
                projected.get('ordering') != 'sparse-qmult-free-design-factor' or
                projected.get('projected_reduction_fallback_reason') != 'census-only-observation-tiled-fallback'):
            raise ValueError(f'{case}: census route or observation-tiled fallback was not recorded')
        if n <= p or p != 2 * atoms or m != atoms or projected.get('tail_rows') != n - p:
            raise ValueError(f'{case}: projected-tail dimensions are inconsistent')
        if not q_nnz or not tail_nnz or tail_nnz > q_nnz or projected.get('raw_nonzeros', 0) <= 0:
            raise ValueError(f'{case}: sparse transform returned invalid nonzero counts')
        if not close(projected['raw_density'], projected['raw_nonzeros'] / (n * m)):
            raise ValueError(f'{case}: raw derivative density does not match its shape')
        if not close(projected['q_transformed_density'], q_nnz / (n * m)):
            raise ValueError(f'{case}: Q-transformed density does not match its shape')
        if not close(projected['tail_density'], tail_nnz / ((n - p) * m)):
            raise ValueError(f'{case}: tail density does not match its shape')
        rows.append({
            'case': case,
            'status': run['status'],
            'n': n,
            'p': p,
            'm': m,
            'D_nnz': projected['raw_nonzeros'],
            'D_density': projected['raw_density'],
            'QTD_nnz': q_nnz,
            'QTD_density': projected['q_transformed_density'],
            'tail_rows': n - p,
            'tail_nnz': tail_nnz,
            'tail_density': projected['tail_density'],
            'Q_transform_seconds': projected['q_transform_seconds'],
            'QTD_storage_bytes': projected['q_transformed_storage_bytes'],
            'tail_storage_bytes': projected['tail_storage_bytes'],
            'process_wall_seconds': run['process_wall_seconds'],
            'peak_rss_bytes': run['peak_rss_bytes'],
            'assessment_seconds': details.get('assessment_seconds'),
        })
    return rows


def write_analysis():
    rows = analyze_campaign()
    with (CAMPAIGN / 'projected-tail-summary.csv').open('w', newline='') as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        'schema_version': 1,
        'campaign': 'joint-projected-tail-r1',
        'branch': 'develop',
        'base_commit': 'f955f4643b1b375489ecd13bc7c2edfc95409096',
        'selection': {
            'cases': list(CASES),
            'reason': 'Measure sparse Q-transformed derivative and projected tail before constructing a tail QR factor.',
        },
        'configuration': {
            'profile': 'solve',
            'measurement_scope': 'joint_search_and_returned_state_assessment',
            'backend': 'SPQR',
            'assessment_reduction': 'compact-stack-qr',
            'projected_reduction': 'projected-tail-census',
            'preconditioner': 'Schwarz',
            'eigen_threads': 1,
            'rss_limit_bytes': RSS_LIMIT_BYTES,
            'timeout_seconds': TIME_LIMIT_SECONDS,
        },
        'transform': {
            'factor_source': 'existing free-design factor reused by assessment preparation',
            'representation': 'SuiteSparseQR sparse Q transpose output inspected in place; tail storage estimated from compressed Eigen value, index, and column-pointer bytes',
            'conceptual_dense_transform_materialized': False,
            'tail_matrix_materialized': False,
            'tail_qr_run': False,
            'observation_tiled_qr_retained': True,
        },
        'decision': {
            'production_promotion': False,
            'reason': 'Census only; numerical Projected Tail QR validation and candidate resource measurements are required before any decision.',
        },
    }
    (CAMPAIGN / 'campaign-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    analysis = {
        'schema_version': 1,
        'campaign': 'joint-projected-tail-r1',
        'representation': 'sparse-qmult',
        'measurement_count': len(rows),
        'production_promotion': False,
        'cases': rows,
    }
    (CAMPAIGN / 'projected-tail-analysis.json').write_text(json.dumps(analysis, indent=2) + '\n')
    return rows


class ProjectedTailAnalysisTest(unittest.TestCase):
    def test_all_census_cases_complete_under_the_original_limits(self):
        rows = analyze_campaign()
        self.assertEqual([row['case'] for row in rows], list(CASES))
        for row in rows:
            with self.subTest(case=row['case']):
                self.assertLess(row['process_wall_seconds'], TIME_LIMIT_SECONDS)
                self.assertLess(row['peak_rss_bytes'], RSS_LIMIT_BYTES)
                self.assertGreater(row['tail_nnz'], 0)
                self.assertGreater(row['Q_transform_seconds'], 0.)

    def test_summary_artifacts_match_recomputed_census(self):
        rows = analyze_campaign()
        summary = json.loads((CAMPAIGN / 'projected-tail-analysis.json').read_text())
        self.assertFalse(summary['production_promotion'])
        self.assertEqual(summary['cases'], rows)
        with (CAMPAIGN / 'projected-tail-summary.csv').open(newline='') as stream:
            csv_rows = list(csv.DictReader(stream))
        self.assertEqual([row['case'] for row in csv_rows], list(CASES))
        manifest = json.loads((CAMPAIGN / 'campaign-manifest.json').read_text())
        self.assertFalse(manifest['transform']['conceptual_dense_transform_materialized'])
        self.assertFalse(manifest['transform']['tail_matrix_materialized'])
        self.assertFalse(manifest['transform']['tail_qr_run'])
        self.assertTrue(manifest['transform']['observation_tiled_qr_retained'])


if __name__ == '__main__':
    write_analysis()
    unittest.main()
