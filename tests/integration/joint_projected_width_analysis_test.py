#!/usr/bin/env python3
"""Validate projected-width acceptance campaign reports and decision gates."""

import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / 'docs/developer/figures/joint-projected-width-r1'
CURRENT_512 = ROOT / 'docs/developer/figures/joint-compact-acceptance-r2/post-promotion-summary_runs'
CASES = ('chain-128', 'cube-128', 'chain-256', 'cube-256')


def load_run(path):
    report = json.loads(path.read_text())
    run = report['execution']['runs'][0]
    details = run.get('result', {}).get('details', {})
    assessment = details.get('assessment_telemetry', {})
    micro = details.get('derivative_reduction_micro_attribution', {})
    projected = micro.get('projected_reduction', {})
    return report, run, details, assessment, projected


def close(left, right, tolerance=1e-10):
    return abs(left - right) <= tolerance * (1. + abs(left))


def analyze_campaign():
    rows = []
    for case in CASES:
        current_path = CAMPAIGN / f'{case}-current.json'
        candidate_path = CAMPAIGN / f'{case}-candidate.json'
        current, current_run, current_details, current_assessment, _ = load_run(current_path)
        candidate, candidate_run, candidate_details, candidate_assessment, candidate_projected = load_run(candidate_path)
        current_numbers = current_run['numerics']
        candidate_numbers = candidate_run['numerics']
        parity = (
            current_run['status'] == candidate_run['status'] == 'completed' and
            current_numbers['rank'] == candidate_numbers['rank'] and
            current_numbers['runtime_convergence'] == candidate_numbers['runtime_convergence'] and
            close(current_numbers['objective'], candidate_numbers['objective']) and
            current_details['stop_reason'] == candidate_details['stop_reason'] and
            current_details['accepted_updates'] == candidate_details['accepted_updates'] and
            current_details['profile_evaluations'] == candidate_details['profile_evaluations'] and
            current_details['endpoint_valid'] == candidate_details['endpoint_valid'] and
            current_details['assessment_execution'] == candidate_details['assessment_execution'] == 'completed'
        )
        if not parity:
            raise ValueError(f'{case}: current/candidate solve outcomes differ')
        if (candidate_projected.get('projected_reduction_kind') != 'structured-compact-qr' or
                candidate_projected.get('projected_reduction_attempts') != 1 or
                candidate_projected.get('projected_reduction_accepted') != 1 or
                candidate_projected.get('projected_reduction_fallbacks') != 0):
            raise ValueError(f'{case}: structured route was not accepted cleanly')
        current_stages = current_assessment['completed_stage_seconds']
        candidate_stages = candidate_assessment['completed_stage_seconds']
        current_projected = current_stages['derivative-projected-qr']
        candidate_projected_seconds = candidate_projected['projected_reduction_seconds']
        rows.append({
            'case': case,
            'current_status': current_run['status'],
            'candidate_status': candidate_run['status'],
            'current_assessment_seconds': current_details['assessment_seconds'],
            'candidate_assessment_seconds': candidate_details['assessment_seconds'],
            'assessment_speedup': current_details['assessment_seconds'] / candidate_details['assessment_seconds'],
            'current_projected_qr_seconds': current_projected,
            'candidate_projected_reduction_seconds': candidate_projected_seconds,
            'projected_speedup': current_projected / candidate_projected_seconds,
            'current_reference_seconds': current_stages.get('reference-evaluation'),
            'candidate_reference_seconds': candidate_stages.get('reference-evaluation'),
            'current_design_spectrum_seconds': current_stages.get('design-spectrum'),
            'candidate_design_spectrum_seconds': candidate_stages.get('design-spectrum'),
            'current_derivative_preparation_seconds': current_stages.get('derivative-preparation'),
            'candidate_derivative_preparation_seconds': candidate_stages.get('derivative-preparation'),
            'current_derivative_reduction_seconds': current_stages.get('derivative-reduction'),
            'candidate_derivative_reduction_seconds': candidate_stages.get('derivative-reduction'),
            'current_rows_seconds': current_stages.get('derivative-rows'),
            'candidate_rows_seconds': candidate_stages.get('derivative-rows'),
            'current_compact_jacobian_qr_seconds': current_stages.get('derivative-compact-jacobian-qr'),
            'candidate_compact_jacobian_qr_seconds': candidate_stages.get('derivative-compact-jacobian-qr'),
            'current_width_spectrum_seconds': current_stages.get('projected-width-spectrum'),
            'candidate_width_spectrum_seconds': candidate_stages.get('projected-width-spectrum'),
            'current_normalized_width_spectrum_seconds': current_stages.get('normalized-width-spectrum'),
            'candidate_normalized_width_spectrum_seconds': candidate_stages.get('normalized-width-spectrum'),
            'current_correction_seconds': current_stages.get('correction-jacobian-spectrum'),
            'candidate_correction_seconds': candidate_stages.get('correction-jacobian-spectrum'),
            'current_process_wall_seconds': current_run['process_wall_seconds'],
            'candidate_process_wall_seconds': candidate_run['process_wall_seconds'],
            'current_peak_rss_bytes': current_run['peak_rss_bytes'],
            'candidate_peak_rss_bytes': candidate_run['peak_rss_bytes'],
            'parity': parity,
            'candidate_fallbacks': candidate_projected['projected_reduction_fallbacks'],
            'candidate_active_stage': None,
            'candidate_sparse_rows': candidate_projected.get('sparse_rows'),
            'candidate_sparse_columns': candidate_projected.get('sparse_columns'),
            'candidate_sparse_nnz': candidate_projected.get('sparse_nnz'),
            'candidate_symbolic_seconds': candidate_projected.get('symbolic_seconds'),
            'candidate_numeric_seconds': candidate_projected.get('numeric_seconds'),
            'candidate_ordering': candidate_projected.get('ordering'),
            'candidate_maximum_dense_bytes': candidate_projected.get('maximum_dense_bytes'),
        })

    for case in ('chain-512', 'cube-512'):
        current_path = next(CURRENT_512.glob(f'{case}-*-compact-stack-qr-solve-r1-w0.json'))
        candidate_path = CAMPAIGN / f'{case}-candidate.json'
        current, current_run, current_details, current_assessment, _ = load_run(current_path)
        candidate, candidate_run, candidate_details, candidate_assessment, candidate_projected = load_run(candidate_path)
        if current_run['status'] != 'completed' or current_details.get('assessment_execution') != 'completed':
            raise ValueError(f'{case}: current production reference did not complete')
        if (current['metadata']['solver_policy']['preconditioner'] != 'schwarz' or
                current['metadata']['solver_policy']['operator_rank_mode'] != 'auto' or
                current['metadata']['solver_policy']['schwarz_core_atoms'] != 128 or
                current['metadata']['solver_policy']['schwarz_overlap_hops'] != 1 or
                current['metadata']['solver_policy']['schwarz_max_block_atoms'] != 512 or
                current['metadata']['assessment_reduction'] != 'compact-stack-qr'):
            raise ValueError(f'{case}: current reference settings do not match the campaign')
        if (candidate_run['status'] != 'rss_limit' or
                candidate_details.get('active_assessment_stage') != 'derivative-projected-qr' or
                candidate_assessment.get('active_assessment_stage') != 'derivative-projected-qr' or
                candidate_projected.get('projected_reduction_attempts') != 1 or
                candidate_projected.get('projected_reduction_accepted') != 0):
            raise ValueError(f'{case}: candidate did not hit the expected projected-QR RSS frontier')
        current_stages = current_assessment['completed_stage_seconds']
        candidate_stages = candidate_assessment.get('completed_stage_seconds', {})
        rows.append({
            'case': case,
            'current_status': current_run['status'],
            'candidate_status': candidate_run['status'],
            'current_assessment_seconds': current_details['assessment_seconds'],
            'candidate_assessment_seconds': None,
            'assessment_speedup': None,
            'current_projected_qr_seconds': current_stages['derivative-projected-qr'],
            'candidate_projected_reduction_seconds': None,
            'projected_speedup': None,
            'current_reference_seconds': current_stages.get('reference-evaluation'),
            'candidate_reference_seconds': candidate_stages.get('reference-evaluation'),
            'current_design_spectrum_seconds': current_stages.get('design-spectrum'),
            'candidate_design_spectrum_seconds': candidate_stages.get('design-spectrum'),
            'current_derivative_preparation_seconds': current_stages.get('derivative-preparation'),
            'candidate_derivative_preparation_seconds': candidate_stages.get('derivative-preparation'),
            'current_derivative_reduction_seconds': current_stages.get('derivative-reduction'),
            'candidate_derivative_reduction_seconds': candidate_stages.get('derivative-reduction'),
            'current_rows_seconds': current_stages.get('derivative-rows'),
            'candidate_rows_seconds': candidate_stages.get('derivative-rows'),
            'current_compact_jacobian_qr_seconds': current_stages.get('derivative-compact-jacobian-qr'),
            'candidate_compact_jacobian_qr_seconds': candidate_stages.get('derivative-compact-jacobian-qr'),
            'current_width_spectrum_seconds': current_stages.get('projected-width-spectrum'),
            'candidate_width_spectrum_seconds': candidate_stages.get('projected-width-spectrum'),
            'current_normalized_width_spectrum_seconds': current_stages.get('normalized-width-spectrum'),
            'candidate_normalized_width_spectrum_seconds': candidate_stages.get('normalized-width-spectrum'),
            'current_correction_seconds': current_stages.get('correction-jacobian-spectrum'),
            'candidate_correction_seconds': candidate_stages.get('correction-jacobian-spectrum'),
            'current_process_wall_seconds': current_run['process_wall_seconds'],
            'candidate_process_wall_seconds': candidate_run['process_wall_seconds'],
            'current_peak_rss_bytes': current_run['peak_rss_bytes'],
            'candidate_peak_rss_bytes': candidate_run['peak_rss_bytes'],
            'parity': None,
            'candidate_fallbacks': candidate_projected['projected_reduction_fallbacks'],
            'candidate_active_stage': candidate_details['active_assessment_stage'],
            'candidate_sparse_rows': candidate_projected.get('sparse_rows'),
            'candidate_sparse_columns': candidate_projected.get('sparse_columns'),
            'candidate_sparse_nnz': candidate_projected.get('sparse_nnz'),
            'candidate_symbolic_seconds': candidate_projected.get('symbolic_seconds'),
            'candidate_numeric_seconds': candidate_projected.get('numeric_seconds'),
            'candidate_ordering': candidate_projected.get('ordering'),
            'candidate_maximum_dense_bytes': candidate_projected.get('maximum_dense_bytes'),
        })
    return rows


def write_analysis():
    rows = analyze_campaign()
    with (CAMPAIGN / 'projected-width-summary.csv').open('w', newline='') as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    analysis = {
        'schema_version': 1,
        'campaign': 'joint-projected-width-r1',
        'acceptance_cases': rows[:4],
        'frontier_cases': rows[4:],
        'production_promotion': False,
        'promotion_blocker': 'candidate exceeded the fixed 4 GiB RSS envelope at derivative-projected-qr on both 512-atom topologies',
    }
    (CAMPAIGN / 'projected-width-analysis.json').write_text(json.dumps(analysis, indent=2) + '\n')
    return rows


class ProjectedWidthAnalysisTest(unittest.TestCase):
    def test_128_256_assessment_and_projected_route_parity(self):
        rows = analyze_campaign()[:4]
        self.assertEqual([row['case'] for row in rows], list(CASES))
        for row in rows:
            with self.subTest(case=row['case']):
                self.assertTrue(row['parity'])
                self.assertGreater(row['assessment_speedup'], 1.)
                self.assertGreater(row['projected_speedup'], 1.)
                self.assertGreater(row['candidate_peak_rss_bytes'], row['current_peak_rss_bytes'])
                self.assertEqual(row['candidate_fallbacks'], 0)

    def test_512_frontier_is_rss_limited_during_projected_qr(self):
        rows = analyze_campaign()[4:]
        self.assertEqual([row['case'] for row in rows], ['chain-512', 'cube-512'])
        for row in rows:
            with self.subTest(case=row['case']):
                self.assertEqual(row['candidate_status'], 'rss_limit')
                self.assertEqual(row['candidate_active_stage'], 'derivative-projected-qr')
                self.assertEqual(row['candidate_fallbacks'], 0)
                self.assertGreater(row['candidate_peak_rss_bytes'], 4 * 1024**3)

    def test_summary_artifacts_match_recomputed_analysis(self):
        rows = analyze_campaign()
        summary = json.loads((CAMPAIGN / 'projected-width-analysis.json').read_text())
        self.assertFalse(summary['production_promotion'])
        with (CAMPAIGN / 'projected-width-summary.csv').open(newline='') as stream:
            csv_rows = list(csv.DictReader(stream))
        self.assertEqual(len(csv_rows), len(rows))
        self.assertEqual([row['case'] for row in csv_rows], [row['case'] for row in rows])
        self.assertEqual(summary['promotion_blocker'],
            'candidate exceeded the fixed 4 GiB RSS envelope at derivative-projected-qr on both 512-atom topologies')


if __name__ == '__main__':
    write_analysis()
    unittest.main()
