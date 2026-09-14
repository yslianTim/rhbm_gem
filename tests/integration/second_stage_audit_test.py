#!/usr/bin/env python3
import copy
import importlib.util
import json
from pathlib import Path
import unittest
import subprocess
import sys
import tempfile

SOURCE = Path(__file__).resolve().parents[2] / 'resources/tools/developer/second_stage_audit.py'
spec = importlib.util.spec_from_file_location('audit', SOURCE)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def batch():
    return {'stages': {}, 'anomalies': {'total': 0, 'shown': 0, 'omitted': 0, 'categories': {}, 'details': []}}


def records():
    skipped = {'executed': False, 'result': 'skipped', 'reason': 'no-ordinary-components'}
    iteration = dict(kind='iteration', attempt=1, accepted_clusters=0, rejected_clusters=2,
                     selection_audit={'ordinary': skipped, 'after_rescue': copy.deepcopy(skipped)},
                     convergence={'reference': 'iteration_previous', 'status': 'not_evaluated'}, **batch())
    terminal = dict(kind='terminal', attempt=1, stop_reason='all-rejected-backtracking-exhausted',
                    final_state_source='best-audit', final_polish={'attempted': True, 'objective_accepted': True,
                    'operator_certified': False, 'applied': False}, **batch())
    return [{'kind': 'start', 'settings': {}}, iteration, terminal]


def log(rows, schema=1):
    return '\n'.join(f'[Debug] Second-stage audit: schema={schema}, payload={json.dumps(row)}' for row in rows)


class AuditParserTest(unittest.TestCase):
    def test_actual_statuses_and_incomplete_run(self):
        parsed = audit.parse(log(records()))
        self.assertTrue(parsed['complete'])
        self.assertFalse(parsed['terminal']['final_polish']['applied'])
        self.assertIn('not_evaluated', audit.report(parsed))
        self.assertFalse(audit.parse(log(records()[:-1]))['complete'])

    def test_salvage_empty_and_unavailable_are_not_success(self):
        rows = records()
        rows[1]['selection_audit']['ordinary'] = {'executed': True, 'result': 'empty_after_salvage', 'reason': 'no-selection-remains'}
        rows[1]['selection_audit']['after_rescue'] = {'executed': False, 'result': 'unavailable', 'reason': 'previous-objective-unavailable'}
        self.assertIn('empty_after_salvage', audit.report(audit.parse(log(rows))))
        rows[1]['selection_audit']['after_rescue']['result'] = 'passed'
        with self.assertRaisesRegex(ValueError, 'unexecuted'):
            audit.parse(log(rows))

    def test_report_keeps_missing_candidate_and_operator_references(self):
        rows = records()
        rows[1]['score'] = {'reference': 'iteration_previous', 'candidate': {'value': None, 'reason': 'nonfinite'},
                            'best': {'value': {'total': 123.0}, 'reason': None}}
        text = audit.report(audit.parse(log(rows)))
        self.assertIn('unavailable: nonfinite', text)
        self.assertIn('iteration_previous', text)
        self.assertNotIn('123', text)

    def test_schema_malformed_json_and_nonfinite_rejected(self):
        with self.assertRaisesRegex(ValueError, 'unsupported'):
            audit.parse(log(records(), 2))
        with self.assertRaises(ValueError):
            audit.parse('Second-stage audit: schema=1, payload={broken')
        rows = records(); rows[0]['value'] = float('inf')
        with self.assertRaisesRegex(ValueError, 'numeric constant'):
            audit.parse(log(rows))
        rows[0]['value'] = None; rows[0]['reason'] = 'nonfinite'
        self.assertTrue(audit.parse(log(rows))['complete'])

    def test_bounded_counts_and_terminal_semantics(self):
        rows = records(); rows[1]['anomalies']['details'] = [{}] * 6
        with self.assertRaisesRegex(ValueError, 'limit'):
            audit.parse(log(rows))
        rows = records(); rows[1]['anomalies']['total'] = 2
        with self.assertRaisesRegex(ValueError, 'mismatch'):
            audit.parse(log(rows))
        rows = records(); rows[-1]['final_polish']['applied'] = True
        with self.assertRaisesRegex(ValueError, 'Applied polish'):
            audit.parse(log(rows))

    def test_no_legacy_parser_or_cross_reference_substitution(self):
        with self.assertRaisesRegex(ValueError, 'No current'):
            audit.parse('Second-stage phase audit: schema=2, payload={}')
        rows = records(); rows[1]['convergence']['reference'] = 'final_candidate'
        with self.assertRaisesRegex(ValueError, 'reference'):
            audit.parse(log(rows))
        self.assertEqual(len(audit.parse('[Info] ordinary progress\n' + log(records()))['iterations']), 1)

    def test_cli_writes_only_current_report_and_json(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'run.log'; source.write_text(log(records()))
            subprocess.run([sys.executable, str(SOURCE), str(source), '--output-dir', str(root / 'out')], check=True)
            self.assertEqual(sorted(p.name for p in (root / 'out').iterdir()), ['audit.json', 'report.md'])
            self.assertTrue(json.loads((root / 'out/audit.json').read_text())['complete'])

    def test_documented_current_commands_and_retired_sources(self):
        root = SOURCE.parents[3]
        guide = (root / 'docs/developer/second-stage-audit.md').read_text()
        self.assertIn('second_stage_audit.py', guide)
        self.assertIn('RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT', guide)
        self.assertIn('stopped by user decision', guide)
        for name in ['PhaseAudit', 'TrustModelAudit', 'ClusterHistoryObserver']:
            self.assertFalse((root / f'src/core/detail/second_stage/observation/{name}.cpp').exists())
            self.assertNotIn(name + '.cpp', (root / 'src/CMakeLists.txt').read_text())
            self.assertNotIn(name + '.cpp', (root / 'tests/CMakeLists.txt').read_text())
        self.assertFalse((root / 'resources/tools/developer/second_stage_phase_audit.py').exists())
        self.assertIn('Second-Stage Local Fitting Summary', (root / 'src/core/detail/second_stage/observation/SecondStageLogging.cpp').read_text())
        self.assertIn('Local-fitting atom cutoff: atoms=', (root / 'src/core/detail/second_stage/observation/SecondStageLogging.cpp').read_text())

    def test_numerical_sources_do_not_access_audit_payloads(self):
        source_dir = SOURCE.parents[3] / 'src/core/detail/second_stage'
        for source in source_dir.glob('*.cpp'):
            with self.subTest(source=source.name):
                text = source.read_text()
                self.assertNotRegex(text, r'\b(?:AuditEvent|AuditBatch|AuditCategory|AuditStage|SecondStageAuditData)\b')
                self.assertNotRegex(text, r'(?:\.|->)\s*(?:Audit|Record)\s*\(')
                self.assertNotIn('RecordJointMemberRejection', text)



if __name__ == '__main__':
    unittest.main()
