import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('phase_audit', Path(__file__).resolve().parents[2] / 'resources/tools/developer/second_stage_phase_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def event(identifier='base', parent='', total=1.0, delta=None, stage='baseline'):
    return dict(attempt=6, domain_id=2, candidate_id=identifier, parent_id=parent,
                stage=stage, key=[], factor=1, disposition='observed', reason='', final_retained=False,
                objective=dict(fit=total, tail_weighted=0, offset=0, total=total),
                delta_parent=delta, delta_baseline=delta, operator_reproduced=True,
                operator=dict(status='available', p99=[.01, .02, .03]))


def log(events):
    return '\n'.join('[Debug] Second-stage phase audit: schema=2, payload=' + json.dumps(e) for e in events)


class PhaseAuditTest(unittest.TestCase):
    def test_precision_and_parent(self):
        base = event(total=1.2345678901234567)
        child = event('child', 'base', .12345678901234567, .12345678901234567-base['objective']['total'])
        parsed, _ = audit.parse(log([child, base]))
        self.assertEqual(parsed, [child, base])

    def test_cross_domain_attempt_missing_and_duplicate_rejected(self):
        for field, value in [('domain_id', 3), ('attempt', 7), ('parent_id', 'missing')]:
            child = event('child', 'base', .9, -.1)
            child[field] = value
            with self.assertRaises(ValueError):
                audit.parse(log([event(), child]))
        with self.assertRaises(ValueError):
            audit.parse(log([event(), event()]))

    def test_schema_and_delta_validated(self):
        for schema in (1, 3):
            for kind in ('audit:', 'audit counters:'):
                with self.subTest(schema=schema, kind=kind), self.assertRaisesRegex(ValueError, 'Unsupported phase audit schema'):
                    audit.parse(log([event()]).replace('schema=2', f'schema={schema}').replace('audit:', kind))
        with self.assertRaises(ValueError):
            audit.parse(log([event(), event('child', 'base', .9, -.2)]))

    def test_unavailable_or_not_reproduced_cannot_support_consistency(self):
        child = event('child', 'base', .9, -.1)
        child['operator']['p99'] = [.02, .03, .04]
        self.assertTrue(audit.summarize([event(), child])[0]['opposite_changes'])
        child['operator_reproduced'] = False
        self.assertFalse(audit.summarize([event(), child])[0]['opposite_changes'])
        child['operator_reproduced'] = True
        child['operator']['status'] = 'unavailable'
        self.assertFalse(audit.summarize([event(), child])[0]['opposite_changes'])

    def test_rejected_global_correction_and_assembly(self):
        polish = event('polish', 'base', .9, -.1, 'local-polish')
        polish['disposition'] = 'accepted'
        assembly = event('assembly', 'base', 1.1, .1, 'assembly-after-polish')
        correction = event('correction', 'assembly', 1.05, -.05, 'boundary-correction')
        correction.update(disposition='rejected', reason='member-best', delta_baseline=.05)
        summary = audit.summarize([event(), polish, assembly, correction])[0]
        self.assertTrue(summary['local_improvement_assembly_worsens'])
        self.assertEqual(summary['gate_blocked_corrections'][0]['reason'], 'member-best')

    def test_empty_disabled_log(self):
        self.assertEqual(audit.parse('ordinary production log'), ([], []))

    def test_cli_writes_current_reports_in_attempt_order(self):
        events = []
        counters = []
        for attempt in (8, 4, 9, 5):
            base = event(identifier=f'{attempt}/base')
            base['attempt'] = attempt
            events.append(base)
            counters.append(dict(attempt=attempt, objective_evaluations=1,
                                 objective_sample_evaluations=3, operator_evaluations=1,
                                 failures=0, elapsed_ms=1.5))
        text = log(events) + '\n' + '\n'.join(
            'Second-stage phase audit counters: schema=2, payload=' + json.dumps(c)
            for c in counters)
        parsed_events, parsed_counters = audit.parse(text)
        self.assertEqual(parsed_events, events)
        self.assertEqual(parsed_counters, counters)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'run.log'
            output = root / 'analysis'
            source.write_text(text)
            result = subprocess.run([sys.executable, str(Path(audit.__file__)), str(source),
                                     '--output-dir', str(output)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual({p.name for p in output.iterdir()},
                             {'candidates.json', 'attempts.json', 'counters.json', 'report.md'})
            self.assertEqual(json.loads((output / 'candidates.json').read_text()), events)
            self.assertEqual(json.loads((output / 'counters.json').read_text()), counters)
            summaries = json.loads((output / 'attempts.json').read_text())
            self.assertEqual([s['attempt'] for s in summaries], [4, 5, 8, 9])
            self.assertEqual(set(summaries[0]), {'attempt', 'domain_id', 'operator_reproduced',
                             'stages', 'opposite_changes', 'gate_blocked_corrections',
                             'local_improvement_assembly_worsens'})
            report = (output / 'report.md').read_text()
            self.assertEqual(sorted([4, 5, 8, 9], key=lambda n: report.index(f'## Attempt {n}')),
                             [4, 5, 8, 9])


if __name__ == '__main__':
    unittest.main()
