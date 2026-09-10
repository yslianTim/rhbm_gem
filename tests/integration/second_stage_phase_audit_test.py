import importlib.util
import json
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('phase_audit', Path(__file__).resolve().parents[2] / 'tools/second_stage_phase_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def event(identifier='base', parent='', total=1.0, delta=None, stage='baseline'):
    return dict(attempt=6, domain_id=2, candidate_id=identifier, parent_id=parent,
                stage=stage, key=[], factor=1, disposition='observed', reason='', final_retained=False,
                objective=dict(fit=total, tail_weighted=0, offset=0, total=total),
                delta_parent=delta, delta_baseline=delta, operator_reproduced=True,
                operator=dict(status='available', p99=[.01, .02, .03]), direction_samples=[])


def log(events):
    return '\n'.join('[Debug] Second-stage phase audit: schema=1, payload=' + json.dumps(e) for e in events)


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
        with self.assertRaises(ValueError):
            audit.parse(log([event()]).replace('schema=1', 'schema=2'))
        with self.assertRaises(ValueError):
            audit.parse(log([event(), event('child', 'base', .9, -.2)]))

    def test_direction_sampling_is_not_proof(self):
        child = event('child', 'base', 1.2, .2, 'production-proposal')
        child['direction_samples'] = [dict(alpha=1, delta_parent=.2), dict(alpha=.5, delta_parent=-.1)]
        summary = audit.summarize([event(), child])[0]
        self.assertEqual(summary['oversized_steps'], ['child'])
        child['direction_samples'][1]['delta_parent'] = .01
        summary = audit.summarize([event(), child])[0]
        self.assertEqual(summary['no_sampled_descent'], ['child'])
        self.assertIn('not a mathematical proof', audit.report([event(), child], [], [summary]))

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


if __name__ == '__main__':
    unittest.main()
