import copy
import unittest
from joint_fixed_validation import improvement, total_bounds
from joint_validation_checks import fixed_action_parity as parity, vector_close
from joint_validation_checks import overlap_eligible
from joint_validation_test_data import fixed_result as fixed_result


class FixedComparison(unittest.TestCase):
    def test_same_state_and_step_controls(self):
        a=fixed_result(); b=copy.deepcopy(a); b.update(mode='normal',normal_q_actions=2)
        self.assertTrue(parity(a,b)['passed'])
        for section,key,value in [('state_control','beta',[2.,.1]),('state_control','objective',.1),
                                  ('state_control','residual',[float('nan'),.2]),('work','reference_solves',1)]:
            broken=copy.deepcopy(b); broken[section][key]=value
            self.assertFalse(parity(a,broken)['passed'])
        b['steps'][0]['step']=[.3]; self.assertFalse(parity(a,b)['passed'])

    def test_missing_nonfinite_and_unavailable_evidence(self):
        a=fixed_result()
        for key in ('rank','apply','adjoint','normal','operator_gradient','steps','state_control'):
            b=copy.deepcopy(a); del b[key]
            self.assertFalse(parity(a,b)['passed'],key)
        b=copy.deepcopy(a); b['steps'][0]['true_residual']=float('nan')
        self.assertFalse(parity(a,b)['passed'])
        self.assertFalse(vector_close([],[]))
        self.assertFalse(vector_close([1],[float('inf')]))

    def test_no_median_or_improvement_from_incomplete_or_failed_runs(self):
        self.assertNotIn('baseline_median',improvement([2,3],[1,1],True))
        self.assertFalse(improvement([2,3,4],[1,1,1],False)['observed'])
        self.assertFalse(improvement([2,3,4],[3,1,1],True)['observed'])
        self.assertTrue(improvement([2,3,4],[1,1,1],True)['observed'])

    def test_old_totals_are_bounded_without_fabricating_gradient_time(self):
        r=dict(steps=[dict(total_seconds=2.)],search_work=dict(operator_adjoint_seconds=.4))
        self.assertEqual(total_bounds(r,0),(2.,2.4))
        r['total_includes_gradient']=True
        self.assertEqual(total_bounds(r,0),(2.,2.))
        self.assertEqual(total_bounds({},0),(None,None))

    def test_overlap_uses_memberships_not_names(self):
        self.assertIsNone(overlap_eligible({}))
        self.assertFalse(overlap_eligible(dict(partition=dict(blocks=1,atom_memberships=[1,1]))))
        self.assertTrue(overlap_eligible(dict(partition=dict(blocks=4,atom_memberships=[1,2,1]))))
        self.assertIsNone(overlap_eligible(dict(partition=dict(blocks=2,atom_memberships=[3]))))


class FixedCli(unittest.TestCase):
    def test_compare_exit_codes_and_receipt_immutability(self):
        import json
        import subprocess
        import sys
        import tempfile
        from pathlib import Path
        from joint_fixed_validation import CASES,MODES
        for expected,missing,bad in ((0,False,False),(3,True,False),(1,False,True),(1,True,True)):
            with self.subTest(expected=expected,missing=missing,bad=bad),tempfile.TemporaryDirectory() as folder:
                root=Path(folder); report=dict(states={c:dict(sha256='state') for c in CASES},runs={})
                for c in CASES:
                    for b in ('eigen','spqr'):
                        for mode in MODES:
                            value=fixed_result()
                            if mode=='C': value.update(mode='normal',normal_q_actions=2)
                            if bad and mode=='B': value['state_control']['objective']=10
                            report['runs'][f'{c}/{b}/{mode}']=[dict(process=dict(status='completed'),state_sha256='state',result=dict(value,stage='complete')) for _ in range(3)]
                if missing: report['runs'].pop(f'{CASES[-1]}/spqr/C')
                receipt=root/'campaign.json'; receipt.write_text(json.dumps(report)); original=receipt.read_bytes()
                process=subprocess.run([sys.executable,str(Path(__file__).with_name('joint_fixed_validation.py')),'--compare','--work-dir',folder],capture_output=True,text=True)
                self.assertEqual(process.returncode,expected,process.stderr)
                self.assertEqual(receipt.read_bytes(),original)
                self.assertEqual(json.loads((root/'comparison.json').read_text())['exit_code'],expected)

if __name__=='__main__': unittest.main()
