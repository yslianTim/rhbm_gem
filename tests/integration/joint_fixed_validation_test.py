import copy
import unittest
from joint_fixed_validation import parity, improvement, vector_close, total_bounds
from joint_search_validation import overlap_eligible


def result():
    return dict(valid=True,reason='full-profile-operator',mode='composed',normal_q_actions=6,
        state_control=dict(eta=[.2],beta=[1.,.1],free_columns=[0,1],scale=2.,rank_rows=10,
                           residual=[.1,.2],gradient=[.03],objective=.025),
        rank=[dict(valid=True,rank=2,rows=2,columns=2,relative_threshold=1e-12,absolute_override=-1,
                   threshold=1e-12,singular_values=[1.,.5])],
        apply=[.1,.2],adjoint=[.3],normal=[.4],operator_gradient=[.03],
        work=dict(reference_solves=0,derivative_preparations=0),
        steps=[dict(kind=k,valid=True,step=[.2],predicted=.001,true_residual=1e-11) for k in ('identity','diagonal','schwarz')])


class FixedComparison(unittest.TestCase):
    def test_same_state_and_step_controls(self):
        a=result(); b=copy.deepcopy(a); b.update(mode='normal',normal_q_actions=2)
        self.assertTrue(parity(a,b)['passed'])
        for section,key,value in [('state_control','beta',[2.,.1]),('state_control','objective',.1),
                                  ('state_control','residual',[float('nan'),.2]),('work','reference_solves',1)]:
            broken=copy.deepcopy(b); broken[section][key]=value
            self.assertFalse(parity(a,broken)['passed'])
        b['steps'][0]['step']=[.3]; self.assertFalse(parity(a,b)['passed'])

    def test_missing_nonfinite_and_unavailable_evidence(self):
        a=result()
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


if __name__=='__main__': unittest.main()
