import copy
import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np
import atom_block_grid_composite as composite
import atom_centered_voxel_union as union
import unique_stencil_grid_test as previous


class AtomBlockGridCompositeRunnerTest(unittest.TestCase):
    def fixture(self):
        fit, context = previous.UniqueStencilGridRunnerTest().fixture()
        context['atoms'] = [{'alpha': .1}]
        fit.update(experiment=composite.NAME, alpha_source='checkpoint', svd_preconditioner='independent-tsqr-8192',
                   membership_count=30, blocks=[{'owner': 0, 'alpha': .1, 'rows': 30, 'lambda': 1}],
                   selected_seed=0, variances=[.01])
        phase = {'iterations': 0, 'stop': 'budget-exhausted', 'stationarity': .01,
                 'beta': fit['beta'], 'variances': [.01]}
        branch = {'primary': phase, 'reference': phase, 'trace': [phase], 'reference_trace': [],
                  'qualified': False, 'initial_variances': [.01]}
        fit['branches'] = [dict(copy.deepcopy(branch), seed=seed) for seed in ('checkpoint', 'constrained-ls')]
        return fit, context, [(np.arange(30), np.zeros(30))]

    def test_scope_and_missing_cases(self):
        index = {'selected_state_ids': ['baseline-best-28'], 'alpha_source': 'checkpoint'}
        self.assertEqual(composite.expected_cases(index), ['baseline-best-28-checkpoint-alpha'])
        for changed in (dict(index, selected_state_ids=[]), dict(index, selected_state_ids=['a','b']), dict(index, alpha_source='common')):
            with self.assertRaises(RuntimeError): composite.expected_cases(changed)
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)
            for name, value in [('state-index',index),('completion',{'complete':False,'cases':[]}),('fit-index',{'cases':[]})]:
                composite.experiment.write(p/(name+'.json'),value)
            with self.assertRaisesRegex(RuntimeError,'Missing requested cases'): composite.summarize(p)

    def test_alpha_b_branches_and_budgets_are_required(self):
        fit, context, blocks = self.fixture(); composite.validate_fit(fit,context,blocks,30)
        for mutation, pattern in ((lambda f:f['blocks'][0].update(alpha=.2),'alpha'),
                                  (lambda f:f['abc'][0].__setitem__(1,.5),'Width'),
                                  (lambda f:f['branches'].pop(),'branch'),
                                  (lambda f:f['branches'][0]['primary'].update(iterations=101),'budget'),
                                  (lambda f:f.update(qualified=True),'Unverified')):
            changed=copy.deepcopy(fit); mutation(changed)
            with self.assertRaisesRegex(RuntimeError,pattern): composite.validate_fit(changed,context,blocks,30)

    def test_qualified_requires_independent_certificate_and_complete_trace(self):
        fit, context, blocks=self.fixture()
        fit['branches'][0]['qualified']=True
        with self.assertRaisesRegex(RuntimeError,'qualified evidence'): composite.validate_fit(fit,context,blocks,30)
        fit['branches'][0]['qualified']=False; fit['branches'][0]['trace']=[]
        with self.assertRaisesRegex(RuntimeError,'trajectory'): composite.validate_fit(fit,context,blocks,30)

    def test_complete_memberships_and_duplicate_rejection(self):
        data={'grid_size':[5]*3,'generation_origin':[-1.]*3,'generation_spacing':[.5,.75,1.],
              'sampling_origin':[-1.1]*3,'sampling_spacing':[.6]*3,'radius':1.5,
              'atoms':[{'position':[0,0,0]},{'position':[.5,0,0]}]}
        indices,coverage,_=union.sphere_membership(data)
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp); rows=[]
            for b,(ids,_) in enumerate(composite.geometry_blocks(data)):
                rows += [(b,int(np.searchsorted(indices,i)),int(i)) for i in ids]
            def write(entries):
                with (p/'memberships.csv').open('w') as stream:
                    writer=csv.writer(stream); writer.writerow(['block','row','index']); writer.writerows(entries)
                data.update(membership_count=len(rows),membership_table_sha256=composite.fold.sha256_file(p/'memberships.csv'))
            write(rows); blocks=composite.validate_memberships(p,data)
            self.assertEqual(sum(len(r) for r,_ in blocks),int(coverage.sum()))
            write([rows[0]]+rows)
            with self.assertRaisesRegex(RuntimeError,'membership'): composite.validate_memberships(p,data)
            write(rows[:-1])
            with self.assertRaisesRegex(RuntimeError,'membership'): composite.validate_memberships(p,data)

    def test_overlap_mass_equations_and_zero_alpha_scales(self):
        r=np.array([.2,-.1,.3,.4]); rows=[np.array([0,1,2]),np.array([1,2,3])]
        blocks=[(x,np.zeros(len(x))) for x in rows]
        columns=[(x,np.ones(len(x)),np.arange(1,len(x)+1,dtype=float)) for x in rows]
        e=composite.composite_evidence(r,[.03,.07],[0,.5],blocks,columns,[1,0,1,0])
        self.assertAlmostEqual(sum(b['linear_weight_share'] for b in e['blocks']),1.)
        q=[b['scaled_prefactor'] for b in e['blocks']]
        expected=np.zeros(4)
        for i,rr in enumerate(rows): expected[rr]+=q[i]*np.exp(-.5*[0,.5][i]*r[rr]**2/[.03,.07][i])
        np.testing.assert_allclose(e['omega'],expected,rtol=1e-14)
        zero=composite.composite_evidence(r,[.01,1],[0,0],blocks,columns,[1,0,1,0])
        self.assertGreater(zero['omega'].max()/zero['omega'].min(),50)

    def test_aggregate_weight_roundoff_is_propagated_without_accepting_changed_weights(self):
        residual=np.array([1e-6,-2e-6,3e-6]); variance=2e-12
        actual=residual+3*np.finfo(float).eps
        table=np.column_stack((np.arange(3),np.ones(3),residual,np.exp(-.5*actual**2/variance)))
        blocks=[(np.arange(3),np.zeros(3))]
        self.assertFalse(np.allclose(np.exp(-.5*residual**2/variance),table[:,3],rtol=1e-10,atol=1e-14))
        composite.validate_aggregate_weights(table,[variance],[1.],blocks)
        table[0,3]+=1e-5
        with self.assertRaisesRegex(RuntimeError,'Aggregate weight replay'):
            composite.validate_aggregate_weights(table,[variance],[1.],blocks)

    def test_compare_does_not_hide_scientific_fields(self):
        a={'seconds':1,'beta':[1],'variances':[.1,.2],'membership_count':8,'iterations':100}
        b=dict(a,seconds=2)
        self.assertEqual(union.scientific(a),union.scientific(b))
        self.assertNotEqual(union.scientific(a),union.scientific(dict(b,membership_count=9)))


if __name__=='__main__': unittest.main()
