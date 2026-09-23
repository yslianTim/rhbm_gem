"""Contract tests for the offline experiment runner; no estimator success oracle is changed."""
import sys
import tempfile
import time
import unittest
from pathlib import Path
import numpy as np
import joint_validation as v
from joint_validation_report import parameter_stats, statistical_summary, wilson

class ValidationTest(unittest.TestCase):
    def test_reference_source_mismatch_fails_before_launch(self):
        from joint_reference_validation import build_info
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'generated/Release').mkdir(parents=True)
            (root/'generated/Release/SimulationBuildInfo.hpp').write_text('wrong source')
            with self.assertRaisesRegex(ValueError,'source fingerprint mismatch'):
                build_info(root,'expected')

    def test_reference_missing_export_is_not_a_completed_result(self):
        from joint_reference_validation import exports
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,'Missing completed-command export'):
                exports(Path(tmp))

    def test_reference_censored_or_missing_samples_have_no_speed_measurement(self):
        from joint_reference_validation import command_group
        for rows in ([],[dict(status='time-limit')]*3,[dict(status='not-run-budget')]*2):
            result,states=command_group(rows,Path('.'),True)
            self.assertFalse(result['passed']); self.assertFalse(result['complete'])
            self.assertIsNone(result['median_seconds']); self.assertEqual(states,[])

    def test_reference_objective_is_not_normalized_twice(self):
        from joint_reference_validation import endpoint_comparison
        state=dict(ac=[1.,0.],b=[.5],objective=0.)
        a=dict(assembled_state=state,atom_ids=['a'],components=[],observation_scale=1e6)
        self.assertTrue(endpoint_comparison(a,a)['passed'])
        b={**a,'assembled_state':{**state,'objective':1e-8}}
        self.assertFalse(endpoint_comparison(a,b)['passed'])
        self.assertEqual(endpoint_comparison(a,b)['normalized_objective_difference'],1e-8)

    def test_reference_requires_three_distinct_samples_and_zero_search_references(self):
        from unittest.mock import patch
        import joint_reference_validation as r
        outcome=dict(runtime_convergence='passed',costs=dict(search_reference_seconds=0),components=[dict(reference_evaluations=0)])
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rows=[]
            for k in (1,2,3):
                directory=root/str(k);directory.mkdir()
                for name in r.EXPORTS:(directory/name).write_text('fixture')
                rows.append(dict(status='completed',repetition=k,directory=str(k),case='single-128',total_command_seconds=1.,
                                 export_hashes={n:v.sha(directory/n) for n in r.EXPORTS}))
            with patch.object(r,'exports',return_value=outcome),patch.object(r.v,'read',return_value=dict(b=[.5])):
                self.assertTrue(r.command_group(rows,root,True)[0]['passed'])
                self.assertFalse(r.command_group(rows[:2],root,True)[0]['passed'])
                self.assertFalse(r.command_group([rows[0]]*3,root,True)[0]['passed'])
                outcome['components'][0]['reference_evaluations']=1
                self.assertFalse(r.command_group(rows,root,True)[0]['passed'])

    def test_reference_reaggregation_is_repeatable_and_rejects_corruption(self):
        from unittest.mock import patch
        from joint_reference_validation import summarize
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); path=root/'evidence.txt';path.write_text('original')
            r=dict(finished=True,historical=dict(compact_gate_passed=True),commands={},audits={},files={'evidence.txt':v.sha(path)})
            with patch('joint_reference_validation.historical',return_value=r['historical']):
                self.assertEqual(summarize(r,root),summarize(r,root))
                self.assertFalse(summarize(r,root)['latest_512_passed'])
            path.write_text('changed')
            with self.assertRaisesRegex(ValueError,'Evidence hash mismatch'):
                summarize(r,root)

    def test_compact_svd_checks_weak_spectrum_and_threshold_not_only_residual(self):
        import json
        from joint_compact_validation import svd_parity
        row=dict(valid=True,rank=2,threshold=1e-9,singular_values=[1.,2e-9],solution=[1.,2.])
        self.assertTrue(svd_parity(row,row)['passed'])
        self.assertFalse(svd_parity(row,{**row,'rank':1})['passed'])
        self.assertFalse(svd_parity(row,{**row,'threshold':1e-8})['passed'])
        self.assertFalse(svd_parity(row,{**row,'singular_values':[1.,4e-9]})['passed'])
        self.assertFalse(svd_parity(row,{**row,'solution':[1.,None]})['passed'])
        self.assertEqual(len(svd_parity(row,row)['weak_values']),2)
        json.dumps(svd_parity(row,row),allow_nan=False)

    def test_compact_audit_requires_trust_and_full_derivative(self):
        from joint_compact_validation import audit_parity
        row=dict(valid=True,free_rank=2,kkt_passed=True,active_atoms=[],beta=[1.,2.],objective=.1,relative_residual=.1,b_gradient=[0.])
        spectrum=dict(valid=True,rank=2,threshold=1e-9,singular_values=[1.,.5],solution=[])
        derivative=dict(valid=True,**{k:[1.,2.] for k in ('coefficients','correction','projected','jacobian','response')})
        a=dict(primary=row,reference=row,trust=dict(passed=True),initial_b=[.5],
               derivative_audit=derivative,derivative_reason='full-profile-derivative',svd_records=[spectrum])
        self.assertTrue(audit_parity(a,a)['passed'])
        self.assertFalse(audit_parity(a,{**a,'trust':dict(passed=False)})['passed'])
        self.assertFalse(audit_parity(a,{**a,'derivative_audit':{**derivative,'correction':[]}})['passed'])

    def test_compact_summary_cannot_pass_missing_or_censored_evidence(self):
        from joint_compact_validation import summary
        result=summary(dict(fixed={},audits={},baseline_controls={},commands={}))
        self.assertFalse(result['compact_gate_passed'])
        self.assertFalse(result['complete_512_passed'])
        self.assertIsNone(result['commands']['single-512']['candidate']['median_seconds'])

    def test_sparse_parity_requires_matching_rank_and_available_coefficients(self):
        from joint_sparse_validation import parity
        row=dict(valid=True,free_rank=2,kkt_passed=True,active_atoms=[],beta=[1.,2.],objective=.1,relative_residual=.1,b_gradient=[0.])
        a=dict(primary=row,reference=row)
        self.assertTrue(parity(a,a)['primary']['passed'])
        self.assertFalse(parity(a,dict(primary={**row,'free_rank':1}))['primary']['passed'])
        self.assertFalse(parity(a,dict(primary={**row,'beta':[None,None]}))['primary']['passed'])
        self.assertFalse(parity(a,{})['primary']['available'])

    def test_empty_statistics_and_zero_failure_uncertainty(self):
        self.assertIsNone(parameter_stats([])['bias'])
        self.assertGreater(wilson(0,20)[1],.1)
        self.assertIsNone(wilson(0,0))

    def test_statistical_denominators_exclude_missing_processes_and_states(self):
        def result(state,status,error):
            return dict(state_available=state,runtime_convergence=status,errors=[[error]*3]*2 if state else None,
                        failure_reasons=[],prediction_rmse=None,residual_rmse=None,residual_neighbor_correlation=None)
        records=[]
        for outcome in (result(True,'passed',1),result(True,'failed',3),result(False,'unavailable',0),None):
            records.append(dict(distance=1.2,shift=0,noise='iid',sigma_fraction=.05,modes=['fixed'],
                                results={'fixed':outcome} if outcome is not None else {}))
        row=statistical_summary(records)[0]
        self.assertEqual((row['attempts'],row['recorded'],row['process_incomplete']),(4,3,1))
        self.assertEqual(row['unavailable_rate'],1/3)
        self.assertEqual(row['nonconverged_available_rate'],.5)
        self.assertEqual(row['all_available']['bias'][0][0],2)
        self.assertEqual(row['converged_only']['bias'][0][0],1)
        self.assertEqual(row['nonconverged_only']['bias'][0][0],3)

    def test_matrix_and_paired_seeds(self):
        cases=list(v.conditions())
        self.assertEqual(len(cases),326)
        self.assertEqual(sum(len(c[-1]) for c in cases),450)
        self.assertEqual(len(set(v.SEEDS)),20)
        for c in cases:
            if c[4] is not None:self.assertIn(c[4],v.SEEDS)
            if 'production' in c[-1]:self.assertTrue(c[3] in (0.,.05))

    def test_noise_is_reproducible_and_variance_matched(self):
        for kind in ('iid','correlated'):
            np.testing.assert_array_equal(v.noise_field(v.SEEDS[0],kind),v.noise_field(v.SEEDS[0],kind))
            fields=np.array([v.noise_field(s,kind) for s in v.SEEDS])
            self.assertAlmostEqual(float(np.mean(fields**2)),1.,delta=.06)
            self.assertEqual(fields.shape,(20,25**3))
            self.assertGreater(float(np.std(fields.mean(axis=1))),0) # no per-field demeaning

    def test_shift_changes_data_not_selection(self):
        a=v.clean_map(1.2,0);b=v.clean_map(1.2,.15)
        self.assertGreater(float(np.linalg.norm(a-b)),0)
        self.assertEqual(int(v.TARGET.sum()),515)
        self.assertTrue(np.isfinite(b).all())

    def test_missing_state_stays_missing(self):
        outcome=dict(assembled_state=None,runtime_convergence='unavailable',initialization={'valid':False},components=[])
        m=v.metrics(outcome,None,np.zeros(1),np.zeros(1),np.array([0]))
        self.assertIsNone(m['errors']);self.assertIsNone(m['parameters'])
        self.assertFalse(m['state_available'])

    def test_timeout_and_memory_are_process_failures(self):
        with tempfile.TemporaryDirectory() as tmp:
            command=[sys.executable,'-c','import time; time.sleep(10)']
            timeout=v.monitored(command,Path(tmp)/'timeout',time.monotonic()+5,seconds=.2)
            self.assertEqual(timeout['status'],'time-limit')
            memory=v.monitored(command,Path(tmp)/'memory',time.monotonic()+5,rss_limit=1)
            self.assertEqual(memory['status'],'rss-limit')
            self.assertGreater(memory['sampled_tree_peak_rss_bytes'],1)

    def test_expired_budget_does_not_launch(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=v.monitored(['does-not-exist'],tmp,time.monotonic()-1)
            self.assertEqual(result['status'],'not-run-budget')

if __name__=='__main__':unittest.main()
