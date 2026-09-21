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
