"""Receipt comparison must distinguish parity, unavailable and incomplete work."""
import copy
import tempfile
import unittest
from pathlib import Path
import joint_operator_validation as runner


def state():
    endpoint=dict(valid=True, reason='qualified-inner', feasible=True, free_rank=2,
                  kkt_passed=True, active_atoms=[], beta=[2., .2], objective=.1,
                  relative_residual=.2, b_gradient=[1e-13])
    return dict(primary=endpoint, reference=copy.deepcopy(endpoint))


class Receipts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=runner.baseline_source_hash()

    def compare(self, left, right, status='completed'):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name, data in [('baseline',left),('candidate',right)]:
                rows=[dict(process=dict(status=status), result=data) for _ in range(3)]
                runner.v.write(root/(name+'.json'),dict(fixed={'eigen/test':rows}, controls={},
                    builds={'eigen':dict(source_sha256=self.source)}))
            return runner.compare(root)['eigen/test']

    def test_complete_and_changed_objective(self):
        baseline=state(); actual=copy.deepcopy(baseline)
        self.assertTrue(self.compare(baseline,actual)['passed'])
        actual['primary']['relative_residual']=.3
        self.assertFalse(self.compare(baseline,actual)['passed'])

    def test_limits_never_pass(self):
        self.assertIsNone(self.compare(state(),state(),'time-limit')['passed'])

    def test_unavailable_is_explicit(self):
        data=state()
        for endpoint in data.values(): endpoint.update(valid=False,reason='rank-deficient')
        result=self.compare(data,data)
        self.assertTrue(result['passed'])
        self.assertTrue(result['checks'][0]['roles']['primary']['expected_unavailable'])
        self.assertFalse(result['checks'][0]['roles']['primary']['available'])

    def test_failed_operator_cannot_hide_behind_primary_parity(self):
        data=state(); actual=copy.deepcopy(data); actual['operator']=dict(valid=True,passed=False)
        self.assertFalse(self.compare(data,actual)['passed'])


if __name__=='__main__': unittest.main()
