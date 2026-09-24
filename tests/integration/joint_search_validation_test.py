import copy
import unittest
import tempfile
import os
from pathlib import Path
from unittest.mock import patch
import joint_search_validation as runner
from joint_search_validation import scientific_parity


def result():
    return dict(stage='complete',search_completed=True,search=dict(residual_scale=2),
        returned_assessment=dict(runtime_convergence='passed',runtime_failure='none',runtime_checks=dict(inner=True,local_correction=True),design_spectrum=dict(rank=2)),
        returned_state=dict(objective=.04,beta=[2.,.2],b=[.5],active_atoms=[]))


class SearchComparison(unittest.TestCase):
    def test_stale_binaries_are_rejected_before_a_campaign(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); source=root/'source'; build=root/'build'
            for path in ('src/a.cpp','include/a.hpp','cmake/a.cmake','CMakeLists.txt','tests/experiments/joint_sparse_benchmark.cpp'):
                target=source/path; target.parent.mkdir(parents=True,exist_ok=True); target.write_text('source')
                os.utime(target,ns=(100,100))
            for path in ('src/librhbm_gem.dylib','bin/joint_sparse_benchmark'):
                target=build/path; target.parent.mkdir(parents=True,exist_ok=True); target.write_text('binary')
                os.utime(target,ns=(200,200))
            (build/'CMakeCache.txt').write_text(f'CMAKE_HOME_DIRECTORY:INTERNAL={source}\n')
            runner.require_current_build(build)
            os.utime(source/'src/a.cpp',ns=(300,300))
            with self.assertRaisesRegex(ValueError,'newer than measured binaries'): runner.require_current_build(build)
            os.utime(source/'src/a.cpp',ns=(100,100))
            os.utime(source/'tests/experiments/joint_sparse_benchmark.cpp',ns=(300,300))
            os.utime(build/'bin/joint_sparse_benchmark',ns=(400,400))
            runner.require_current_build(build)

    def test_trajectory_does_not_replace_endpoint_check(self):
        a=result(); b=copy.deepcopy(a); b['search']['accepted_updates']=99
        self.assertTrue(scientific_parity(a,b)['passed'])
        b['returned_state']['b'][0]+=.001
        self.assertFalse(scientific_parity(a,b)['passed'])

    def test_evidence_and_availability(self):
        a=result(); b=copy.deepcopy(a); b['returned_assessment']['runtime_checks']['local_correction']=False
        self.assertFalse(scientific_parity(a,b)['passed'])
        b['returned_assessment']=None; self.assertFalse(scientific_parity(a,b)['passed'])

    def test_normalized_objective_and_search_failure(self):
        a=result(); b=copy.deepcopy(a); b['returned_state']['objective']+=1e-8
        self.assertFalse(scientific_parity(a,b)['passed'])
        b=copy.deepcopy(a); b['search_completed']=False
        self.assertFalse(scientific_parity(a,b)['passed'])

    def test_explicit_unavailability_is_not_a_qualified_endpoint(self):
        a=dict(returned_assessment=None,search=dict(initial=dict(valid=False,reason='rank-deficient')))
        self.assertTrue(scientific_parity(a,a)['passed'])
        b=copy.deepcopy(a); b['search']['initial']['reason']='unknown'
        self.assertFalse(scientific_parity(a,b)['passed'])


    def test_incomplete_campaign_never_promotes(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            for stage in ('baseline','eigen','spqr','large-local'):
                runner.v.write(root/(stage+'.json'),dict(builds={'eigen':dict(source_sha256='frozen')},inputs={},runs={},finished=False))
            with patch.object(runner,'frozen_hash',return_value='frozen'):
                compared=runner.compare(root)
            self.assertFalse(compared['benchmark_promotion_eligible'])
            self.assertFalse(compared['comparisons_complete'])
            self.assertFalse(compared['promote'])


if __name__=='__main__': unittest.main()
