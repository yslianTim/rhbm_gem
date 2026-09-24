import copy
import os
import subprocess
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import joint_bounded_validation as runner
from joint_fixed_validation_test import result as fixed_result
from joint_search_validation_test import result as search_result


def row(value): return dict(process=dict(status='completed'),result=value)


class BoundedValidation(unittest.TestCase):
    def test_schedule_and_status_precedence(self):
        cases=list(runner.search_schedule()); self.assertEqual(len(cases),28)
        for backend in runner.BACKENDS:
            for kind in ('legacy','schwarz'):
                self.assertEqual([r for b,c,k,r in cases if b==backend and c=='single-512' and k==kind],[1,2,3])
        self.assertEqual(runner.status([True])['exit_code'],0)
        self.assertEqual(runner.status([None])['exit_code'],3)
        self.assertEqual(runner.status([None,False])['exit_code'],1)

    def test_no_performance_qualification_for_unconverged_or_missing_evidence(self):
        value=search_result(); self.assertTrue(runner.qualified(row(value)))
        for assessment in (None,{},dict(runtime_convergence='failed')):
            value['returned_assessment']=assessment
            self.assertFalse(runner.qualified(row(value)))
        value=search_result(); value['returned_state']['b']=[float('nan')]
        self.assertFalse(runner.qualified(row(value)))
        self.assertFalse(runner.qualified({}))

    def test_wall_timer_scope_and_nonfinite_inputs(self):
        value=fixed_result(); value['stage']='complete'; value['total_includes_gradient']=True
        for step in value['steps']:
            step.update(fixed_step_wall_seconds=1.,total_seconds=.6,preparation_seconds=.1,gradient_seconds=.1,
                        metric_seconds=.1,partition_seconds=.1,build_seconds=.1,solve_seconds=.1)
        kinds=('identity','diagonal','schwarz'); self.assertTrue(runner.fixed_check(row(value),kinds))
        for key,bad in (('fixed_step_wall_seconds',.1),('gradient_seconds',float('nan')),('fixed_step_wall_seconds',None)):
            broken=copy.deepcopy(value); broken['steps'][0][key]=bad
            self.assertFalse(runner.fixed_check(row(broken),kinds))
        self.assertIsNone(runner.fixed_check({},kinds))

    def test_rank_evidence_and_oracle_identity(self):
        result=dict(stage='complete',input_sha256='same',free_columns=2,rank_compact_extractions=0,
                    work=dict(free_design_svds=0,reference_solves=0,derivative_preparations=0),
                    rank_result=dict(status='full-rank',rank_lower=2,rank_upper=2,minimum_lower=.4,maximum_lower=1.,maximum_upper=2.,threshold_lower=1e-12,threshold_upper=2e-12,reconstruction_error=1e-14,orthogonal_minimum=.999))
        oracle=row(dict(stage='complete',input_sha256='same',free_columns=2,rank_result=dict(valid=True,rank=2,singular_values=[1.,.5],threshold=1e-12)))
        self.assertTrue(runner.rank_check(row(result),[oracle])['passed'])
        oracle['result']['input_sha256']='different'
        self.assertFalse(runner.rank_check(row(result),[oracle])['passed'])
        result['rank_result']['status']='unavailable'
        self.assertIsNone(runner.rank_check(row(result),[oracle])['passed'])
        result['rank_compact_extractions']=1
        self.assertFalse(runner.rank_check(row(result),[])['passed'])

    def test_campaign_cancel_preserves_partial_receipt_and_refuses_resume(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)/'campaign'; builds={b:Path(folder)/b for b in runner.BACKENDS}
            with patch.object(runner,'validate_builds',return_value={b:dict(source_sha256='same') for b in runner.BACKENDS}),patch.object(runner.v,'process_tree_rss',return_value=0),patch.object(runner.v,'monitored',side_effect=KeyboardInterrupt):
                with self.assertRaises(KeyboardInterrupt): runner.campaign(root,Path(folder)/'inputs',builds)
            receipt=json.loads((root/'campaign.json').read_text())
            self.assertFalse(receipt['finished']); self.assertEqual(receipt['interruption'],'user-cancelled')
            self.assertEqual(receipt['runs']['freeze/chain-8']['process']['status'],'user-cancelled')
            self.assertEqual(receipt['runs']['search/eigen/single-512/legacy/1']['process']['status'],'not-run-user-cancelled')
            with self.assertRaises(ValueError): runner.campaign(root,Path(folder),builds)

    def test_resource_stop_ends_repetitions_and_larger_rank_cases(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)/'campaign'; inputs=Path(folder)/'inputs'
            for case in ('single-128','single-512','heterogeneous-168'):
                path=inputs/case; path.mkdir(parents=True)
                for name in ('input.cif','input.map','widths.json'): (path/name).write_text('{}')
            calls=[]
            def monitor(command,directory,deadline,**limits):
                calls.append((command,deadline,limits))
                if '--fixed' in command and command[command.index('--fixed')+1]=='freeze':
                    Path(command[command.index('--fixed')-1]).write_text('{}')
                    return dict(status='completed')
                return dict(status='time-limit')
            builds={b:Path(folder)/b for b in runner.BACKENDS}
            with patch.object(runner,'validate_builds',return_value={b:dict(source_sha256='same') for b in runner.BACKENDS}),patch.object(runner.v,'process_tree_rss',return_value=0),patch.object(runner.v,'monitored',side_effect=monitor),patch.object(runner,'unpack',return_value=inputs/'heterogeneous-168'):
                self.assertEqual(runner.campaign(root,inputs,builds),3)
            report=json.loads((root/'campaign.json').read_text())
            self.assertEqual(report['runs']['search/eigen/single-512/legacy/2']['process']['status'],'not-run-group-stopped')
            self.assertEqual(report['runs']['rank/spqr/chain-2000/prototype']['process']['status'],'not-run-topology-stopped')
            self.assertEqual(sum(runner.STAGE_SECONDS.values()),3600)
            self.assertTrue(all(c[2]['seconds']==(600 if '--search' in c[0] else 180) for c in calls))

    def test_mismatched_source_fingerprints_rejected_before_running(self):
        builds={b:Path(b) for b in runner.BACKENDS}
        with patch.object(runner,'require_current_build'),patch.object(runner,'fingerprint',side_effect=[dict(source_sha256='one'),dict(source_sha256='two')]):
            with self.assertRaisesRegex(ValueError,'source fingerprints'): runner.validate_builds(builds)


@unittest.skipUnless(os.environ.get('JOINT_BOUNDED_TEST_EXECUTABLE'),'driver supplied by CTest')
class DiagnosticDriver(unittest.TestCase):
    def test_real_fixed_timer_and_rank_modes(self):
        executable=os.environ['JOINT_BOUNDED_TEST_EXECUTABLE']
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); state=root/'state.json'
            def run(phase,output,*options):
                process=subprocess.run([executable,'synthetic','chain','2',phase,str(output),*map(str,options)],env=runner.v.ENV,capture_output=True,text=True,timeout=30)
                self.assertEqual(process.returncode,0,process.stderr)
                return json.loads(output.read_text())
            run('fixed',state,'--fixed','freeze')
            measured=run('fixed',root/'fixed.json','--fixed','normal','--state',state,'--resources')
            self.assertTrue(runner.fixed_check(row(measured),('identity','diagonal','schwarz')))
            selected=run('fixed',root/'selected.json','--fixed','normal','--state',state,'--fixed-preconditioner','schwarz','--resources')
            self.assertTrue(runner.fixed_check(row(selected),('schwarz',)))
            prototype=run('rank',root/'prototype.json','--resources')
            oracle=run('rank-oracle',root/'oracle.json','--resources')
            if os.environ.get('JOINT_BOUNDED_TEST_BACKEND')=='SPQR':
                self.assertEqual(prototype['rank_result']['status'],'full-rank')
                self.assertTrue(runner.rank_check(row(prototype),[row(oracle)])['passed'])
            else: self.assertEqual(prototype['rank_result']['reason'],'rank-backend-unavailable')


if __name__=='__main__': unittest.main()
