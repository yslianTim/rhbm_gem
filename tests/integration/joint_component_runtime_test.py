import json
from pathlib import Path
import tempfile
import unittest
from joint_component_runtime import frozen_parity


class FrozenRuntimeParityTest(unittest.TestCase):
    def test_missing_reference_is_not_a_passing_empty_comparison(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(RuntimeError):
                frozen_parity(root/'old', root/'new', root/'audits')

    def test_preserves_old_scopes_and_detects_changed_state(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for folder in ('old/datasets', 'new/datasets', 'audits', 'new/audits'):
                target = root/folder/'baseline'/'fit.json'
                target.parent.mkdir(parents=True)
                target.write_text(json.dumps(dict(beta=[1, 2], search_seconds=1)))
            extra = root/'new/audits/baseline/local-components/fit.json'
            extra.parent.mkdir(); extra.write_text('{}')
            args = root/'old', root/'new', root/'audits'
            self.assertTrue(frozen_parity(*args)['passed'])
            (root/'new/datasets/baseline/fit.json').write_text(json.dumps(dict(beta=[1, 3], search_seconds=2)))
            self.assertEqual(frozen_parity(*args)['differences'], ['search/baseline/fit.json'])


if __name__ == '__main__': unittest.main()
