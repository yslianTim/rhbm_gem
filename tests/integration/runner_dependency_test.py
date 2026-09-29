"""Prevent historical runner-to-runner dependencies from returning."""
import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
INTEGRATION = ROOT / 'tests/integration'
RUNNERS = {
    'endpoint_refinement',
    'failed_only_refinement',
    'fold_168_regression',
    'joint_bounded_validation',
    'joint_compact_validation',
    'joint_component_audit',
    'joint_component_runtime',
    'joint_fixed_validation',
    'joint_operator_validation',
    'joint_partial_selection',
    'joint_postprocessing_benchmark',
    'joint_reference_validation',
    'joint_search_validation',
    'joint_sparse_validation',
    'joint_validation',
    'joint_validation_profile',
    'joint_validation_report',
    'mdpde_experiment',
    'simulation_contract',
}


def imported_modules(tree):
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            yield from (alias.name.split('.')[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            yield node.module.split('.')[0]


def literal_strings(node):
    return [item.value for item in ast.walk(node)
            if isinstance(item, ast.Constant) and isinstance(item.value, str)]


def runner_command_reference(tree):
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = node.func.attr if isinstance(node.func, ast.Attribute) else (
            node.func.id if isinstance(node.func, ast.Name) else '')
        if name not in {'run', 'Popen', 'check_call', 'check_output', 'system', 'popen',
                        'import_module', '__import__', 'run_module', 'run_path'}:
            continue
        for value in literal_strings(node):
            for runner in RUNNERS:
                if value == runner or value.endswith('/' + runner + '.py') or value == runner + '.py':
                    yield runner


class RunnerDependencyTest(unittest.TestCase):
    def test_runner_modules_do_not_import_other_runners(self):
        for runner in sorted(RUNNERS):
            path = INTEGRATION / f'{runner}.py'
            if not path.is_file():
                continue
            tree = ast.parse(path.read_text(), filename=str(path))
            dependencies = set(imported_modules(tree)) & RUNNERS - {runner}
            self.assertFalse(dependencies, f'{path.name} imports runner(s): {sorted(dependencies)}')
            launched = set(runner_command_reference(tree)) - {runner}
            self.assertFalse(launched, f'{path.name} launches/imports runner(s): {sorted(launched)}')

    def test_support_modules_do_not_import_or_launch_runners(self):
        for path in sorted(INTEGRATION.glob('*.py')):
            if path.stem in RUNNERS or path.stem.endswith('_test') or path.name == Path(__file__).name:
                continue
            tree = ast.parse(path.read_text(), filename=str(path))
            dependencies = set(imported_modules(tree)) & RUNNERS
            self.assertFalse(dependencies, f'{path.name} imports runner(s): {sorted(dependencies)}')
            launched = set(runner_command_reference(tree))
            self.assertFalse(launched, f'{path.name} launches/imports runner(s): {sorted(launched)}')

    def test_integration_tests_do_not_import_each_other(self):
        test_modules = {path.stem for path in INTEGRATION.glob('*_test.py')}
        for path in sorted(INTEGRATION.glob('*_test.py')):
            tree = ast.parse(path.read_text(), filename=str(path))
            dependencies = set(imported_modules(tree)) & test_modules - {path.stem}
            self.assertFalse(dependencies, f'{path.name} imports test module(s): {sorted(dependencies)}')


if __name__ == '__main__':
    unittest.main()
