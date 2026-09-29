"""Check historical artifact provenance and active-tree ownership."""
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = ROOT / 'docs/developer/figures/experiment-retirement-baseline/artifact-manifest.json'
ALLOWED_STATUSES = {
    'present_active',
    'present_historical',
    'removed_from_worktree',
    'missing',
    'external_archive',
}
LARGE_ARTIFACT_BYTES = 5 * 1024 * 1024


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def manifest_paths(value):
    paths = set()
    if isinstance(value, dict):
        for key, item in value.items():
            if key == 'path' and isinstance(item, str):
                paths.add(item)
            paths.update(manifest_paths(item))
    elif isinstance(value, list):
        for item in value:
            paths.update(manifest_paths(item))
    return paths


def active_text_files():
    roots = [ROOT / name for name in ('src', 'include', 'tests', 'resources', 'cmake')]
    suffixes = {'.c', '.cc', '.cpp', '.h', '.hpp', '.py', '.cmake', '.md', '.txt', '.json'}
    files = {ROOT / 'CMakeLists.txt', ROOT / 'README.md'}
    for root in roots:
        if root.exists():
            files.update(path for path in root.rglob('*') if path.is_file() and path.suffix in suffixes)
    docs = ROOT / 'docs/developer'
    for path in docs.rglob('*.md'):
        relative_parts = path.relative_to(docs).parts
        if 'figures' not in relative_parts and path.name != 'experiment-retirement-baseline.md':
            files.add(path)
    return sorted(path for path in files if path.is_file())


def code_text_files():
    roots = [ROOT / name for name in ('src', 'include', 'tests', 'resources', 'cmake')]
    suffixes = {'.c', '.cc', '.cpp', '.h', '.hpp', '.py', '.cmake', '.md', '.txt', '.json'}
    files = {ROOT / 'CMakeLists.txt'}
    for root in roots:
        if root.exists():
            files.update(path for path in root.rglob('*') if path.is_file() and path.suffix in suffixes)
    return sorted(path for path in files if path.is_file())


class HistoricalArtifactManifestTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding='utf-8'))
        if cls.manifest['schema_version'] != 1:
            raise AssertionError('unexpected artifact manifest schema version')
        cls.update = cls.manifest['retirement_update']['pr6_update']
        cls.artifacts = cls.update['artifacts']
        if cls.update['entry_count'] != len(cls.artifacts):
            raise AssertionError('artifact entry count does not match the manifest')

    def test_removed_artifacts_have_unique_provenance_and_are_absent(self):
        ids = [record['id'] for record in self.artifacts]
        paths = [record['path'] for record in self.artifacts]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(len(paths), len(set(paths)))
        self.assertEqual(set(self.update['status_values']), {record['status'] for record in self.artifacts})
        self.assertTrue(set(self.update['status_values']).issubset(ALLOWED_STATUSES))

        for record in self.artifacts:
            with self.subTest(artifact=record['id']):
                self.assertIn(record['status'], ALLOWED_STATUSES)
                self.assertIn(record['classification'], {
                    'ACTIVE_FIXTURE',
                    'ACTIVE_REFERENCE',
                    'ACTIVE_BENCHMARK_INPUT',
                    'ACTIVE_RESEARCH_OUTPUT',
                    'HISTORICAL_EVIDENCE',
                    'DUPLICATE_HISTORICAL_EVIDENCE',
                    'OBSOLETE_RAW_LOG',
                    'MISSING_HISTORICAL_ARTIFACT',
                })
                self.assertFalse(record['active_dependency'])
                self.assertEqual(record['active_callers'], [])
                self.assertTrue(record['historical_evidence'])
                self.assertRegex(record['sha256'], r'^[0-9a-f]{64}$')
                self.assertGreater(record['size_bytes'], 0)
                self.assertRegex(record['source_commit'], r'^[0-9a-f]{40}$')
                self.assertRegex(record['last_present_commit'], r'^[0-9a-f]{40}$')
                artifact_path = ROOT / record['path']
                self.assertEqual(record['status'], 'removed_from_worktree')
                self.assertFalse(artifact_path.exists(), record['path'])

    def test_removed_artifacts_are_not_current_runtime_or_documentation_inputs(self):
        files = active_text_files()
        for record in self.artifacts:
            artifact_path = record['path']
            name = Path(artifact_path).name
            stem = name
            for suffix in ('.tar.gz', '.tar.xz'):
                if stem.endswith(suffix):
                    stem = stem[:-len(suffix)]
            tokens = {artifact_path, name}
            if '-' in stem or len(stem) >= 12:
                tokens.add(stem)
            for source in files:
                try:
                    content = source.read_text(encoding='utf-8')
                except UnicodeDecodeError:
                    continue
                for token in tokens:
                    if token in content:
                        self.fail(f'{source.relative_to(ROOT)} references removed artifact {token}')

    def test_fixture_inventory_hashes_and_sizes_match(self):
        for record in self.manifest['fixture_inventory']['files']:
            with self.subTest(fixture=record['path']):
                fixture = ROOT / record['path']
                self.assertTrue(fixture.is_file(), record['path'])
                self.assertEqual(fixture.stat().st_size, record['bytes'])
                self.assertEqual(sha256_file(fixture), record['sha256'])

    def test_previously_missing_artifacts_remain_missing_and_have_no_code_readers(self):
        files = code_text_files()
        for record in self.manifest['missing_current_tree']:
            reference = record['reported_reference']
            with self.subTest(artifact=reference):
                self.assertEqual(record['status'], 'missing')
                if record['directory_specified']:
                    self.assertFalse((ROOT / reference).exists(), reference)
                name = Path(reference).name
                for source in files:
                    try:
                        content = source.read_text(encoding='utf-8')
                    except UnicodeDecodeError:
                        continue
                    self.assertNotIn(reference, content, source.relative_to(ROOT))
                    self.assertNotIn(name, content, source.relative_to(ROOT))

    def test_large_figure_files_have_manifest_ownership(self):
        known_paths = manifest_paths(self.manifest)
        figures = ROOT / 'docs/developer/figures'
        for artifact in figures.rglob('*'):
            if artifact.is_file() and artifact.stat().st_size >= LARGE_ARTIFACT_BYTES:
                with self.subTest(artifact=str(artifact.relative_to(ROOT))):
                    self.assertIn(artifact.relative_to(ROOT).as_posix(), known_paths)


if __name__ == '__main__':
    unittest.main()
