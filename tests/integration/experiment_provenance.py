"""Source and build fingerprint calculations, without campaign expectations."""
import hashlib
import io
import platform
import subprocess
import tarfile
from pathlib import Path

from experiment_io import ROOT, digest, sha, sha256_file


def source_hash(root):
    files = sorted(p for folder in ('src', 'include', 'cmake') for p in (root / folder).rglob('*')
                   if p.is_file() and p.suffix in ('.cpp', '.hpp', '.h', '.cmake', '.txt'))
    files += [root / 'CMakeLists.txt']
    return digest({str(p.relative_to(root)): sha(p) for p in files})


def git_source_hash(commit, root=ROOT):
    data = subprocess.check_output(['git', 'archive', commit, 'src', 'include', 'cmake', 'CMakeLists.txt'], cwd=root)
    with tarfile.open(fileobj=io.BytesIO(data)) as archive:
        return digest({m.name: hashlib.sha256(archive.extractfile(m).read()).hexdigest() for m in archive
                       if m.isfile() and Path(m.name).suffix in ('.cpp', '.hpp', '.h', '.cmake', '.txt')})


def build_fingerprint(root):
    cache = (root / 'CMakeCache.txt').read_text()
    source = Path(next(line.split('=', 1)[1] for line in cache.splitlines()
                       if line.startswith('CMAKE_HOME_DIRECTORY:')))
    binaries = [root / 'bin/joint_fixed_neighbor_experiment', root / 'bin/joint_component_runtime',
                root / 'bin/joint_offline_diagnostic', root / 'bin/joint_statistical_experiment',
                root / 'bin/RHBM-GEM', *sorted((root / 'src').glob('librhbm_gem.*'))]
    return dict(source_root=str(source), source_sha256=source_hash(source), cache=cache,
                binaries={str(p.relative_to(root)): sha(p) for p in binaries if p.is_file()},
                configuration=(root / 'generated/Release/SimulationConfiguration-CXX.txt').read_text(),
                linked_libraries=subprocess.check_output(
                    ['otool', '-L', str(root / 'src/librhbm_gem.dylib')] if platform.system() == 'Darwin'
                    else ['ldd', str(root / 'src/librhbm_gem.so')], text=True))


def require_current_build(build):
    cache = (build / 'CMakeCache.txt').read_text()
    source = Path(next(line.split('=', 1)[1] for line in cache.splitlines() if line.startswith('CMAKE_HOME_DIRECTORY:')))
    libraries = list((build / 'src').glob('librhbm_gem.*'))
    outputs = [build / 'bin/joint_fixed_neighbor_experiment', *libraries]
    if not libraries or any(not p.is_file() for p in outputs): raise ValueError(f'Build is incomplete: {build}')
    sources = [p for folder in ('src', 'include', 'cmake') for p in (source / folder).rglob('*')
               if p.is_file() and p.suffix in ('.cpp', '.hpp', '.h', '.cmake', '.txt')]
    sources += [source / 'CMakeLists.txt']
    core_time = max(p.stat().st_mtime_ns for p in sources)
    driver_sources = [source / 'tests/experiments/joint_fixed_neighbor.cpp',
                      *list((source / 'tests/support').glob('Joint*.hpp')),
                      *list((source / 'tests/support').glob('Joint*.cpp'))]
    driver_time = max(p.stat().st_mtime_ns for p in driver_sources if p.is_file())
    if core_time > min(p.stat().st_mtime_ns for p in libraries) or max(core_time, driver_time) > (build / 'bin/joint_fixed_neighbor_experiment').stat().st_mtime_ns:
        raise ValueError(f'Sources are newer than measured binaries; rebuild before starting a campaign: {build}')


def source_and_build_provenance(executable, source):
    build = executable.parent.parent
    files = [executable, build / 'src/librhbm_gem.dylib', build / 'CMakeCache.txt']
    links = {str(p): subprocess.check_output(['otool', '-L', str(p)], text=True)
             for p in files[:2] if p.is_file()}
    dependencies = {Path(line.strip().split(' (')[0]) for output in links.values()
                    for line in output.splitlines()[1:] if line.strip().startswith('/')}
    eigen_metadata = []
    if files[2].is_file():
        for line in files[2].read_text().splitlines():
            if line.startswith('Eigen3_DIR:PATH='):
                eigen_metadata = sorted(Path(line.split('=', 1)[1]).glob('*.cmake'))
    source_files = sorted(p for folder in ('src', 'include', 'tests', 'cmake') for p in (source / folder).rglob('*')
                          if p.is_file() and p.suffix in {'.cpp', '.hpp', '.py', '.cmake', '.txt', '.json'})
    source_files += [source / 'CMakeLists.txt']
    return {'artifacts': {str(p): sha256_file(p) for p in files if p.is_file()},
            'source_files': {str(p.relative_to(source)): sha256_file(p) for p in source_files},
            'linked_libraries': links,
            'dependency_hashes': {str(p): sha256_file(p) for p in sorted(dependencies) + eigen_metadata if p.is_file()}}
