"""Shared I/O and comparisons for the runtime regression and offline tools."""
import hashlib
import json
import math
import os
from pathlib import Path
import tarfile

NUMERIC_ENV = {name: '1' for name in
               ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')}
os.environ.update(NUMERIC_ENV)


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def scientific(value):
    if isinstance(value, dict):
        return {key: scientific(item) for key, item in value.items()
                if not key.endswith('seconds') and key != 'process_peak_rss_bytes'}
    if isinstance(value, list):
        return [scientific(item) for item in value]
    return value


def differences(expected, actual, path=''):
    """Exact decisions/availability; numerical values use scaled 1e-10 tolerance."""
    if isinstance(expected, dict) and isinstance(actual, dict):
        if expected.keys() != actual.keys():
            return [path+': keys differ']
        return [item for key in expected for item in differences(expected[key], actual[key], path+'/'+str(key))]
    if isinstance(expected, list) and isinstance(actual, list):
        if len(expected) != len(actual):
            return [path+': lengths differ']
        return [item for i, (a, b) in enumerate(zip(expected, actual)) for item in differences(a, b, path+'/'+str(i))]
    if type(expected) is bool or expected is None or isinstance(expected, str) or type(expected) is int:
        return [] if type(expected) is type(actual) and expected == actual else [path]
    if isinstance(expected, float) and type(actual) in (int, float):
        return [] if math.isfinite(actual) and abs(expected-actual) <= 1e-10*(1+abs(expected)) else [path]
    return [] if expected == actual else [path]


def unpack(catalog_path, dataset, destination):
    catalog_path = Path(catalog_path)
    catalog = read(catalog_path)
    entry = catalog['datasets'][dataset]
    archive = catalog_path.parent/entry['archive']
    require(sha(archive.read_bytes()) == entry['sha256'], 'Fixture archive hash mismatch: '+dataset)
    target = Path(destination)/dataset
    with tarfile.open(archive) as bundle:
        expected = {dataset+'/'+name for name in entry['files']}
        require(set(bundle.getnames()) == expected, 'Unexpected fixture members: '+dataset)
        for name, digest in entry['files'].items():
            member = bundle.getmember(dataset+'/'+name)
            require(member.isfile() and Path(name).name == name, 'Invalid fixture member')
            data = bundle.extractfile(member).read()
            require(sha(data) == digest, 'Fixture member hash mismatch: '+name)
            target.mkdir(parents=True, exist_ok=True)
            file = target/name
            if not file.exists() or sha(file.read_bytes()) != digest:
                file.write_bytes(data)
    return target
