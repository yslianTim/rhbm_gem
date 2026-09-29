"""Shared file, JSON, and hashing helpers for developer experiments."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def write_json(path, value):
    Path(path).write_text(
        json.dumps(value, indent=2, sort_keys=False, allow_nan=False) + "\n",
        encoding="utf-8")


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def digest(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
