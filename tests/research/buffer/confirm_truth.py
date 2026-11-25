"""Prepare reproducible complete synthetic maps from frozen numerical inputs."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct

import numpy as np

from compare import map_header, read_map
from geometry import dump


def fixtures(source, destination):
    frozen = json.loads(source.read_text())
    cases = [x for x in frozen if x["name"] in ("narrow", "coarse_grid", "wide_signed_charge")]
    train = np.random.SeedSequence(851251).spawn(len(cases))
    independent = np.random.SeedSequence(851269).spawn(len(cases))
    for case, train_seed, independent_seed in zip(cases, train, independent, strict=True):
        noise = np.asarray(case["noise"])
        rows = noise[:, 0].astype(np.int64)
        xyz = np.array(np.unravel_index(rows, tuple(case["dimensions"]), order="F")).T
        panels, indices = np.unique(xyz // 4, axis=0, return_inverse=True)
        rms = float(np.sqrt(np.mean(noise[:, 1:] ** 2)))
        repeats = 24
        for key, seed in (("noise", train_seed), ("validation_noise", independent_seed)):
            values = np.random.default_rng(seed).normal(0, rms, size=(len(panels), repeats))[indices]
            case[key] = [[int(row), *v.tolist()] for row, v in zip(rows, values, strict=True)]
        case["seed"] = int(train_seed.entropy)
        case["noise_model"] = {"kind": "block Gaussian", "integer_panel_side": 4, "rms": rms,
                               "train_spawn_key": train_seed.spawn_key, "validation_spawn_key": independent_seed.spawn_key,
                               "validation_seed": int(independent_seed.entropy)}
    dump(destination, cases)


def sha256(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def workflow_map(source, row_file, seed, destination):
    header = map_header(source)
    rows = np.loadtxt(row_file, delimiter="\t", skiprows=1, dtype=np.int64)[:, 0]
    signal = read_map(source, rows)
    sigma = .005 * float(np.sqrt(np.mean(signal ** 2)))
    shutil.copyfile(source, destination)
    m = np.memmap(destination, dtype="<f4", mode="r+", offset=1024 + header["extra"],
                  shape=(int(np.prod(header["shape"])),))
    rng = np.random.default_rng(seed)
    total, sum_square, lo, hi = 0., 0., np.inf, -np.inf
    for start in range(0, len(m), 1 << 18):
        block = m[start:start + (1 << 18)]
        block[:] = block + rng.normal(0, sigma, size=len(block))
        values = block.astype(float)
        total += float(values.sum())
        sum_square += float(values @ values)
        lo, hi = min(lo, float(values.min())), max(hi, float(values.max()))
    m.flush()
    mean = total / len(m)
    rms = max(0, sum_square / len(m) - mean**2) ** .5
    with destination.open("r+b") as f:
        f.seek(76)
        f.write(struct.pack("<3f", lo, hi, mean))
        f.seek(216)
        f.write(struct.pack("<f", rms))
    dump(destination.with_suffix(".json"), {"source": str(source), "source_sha256": sha256(source),
        "sha256": sha256(destination), "seed": seed, "model": "independent voxel Gaussian",
        "sigma": sigma, "signal_scope": "frozen Original Domain rows", "normalization": "raw simulation units"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    small = modes.add_parser("fixtures")
    small.add_argument("source", type=Path)
    small.add_argument("destination", type=Path)
    large = modes.add_parser("workflow")
    large.add_argument("source", type=Path)
    large.add_argument("row_file", type=Path)
    large.add_argument("seed", type=int)
    large.add_argument("destination", type=Path)
    args = parser.parse_args()
    if args.mode == "fixtures":
        fixtures(args.source, args.destination)
    else:
        workflow_map(args.source, args.row_file, args.seed, args.destination)
