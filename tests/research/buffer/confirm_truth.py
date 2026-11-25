"""Prepare reproducible complete synthetic maps from frozen numerical inputs."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct

import numpy as np

from compare import map_header, read_map
from bundle_selection import choose
from geometry import Geometry, dump


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


def study_case(base, name, changes):
    """Change named true model entries without changing a draw or geometry."""
    truth = [list(p) for p in base["truth"]]
    for atom, columns in changes.items():
        for column, value in columns.items():
            truth[atom-1][column] = value
    if not all(.25 * base["spacing"] <= p[1] <= base["upper"] for p in truth):
        raise ValueError(f"{name}: simulation truth outside the existing B1 interval")
    return base | {"name": name, "truth": truth, "width_diagnostics": True}


def shared_panels(grids, train_seed, validation_seed, repeats):
    nominal = grids["nominal"]
    sigma = nominal["noise_model"]["rms"]
    origin = np.asarray(nominal["positions"][0])
    panel_size = 4 * nominal["spacing"]
    panels = {}
    for name, fixture in grids.items():
        rows = np.array([row[0] for row in fixture["noise"]],dtype=np.int64)
        local = np.array(np.unravel_index(rows,tuple(fixture["dimensions"]),order="F")).T
        xyz = np.asarray(fixture["origin"]) + local * fixture["spacing"]
        indices = np.floor((xyz-origin)/panel_size).astype(np.int64)
        panels[name] = rows, [tuple(i) for i in indices]
    keys = sorted({p for _, ids in panels.values() for p in ids})
    lookup = {p:i for i,p in enumerate(keys)}
    draws = {key:np.random.default_rng(seed).normal(0,sigma,size=(len(keys),repeats))
        for key,seed in (("noise",train_seed),("validation_noise",validation_seed))}
    generated = {}
    for name,(rows,ids) in panels.items():
        slots = [lookup[i] for i in ids]
        data = {key:[[int(row),*sample.tolist()] for row,sample in zip(rows,d[slots],strict=True)]
            for key,d in draws.items()}
        data["seed"] = train_seed
        data["noise_model"] = {"kind":"shared physical Gaussian panels","rms":sigma,"panel_size_A":panel_size,
            "phase_origin":origin.tolist(),"validation_seed":validation_seed,"repeats":repeats}
        generated[name] = grids[name] | data
    return generated


def width_inputs(source, root, geometry_only):
    frozen = {f["name"]:f for f in json.loads(source.read_text())}
    grids = {"nominal":frozen["narrow"],"coarse":frozen["coarse_grid"]}
    geometry = root/"native-geometry"
    if geometry_only:
        dump(root/"geometry-fixtures.json", [x | {"geometry_directory":str(geometry/g)} for g,x in grids.items()])
        return
    selectors = {g:choose(Geometry(geometry/g),.5,1,"G2b") for g in grids}
    grids = {g:x | {"buffer":selectors[g]["buffer"]} for g,x in grids.items()}
    nom = grids["nominal"]
    narrow = nom["truth"][0][1]
    by_stage = {k:[] for k in ("W2","W3","W4","W5","confirm")}

    def add(stage, name, base, tau, amplitude=None, charge=None, other=None):
        change = {1:{1:tau}}
        if amplitude is not None: change[1][0] = amplitude
        if charge is not None: change[1][2] = charge
        if other is not None: change[other] = {1:2.85}
        by_stage[stage].append(study_case(base,f"{stage}-{name}",change))

    for tau in (narrow,.25,.50,.75,1.,1.25,1.50,1.75,2.,2.25,2.85):
        add("W2",f"{tau:.4f}",nom,tau)
    signed_c = frozen["wide_signed_charge"]["truth"][0][2]
    wide_a = frozen["wide_signed_charge"]["truth"][0][0]
    for tau in (.75,1.25,1.50,1.75,2.,2.85):
        add("W3",f"{tau:.2f}-C45",nom,tau,charge=signed_c)
    for tau in (2.,2.85):
        for charge in (nom["truth"][0][2],signed_c):
            add("W3",f"{tau:.2f}-A4000-C{abs(charge):g}",nom,tau,amplitude=wide_a,charge=charge)
        add("W3",f"{tau:.2f}-matched-center",nom,tau,amplitude=nom["truth"][0][0]*(tau/.75)**3)

    spacing = shared_panels(grids,512983,512989,24)
    for grid,base in spacing.items():
        for tau in (.24,.50,1.,1.50,2.):
            add("W4",f"{grid}-absolute-{tau:.2f}",base,tau)
        if grid=="nominal":
            add("W4","nominal-relative-0.3",base,narrow)
        if grid=="coarse" and base["buffer"]!=[3]:
            for tau in (.24,.50,1.,1.50,2.):
                add("W4",f"coarse-forced3-{tau:.2f}",base | {"buffer":[3]},tau)

    for tau in (narrow,.75):
        add("W5",f"M2-neighbor-{tau:.4f}",nom,tau,other=2)
        add("W5",f"M3-buffer-{tau:.4f}",nom,tau,other=3)
    held = shared_panels(grids,513203,513211,48)
    add("confirm","boundary-C3",held["nominal"],2.)
    add("confirm","boundary-C45",held["nominal"],2.,charge=signed_c)
    add("confirm","nominal-relative-0.3",held["nominal"],narrow)
    add("confirm","coarse-relative-0.3",held["coarse"],.24)
    if held["coarse"]["buffer"]!=[3]:
        add("confirm","coarse-forced3-0.24",held["coarse"] | {"buffer":[3]},.24)
    add("confirm","broad-buffer-narrow",held["nominal"],narrow,other=3)

    manifest = {"selector":selectors,"out_of_range":[{"tau":narrow,"spacing":.8,"lower":.2}],"inputs":{}}
    for stage,cases in by_stage.items():
        p = root/f"{stage}-fixtures.json"
        dump(p,cases)
        manifest["inputs"][stage] = {"sha256":sha256(p),"cases":[
            {k:x[k] for k in ("name","buffer","truth","spacing","seed","noise_model","upper")} for x in cases]}
    dump(root/"width-matrix-manifest.json",manifest)


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
    width = modes.add_parser("widths")
    width.add_argument("source", type=Path)
    width.add_argument("output", type=lambda p:Path(p).resolve())
    width.add_argument("--geometry-only", action="store_true")
    args = parser.parse_args()
    if args.mode == "fixtures":
        fixtures(args.source, args.destination)
    elif args.mode == "widths":
        width_inputs(args.source,args.output,args.geometry_only)
    else:
        workflow_map(args.source, args.row_file, args.seed, args.destination)
