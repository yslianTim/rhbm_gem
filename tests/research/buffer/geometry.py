"""Audit production support exports and summarize fixed-domain recovery geometry."""

import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import sqlite3

import numpy as np


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def baseline(snapshot, dataset):
    with sqlite3.connect((snapshot / dataset["baseline"]).as_uri() + "?mode=ro", uri=True) as db:
        result = json.loads(db.execute("select result_json from model_joint_result").fetchone()[0])
        identities = db.execute(
            "select serial_id, chain_id, component_id, atom_id, element "
            "from model_atom where is_selected=1"
        ).fetchall()
    return result, identities


class Geometry:
    def __init__(self, path):
        metadata = json.loads((path / "geometry.json").read_text())
        self.atoms = {a["id"]: a for a in metadata.pop("atoms")}
        self.meta = metadata
        self.ids = list(self.atoms)
        self.output = [i for i, a in self.atoms.items() if a["output"]]
        self.slots = dict(zip(self.ids, np.cumsum([0] + [a["support"] for a in self.atoms.values()])))
        self.support = np.memmap(
            path / "support.bin", dtype=[("row", "=u4"), ("square", "=f8")], mode="r"
        )
        if len(self.support) != sum(a["support"] for a in self.atoms.values()):
            raise ValueError("Incomplete support export")
        rows = np.loadtxt(path / "rows.tsv", delimiter="\t", skiprows=1, dtype=np.int64)
        self.voxels = rows[:, 0]
        self.retained = rows[:, 1].astype(bool)
        if len(self.voxels) != self.meta["original_rows"]:
            raise ValueError("Incomplete row export")
        halo = defaultdict(list)
        for i in self.ids:
            if not self.atoms[i]["output"]:
                for row in self.atom_support(i)["row"].tolist():
                    halo[row].append(i)
        self.blockers = {row: tuple(sorted(ids)) for row, ids in halo.items()}
        if set(self.blockers) != set(np.flatnonzero(~self.retained)):
            raise ValueError("Native fitted mask differs from complete Halo incidence")

    def atom_support(self, atom):
        start = self.slots[atom]
        return self.support[start:start + self.atoms[atom]["support"]]

    def rows_retained(self, rows, chosen):
        return np.fromiter(
            (all(b in chosen for b in self.blockers.get(int(r), ())) for r in rows),
            dtype=bool, count=len(rows),
        )


def equal_radius(square, closest):
    tolerance = 8 * np.finfo(float).eps * max(closest, np.finfo(float).tiny)
    return square <= closest + tolerance


def spread(square):
    if not len(square):
        return None
    r = np.sqrt(square)
    return {"quantiles": np.quantile(r, [0, .25, .5, .75, 1]).tolist(), "sd": float(np.std(r))}


def gaussian(square, kept, tau):
    def weights(values):
        d = values-values.min()
        scaled = np.zeros(len(d))
        positive = d>0
        with np.errstate(over="ignore", under="ignore"):
            scaled[positive] = np.exp(np.log(d[positive])-2*np.log(tau))
        return scaled, np.exp(-scaled)

    def weighted_spread(offset, w):
        use = w>0
        x, v = offset[use], w[use]
        return float(np.sqrt(np.average((x-np.average(x,weights=v))**2,weights=v)))

    shift, w = weights(square)
    fit_shift, fitted = weights(square[kept])
    return {"retention":float(w[kept].sum()/w.sum()),
            "original_effective":float(w.sum()**2/(w@w)),
            "fitting_effective":float(fitted.sum()**2/(fitted@fitted)),
            "original_square_sd_over_tau2":weighted_spread(shift,w),
            "fitting_square_sd_over_tau2":weighted_spread(fit_shift,fitted)}


def audit(out, tag, dataset, scales):
    snapshot = out / "snapshot"
    g = Geometry(out / "geometry" / tag)
    saved, identity = baseline(snapshot, dataset)
    saved_ids = [int(saved["atom_ids"][i]) for i in saved["fitting_to_original_contributors"]]
    assert g.output == saved_ids == sorted(t[0] for t in identity)
    assert {a[1] for a in identity} == set(dataset["chains"].split(","))
    for atom, chain, residue, name, element in identity:
        a = g.atoms[atom]
        assert (chain, residue, name, element) == (a["chain"], a["residue"], a["atom"], a["element"])
    np.testing.assert_array_equal(g.voxels, np.asarray(saved["row_ids"], dtype=np.int64))
    np.testing.assert_array_equal(g.voxels[g.retained], np.asarray(saved["fitting_row_ids"], dtype=np.int64))
    counts = defaultdict(Counter)
    breakdown = defaultdict(Counter)
    q = next(e["value"] for e in saved["evidence"] if e["name"] == "global-width-calibration")
    seeds = {int(a["id"]): a for a in saved["initialization"]["atoms"]}
    state = saved["assembled_state"]
    points = {int(saved["atom_ids"][p["atom"]]):p for p in saved["target_evidence"]["point_atoms"]}
    h = g.meta["isotropic_spacing"]
    radii = {"0.5A": .5}
    if scales:
        radii.update({k: h*s if h is not None else None for k, s in
                      [("0.5h", .5), ("h", 1.), ("half-diagonal", np.sqrt(3)/2)]})
    bundle_index = {}
    phase = "scales" if scales else "population"
    with gzip.open(out/f"{tag}-{phase}.jsonl.gz", "wt") as f:
        for index, atom in enumerate(g.output):
            support = g.atom_support(atom)
            rows, squares = support["row"], support["square"]
            kept = g.retained[rows]
            nearest = equal_radius(squares, squares.min())
            blocked_near = [g.blockers[int(r)] for r in rows[nearest & ~kept]]
            record = dict(g.atoms[atom], original=len(rows), fitting=int(kept.sum()),
                          excluded=int((~kept).sum()), original_nearest=float(np.sqrt(squares.min())),
                          fitting_nearest=float(np.sqrt(squares[kept].min())),
                          original_radial=spread(squares), fitting_radial=spread(squares[kept]),
                          excluded_radial=spread(squares[~kept]),
                          blocking_sizes=dict(Counter(len(g.blockers[int(r)]) for r in rows[~kept])),
                          blocked_incidence={int(g.voxels[r]):g.blockers[int(r)] for r in rows[~kept]} if not scales else {},
                          direct_nearest=min(blocked_near,key=lambda b:(len(b),b)) if blocked_near else None,
                          endpoint=[state["ac"][2*index], state["tau"][index], state["ac"][2*index+1]],
                          evidence=points[atom], cores={})
            for key, radius in radii.items():
                if radius is None:
                    counts[key]["unsupported"] += 1
                    record["cores"][key] = None
                    continue
                mask = squares <= radius**2
                n, b = int(mask.sum()), int((mask & kept).sum())
                lost = sorted(set(g.blockers[int(r)] for r in rows[mask & ~kept]),key=lambda v:(len(v),v))
                core = record["cores"][key] = {"original":n,"fitting":b,"retention":b/n if n else None,
                                              "radial":spread(squares[mask & kept]),"bundles":lost}
                totals = counts[key]
                totals["no_original_core"] += n == 0
                totals["core_lost"] += n > 0 and b == 0
                totals["core_partial"] += 0 < b < n
                totals["radially_concentrated"] += bool(b >= 2 and core["radial"]["quantiles"][-1]
                                                       - core["radial"]["quantiles"][0] <= .05)
                totals["nearest_lost"] += bool(np.any(nearest & ~kept))
                if lost:
                    totals["single_can_restore_any"] += len(lost[0]) == 1
                    totals["needs_multiple_for_any"] += len(lost[0]) > 1
                    totals["union_for_all_requires_multiple"] += len(set().union(*map(set,lost))) > 1
                if n > 0 and b == 0:
                    for field in ("chain", "element", "atom"):
                        breakdown[f"{key}:{field}"][str(record[field])] += 1
                for bundle in lost:
                    if bundle in bundle_index:
                        continue
                    chosen = set(bundle)
                    entries = [g.atom_support(i)["row"] for i in bundle]
                    fitting_counts = [int(g.rows_retained(e, chosen).sum()) for e in entries]
                    bundle_index[bundle] = {
                        "buffer_support": dict(zip(bundle, fitting_counts)),
                        "new_buffer_fit_incidence":sum(fitting_counts),
                    }
            record["nearest_change"] = record["fitting_nearest"] - record["original_nearest"]
            if scales:
                seed = seeds[atom]
                ref = seed["mdpde"][1] if seed["initial_stage_fit_qualified"] and seed["seed_source"] != "median-fallback" else None
                record["gaussian"] = {key: gaussian(squares, kept, tau) if tau and tau > 0 else None
                                      for key,tau in [("floor",h/4 if h else None),("initial",ref),("q50",q)]}
            f.write(json.dumps(record, allow_nan=False) + "\n")
    return {"output_count": len(g.output), "rows": g.meta["original_rows"],
            "scales": radii, "counts": counts, "core_loss_breakdown": breakdown,
            "blocking_size_all_voxels": Counter(map(len, g.blockers.values())),
            "bundles": [{"ids": key, **value} for key, value in sorted(bundle_index.items())],
            "geometry": g.meta, "baseline_observation_scale": saved["observation_scale"]}


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("out", type=lambda x: Path(x).resolve())
    cli.add_argument("--scales", action="store_true")
    args = cli.parse_args()
    sources = json.loads((args.out/"snapshot/SHA256.json").read_text())["datasets"]
    summary = {}
    for tag, data in sources.items():
        summary[tag] = audit(args.out, tag, data, args.scales)
        print(tag, summary[tag]["output_count"], summary[tag]["counts"]["0.5A"], flush=True)
    assert sum(x["output_count"] for x in summary.values()) == 25558
    dump(args.out/("scales.json" if args.scales else "population.json"), summary)
