"""Independent snapshot and structural checks for the runtime fixture catalog."""
from pathlib import Path
from joint_runtime_support import read, require, sha
import numpy as np


def digest(path):
    return sha(Path(path).read_bytes())


def load(path):
    path = Path(path)
    snapshot = read(path/"snapshot.json")
    for name in ("voxels", "contributors"):
        require(digest(path/(name+".csv")) == snapshot[name+"_sha256"], "Snapshot hash mismatch.")
    table = np.atleast_1d(np.genfromtxt(path/"contributors.csv", delimiter=",", names=True,
                                     dtype=[("row", int), ("atom", int), ("square", float)]))
    voxels = np.atleast_2d(np.loadtxt(path/"voxels.csv", delimiter=",", skiprows=1))
    ids = [str(a["serial_id"]) for a in read(path/"dataset.json")["atoms"]]
    require(len(table) == snapshot["memberships"] and len(voxels) == snapshot["rows"] and
            len(ids) == snapshot["atoms"] and len(set(ids)) == len(ids), "Snapshot population mismatch.")
    require(np.all((table["row"] >= 0) & (table["row"] < len(voxels))) and
            np.all((table["atom"] >= 0) & (table["atom"] < len(ids))) and
            np.all(np.isfinite(table["square"]) & (table["square"] >= 0) & (table["square"] <= 6.25)),
            "Invalid contributor entries.")
    pairs = list(zip(table["row"], table["atom"]))
    require(pairs == sorted(set(pairs)), "Contributors are not unique ordered CSR entries.")
    offsets = np.r_[0, np.cumsum(np.bincount(table["row"], minlength=len(voxels)))]
    require(offsets.tolist() == snapshot["row_offsets"], "Changed CSR offsets.")
    return {"table": table, "y64": voxels[:, 7], "y32": voxels[:, 8], "ids": ids,
            "hash": digest(path/"snapshot.json")}


def verify_census(data, census):
    """Reconstruct connectivity independently; no coefficient values participate."""
    ids, table = data["ids"], data["table"]
    parent = list(range(len(ids)))

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a

    first = {}
    for row, atom, _ in table:
        if row in first:
            parent[root(int(atom))] = root(first[row])
        else:
            first[row] = int(atom)
    groups = {}
    for a in range(len(ids)):
        groups.setdefault(root(a), []).append(a)
    components = sorted(groups.values(), key=lambda atoms: min(ids[a] for a in atoms))
    require(len(components) == census["component_count"], "Wrong component count.")
    atom_component = [-1]*len(ids); row_component = [-1]*len(data["y64"])
    for c, (atoms, record) in enumerate(zip(components, census["components"])):
        rows = sorted(row for row, atom in first.items() if root(atom) == root(atoms[0]))
        require(record["id"] == min(ids[a] for a in atoms) and record["atoms"] == atoms and
                record["rows"] == rows, "Wrong component identity or membership.")
        atom_map = [-1]*len(ids); row_map = [-1]*len(data["y64"])
        for k, a in enumerate(atoms): atom_map[a] = k
        for k, r in enumerate(rows): row_map[r] = k
        require(record["parent_atom_to_local"] == atom_map and record["parent_row_to_local"] == row_map,
                "Wrong inverse local mapping.")
        count = int(np.count_nonzero(np.isin(table["atom"], atoms)))
        require(record["memberships"] == count, "Lost component memberships.")
        for a in atoms: atom_component[a] = c
        for r in rows: row_component[r] = c
    require(atom_component == census["atom_component"] and row_component == census["row_component"],
            "Wrong inverse component mapping.")
    require(census["constant_rows"] == [r for r, c in enumerate(row_component) if c < 0], "Lost constant rows.")
    require(census["unobserved_atoms"] == [a for a in range(len(ids)) if a not in table["atom"]], "Lost unobserved atoms.")
    require(census["snapshot_sha256"] == data["hash"], "Census snapshot mismatch.")
    return True
