"""Deterministic research-only selection on the native Original Domain."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from geometry import Geometry, dump, equal_radius


def recoverable(blockers, chosen):
    return all(i in chosen for i in blockers)


def greedy(bundles, benefits, budget, supported):
    """Greedy complete-bundle marginal gain per added atom, with canonical ties."""
    chosen = frozenset()
    actions = []
    stop = "no-positive-feasible-gain"
    for _ in range(budget):
        ranked = []
        for bundle in bundles:
            added = set(bundle) - chosen
            if not added or len(chosen) + len(added) > budget:
                continue
            proposed = chosen | added
            gain = math.fsum(weight for blocking, weight in benefits.items()
                       if recoverable(blocking, proposed) and not recoverable(blocking, chosen))
            if gain > 0:
                ranked.append((-gain / len(added), tuple(sorted(added)), proposed))
        for score, _, proposed in sorted(ranked):
            if supported(proposed):
                actions.append({"added": sorted(proposed-chosen), "gain_per_atom": -score})
                chosen = proposed
                break
        else:
            break
        if len(chosen) == budget:
            stop = "atom-budget"
            break
    return chosen, actions, stop


def scopes(geometry, radius, focus, strategy):
    rows = {}
    for atom in geometry.output:
        if focus is not None and atom not in focus:
            continue
        support = geometry.atom_support(atom)
        native_rows, square = support["row"], support["square"]
        core = square <= radius**2
        eligible = equal_radius(square, square.min()) if strategy == "G2a" else core
        if np.any(eligible & ~geometry.retained[native_rows]):
            rows[atom] = (native_rows, square, core)
    return rows


def choose(g, radius, budget, strategy, focus=None):
    tracked = scopes(g, radius, focus, strategy)
    benefits = {}
    candidates = set()
    near_halos = {}
    for atom, (rows, square, core) in tracked.items():
        mask = equal_radius(square, square.min()) if strategy == "G2a" else core
        if strategy == "G3":
            h = g.meta["isotropic_spacing"]
            if h is None:
                continue
            mask = np.ones(len(square), dtype=bool)
            weights = np.exp(-(square-square.min())/(h/4)**2)
        else:
            weights = np.ones(len(square))
        weights = weights * mask
        for row,w in zip(rows[mask & ~g.retained[rows]],weights[mask & ~g.retained[rows]]):
            blocking = g.blockers[int(row)]
            candidates.add(blocking)
            benefits[blocking] = benefits.get(blocking, 0.) + w/weights.sum()
            if strategy == "G1":
                for h in blocking:
                    d = float(np.linalg.norm(np.array(g.atoms[atom]["position"])-g.atoms[h]["position"]))
                    near_halos[h] = min(d, near_halos.get(h, np.inf))

    def supported(chosen):
        for atom in chosen:
            support = g.atom_support(atom)
            kept = g.rows_retained(support["row"], chosen)
            x = support["square"][kept]
            if len(x) < 3 or np.all(equal_radius(x, x.min())):
                return False
        return True

    if strategy == "G0":
        selection, history, stop = frozenset(), [], "baseline"
    elif strategy == "G1":
        selection, history, stop = frozenset(), [], "no-nearby-supported-atom"
        for distance, halo in sorted((d,h) for h,d in near_halos.items()):
            if distance > 3 or len(selection) == budget:
                break
            updated = selection | {halo}
            if supported(updated):
                selection = updated
                history.append({"added":[halo], "distance":distance})
        if len(selection) == budget:
            stop = "atom-budget"
    else:
        selection, history, stop = greedy(candidates, benefits, budget, supported)
    affected = {}
    recovered = {row for row,blocking in g.blockers.items() if recoverable(blocking, selection)}
    for atom, (rows, square, core) in tracked.items():
        new = np.fromiter((int(r) in recovered for r in rows), bool, count=len(rows))
        if new.any():
            retained = g.retained[rows] | new
            affected[str(atom)] = {
                "core_recovered":int((core&new).sum()), "total_recovered":int(new.sum()),
                "nearest_A":float(np.sqrt(square[retained].min())),
                "core_sd_A":float(np.sqrt(square[retained&core]).std()) if (retained&core).any() else None,
            }
    supports = {
        i:int(g.rows_retained(g.atom_support(i)["row"],selection).sum()) for i in sorted(selection)
    }
    # This is a geometry-only count. The native builder/solver do resource preflight.
    new_output_incidence = sum(sum(int(r) in recovered for r in g.atom_support(i)["row"])
                               for i in g.output) if recovered else 0
    return {"buffer":sorted(selection), "steps":history, "stop":stop, "recovered_rows":len(recovered),
            "buffer_support":supports, "affected":affected, "added_support_incidence":sum(supports.values())+new_output_incidence}


def main(out):
    result = {}
    for path in sorted((out/"geometry").iterdir()):
        g = Geometry(path)
        for label, focus in ([("population",None),("2770",[2770])] if path.name=="9evx-unsharp" else [("population",None)]):
            result[f"{path.name}:{label}"] = {
                strategy:{str(k):choose(g,.5,k,strategy,focus) for k in (1,2,4)}
                for strategy in ("G0","G1","G2a","G2b","G3")
            }
        print(path.name, result[f"{path.name}:population"]["G2b"]["4"], flush=True)
    dump(out/"selection.json",result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=lambda p:Path(p).resolve())
    main(p.parse_args().out)
