"""Compare certified research endpoints on fixed voxel identities."""

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

from geometry import Geometry, baseline, dump, gaussian


def map_header(path):
    with path.open("rb") as file:
        h = file.read(1024)
    i = np.frombuffer(h, dtype="<i4")
    f = np.frombuffer(h, dtype="<f4")
    return {"shape":i[:3].tolist(),"mode":int(i[3]),"start":i[4:7].tolist(),
            "sampling":i[7:10].tolist(),"cell":f[10:16].tolist(),"axes":i[16:19].tolist(),
            "origin":f[49:52].tolist(),"extra":int(i[23])}


def read_map(path, ids):
    h = map_header(path)
    if h["mode"]!=2 or h["axes"]!=[1,2,3]:
        raise ValueError(f"Unsupported shared grid layout {h}")
    return np.memmap(path, dtype="<f4", mode="r", offset=1024+h["extra"])[ids].astype(float)


def kernel_columns(square, width):
    distance = np.sqrt(square)
    g = np.exp(-.5*square/(width*width))/(2*np.pi*width*width)**1.5
    at_zero = math.sqrt(2/math.pi)/width
    k = np.fromiter((math.erf(r/(width*math.sqrt(2)))/r if r>=1e-5 else at_zero
                     for r in distance), float, len(distance))
    return g, k


def field(g, points):
    result=np.zeros(len(g.voxels))
    for atom,(a,t,c) in points.items():
        incidence=g.atom_support(atom)
        gaussian, charge = kernel_columns(incidence["square"], t)
        result[incidence["row"]]+=a*gaussian+c*charge
    return result


def read_parameters(path):
    with (path/"parameters.tsv").open() as file:
        return {int(r["serial"]):r for r in csv.DictReader(file,delimiter="\t")}


def endpoint(parameters):
    return np.array([[float(p[key]) for key in ("A","tau","C")] for p in parameters])


def archived(root, tag, g, ds):
    s, rows = baseline(root/"snapshot", ds)
    keys = [int(s["atom_ids"][i]) for i in s["fitting_to_original_contributors"]]
    assert keys==g.output
    fit=s["assembled_state"]
    points={i:(fit["ac"][2*k],fit["tau"][k],fit["ac"][2*k+1]) for k,i in enumerate(keys)}
    halo={int(p["atom_id"]):(p["amplitude"],p["width"],p["charge"]) for p in s["halo_median_reconstruction"]["halos"]}
    signal=field(g,points)
    full=signal+field(g,halo)
    # Consume the frozen, officially generated all-target stage results read-only.
    sys.path.insert(0,str(root/"snapshot/scripts"))
    from domain import Saved
    saved=Saved(root/"snapshot"/ds["baseline"])
    stages,groups=saved.replay_stage(())
    local={str(i):{"reason":x["post_joint_local"]["reason"],
                   "available":x["post_joint_local"]["information"] is not None,
                   "qscore":x["post_joint_qscore"]["scores"]}
           for i,x in stages.items() if i in points}
    group_points={str(row["key"]):row["summary"]["inference"]["mdpde"] for row in groups
                  if row.get("summary") and row["summary"].get("inference")}
    return {"source":"frozen-native-G0-with-reconstructed-projection","kept":g.retained,"fit":signal,"total":full,
            "point":np.array(list(points.values())), "point_ids":keys, "fit_points":points,
            "norm":s["metadata"]["map_normalization"]["divisor"],"local":local,
            "groups":len(groups),"group_points":group_points,"state":s}


def measured(path,g):
    s=json.loads((path/"summary.json").read_text())
    p=np.loadtxt(path/"prediction.tsv",delimiter="\t",skiprows=1)
    np.testing.assert_array_equal(p[:,0],g.voxels)
    parameters=read_parameters(path)
    assert [i for i,p in parameters.items() if p["output"]=="1"]==g.output
    assert set(parameters)-set(g.output)==set(s["buffer"])
    group_file=path/"groups.json"
    groups=json.loads(group_file.read_text()) if group_file.exists() else []
    return {"source":"fresh-research-native-solver","kept":p[:,1].astype(bool),"fit":p[:,2],"total":p[:,3],
            "alternative":p[:,4],"point":endpoint(parameters[i] for i in g.output),"point_ids":g.output,
            "fit_points":dict(zip(parameters,endpoint(parameters.values()))),
            "norm":s["map_divisor"],"raw":s,"local":s.get("post_joint",{}),
            "parameters":parameters, "groups":len(groups) if group_file.exists() else None,
            "group_points":{str(row["key"]):row["mdpde"] for row in groups if row.get("mdpde")}}


def describe(array):
    return np.quantile(array,[0,.05,.5,.95,.99,1],axis=0).tolist() if len(array) else None


def error(e):
    return {"count":len(e),"RMSE":float(np.sqrt(np.mean(e*e))),"max_abs":float(np.max(np.abs(e)))} if len(e) else {"count":0}


def stage(m):
    local=m["local"]
    reasons=Counter(x.get("reason") or "available" for x in local.values())
    return {"counts":reasons,"qscore":describe(np.array([x["qscore"] for x in local.values()
                                                         if x.get("qscore") is not None])),
            "groups":m.get("groups"),"recorded_outputs":len(local),
            "peeling_available":sum(x.get("peeling_available",0) for x in local.values()) if any("peeling_available" in x for x in local.values()) else None,
            "peeling_samples":sum(x.get("peeling_samples",0) for x in local.values()) if any("peeling_samples" in x for x in local.values()) else None}


def feature_stats(path,g):
    if not (path/"features.tsv").exists(): return None
    rows=np.loadtxt(path/"features.tsv",delimiter="\t")
    assert set(rows[:,0]) <= set(g.output)
    a=rows[:,5]/max(abs(rows[:,5]))
    z=(a-np.mean(a))**2
    ids=np.argsort(z)[-5:][::-1]
    standardized=path/"standardized.tsv"
    return {"raw_missing":(~np.isfinite(rows[:,1:])).sum(axis=0).tolist(),
            "native_standardized_quantiles":describe(np.loadtxt(standardized)[:,1:]) if standardized.exists() else None,
            "rows":len(rows),"amplitude_variance_share":[[int(rows[i,0]),float(z[i]/z.sum())] for i in ids],
            "umap":json.loads((path/"umap.json").read_text())}


def embedding(path,g):
    direct=path/"umap.tsv"
    if direct.is_file(): return np.loadtxt(direct)
    csvs=list(path.glob("umap_embedding_*.csv"))
    for data in csvs:
        with data.open() as file:
            reader=csv.DictReader(file)
            columns=("serial id","umap x","umap y","second umap x","second umap y")
            if not set(columns)<=set(reader.fieldnames): continue
            values=np.array([[float(row[k]) for k in columns] for row in reader])
            assert set(values[:,0])<=set(g.output)
            return values
    return None


def knn_overlap(first, second, k=15):
    ids,ia,ib=np.intersect1d(first[:,0],second[:,0],return_indices=True)
    if len(ids)<=k: return None
    overlap={}
    for stage,columns in (("first",slice(1,3)),("second",slice(3,5))):
        a,b=first[ia,columns],second[ib,columns]
        fractions=[]
        for begin in range(0,len(ids),128):
            stop=min(len(ids),begin+128)
            nbors=[]
            for x in (a,b):
                d=np.sum((x[begin:stop,None,:]-x[None,:,:])**2,axis=2)
                d[np.arange(stop-begin),np.arange(begin,stop)]=np.inf
                nbors.append(np.argpartition(d,k-1,axis=1)[:,:k])
            fractions.extend(len(set(u)&set(v))/k for u,v in zip(*nbors))
        overlap[stage]=describe(np.asarray(fractions))
    return {"common_ids":len(ids),"k":k,"overlap_quantiles":overlap}


def radial_diagnostics(g, atom, selected, point, fit_points):
    support=g.atom_support(atom)
    square=support["square"]
    kept=selected[support["row"]]
    a,t,c=point
    x=square[kept]
    exp=np.exp(-.5*x/(t*t))
    at_zero=math.sqrt(2/math.pi)/t
    gcol,kcol=kernel_columns(x,t)
    deriv=a*gcol*(x/(t*t)-3)-c*at_zero*exp
    basis=np.column_stack((gcol,kcol))
    norm=np.linalg.norm
    projected=deriv-basis@np.linalg.lstsq(basis,deriv,rcond=None)[0]
    total=norm(a*gcol+c*kcol)
    # This is a focal diagnostic over the native retained rows. Hold the
    # co-fitted neighbor widths fixed and include their free A/C columns.
    rows = support["row"][kept]
    row_slots = np.full(len(g.voxels), -1, dtype=np.intp)
    row_slots[rows] = np.arange(len(rows))
    nuisance = []
    for other, (amplitude, width, _) in fit_points.items():
        if other == atom:
            continue
        shared = g.atom_support(other)
        slots = row_slots[shared["row"]]
        overlap = slots >= 0
        if not np.any(overlap):
            continue
        for column in kernel_columns(shared["square"][overlap],width)[0 if amplitude>0 else 1:]:
            length = norm(column)
            if length > 0:
                mapped = np.zeros(len(rows))
                mapped[slots[overlap]] = column / length
                nuisance.append(mapped)
    conditional = basis
    if nuisance:
        controls = np.column_stack(nuisance)
        conditional = basis - controls @ np.linalg.lstsq(controls,basis,rcond=None)[0]
    norms = norm(conditional,axis=0)
    correlation = float(conditional[:,0] @ conditional[:,1] / np.prod(norms)) if np.all(norms>0) else None
    return {"gaussian":gaussian(square,kept,t),
            "projected_log_width_norm_own_AC_only":float(norm(projected)),
            "AC_cosine_conditioned_on_overlapping_fitting_AC_fixed_widths":correlation,
            "overlapping_fitting_nuisance_columns":len(nuisance),
            "own_gaussian_charge_cosine":float(gcol@kcol/(norm(gcol)*norm(kcol))) if norm(gcol)>0 else None,
            "own_AC_cancellation":float((norm(a*gcol)+norm(c*kcol))/total) if total>0 else None}


def compare(g, left, right, reference, focus):
    base=left["kept"]
    cand=right["kept"]
    assert np.all(~base|cand), "Fixed-domain nesting failed"
    core=np.zeros(len(base),dtype=bool)
    for i in g.output:
        s=g.atom_support(i)
        core[s["row"][s["square"]<=.25]]=True
    domains={"common":base&cand,"recovered":cand&~base,"original":np.ones(len(base),bool),
             "core":core,"recovered_core":core&cand&~base}
    if focus in g.output:
        s=g.atom_support(focus)
        domains["focal_neighborhood"]=np.isin(np.arange(len(base)),s["row"])
        f=np.zeros(len(base),dtype=bool)
        f[s["row"][s["square"]<=.25]]=True
        domains["focal_core"]=f
    # Report in raw map units so both training normalization factors remain explicit.
    errs=[m["total"]*m["norm"]-reference for m in (left,right)]
    target_errors={}
    for k,mask in domains.items():
        target_errors[k] = [error(e[mask]) for e in errs]
    per_target, per_target_recovered = [], []
    for i in g.output:
        rows=g.atom_support(i)["row"]
        per_target.append([error(e[rows])["RMSE"] for e in errs])
        recovered=rows[domains["recovered"][rows]]
        if len(recovered):
            per_target_recovered.append([error(e[recovered])["RMSE"] for e in errs])
    p=np.array(per_target)
    tracked={}
    for i in (focus, 2770, 8607, 7995, 5159, 3583, 3242, 32994, 372, 735):
        if i not in g.output: continue
        k=g.output.index(i)
        tracked[i]={"A_tau_C":[m["point"][k,:].tolist() for m in (left,right)],
                    "raw_RMSE":p[k,:].tolist(),
                    "local":[m["local"].get(str(i)) for m in (left,right)],
                    "shape":[radial_diagnostics(g,i,m["kept"],m["point"][k,:],m["fit_points"]) for m in (left,right)]}
    rec=np.array(per_target_recovered).reshape(-1,2)
    joint_scores=[]
    transitions=Counter()
    for atom in sorted(set(left["local"])&set(right["local"])):
        l,r=left["local"][atom],right["local"][atom]
        transitions[(l.get("reason") or "available")+" → "+(r.get("reason") or "available")]+=1
        if l.get("qscore") is not None and r.get("qscore") is not None:
            joint_scores.append(np.array(r["qscore"])-l["qscore"])
    groups=[left.get("group_points",{}),right.get("group_points",{})]
    group_deltas=[np.array(groups[1][k])-groups[0][k] for k in sorted(set(groups[0])&set(groups[1]))]
    return {"domains_raw":target_errors,"paired_qscore_delta":describe(np.array(joint_scores)),
            "paired_local_reason_transitions":transitions,
            "paired_group_mdpde_delta":describe(np.array(group_deltas)),
            "fit_only_common_raw":[error((m["fit"]*m["norm"]-reference)[domains["common"]]) for m in (left,right)],
            "external_halo_signal_original_raw":[error((m["total"]-m["fit"])*m["norm"]) for m in (left,right)], "paired_target_delta_RMSE":describe(p[:,1]-p[:,0]),
            "paired_target_delta_recovered_RMSE":describe(rec[:,1]-rec[:,0]),
            "target_recovered_count":len(rec),
            "parameters":[describe(m["point"]) for m in (left,right)],
            "zero_A":[int((m["point"][:,0]==0).sum()) for m in (left,right)],
            "parameter_delta":describe(right["point"]-left["point"]),
            "halo_sensitivity_RMS":error((right["alternative"]-right["total"])*right["norm"]) if "alternative" in right else None,
            "stages":[stage(m) for m in (left,right)],"tracked":tracked}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root",type=lambda p:Path(p).resolve())
    p.add_argument("tag")
    p.add_argument("first")
    p.add_argument("second")
    p.add_argument("--reference",type=Path)
    p.add_argument("--focus",type=int,default=2770)
    a=p.parse_args()
    g=Geometry(a.root/"geometry"/a.tag)
    ds=json.loads((a.root/"snapshot/SHA256.json").read_text())["datasets"][a.tag]
    def run(name):
        return a.root/"runs"/f"{a.tag}-{name}"
    def load(name):
        return archived(a.root,a.tag,g,ds) if name=="archive" else measured(run(name),g)
    l,r=load(a.first),load(a.second)
    reference=a.reference or a.root/"snapshot"/ds["map"]
    result=compare(g,l,r,read_map(reference,g.voxels),a.focus)
    result["sources"]=[l["source"],r["source"]]
    result["regression"]=bool(np.array_equal(l["point"],r["point"]))
    result["max_parameter_abs_delta"]=np.max(abs(r["point"]-l["point"]),axis=0).tolist()
    result["features"]=[feature_stats(run(a.first),g),feature_stats(run(a.second),g)]
    embeddings=[embedding(run(n),g) for n in (a.first,a.second)]
    if all(e is not None for e in embeddings):
        result["knn"]=knn_overlap(*embeddings)
    dump(a.root/f"compare-{a.tag}-{a.first}-{a.second}.json",result)
