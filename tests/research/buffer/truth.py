"""Prepare frozen D4 draws and analyze native JCE paired truth estimates."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

from geometry import dump


def prepare(root):
    sys.path.insert(0,str(root/"snapshot/scripts"))
    from multiatom import grid_for, SCENARIOS, IDS, N, Saved, S0, bounds
    s=Saved(S0)
    seed=2907
    fixtures=[]
    for name,(h,ratio,a,halo_a) in SCENARIOS.items():
        if h==.574: h=s.h[0]
        positions=np.array([s.xyz[i] for i in IDS])
        xyz,_,index=grid_for(s,positions,h)
        phase=((positions[0]-s.o)/s.h)%1
        lower=index.min(axis=0)
        origin=(positions[0]-phase*h)+lower*h
        dims=index.max(axis=0)-lower+1
        x,y,z=(index-lower).T
        rows=x+dims[0]*(y+dims[1]*z)
        nearest=np.rint((xyz-s.o)/s.h).astype(int)
        x,y,z=nearest.T
        h1,h2=s.halves(x+s.n[0]*(y+s.n[1]*z))
        panels,inverse=np.unique(index//4,axis=0,return_inverse=True)
        rng=np.random.default_rng(seed)
        delta=(h1-h2)/2
        noise=np.array([rng.choice([-1.,1.],len(panels))[inverse]*delta for _ in range(N)])
        truth=[[a,ratio*h if ratio else 2.85,-3. if ratio else -45.],
               [165.,.72,-2.6],[halo_a,.79,-3.6],[170.,.64,-1.8]]
        fixtures.append({"name":name,"seed":seed,"spacing":h,"origin":origin.tolist(),
                         "dimensions":dims.tolist(),"positions":positions.tolist(),"truth":truth,
                         "upper":bounds(s)[1],"noise":[[int(i),*v] for i,v in zip(rows,noise.T)]})
    dump(root/"truth-fixtures.json",fixtures)


def metrics(e):
    if not len(e): return None
    absolute=np.abs(e)
    return {"bias":e.mean(axis=0).tolist(),"rmse":np.sqrt(np.mean(e**2,axis=0)).tolist(),
            "median_absolute":np.median(absolute,axis=0).tolist(),
            "p95_absolute":np.quantile(absolute,.95,axis=0).tolist(),
            "p99_absolute":np.quantile(absolute,.99,axis=0).tolist(),
            "maximum_absolute":absolute.max(axis=0).tolist()}


def summarize(root):
    result={}
    for line in (root/"truth-results.jsonl").read_text().splitlines():
        d=json.loads(line)
        truth=np.array(d["truth"])[:2]
        buffers={name:{v["replicate"]:v for v in d["fits"] if v["strategy"]==name}
                 for name in ("G0","single")}
        pair=[rep for rep in buffers["G0"] if all(x[rep]["qualified"] for x in buffers.values())]
        r={}
        for name,runs in buffers.items():
            valid=[v for v in runs.values() if v["qualified"]]
            errors=np.array([np.array(v["points"])[:2]-truth for v in valid])
            paired=np.array([np.array(runs[rep]["points"])[:2]-truth for rep in pair])
            pred={k:{} for k in ("common","original","core","recovered")}
            for v in valid:
                rep=v["replicate"]
                for key, value in v["prediction_RMSE"].items():
                    if value is not None:
                        pred[key][rep]=value
                outer=buffers["single"][rep].get("recovered_rows")
                if outer:
                    values=np.array(v["residuals"])[outer]
                    pred["recovered"][rep]=float(np.sqrt(np.mean(values**2)))
            all_predictions, paired_predictions = {}, {}
            for domain, values in pred.items():
                all_predictions[domain] = float(np.median(list(values.values()))) if values else None
                matched = [values[i] for i in pair if i in values]
                paired_predictions[domain] = np.quantile(matched,[.5,.95,.99,1]).tolist() if matched else None
            r[name]={"attempted":len(runs),"qualified":len(valid),
                     "conditional_errors":metrics(errors),"paired_errors":metrics(paired),
                     "prediction_median":all_predictions,"prediction_paired_quantiles":paired_predictions,
                     "profile_evaluations":sum(v["evaluations"] for v in runs.values())}
        if pair:
            delta=np.array([abs(np.array(buffers["single"][rep]["points"])[:2]-truth)-
                            abs(np.array(buffers["G0"][rep]["points"])[:2]-truth) for rep in pair])
            rng=np.random.default_rng(2907)
            boots=rng.choice(len(pair),size=(2000,len(pair)),replace=True)
            r["paired_absolute_error_delta"]={"mean":delta.mean(axis=0).tolist(),
                "paired_replicate_bootstrap_95":np.quantile(delta[boots].mean(axis=1),[.025,.975],axis=0).tolist()}
        r["qualified_pairs"]=len(pair)
        result[d["scenario"]]=r
    dump(root/"truth-summary.json",result)


def workflow_truth(root):
    from compare import Geometry, archived, measured
    tag="6z6u-sim"
    ds=json.loads((root/"snapshot/SHA256.json").read_text())["datasets"][tag]
    g=Geometry(root/"geometry"/tag)
    path=root/"runs/6z6u-sim-g2"
    status=json.loads((path/"exit.json").read_text())
    manifest=json.loads(Path(str(root/"snapshot"/ds["map"])+".simulation.json").read_text())
    assert manifest["output"]["map_sha256"]==ds["map_sha256"]
    reference={r["serial_id"]:[r["element"],r["effective_gaussian_width"],r["charge_used"]]
               for r in manifest["atoms"]}
    truth=np.array([reference[i] for i in g.output])
    base=archived(root,tag,g,ds)
    result={"population":len(g.output),"baseline":"certified frozen G0 workflow",
            "G0":metrics(base["point"]-truth),"candidate_exit":status["exit"]}
    if status["exit"]==0:
        candidate=measured(path,g)
        errors=candidate["point"]-truth
        result["G2b"]=metrics(errors)
        result["pair_absolute_delta_mean"]=(abs(errors)-abs(base["point"]-truth)).mean(axis=0).tolist()
    dump(root/"workflow-truth.json",result)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root",type=lambda s:Path(s).resolve())
    p.add_argument("phase",choices=("prepare","summarize","workflow"))
    a=p.parse_args()
    {"prepare":prepare,"summarize":summarize,"workflow":workflow_truth}[a.phase](a.root)
