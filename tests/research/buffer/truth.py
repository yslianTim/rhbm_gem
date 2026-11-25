"""Analyze frozen D4 draws and native JCE paired truth estimates."""

import argparse
import json
from pathlib import Path

import numpy as np

from geometry import dump

WIDTH_BINS = ("≤0.5", "0.5–1.0", "1.0–1.5", "1.5–2.0", ">2.0")


def width_bin(true_width):
    return WIDTH_BINS[np.searchsorted([.5, 1., 1.5, 2.], true_width, side="left")]


def paired_statistics(left, right):
    if not len(left):
        return None
    delta = abs(right) - abs(left)
    boots = np.random.default_rng(92173).integers(len(left), size=(2000, len(left)))
    rms_delta = np.sqrt(np.mean(right[boots]**2, axis=1)) - np.sqrt(np.mean(left[boots]**2, axis=1))
    return {"mean_absolute_error_delta": delta.mean(axis=0).tolist(),
            "median_absolute_error_delta": np.median(delta, axis=0).tolist(),
            "rmse_delta": (np.sqrt(np.mean(right**2,axis=0))-np.sqrt(np.mean(left**2,axis=0))).tolist(),
            "mean_absolute_error_delta_ci": np.quantile(delta[boots].mean(axis=1), [.025,.975], axis=0).tolist(),
            "rmse_delta_ci": np.quantile(rms_delta,[.025,.975],axis=0).tolist()}


def qualification_pairs(arms):
    baseline, candidate = arms["G0"], arms["single"]
    if baseline.keys() != candidate.keys():
        raise ValueError("Qualification comparison must contain both arms of every replicate.")
    outcomes = {k: [] for k in ("n11", "n10", "n01", "n00")}
    for rep in sorted(baseline):
        key = f"n{int(baseline[rep]['qualified'])}{int(candidate[rep]['qualified'])}"
        outcomes[key].append(rep)
    return outcomes


def fitted_points(fit):
    if not fit["qualified"]:
        return {}
    points = fit["points"]
    ids = fit.get("fitting_ids", range(1, len(points)+1))
    return {int(i): np.asarray(point) for i,point in zip(ids,points,strict=True)}


def prediction_trials(fit, recovered):
    if not fit["qualified"]:
        return {}
    scores={"truth":{k:v for k,v in fit.get("prediction_RMSE",{}).items() if v is not None}}
    if recovered:
        scores["truth"]["recovered"]=float(np.sqrt(np.mean(np.asarray(fit["residuals"])[recovered]**2)))
    for key in ("training","independent"):
        if (residual := fit.get(f"oracle_{key}_residuals")) is None:
            continue
        values=np.asarray(residual)
        reduced={"original":float(np.sqrt(np.mean(values**2)))}
        for domain, rows in (fit["oracle_domain_rows"] | {"recovered":recovered}).items():
            if rows:
                reduced[domain]=float(np.sqrt(np.mean(values[rows]**2)))
        scores[key]=reduced
    return scores


def prediction_summary(arms, paired):
    values={name:{r:prediction_trials(f,arms["single"][r].get("recovered_rows"))
        for r,f in fits.items()} for name,fits in arms.items()}
    combined={}
    for key in ("truth","training","independent"):
        domains=sorted({d for arm in values.values() for p in arm.values() for d in p.get(key,{})})
        if not domains:
            continue
        combined[key]={}
        for domain in domains:
            samples={a:{r:p[key][domain] for r,p in by_rep.items() if domain in p.get(key,{})}
                for a,by_rep in values.items()}
            ids=[r for r in paired if all(r in group for group in samples.values())]
            l,r=(np.array([samples[a][i] for i in ids]) for a in ("G0","single"))
            delta=r-l
            ci=None
            if len(delta):
                boot=np.random.default_rng(92173).integers(len(delta),size=(2000,len(delta)))
                ci=np.quantile(delta[boot].mean(axis=1),[.025,.975]).tolist()
            combined[key][domain]={"conditional_median":{a:float(np.median(list(v.values()))) if v else None
                for a,v in samples.items()},"matched_median": [float(np.median(x)) for x in (l,r)] if len(ids) else None,
                "matched_count":len(ids),"mean_delta":float(delta.mean()) if len(delta) else None,"mean_delta_ci":ci}
    return combined


def fit_diagnostics(arms, paired):
    result={}
    for name,fits in arms.items():
        failures={}
        atom_geometry={}
        for rep, fit in fits.items():
            if not fit["qualified"]:
                failures[rep]={"components":[{k:c[k] for k in ("component","stop","replay_reason","limiting_atom",
                    "primary_projected_width_kkt","independent_projected_width_kkt","primary_ac_kkt","independent_ac_kkt")}
                    for c in fit.get("width_diagnostics",())],"diagnostic_error":fit.get("diagnostic_error"),
                    "sweep_tail":fit.get("recent_sweeps")}
            if rep not in paired:
                continue
            for c in fit.get("width_diagnostics",()):
                for atom in c["atoms"]:
                    if "fit_diagnostics" in atom:
                        atom_geometry.setdefault(atom["id"],[]).append(atom["fit_diagnostics"] | {
                            k:float(atom[k]) for k in ("amplitude_lower_face","width_lower_face","width_upper_face")})
        result[name]={"failures":failures,"paired_diag_medians":{a:{k:float(np.median(col)) if col else None
            for k in rows[0] for col in [[r[k] for r in rows if r.get(k) is not None]]} for a,rows in atom_geometry.items()},
            "max_construction_estimate":max((f.get("construction_estimate",0) for f in fits.values()),default=0),
            "profile_evaluations":sum(f.get("evaluations",0) for f in fits.values() if not f.get("reused_identical_G0")),
            "native_solve_wall_seconds":sum(f.get("wall_seconds",0) for f in fits.values() if not f.get("reused_identical_G0"))}
    return result


def stratify_scenario(data):
    arms = {name:{v["replicate"]:v for v in data["fits"] if v["strategy"]==name}
            for name in ("G0","single")}
    transitions = qualification_pairs(arms)
    by_arm = {name:{rep:fitted_points(v) for rep,v in group.items()} for name,group in arms.items()}
    qualified = {name:sum(v["qualified"] for v in group.values()) for name,group in arms.items()}
    counts = np.array([int(arms['single'][rep]['qualified'])-int(v['qualified']) for rep,v in sorted(arms['G0'].items())])
    boot = np.random.default_rng(92173).integers(len(counts), size=(2000,len(counts)))
    result = {"truth":data["truth"],"population":"B" if any(a[1]>2 for a in data["truth"]) else "A",
              "attempted":len(counts),"qualified":qualified,"qualification_delta":float(counts.mean()),
              "qualification_delta_ci":np.quantile(counts[boot].mean(axis=1),[.025,.975]).tolist(),
              "transitions":transitions,"atoms":{},"bins":dict.fromkeys(WIDTH_BINS,0),
              "prediction":prediction_summary(arms,transitions['n11']),"fit":fit_diagnostics(arms,transitions['n11'])}
    for atom, point in enumerate(data["truth"],1):
        group = {name:{rep:p[atom]-point for rep,p in fits.items() if atom in p} for name,fits in by_arm.items()}
        matched = [rep for rep in transitions['n11'] if all(rep in g for g in group.values())]
        paired = [np.array([group[name][rep] for rep in matched]) for name in ("G0","single")]
        role = ("output", "output", "buffer", "external_halo")[atom-1]
        bin_name = width_bin(point[1])
        if role == 'output':
            result['bins'][bin_name] += 1
        result["atoms"][str(atom)] = {"role":role,"width_bin":bin_name,"in_interest":role=="output" and point[1]<=2.,
            "conditional":{k:metrics(np.array(list(v.values()))) for k,v in group.items()},
            "coverage":{k:len(v) for k,v in group.items()},"paired_count":len(matched),
            "paired":dict(zip(("G0","single"),map(metrics,paired),strict=True)),"effect":paired_statistics(*paired)}
    return result


def stratify_matrix(path, destination):
    dump(destination, {data["scenario"]:stratify_scenario(data)
        for data in map(json.loads,path.read_text().splitlines())})


def metrics(e):
    if not len(e): return None
    absolute=np.abs(e)
    return {"bias":e.mean(axis=0).tolist(),"rmse":np.sqrt(np.mean(e**2,axis=0)).tolist(),
            "parameter_sd":e.std(axis=0).tolist(),
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
            synthetic={kind:{k:{} for k in pred} for kind in ("training","independent")}
            for v in valid:
                rep=v["replicate"]
                for key, value in v["prediction_RMSE"].items():
                    if value is not None:
                        pred[key][rep]=value
                outer=buffers["single"][rep].get("recovered_rows")
                if outer:
                    values=np.array(v["residuals"])[outer]
                    pred["recovered"][rep]=float(np.sqrt(np.mean(values**2)))
                for kind, domains in synthetic.items():
                    if f"oracle_{kind}_residuals" not in v:
                        continue
                    values=np.asarray(v[f"oracle_{kind}_residuals"])
                    indices=v["oracle_domain_rows"] | {"original":list(range(len(values))),"recovered":outer}
                    for domain, rows in indices.items():
                        if rows:
                            domains[domain][rep]=float(np.sqrt(np.mean(values[rows]**2)))
            all_predictions, paired_predictions = {}, {}
            for domain, values in pred.items():
                all_predictions[domain] = float(np.median(list(values.values()))) if values else None
                matched = [values[i] for i in pair if i in values]
                paired_predictions[domain] = np.quantile(matched,[.5,.95,.99,1]).tolist() if matched else None
            r[name]={"attempted":len(runs),"qualified":len(valid),
                     "conditional_errors":metrics(errors),"paired_errors":metrics(paired),
                     "prediction_median":all_predictions,"prediction_paired_quantiles":paired_predictions,
                     "bound_activity":np.sum([v.get("bound_activity",[0,0,0]) for v in valid],axis=0).tolist() if valid else None,
                     "profile_evaluations":sum(v["evaluations"] for v in runs.values())}
            for kind, domains in synthetic.items():
                if any(domains.values()):
                    r[name][f"oracle_{kind}_RMSE"]={domain:{"conditional_median":float(np.median(list(values.values()))),
                        "paired_median":float(np.median([values[i] for i in pair if i in values]))
                        if any(i in values for i in pair) else None} for domain,values in domains.items() if values}
        if pair:
            delta=np.array([abs(np.array(buffers["single"][rep]["points"])[:2]-truth)-
                            abs(np.array(buffers["G0"][rep]["points"])[:2]-truth) for rep in pair])
            rng=np.random.default_rng(2907)
            boots=rng.choice(len(pair),size=(2000,len(pair)),replace=True)
            r["paired_absolute_error_delta"]={"mean":delta.mean(axis=0).tolist(),
                "paired_replicate_bootstrap_95":np.quantile(delta[boots].mean(axis=1),[.025,.975],axis=0).tolist()}
        r["qualified_pairs"]=len(pair)
        r["qualification_transitions"]={f"n{a}{b}":[i for i in buffers["G0"]
            if (buffers["G0"][i]["qualified"],buffers["single"][i]["qualified"])==(bool(a),bool(b))]
            for a,b in ((1,1),(1,0),(0,1),(0,0))}
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
    p.add_argument("phase",choices=("summarize","workflow"))
    a=p.parse_args()
    {"summarize":summarize,"workflow":workflow_truth}[a.phase](a.root)
