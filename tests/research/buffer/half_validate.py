"""Train the native solver from each half map, then score the other half."""

import argparse
import json
from pathlib import Path

import numpy as np

from compare import Geometry, compare, dump, endpoint, map_header, measured, read_map, read_parameters
from run import fit


HALF_DATA = {"9evx-unsharp":"50019","6cvm":"7770","6z6u-exp":"11103"}


def training_evidence(root, tag):
    checks = ("conditional-ac-optimality", "width-stationarity", "endpoint-verification",
              "independent-replay", "global-width-calibration", "global-width-lower-bound",
              "active-width-lower-bound", "active-width-upper-bound")
    fields = ("output_targets", "fitting_targets", "buffer", "map_divisor", "observation_scale",
              "runtime", "target_runtime", "construction_estimated_peak_bytes",
              "peak_rss_platform_units", "fitting_incidence", "search_seconds")
    result = {}
    for arm in ("g0", "g2"):
        result[arm] = {}
        endpoints = []
        for fold in ("h1", "h2"):
            run = root / "runs" / f"{tag}-{fold}-{arm}"
            attempt = json.loads((run / "exit.json").read_text())
            metadata = run / "summary.json"
            native = json.loads(metadata.read_text()) if metadata.is_file() else {}
            result[arm][fold] = {
                **attempt, **{key: native.get(key) for key in fields},
                "checks": {e["name"]: e for e in native.get("evidence", ()) if e["name"] in checks},
                "components": [{k: c[k] for k in ("id", "sweeps", "evaluations", "stop", "state")}
                               for c in native.get("components", ())],
                "run_log": str(run / "run.log"),
            }
            if attempt["exit"] == 0:
                points = read_parameters(run)
                result[arm][fold]["qualified_output_count"] = sum(p["output"] == "1" for p in points.values())
                endpoints.append(points)
        if len(endpoints) == 2:
            first, second = endpoints
            assert list(first) == list(second)
            output = [i for i, p in first.items() if p["output"] == "1"]
            delta = endpoint(second[i] for i in output) - endpoint(first[i] for i in output)
            result[arm]["paired_half_output_abs_delta_quantiles"] = np.quantile(
                abs(delta), [.5, .95, .99, 1], axis=0).tolist()
            result[arm]["paired_half_buffer_abs_delta"] = {
                i: abs(endpoint([second[i]]) - endpoint([first[i]]))[0].tolist()
                for i, p in first.items() if p["output"] == "0"
            }
    dump(root / f"half-runtime-{tag}.json", result)


def exchange(root, binary, tag):
    metadata=json.loads((root/"snapshot/SHA256.json").read_text())
    ds=metadata["datasets"][tag]
    data=root/"snapshot/data"
    accession=HALF_DATA[tag]
    half=[data/f"emd_{accession}_half_map_{i}.map" for i in (1,2)]
    reference=root/"snapshot"/ds["map"]
    headers=[map_header(p) for p in (reference, *half)]
    fields=("shape","mode","start","sampling","cell","axes","origin")
    assert all(all(h[f]==headers[0][f] for f in fields) for h in headers[1:])
    g=Geometry(root/"geometry"/tag)
    # Selection depends only on the same model positions and the checked grid.
    selection=json.loads((root/"selection.json").read_text())
    cand=selection[f"{tag}:population"]["G2b"]["1"]["buffer"]
    result={}
    for fold in range(2):
        prefix=f"h{fold+1}"
        # The runner derives map divisor, parent scale, InitialStage widths, output-only
        # width priors, fitted nuisance parameters and medians from this training path.
        records={label:fit(root,binary,tag,f"{prefix}-{label}",buffer,
                           map_override=half[fold],stage="prediction")
                 for label,buffer in (("g0",[]),("g2",cand))}
        if any(s["exit"] for s in records.values()):
            result[prefix]={"status":"training-failed","arms":records}
            continue
        l,r=[measured(root/"runs"/f"{tag}-{prefix}-{arm}",g) for arm in ("g0","g2")]
        assert l["norm"]==r["norm"]
        evaluation=read_map(half[1-fold],g.voxels)
        result[prefix]={"status":"evaluated","training_map":half[fold].name,
                        "test_map":half[1-fold].name,"training_divisor":l["norm"],
                        "prediction":compare(g,l,r,evaluation,2770 if tag.startswith("9evx")
                                             else 32994 if tag=="6cvm" else 372)}
    dump(root/f"half-validation-{tag}.json",{"dataset":tag,"header":headers,
        "filtering":"frozen matched half maps for unsharpened experimental volume",
        "half_map_digests":{p.name:metadata["files"]['data/'+p.name] for p in half},
        "fixed_structure":True,"buffer":cand,"folds":result})
    training_evidence(root, tag)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root",type=lambda p:Path(p).resolve())
    p.add_argument("binary",type=Path)
    p.add_argument("datasets",choices=tuple(HALF_DATA),nargs="+")
    a=p.parse_args()
    for name in a.datasets:
        exchange(a.root,a.binary,name)
