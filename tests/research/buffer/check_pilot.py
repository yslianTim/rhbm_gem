"""Enforce the fresh empty-buffer regression before cross-dataset fitting."""

import json
from pathlib import Path
import sys
import numpy as np
from compare import baseline, read_parameters, dump

root=Path(sys.argv[1]).resolve()
runs=root/"runs"
s0=runs/"9evx-unsharp-pilot-g0"
g2=runs/"9evx-unsharp-pilot-g2"
sources=json.loads((root/"snapshot/SHA256.json").read_text())["datasets"]
old,_=baseline(root/"snapshot",sources["9evx-unsharp"])
sums=[json.loads((p/"summary.json").read_text()) for p in (s0,g2)]
exits=[json.loads((p/"exit.json").read_text())["exit"] for p in (s0,g2)]
assert exits==[0,0],exits
for s in sums:
    checks={e["name"]:e["status"] for e in s["evidence"]}
    assert all(checks[k]=="passed" for k in ("width-stationarity","endpoint-verification","independent-replay"))
p=read_parameters(s0)
p2=read_parameters(g2)
ids=[int(old["atom_ids"][i]) for i in old["fitting_to_original_contributors"]]
assert list(p)==[i for i,v in p2.items() if v["output"]=="1"]==ids
v=old["assembled_state"]
np.testing.assert_array_equal([[float(p[i][k]) for k in ("A","tau","C")] for i in ids],
                              np.column_stack((v["ac"][::2],v["tau"],v["ac"][1::2])))
np.testing.assert_array_equal([float(p[i]["init_tau"]) for i in ids],old["initialization"]["tau"])
assert sums[0]["objective"]==old["objective"]
assert all(s["map_divisor"]==old["metadata"]["map_normalization"]["divisor"] and
           s["observation_scale"]==old["observation_scale"] for s in sums)
seed_path=g2/"initialization.json"
check_seed=seed_path.exists()
if check_seed:
    fresh={v["id"]:v for v in json.loads(seed_path.read_text())}
    for seed in old["initialization"]["atoms"]:
        new=fresh[seed["id"]]
        assert new["initial_mdpde"]==seed["mdpde"] and new["seed_source"]==seed["seed_source"],seed["id"]
        assert new["used_tau"]==seed["used_tau"],seed["id"]
assert sums[1]["buffer"]==[83290] and sums[1]["post_joint"]["2770"]["available"]
g=np.genfromtxt(g2/"focal-grid.tsv",names=True,delimiter="\t")
assert int(((g["squared_distance"]<=.25)&(g["fitting"]==1)).sum())==4
previous=root/"e2-matrix.json"
if previous.is_file():
    frozen_e2=json.loads(previous.read_text())["b1-fixed"]["points"]["2770"]
    np.testing.assert_allclose([float(p2[2770][x]) for x in ("A","tau","C")],
                               [float(frozen_e2[x]) for x in ("A","tau","C")],rtol=1e-6,atol=0)
dump(root/"pilot-regression.json",{"verified":True,"g0_exact":True,"initialization_mdpde_checked":check_seed,
                                   "buffer":sums[1]["buffer"],"output_count":len(ids),
                                   "near_recovered":4,"deduped_same_identity":["G1-focal","G2a","G2b"]})
