"""Launch one production-solver research fit at a time with explicit provenance."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from geometry import dump


def fit(root, binary, tag, arm, buffer):
    snapshot = root/"snapshot"
    ds = json.loads((snapshot/"SHA256.json").read_text())["datasets"][tag]
    result = root/"runs"/f"{tag}-{arm}"
    result.mkdir(parents=True,exist_ok=True)
    stage = "full-simulation" if tag == "6z6u-sim" else "full"
    m = snapshot/ds["map"]
    executable = binary.resolve()
    argv = [str(executable),str(snapshot/ds["model"]),str(m),ds["chains"],
            ",".join(map(str,buffer)) or "-", "fixed","output",
            "2770" if tag.startswith("9evx") else ("3242" if tag=="6cvm" else "372"),
            str(result),stage]
    caps = dict(OMP_NUM_THREADS="1",OMP_THREAD_LIMIT="4",OPENBLAS_NUM_THREADS="1",
                VECLIB_MAXIMUM_THREADS="1",MKL_NUM_THREADS="1")
    with m.open("rb") as f:
        digest = hashlib.file_digest(f,"sha256").hexdigest()
    if digest!=ds["map_sha256"]:
        raise ValueError(f"{tag} map hash differs from the frozen source")
    dump(result/"provenance.json", {"argv":argv,"thread_caps":caps,"map_sha256":digest,
        "model_sha256":ds["model_sha256"]})
    started=time.monotonic()
    env=os.environ|caps
    with (result/"run.log").open("w") as log:
        child=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,env=env,check=False)
    record={"exit":child.returncode,"wall_seconds":time.monotonic()-started}
    dump(result/"exit.json",record)
    print(tag,arm,record,flush=True)
    return record


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root",type=lambda p:Path(p).resolve())
    p.add_argument("binary",type=Path)
    p.add_argument("phase",choices=("pilot","matrix"))
    args=p.parse_args()
    chosen=json.loads((args.root/"selection.json").read_text())
    if args.phase=="pilot":
        for arm,buffer in (("pilot-g0",[]),("pilot-g2",chosen["9evx-unsharp:2770"]["G2b"]["1"]["buffer"])):
            fit(args.root,args.binary,"9evx-unsharp",arm,buffer)
    else:
        for tag in ("9evx-sharp","6cvm","6z6u-exp","6z6u-sim"):
            fit(args.root,args.binary,tag,"g2",chosen[f"{tag}:population"]["G2b"]["1"]["buffer"])
