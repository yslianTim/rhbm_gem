"""Run the native CLI on a disposable copy of a frozen baseline database."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from geometry import dump


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", type=lambda p: Path(p).resolve())
    p.add_argument("binary", type=lambda p: Path(p).resolve())
    p.add_argument("tag")
    a=p.parse_args()
    source=a.root/"snapshot"
    d=json.loads((source/"SHA256.json").read_text())["datasets"][a.tag]
    out=a.root/"runs"/f"{a.tag}-archive"
    out.mkdir(parents=True,exist_ok=True)
    # Verify the archived SQLite by the frozen manifest outside this command.
    # The CLI is given a disposable file even on platforms that open SQLite read-write.
    db=out/"disposable.sqlite"
    shutil.copyfile(source/d["baseline"],db)
    argv=[str(a.binary),"umap_embedding","-d",str(db),"-o",str(out),"-k",a.tag,"-j","4"]
    try:
        with (out/"umap.log").open("w") as file:
            r=subprocess.run(argv,stdout=file,stderr=subprocess.STDOUT,check=False)
        dump(out/"umap-command.json",{"argv":argv,"exit":r.returncode})
    finally:
        db.unlink()
    raise SystemExit(r.returncode)
