"""Replay the native UMAP feature transforms on the recorded observations."""

import argparse
import csv
import json
from pathlib import Path
import re
import subprocess

import numpy as np

from compare import describe, dump


def inputs(run, csv_path):
    direct = run / "features.tsv"
    if direct.is_file():
        features = np.loadtxt(direct)
        embed = np.loadtxt(run / "umap.tsv")
        points = {int(row[0]): row[1:3] for row in embed}
        kept = [[*row, *points[int(row[0])], 0] for row in features if int(row[0]) in points]
        return np.asarray(kept), len(features), (~np.isfinite(features[:, 1:])).sum(axis=0).tolist()
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        columns = ("serial id", *reader.fieldnames[3:10], "umap x", "umap y", "bin number")
        rows = np.asarray([[float(row[k]) if row[k] else np.nan for k in columns] for row in reader])
    return rows, len(rows), None


def audit(root, binary, tag, arm):
    run = root / "runs" / f"{tag}-{arm}"
    ds = json.loads((root / "snapshot/SHA256.json").read_text())["datasets"][tag]
    csv_path = run / f"umap_embedding_{tag}.csv"
    if tag == "9evx-unsharp" and arm == "archive" and not csv_path.exists():
        csv_path = root / "snapshot/baselines/umap-b1" / csv_path.name
        run.mkdir(exist_ok=True)
    rows, input_count, original_missing = inputs(run, csv_path)
    audit_input = run / "audited-input.tsv"
    np.savetxt(audit_input, rows, fmt=["%d", *["%.17g"] * 9, "%d"])
    command = [str(binary.resolve()), "umap-audit", str(root / "snapshot" / ds["model"]),
               str(audit_input), str(run)]
    log = run / "audited-features.log"
    with log.open("w") as out:
        subprocess.run(command, stdout=out, stderr=subprocess.STDOUT, check=True)
    summary = json.loads((run / "audited-features.json").read_text())
    matrices = {stage: np.loadtxt(run / f"audited-{stage}-features.tsv") for stage in ("first", "second")}
    for stage, standardized in matrices.items():
        np.testing.assert_array_equal(standardized[:, 0], rows[:, 0])
        assert np.isfinite(standardized).all()
        summary[f"{stage}_quantiles"] = describe(standardized[:, 1:])
    messages = log.read_text()
    imputed = dict.fromkeys(summary["second"], 0)
    for name, count in re.findall(r"Second UMAP feature '([^']+)' filled (\d+) missing values", messages):
        imputed[name] = int(count)
    for name in re.findall(r"Second UMAP feature '([^']+)' has no finite values", messages):
        imputed[name] = len(rows)
    amplitude_column = summary["second"].index("amplitude joint") + 1
    square = matrices["second"][:, amplitude_column] ** 2
    total_variance = square.sum()
    summary.update({"input_count": input_count, "rows": len(rows), "preselection_raw_missing": original_missing,
                    "retained_raw_missing": (~np.isfinite(rows[:, 1:8])).sum(axis=0).tolist(),
                    "second_imputed": imputed, "argv": command,
                    "variance_share": [[int(rows[i, 0]), float(square[i] / total_variance)]
                                       for i in np.argsort(square)[-5:][::-1]] if total_variance else []})
    dump(run / "audited-features.json", summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=lambda p: Path(p).resolve())
    parser.add_argument("binary", type=Path)
    args = parser.parse_args()
    audit(args.root, args.binary, "9evx-unsharp", "pilot-g0")
    audit(args.root, args.binary, "9evx-unsharp", "pilot-g2")
    for tag in ("9evx-unsharp", "9evx-sharp", "6cvm", "6z6u-exp", "6z6u-sim"):
        audit(args.root, args.binary, tag, "archive")
        if tag != "9evx-unsharp":
            audit(args.root, args.binary, tag, "g2")
