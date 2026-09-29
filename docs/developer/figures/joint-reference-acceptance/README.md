# Endpoint-reference acceptance evidence

> Historical artifact bundle. Its campaign runner was retired; current solve
> measurements use [`../../joint-benchmark.md`](../../joint-benchmark.md), and
> reference certification remains owned by the permanent offline tests.

[Acceptance report](../../joint-component-reference-acceptance.md) ·
[Summary](summary.json) · [Verification](verification.json) ·
[Archive hashes and contents](archives.json) · [Campaign log](campaign.txt)

The six archives preserve the complete historical compact campaign and the
incremental endpoint-reference experiment, including raw samples, compact
matrices and RHS, full numerical audits, immutable inputs and initialization,
SQLite databases, JSON/CSV exports, resource receipts, source proofs, build
settings, recorded executable fingerprints, and regression logs.

- [Historical inputs](historical-inputs.tar.gz)
- [Historical audits and captured matrices](historical-audits.tar.gz)
- [Historical results, replay and resource stops](historical-results.tar.gz)
- [Incremental inputs](incremental-inputs.tar.gz)
- [Incremental fixed-state audits](incremental-audits.tar.gz)
- [Incremental commands, initialization, receipts and verification](incremental-results.tar.gz)

Combined compressed size is about 226 MiB; every individual archive is below
100 MiB. Archives contain disjoint paths and must all be extracted into the same
fresh directory. This retains complete data without relying on ignored local
build directories. Historical binaries are represented by their recorded hashes
and source/configuration proofs; current binaries are not substitutes for them.

Current search measurements use the unified solve profile, for example:

```sh
python3 tests/integration/joint_benchmark.py \
  --profile solve --case baseline:first-stage-double \
  --build-dir build/joint-reference-spqr --output build/reference-solve.json
```

Historical archive reaggregation requires the retired source from Git history
(last changed at `7f84f931215f7b5a52948f7ecdbfc15fad79d7c5`) plus Python and NumPy:

```sh
python3 - <<'PY'
import hashlib, json, tarfile
from pathlib import Path
bundle = Path('docs/developer/figures/joint-reference-acceptance')
work = Path('build/joint-reference-relocated')
work.mkdir(parents=True, exist_ok=False)
for name, record in json.loads((bundle/'archives.json').read_text()).items():
    archive = bundle/name
    with archive.open('rb') as stream:
        assert hashlib.file_digest(stream, 'sha256').hexdigest() == record['sha256']
    assert archive.stat().st_size == record['bytes']
    with tarfile.open(archive) as stream:
        stream.extractall(work, filter='data')
PY
python3 tests/integration/joint_reference_validation.py \
  --work-dir build/joint-reference-relocated --report-only
```

Reaggregation verifies the file manifest, recomputes historical compact gates,
checks the corrected normalized objective comparison, verifies completed exports
and initialization, and recomputes all latest numerical and command gates. It
needs neither the original executable paths nor Git and does not run estimators.
The relocated result was checked for exact equality with the saved summary;
see [relocation verification](relocation-verification.txt).
