# Joint component runtime evidence

The numerical baseline is `6f30510c`. These records distinguish additional local
certification from unchanged global certification, and frozen numerical parity
from fresh Map/Model parity.

- `local-summary.json`, `local-repeat.json`: 192 component-local assessments,
  128 regular certificates; 1,314 repeated scientific records agree exactly.
- `isolated-local.json`: eight baseline components in active/weak/zero/duplicate
  composites, both precisions, rerun without any historical fit files present.
- `frozen-monolithic-regression.json`: all 216 legacy/guarded/guarded-log searches;
  432 fit/audit records unchanged, regular counts 38/40/40.
- `physical-summary.json`, `physical-repeat.json`: two physically separated
  12-atom blocks, double/float32 observations and four starts. All eight paired
  branches and sixteen local certificates are regular. All 228 repeated
  scientific records agree.
- `map-summary.json`, `map-repeat.json`, `final-api-parity.json`: ordinary CIF/MRC
  loading, fresh first-stage initialization and four starts. All four paired
  branches and eight local certificates are regular. All 114 scientific records
  agree across the two fresh runs and the final public API verification.
- `geometry.json`: generated geometry and the actual float32 MRC header geometry;
  the latter supplies the loaded-map problem and its monolithic reference.
- `runtime-engineering.json`: focused C++/Python tests and an installed-library
  consumer built with testing disabled. No offline regular certificate is
  claimed by the runtime.

Compressed archives preserve the complete scientific records for the first run,
with source/executable provenance. Repeat summaries preserve the exact comparison
result; elapsed seconds and process peak RSS are excluded from scientific equality.
The map archive also contains the physical input CIF/MRC and both file-entry runs.

`local-records.tar.gz` contains the original local audits and isolated reruns;
`frozen-regression-records.tar.gz` contains the newly executed 216-case regression;
`physical-records.tar.gz` and `map-records.tar.gz` contain fresh-input experiments.
Initial/source identities are retained in these records. Historical input archives
under `../joint-abc-components/` remain the independent, pre-extraction reference.

`fresh-input-costs.csv` separates initialization, search, operational assessment,
assembly, reported offline audit time and boundary-audit time. Search reference
time is nested within search time. These are observed run costs under concurrent
offline validation, not isolated throughput benchmarks.

The final full matrix adds:

- `component-summary.json`: 128 pairs, all 64 required regular pairs pass;
  original/global counts remain 40 and composite/global counts remain 24.
- `frozen-component-parity.json`: 1,248 search plus 2,096 existing audit
  scientific records agree with the pre-extraction archives. Process completion
  manifests are excluded because the old audits used one process per case and
  the new full runner uses one process per dataset, resetting caches per case.
- `component-repeat.json`: all 5,008 scientific records agree between fresh runs.
- `local-core-parity.json`: all 1,168 local scientific records agree before/after
  extraction.
- `scope-comparison.json`: distinct global-restricted/local qualification counts;
  old certificates remain intact.
- `standalone-local-bundle.json`: a one-file rerun, input SHA-256 binding, and
  rejection of mismatched parent observations or invalid widths.
- `component-search-records.tar.gz` and `component-audit-records.tar.gz`: extract
  both into the same directory to reconstruct the first complete full-matrix run.
  The second run's provenance and costs are retained separately alongside the
  exact repeat comparison.
- `standalone-bundle-records.tar.gz`, `final-api-records.tar.gz` and
  `engineering-logs.tar.gz`: isolated inputs/results, final API checks, and
  focused test/build logs.

`archive-hashes.json` records the compressed artifact sizes and SHA-256 digests.
