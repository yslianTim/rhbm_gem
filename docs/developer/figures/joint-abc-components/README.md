# Exact-component experiment records

The [report](../../joint-abc-components-experiment.md) and
[version 1 contract](../../joint_abc_components_contract.md) define the claims and
acceptance rules. The compact tables use the final endpoint-scoped certificates.
Initial-direction audits remain separately labelled as superseded evidence.

## Reviewable results

- `scientific-validation.json`: sequential gates and final repeat/rerun results.
- `monolithic-regression.json`: 216 unchanged historical branches and certificates.
- `original-census.json`, `census.csv`: structural connectivity and component sizes.
- `same-state-gate.json`, `same-state-parity.csv`: frozen initial/Guarded endpoint comparisons.
- `same-state-at-fit-endpoints.csv`, `endpoint-parity.csv`: both independently searched endpoints and paired parameter/prediction gates.
- `failure-matrix.csv`: search completion, usable states and certification limitations.
- `costs-a.csv`, `costs-b.csv`, `search-cost-totals-a.csv`, `cost-summary.json`: per-scope costs and non-overlapping totals.
- `context-setup-costs-a.csv`, `context-setup-costs-b.csv`: context construction, excluded from search budgets.
- `repeat-comparison.json`, `isolated-reruns.json`: exact scientific reproducibility, including failures.
- `endpoint-audit-correction.json`: why initial directions cannot certify an endpoint's weakest mode.
- `audit-revision.json`: unchanged endpoint states and removed unsupported qualifications.
- `engineering-validation.json`, `final-cpp-controls.json`, `shared-solver-controls.json`: tests, control evidence and production linkage checks.
- `build-environment.json`, `reporting-provenance.json`: compiler/library versions, numeric thread settings and the final reporting source hashes.
- `artifact-index.json`: SHA-256 and byte size of each delivered artifact.

The same-state pass flag means all applicable comparisons pass. It does not turn
an unavailable profile, changing active face or rank-deficient correction into
complete equivalence. Component certificates are supplemental and distinct from
the global assembled certificate.

## Complete records and extraction

Extract the following archives into one directory. Each run archive contains its
own `formal-a/` or `formal-b/` prefix; the archives therefore combine without
renaming files or changing recipe paths.

| Archive | Contents |
|---|---|
| `prerequisite-records.tar.gz` | Stage-one 216-branch regression, failed exploratory stage-two records and validated same-state gate |
| `search-records-a.tar.gz`, `search-records-b.tar.gz` | Frozen observations/contributor CSR, original source fits/certificates, composition recipes, all searches/trials, actual assemblies, endpoint same-state controls, contexts, original source snapshots and provenance |
| `superseded-audits-a.tar.gz`, `superseded-audits-b.tar.gz` | Complete initial-direction audits; excluded from final qualification claims |
| `endpoint-audits-a.tar.gz`, `endpoint-audits-b.tar.gz` | Independently recomputed final audits, scoped fits, global/subset directions, precision references, boundary scans, certificates, summaries and audit source/provenance |
| `isolated-reruns.tar.gz` | Standalone component fits, frozen context bundles, fresh audits and exact comparisons |
| `engineering-logs.tar.gz` | Build, test and lint logs, plus final reporting/test source snapshots |

Every run contains all nine original and seven disjoint-composition datasets.
Its `inputs/` directory is immutable and checked against `input-hashes.json`.
Full basis/Jacobian/prediction arrays can be reconstructed from these observations,
support distances, coefficients, saved directions and implementation sources;
comparisons retain their numerical discrepancies and complete spectral evidence.

The contract documents `run`, `audit`, `summarize`, `compare` and
`rerun-component`. The prior certification's
[record archive](../joint-abc-certification/scientific-records.tar.gz) supplies the
`FROZEN_CERTIFICATION` input for a fresh run. Each run's final audit tree is
`formal-*/final-audits`; select it explicitly with `--audits`, or with
`--left-audits` and `--right-audits` for comparison. Original search source and
final audit source are preserved separately. The audit wrapper changed; the
numerical search kernels did not. All re-audits start from empty case-local
caches and independently recalculate boundary scans.

Resource timings are descriptive measurements from overlapping runs on one
host. Search reference time is nested inside search time; precision time is
nested inside audit time. Assembled search counters summarize child searches,
whose time is reported on the individual component rows. Do not add that
aggregate counter to its children. Process peak RSS is a high-water mark and is
not a component allocation measurement or a measured concurrent-memory total.
