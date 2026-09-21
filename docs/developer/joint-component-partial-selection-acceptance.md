# Joint component partial-selection acceptance

## Accepted implementation

**All four validation lanes passed.** The tested implementation is the uncommitted
working tree based on `bb72080c4878ee2333ad3e084f5e7a9d687abf33`, identified by
production source SHA-256 `e346827a2d7ffe05f4d83bcf06c8b92f416d85825f9d4cf3d8f22a20dee73e34`. No new commit, tag or package release
was created. C++ consumers must rebuild; package version remains 2.0.0.

The [machine-readable receipt](joint-component-partial-selection-acceptance.json)
records changed implementation/test file hashes, the shared production source
fingerprint, per-build configuration/build fingerprints, checks, and all twelve
resource runs. The original [v1 freeze](joint-component-v1-acceptance.md) remains
unchanged. Kernel, 2.5 Angstrom support, objective, numerical thresholds, search
budgets and frozen fixtures are unchanged.

Selected non-hydrogen atoms fix observation rows. Every non-hydrogen contributor
whose support intersects those rows is fitted once, including halo atoms outside
the map. Halo cannot expand rows. Initialization uses the full model copy and
may sample outside the target domain; per-atom exceptions and invalid widths
block only their structural component. Target/halo roles, initialization scope
and reasons are saved in schema-3 outcomes. SQLite remains v17; old joint JSON
is rejected without writes. CSV appends SelectionRole and retains all contributors.

## Validation

| Lane | Result |
| --- | --- |
| Default Release | All 25 CTests passed, including CLI, C++ request, Python request and eight frozen runtime cases. Final joint-group rerun passed after adding full-initializer parity and group/bond preservation coverage. |
| Extended/offline Release | All four CTests passed, retaining derivative/precision checks and the extended heterogeneous-168 case. |
| Testing-disabled install | Library and CLI built with testing/Python/offline/extended disabled; installed consumer compiled, executed and exported schema-3 JSON/CSV with selected-domain metadata. |
| Source package without Git | A source copy outside any Git repository passed the same testing-disabled build/install/consumer check, with the identical production source fingerprint. |

Compact logs: [default](figures/joint-component-partial-selection-acceptance/default-tests.txt),
[final joint group](figures/joint-component-partial-selection-acceptance/joint-tests.txt),
[offline](figures/joint-component-partial-selection-acceptance/offline-tests.txt),
[installed](figures/joint-component-partial-selection-acceptance/installed.txt),
[source package](figures/joint-component-partial-selection-acceptance/source.txt),
[repository guards](figures/joint-component-partial-selection-acceptance/guards.txt).

The host/configuration is recorded per lane in the receipt. Builds reused the
installed dependencies and cached pinned UMAP sources; this is not a vendored or
dependency-free source distribution. Build/test lanes overlapped, so their elapsed
times are diagnostics rather than performance measurements.

### Partial-selection controls

- An independent exhaustive grid/catalogue reference verifies exact rows, atom
  ordering, memberships and squared distances, including all-selected input,
  zero/negative values, nonzero origin, anisotropic spacing, clipping, exact
  cutoff boundaries, outside-map halo and hydrogen exclusion.
- Nonrecursive closure excludes atoms that only neighbor halo, and preserves
  unobserved selected targets. A halo at x=0 bridges targets at x=±2.8 Angstrom;
  it has one parameter vector. Partition and same-domain monolithic results agree.
- A noiseless identifiable target/halo control recovers A/B/C within the existing
  scaled 1e-10 criterion. Fixed initial widths and normalization isolate the joint
  objective from modifications outside its target rows.
- Tail-only, zero-signal, duplicated-support/rank-deficient, insufficient-row,
  unobserved-target and invalid-width cases retain their numerical limitations.
  A NaN outside both target supports but inside one initializer's cubic stencil
  produces a captured local initialization exception; the independent component
  remains usable and no full objective is assembled.
- Full-selection initialization matches the original v1 batch workflow's widths,
  alpha values and sample counts. Partial initialization preserves original atom
  and bond selection, group alpha, second-stage results and halo first-stage history.
- CLI, C++ and Python requests use a real backbone target plus side-chain halo,
  then delete their temporary source inputs before export. Role metadata,
  initialization scope and saved evidence survive. Invalid selection metadata
  rolls back storage; schema-2 reads leave database bytes unchanged. Unknown
  selection in hand-built inputs remains null / not-recorded.

For the fixed-B0 matched/omitted control, complete closure gave normalized
objective 8.139650468e-32 and target A/B/C errors below 3e-15. Omitting the halo
on **the same rows** gave objective 0.07347280397 and target errors
ΔA=1.909360106, ΔB=0.1355602156 Angstrom, ΔC=0.1243884602. This is a synthetic model-compensation
control, not a real-data accuracy claim or a comparison against fitting a larger map.

## Complete-process resource baseline

Each case ran in three independent serial processes **after builds and CTests
finished**. Times include the actual Map/Model builder, model copy, contributor
initialization, search, assessment and assembly. Fixture generation precedes the
timed estimator call but contributes to peak RSS. Additional matched/omitted
controls run after timing/RSS capture. Per-phase costs and all results are in the
receipt; table values are medians.

| Case | Targets / halo | Rows | Memberships | Largest component | Total ms | Peak MiB | Runtime convergence |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all | 2 / 0 | 3,345 | 4,938 | 2 | 22.86 | 24.91 | passed |
| partial | 1 / 1 | 2,469 | 4,062 | 2 | 19.66 | 24.03 | passed |
| bridge | 2 / 1 | 4,874 | 6,086 | 3 | 33.61 | 25.86 | passed |
| weak | 1 / 1 | 2,469 | 2,490 | 2 | 56.43 | 24.12 | failed |

All three runs of all, partial and bridge passed runtime convergence. All three
weak runs retained available states but failed **local-correction**; their
numerical-identifiability check passed at the returned state. This failure is
not relabeled as rank deficiency or successful convergence. A separate narrower
initial-width tail-only stress case produces an unusable zero-column problem.
Runtime offline certificates remain NotRun.

These measurements cover only 2–3-contributor synthetic problems. The builder
still scans the complete eligible catalogue, and per-atom initialization rebuilds
selection in a full model copy. A small selected region therefore does not
establish sublinear catalogue cost or a general memory reduction. No maximum
supported component size, total memory/time ceiling, cross-platform behavior,
real-data accuracy or global optimum is established. Weak halo is never removed,
fixed at zero or assigned a separate convergence criterion.

## Reproduction

Use the exact lane settings and commands in the receipt, adjusting local paths.
The existing v1 build directories can be reused with those settings. Build
`tests_all` and run default CTests; in the offline build run
`ctest -L 'joint:extended|joint:offline' --output-on-failure`. Testing-disabled
and no-Git builds run `lint_install_smoke`. Run `lint_repo` in the checkout.
For a working-tree source package, include new implementation/test files and
exclude `.git` and build artifacts; verify its production fingerprint against
the other lanes. Do not represent an archive of the unchanged base commit as
this implementation.

```sh
python3 tests/integration/joint_partial_selection.py \
  --executable build/joint-v1-default/bin/joint_partial_selection \
  --output build/joint-partial-resources.json
```

The resource runner defaults to three repetitions per case. Existing frozen
numerical fixtures and tolerances are not regenerated.
