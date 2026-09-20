# Immutable views and tiled backend acceptance

Baseline: `085e94b80301556a4ee4402ad8ea98608c2ab97a`. Ownership-only change:
`cab0b53f`; tiled derivative/assessment: `3c08095d`; compact LM and backend
parity: `80b4e527`. The public API, objective, Guarded
acceptance, rank policy, budgets and offline certificate meaning are unchanged.
The [compact evidence](joint-component-tiled-backend.json) records all benchmark
samples, applicable error maxima and per-case scope availability.

## Ownership and numerical validation

The ownership-only stage passed all nine joint CTests, including the original
exact search decisions/counts and eight default plus three extended fixtures.
Tests verify observation pointer identity, support storage identity, shared
partition mappings across 128 components, retained views after problem/input
handles are destroyed, immutable input behavior, and existing Map/CIF/MRC
observation/support contracts. No immutable fixture package/hash changed.

The final tiled code passed 22 default CTests, three offline CTests and the
extended CTest. Testing-disabled Release build/install and the installed consumer
passed; the installed library contains no dense-reference or offline-audit
symbols. Synthetic tests cover nonzero residual correction, 17-row tiles,
8191/8192/8193-row boundaries, compact/dense LM steps and predicted decrease,
and preservation of the full residual norm including its orthogonal tail.
Instrumentation bounds generated rows by the tile size and reduction rows by
tile size plus free-design columns.

Eleven frozen representative cases supply two fixed-state comparisons each:
historical and actual endpoints. Eighteen scopes pass all applicable checks,
two zero-signal scopes have unavailable local corrections, and two duplicate
scopes have no trusted state. These are **not 22 full certificates**. At available
identical states, ranks match; the largest observed discrepancies are:

| Quantity | Observed maximum | Fixed tolerance |
| --- | ---: | ---: |
| Projected derivative, relative Frobenius norm | 4.59e-16 | 1e-8 |
| Full Jacobian, relative Frobenius norm | 4.59e-16 | 1e-8 |
| Spectrum, normalized by largest singular value | 7.01e-11 | 1e-10 |
| Local correction, scaled component difference | 1.53e-16 | 1e-10 |

Converged identifiable endpoints retain scaled 1e-10 parameter comparison and
1e-12 normalized objective comparison. Nonconverged endpoints may differ but
must not worsen normalized objective by more than 1e-12. Guarded replay,
availability, runtime checks and budgets are still checked. No numerical
acceptance tolerance was changed during implementation.

Active-a has a different nonconverged trajectory/endpoint. Zero-signal changes
active face (historical free rank 24, actual free rank 23), while width
identification still fails. A free-face rank cannot be compared between those
different faces: the backend checks ranks independently at each fixed state.
An initial comparator incorrectly required cross-face endpoint rank equality;
it was corrected to the agreed fixed-state contract, with a regression test.
The failed check was rerun successfully. Full-rank local-correction parity at both endpoints
remains unavailable. Weak-1e-4 retains its offline derivative limitation, active-a
retains its failed local-correction/derivative limitations, and duplicate still
has no usable state. Optional two-step and full offline certificate controls pass.

## Cost measurement

Same macOS arm64/AppleClang Release settings, single-thread environment, three
fresh processes per version, no concurrent builds/tests. Both executables use the
same public-API-only benchmark driver. This measures loading the frozen input,
constructing `JointProblem`, and one public fit; it excludes dense test oracles
and does not measure allocation of a full `MapObject`.

Heterogeneous-168/first-stage-float32 is one structural component with 139,551
rows and 168 atoms. One full dense N-by-m double matrix is 178.87 MiB.

| Median metric | Baseline | Tiled |
| --- | ---: | ---: |
| Construction | 0.0432 s | 0.0391 s |
| Search | 22.7939 s | 21.6471 s |
| Assessment | 5.4699 s | 5.5868 s |
| Assembly | 0.2579 s | 0.2600 s |
| Total | 28.7429 s | 27.6702 s |
| Process peak RSS | 1,252.14 MiB | 359.14 MiB |

Peak RSS decreases **71.3%**. Total time decreases 3.7%; assessment is slightly
slower. The improvement claimed here is removal of the full dense derivative/LM
storage, not a guaranteed speedup or cross-platform performance result. Sparse
design/factorization storage and quadratic compact factors remain. Parent-global
assessment still scales with the total number of atoms, even with many small
components.

## Reproduction

Use the existing runtime guide to build with extended tests and offline audits.
Run default, extended and offline CTest lanes; no historical matrices are needed.
The normal fixture runner uses backend numerical parity. `--strict-history`
preserves exact-trajectory comparison for the baseline executable and was also
rerun against the independently built baseline.

The regression runner extracts the hashed fixture into its work directory. With
that dataset directory, run the public benchmark three times automatically:

```sh
python3 tests/integration/joint_component_runtime.py benchmark \
  --executable build/joint-cleanup/bin/joint_component_benchmark \
  --input build/joint-cleanup/joint-extended/fixtures/heterogeneous-168 \
  --output build/joint-cleanup/tiled-benchmark.json
```

For baseline timing, build the same `joint_component_benchmark.cpp` public-only
driver against `085e94b8` in a separate build and library directory. Keeping only
a copy of the executable is insufficient when its runtime search path points to
a library that will be rebuilt. Run baseline and candidate separately. The JSON
contains all three samples and medians, rather than a single favorable run.
