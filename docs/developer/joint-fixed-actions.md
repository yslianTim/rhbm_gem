# Fixed-state factor and normal actions

> Current fixed-state measurements use [`joint-benchmark.md`](joint-benchmark.md).
> Historical acceptance results and provenance are indexed in
> [joint-component-evidence.md](joint-component-evidence.md).

`SearchMethod::LegacyCompact` remains the production default. Fixed-action
optimization affects the internal operator path, not the model, active-set A/C
solver, independent reference, trust-region policy, or uncertainty calculation.

## Rank contract

`EvaluateRank` accepts a `RankRequest`: original `RankPolicy`, decision-column
count, absolute override, and boundary rule. Decision dimensions can exceed a
component compact's dimensions. Operator and tiled derivative preparation and
endpoint X/W/normalized-W/J evidence share this entry point. Normalization and
parent residual scaling remain at the existing callers.

`SvdNative` preserves Eigen's cutoff (including its positive-minimum guard);
`StrictGreater` preserves the explicit `sigma > threshold` spectrum counts.
These intentionally differ at equality. The BDC dispatch, Jacobi retry near the
boundary, solution and right-vector requests are unchanged. Valid deficient
rank is evidence, while decomposition failure is unavailable evidence. Neither
can qualify full rank. No iterative spectrum approximation is introduced.

## Fixed factors and normal action

Only operator-owned SPQR factors use the public one-shot R/E/H/HPinv/HTau API.
The exported sparse R and Householder data belong to that immutable factor;
compact extraction restores the original column order from R. Primary A/C
workspaces still use expert symbolic/numeric factorization and generation checks.
Independent reference factors are not shared with either path.

For normalized free design Z, projection P=ZZ+, raw width derivative D and
owner-form residual contraction T, the operator computes

`J'J v = (D'(I-P)D v + T'(Z'Z)^-1 T v) / s^2`.

`NormalSolve` uses two triangular solves with the factor's column permutation.
No normal matrix or inverse is formed. The two derivative terms lie in
orthogonal subspaces in exact arithmetic; nonzero residual correction is kept.
Finite-precision parity, symmetry, curvature and weak directions are tested.
PCG now uses this action. Its tolerance, refresh policy, damping, metric,
predicted reduction from `Apply(step)`, and failure behavior are unchanged.

## Cost interpretation

Preparation includes design preparation, factor acquisition and rank work.
Rank work includes compact extraction and SVD. SVD retry time includes its
fallback work. Action timers and factor-operation timers are nested; do not
sum them. One-shot SPQR factorization has one combined timer, not guessed
symbolic/numeric components. Eigen's opaque least-squares solve is timed as
`eigen_least_squares_seconds`; its Q and triangular call counts are recorded,
but its combined time is not assigned to the Q-only timer.

R nonzeros are not full factor memory. `exported_factor_bytes` describes only
the exported R/H/permutation arrays, excluding the retained sparse design,
allocator overhead and factorization scratch; unavailable values are null.
Process-tree RSS and OS peaks remain separate measurements. Dense probes are
known allocations, not an allocator trace. A p-by-p compact still exists.

## Current measurement

Use the unified benchmark profile for fixed-action timing:

```sh
python3 tests/integration/joint_benchmark.py \
  --profile fixed --case chain-8 --fixed-action normal \
  --preconditioner schwarz --build-dir build/joint-eigen \
  --output build/fixed-normal-schwarz.json
```

Historical measured results and their limits are summarized in the
[canonical evidence index](joint-component-evidence.md).
