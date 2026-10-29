# Joint component estimator runtime

This guide describes the production target-only Pure Halo-Support Exclusion
estimator. The persisted parameterization is
`target-only-pure-halo-support-exclusion-v1`.

## Production flow

`RunJointComponentWorkflow` builds the Pure Exclusion problem and obtains
qualified FirstStage target widths. A width calibration is available only with
at least 20 qualified target donors and 95% target coverage. The sole global
width bound is `Bmax = 4 × Q0.95(qualified FirstStage target widths)`.
Nonfinite or nonpositive donor widths, median-fallback donors, initialization
failures, and overflowed bounds do not produce a usable calibration. Initial
widths are projected into this finite bound; unavailable calibration prevents
Joint search from starting.

Targets are the only free Joint components. Non-target Halo support is used
only for the Pure Halo-Support Exclusion row selection, with the production
support radius of 2.5 Å. The model is
`f_i(r) = A_i G(r; B_i) + C_i K(r; B_i)`, with `A_i ≥ 0`,
`0 < B_i ≤ Bmax`, and signed `C_i`. C multiplies the charge kernel; it is not
a constant background offset.

Production uses FixedNeighbor blocks with core size 12, forward serial
Gauss-Seidel order, and at most 200 sweeps. A block visit may accept at most
one trusted update. At fixed widths, the local A/C problem uses column
normalization, Eigen SparseQR, a nonnegative-amplitude active set, strict
`1e-10` local KKT certification, and the existing local QR iterative
refinement when needed. C remains signed. Width updates use bounded log-width
Active-set LM with a bounded trust-region subproblem, damping adjustment, and
the weak-width safeguard. The full profile derivative includes its
nonzero-residual correction and cancellation protection. Compact SVD supports
rank and weak-direction checks; it is not an alternate A/C solver.

## Qualification

The global search preserves objective monotonicity and practical-v2
convergence: projected A/C KKT ≤ `1e-7` and bound-aware projected Width KKT ≤
`1e-9`. Confirmation requires at least two completed sweeps and no block
failure on the confirming sweep. Global feasibility, the full sweep replay
audit, recovery, and candidate rollback remain part of search. Endpoint
Certification and Independent Replay verify the accepted state. Component
assembly consistency, Local Conditioning, Parameter Stability, target point
quality, and resource preflight remain explicit evidence; local LM termination
alone does not qualify a component. Failed initialization, unavailable width
calibration, incomplete search, unqualified endpoint, missing component state,
and nullable diagnostics remain representable saved outcomes.

Local Conditioning assesses each block's numerical conditioning. Perturbation-
based Parameter Stability deliberately perturbs local parameters, re-estimates
the block, and records parameter, prediction, and objective sensitivity.
These point-quality assessments do not use sweep-to-sweep parameter movement
and do not replace global convergence evidence. Target point validation
records component convergence, parameter validity, conditioning, stability,
point quality, and reasons. `Unavailable` and `NotAssessed` remain valid states
when an assessment cannot be completed.

## Halo reconstruction and downstream work

Qualified Joint Target points are the donors for geometry-local Halo Median
Reconstruction. Each Halo uses up to its four nearest qualified Target points
by Euclidean distance in Å; equal-distance ordering uses Atom ID. A result
records the donor IDs and the geometry-local kNN-4 contract. No retained Halo
or unavailable donors is a valid result state.

The reconstructed Halo points support Original-domain PostFitPeeling. The
workflow then runs Post-Joint Local MDPDE for qualified selected targets and
Group Fitting from those local results. These downstream stages retain their
own availability and qualification evidence.

Memory preflight and memory-bounded evaluation remain part of problem
construction, FixedNeighbor search, endpoint verification, local conditioning,
parameter stability, Halo reconstruction, and peeling.

## Persistence contract

| Payload | Accepted format |
| --- | --- |
| SQLite database | Schema v20 |
| Stage Result JSON | Version 4 |
| Joint Result JSON | Schema v21 |

SQLite opens a new empty database by initializing schema v20 and opens an
existing database only when its exact schema is v20. Stage JSON accepts only
version 4. Joint JSON accepts only schema 21 and the `strict-v2` or
`practical-v2` convergence identifiers. It stores target-only Pure Exclusion
results. No format performs migrations or compatibility decoding. Unsupported
inputs are rejected without repair or overwrite.

Current documents can store qualified and unqualified results, nullable
diagnostics, missing component states, unavailable width calibration, and
unavailable Halo reconstruction. The writer does not select a schema from
result contents.

The [data I/O architecture](architecture/dataobject-io-architecture.md)
documents the database and JSON boundaries. The
[stage integration guide](joint-stage-integration.md) describes downstream
ownership and invalidation.
