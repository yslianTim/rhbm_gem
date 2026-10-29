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
failure on the current sweep. Maximum log-width change remains a diagnostic and
does not gate certification. Endpoint certification and Independent Replay
verify the accepted state. Component assembly consistency, Local Conditioning,
Parameter Stability, target point quality, and resource preflight remain
explicit evidence; local LM termination alone does not qualify a component.
Failed initialization, unavailable width calibration, incomplete search,
unqualified endpoint, missing component state, and nullable diagnostics remain
representable saved outcomes.

Local Conditioning and Parameter Stability describe the accepted solution's
block behavior. Target point validation records component convergence,
parameter validity, conditioning, stability, point quality, and reasons.
These diagnostics retain `Unavailable` and `NotAssessed` states where work
could not be assessed; they do not replace global convergence evidence.

### practical-v2 evaluation

This evaluation supports practical-v2 as the production default. The
measurements below cover numerical qualification, parameter accuracy, and
resource use.

The practical-v2 policy was compared with a frozen practical-v1 Release build
using one thread and five measured runs after one warmup. On the 3, 30, 300,
and 1,000 target synthetic chains, both versions qualified the endpoint and
passed Independent Replay. KKT-only saved one sweep on 3, 300, and 1,000
targets; median wall time changed by -21.5%, +1.2%, -16.2%, and -16.2%,
respectively. Peak RSS changed by less than 1.2% in either direction. At 300
and 1,000 targets, maximum paired per-target differences were `3.81e-6` for A,
`3.16e-7` for C, and `2.62e-7` for B; the 30-target result was bitwise
identical. The largest objective increase was `1.89e-14`.

Signed cancellation now certifies after 19 sweeps where practical-v1 exhausted
its 500-sweep budget, with log-width change `2.13e-6`. The noise-only fixture
reaches a certified KKT endpoint after two sweeps while parameter-stability
evidence remains `insufficient-history` and `NotAssessed`. This is numerical
stationarity evidence, not evidence that parameters are identifiable or
stable. Rank-deficient and near-singular fixtures still exhaust recovery
without endpoint qualification. No real cryo-EM map workload was available.

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
| Joint Result JSON | Schema v20 |

SQLite opens a new empty database by initializing schema v20 and opens an
existing database only when its exact schema is v20. Stage JSON accepts only
version 4. Joint JSON accepts only schema 20 and stores target-only Pure
Exclusion results. No format performs migrations or compatibility decoding.
Unsupported inputs are rejected without repair or overwrite.

Current documents can store qualified and unqualified results, nullable
diagnostics, missing component states, unavailable width calibration, and
unavailable Halo reconstruction. The writer does not select a schema from
result contents.

The [data I/O architecture](architecture/dataobject-io-architecture.md)
documents the database and JSON boundaries. The
[stage integration guide](joint-stage-integration.md) describes downstream
ownership and invalidation.
