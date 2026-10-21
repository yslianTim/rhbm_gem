# Joint component estimator runtime

This guide describes the production target-only Pure Halo-Support Exclusion
estimator. The persisted parameterization is
`target-only-pure-halo-support-exclusion-v1`.

## Production flow

`RunJointComponentWorkflow` builds the Pure Exclusion problem and obtains
qualified target widths from FirstStage fitting. A width calibration is
available only with at least 20 qualified target donors and 95% target
coverage. The global width bound is
`Bmax = 4 × P95(qualified FirstStage target widths)`.
Calibration failure remains an unavailable diagnostic and the search does not
start with an unbounded substitute.

Targets are the only free Joint components. Non-target Halo support is used
only for the Pure Halo-Support Exclusion row selection, with the production
support radius of 2.5 Å. The target model estimates Gaussian amplitude A,
width B, and signed Coulomb charge C.

Production uses FixedNeighbor blocks with core size 12, forward serial
Gauss-Seidel order, and at most 200 sweeps. A block visit may accept at most
one trusted update. Width updates use bounded Active-set Levenberg–Marquardt;
constrained A/C subproblems use SparseQR with the Eigen backend. Amplitudes
satisfy A ≥ 0 and charge C may have either sign.

## Qualification

The global search preserves objective monotonicity, practical global KKT
convergence, complete-sweep stationarity, and endpoint verification. A
qualified endpoint must also pass Independent Replay. Failed initialization,
unavailable width calibration, incomplete search, unqualified endpoint,
missing component state, and nullable diagnostics remain representable saved
outcomes.

Local Conditioning and Parameter Stability describe the accepted solution's
block behavior. Target point validation records component convergence,
parameter validity, conditioning, stability, point quality, and reasons.
These diagnostics retain `Unavailable` and `NotAssessed` states where work
could not be assessed; they do not replace global convergence evidence.

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
| Joint Result JSON | Schema v19 |

SQLite opens a new empty database by initializing schema v20 and opens an
existing database only when its exact schema is v20. Stage JSON accepts only
version 4. Joint JSON accepts only schema 19 and stores target-only Pure
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
