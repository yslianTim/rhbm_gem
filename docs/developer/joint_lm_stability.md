# Joint bounded LM stability audit

Baseline: develop, 32609b2339736972c030a503778e868dc6647d67; clean tracked tree.
The ignored `output/diagnostics/7a6a_joint_components` evidence was available and read.

## Production responsibilities

| Responsibility | Implementation |
| --- | --- |
| FirstStage donor coverage, bounded initialization | `JointComponentEstimator.cpp`: `InitializeContributors`, `FitJointComponentsInitializedWithProgress`; `WidthCalibration.cpp` |
| Component partition and orchestration | `Components.cpp`: `OrchestrateFixedNeighborComponents`, `SolveFixedNeighborComponentView` |
| FixedNeighbor cores and effective local response | `FixedNeighborBlockCoordinate.cpp`: `SearchFixedNeighborComponent` |
| Constrained A/C profile, compact Width derivative | `Search.cpp`: `Profile::Get`, `linearize`; `Numerics.cpp`, `LinearSolve.cpp`, `TiledDerivative.cpp` |
| Bounded Width LM, damping, active bounds | `InstrumentedLM.hpp`: `minimizeOneStep`; `ActiveSetBoundedLM.cpp`: `SolveBoundedTrustRegionStep`, `SolveAtDamping` |
| Local candidate replay and acceptance | `Search.cpp`: `Profile::Trial`, `Accept`; initial A/C acceptance is separate from accepted Width updates |
| Global objective acceptance and safe state commit | `FixedNeighborBlockCoordinate.cpp`: candidate finite/feasibility/replay checks precede cache mutation; failed/rejected candidates retain accepted state |
| Recovery and sweep stationarity | same file: failure visits/retries, full sweep replay, original-coordinate A/C and projected Width KKT |
| Endpoint and independent verification | `CertifyFixedNeighborEndpoint`: streaming primary assessment, independent replay, trust and numerical qualification, memory preflight |
| Component qualification | `Components.cpp`: evidence/state creation; `JointComponentEstimator.cpp`: `JointComponentResult::RuntimeConvergence` |
| Assembly | `Components.cpp`: `AssembleComponents`: qualified component mapping, preserved coefficients, assembled replay/consistency |
| Downstream | `PotentialFittingWorkflow.cpp`: qualified snapshot gates Halo reconstruction, PostFitPeeling, Local and Group fitting |

## S0 baseline contracts

Controlled regression tests record (without changing production) the five-column
rank/bracket failure, near-unit Width radius collapse and coordinate-origin
discontinuity, and historical LM failure veto. Subsequent stages replace these
characterization assertions with the corrected contracts. Existing recovery,
rollback, nonfinite, endpoint and replay tests remain mandatory.

All build/test commands use at most four threads. Real-data runs use new temporary
SQLite files; diagnostic programs, databases and large output remain untracked.

## S1: increment coordinates and null sensitivity

The existing augmented QR solves in `z_i=S_i p_i`, with
`S_i=max(D_i, ||J_i||_2, sqrt(min_normal))`. The damping block is
`sqrt(lambda) D_i/S_i`; conversion back divides the solution by `S_i`.
Bounds, bound hits/releases, multipliers and KKT are checked directly in `p`.
No transformed bound tolerance is used. Scaling fills the existing augmented
workspace and adds only O(number of Widths) storage, with no dense fallback.

A direction is held at zero for this linearization only when its column norm is
below `64 epsilon` times the largest column norm AND its original gradient is at
most `64 epsilon` (normalized objective units). The full candidate gradient is
checked again, including held directions. This does not certify identifiability;
a later A/C or neighbor change recomputes the criterion. Small sensitivity with
nonstationary gradient is explicitly tested and remains movable.

The controlled five-Width case now takes one QR and zero damping brackets rather
than 81 failed QRs. Genuine rank deficiency still uses the existing bounded
positive-damping QR search. No SVD or additional production solver is needed.
S1 validation: Release build; 69 tests across ActiveSetBoundedLM, Search, Runtime,
and FixedNeighbor (including rollback/recovery) passed.

## S2: origin-independent trust radius

The initial radius is `factor * ||D||_2`, the weighted norm of a representative
unit log-width increment (production factor 0.1). It is independent of eta's
origin and limits the initial aggregate model displacement; bounds, replay,
weak-direction safeguards and actual objective reduction still control updates.
The trust-ratio regression explicitly selects a radius that reaches its active
bound rather than depending on the absolute starting coordinates.

Tiny steps use the original log-width increment infinity norm, not `||D p||`.
Production tiny-step success requires the original projected Width KKT plus valid,
feasible A/C. Nonstationary tiny steps return a numerical stop for block recovery.
A cosine-based early stop also cannot bypass production projected KKT.
The radius regression covers B=1, B=1.000000000626 and shifted reference origins,
finite evaluated candidates, objective improvement and final bound-aware KKT.
S2 validation: Release build and all 75 focused tests passed.

## S3: terminal qualification and recovery

`SearchResult::width_failed` separates LM subproblem/nonstationary tiny-step
failure from initial A/C profile acceptance. FixedNeighbor routes these failures
through the existing retry/rollback path even when the initial profile is trusted.
Failure logs retain initial profile trust, Width stop reason and solver history;
progress failed-block counts now include these block visits.

`accepted-state-integrity` is required component evidence: all accepted mutations
and sweep replays must remain intact, and no block recovery may remain pending.
Its value is the unresolved block count; its reason retains recovery attempts and
successes. Fatal stops cannot create qualified component state. Historical
`active-set-lm-subproblems` failures remain visible (including Failed status) but
are no longer themselves a convergence veto. Complete sweep stationarity, both
KKT gates, endpoint verification, independent replay, feasibility and terminal
integrity remain required. Missing integrity evidence rejects live assembly.

The existing extensible evidence serialization carries the new check without
changing JSON structure or introducing compatibility paths. The diagnostic test
facade now forwards actual terminal integrity instead of omitting it. Regression
coverage includes historical failure with qualified assembly, failed/missing
terminal evidence, endpoint/replay/KKT rejection, nonfinite state, successful and
exhausted recovery, and trusted initial A/C followed by an injected Width failure.

S3 validation: 77 focused tests and all four Joint CTest groups (contract,
workflow, numerical, persistence) passed. The ParameterStability replay-failure
regression now prescribes a nonstationary Width directly; it no longer relies on
the old solver stopping at a particular block under loose search criteria.
Production diagnostics also retain damping solve, bracket and QR totals in the
existing LM evidence reason.
