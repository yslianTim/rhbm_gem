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
