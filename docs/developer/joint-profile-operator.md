# Retired ProfileJacobianOperator investigation

Status: **Closed / Removed**

`ProfileJacobianOperator`, `OperatorSearch`, and the associated
preconditioner/Schwarz implementation were historical research infrastructure.
They are no longer part of the current Joint solver, benchmark surface, or
public API.

The active numerical path is FixedNeighbor. Each structural block invokes the
existing local Profile LM `SearchProfile()` operation, then applies
the existing trust, global replay, stationarity, endpoint-certification and
`RuntimeConvergence` contracts. `LinearSolve` owns the EIGEN sparse solve;
`StructuralPartition`, `CompactSvd` and `TiledDerivative` provide the current
partitioning and local numerical work. `ResourceWork` owns resource
attribution; `NumericsWork` owns numerical profiling counters.

Historical fixed-state parity, factor-residency and preparation conclusions
are summarized in [joint-component-evidence.md](joint-component-evidence.md).
They do not define a current operator class or a route that can be selected at
runtime. Closed campaign implementations and outputs remain in Git history.
