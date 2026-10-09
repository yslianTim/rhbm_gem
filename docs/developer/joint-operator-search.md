# Retired OperatorPcg investigation

Status: **Closed / Removed**

This is a historical decision record, not an implementation guide. The
OperatorPcg route was investigated as an alternative global Joint solver and
was not retained after the FixedNeighbor qualification.

## Decision

OperatorPcg introduced a large operator, rank, preconditioner and Schwarz
surface around the Joint estimator. The associated resource-scaling evidence
did not justify keeping a second production route, while qualified FixedNeighbor
met the existing stationarity, endpoint-certification and `RuntimeConvergence`
contracts with a fixed core of 12 atoms, Forward order, 30 maximum sweeps, and
one trusted local update per block visit. FixedNeighbor is therefore the current
production/default estimator.

The current implementation no longer contains `OperatorPcg`,
`SearchOperatorProfile`, `ProfileJacobianOperator`, the OperatorPcg PCG
implementation, or its preconditioner/Schwarz implementation. There is no
current CLI, Python, benchmark, or route-selection option for `operator-pcg`.
The old option is rejected rather than silently mapped to FixedNeighbor.

## Current relationship to LegacyCompact

The former top-level `LegacyCompact` route was retired after FixedNeighbor
promotion. The existing `LegacyCompact` profile solver is still active inside every
FixedNeighbor block visit: it performs the local width search, after which at
most one trusted update is replayed globally. Removing OperatorPcg did not
remove this local numerical primitive.

## Historical compatibility

`JointSolverProvenance` retains the old operator, rank, preconditioner and
Schwarz fields as legacy provenance only. Historical v2 records with
`search_method = operator-pcg` can still decode and round-trip those values;
current writers do not emit them as solver controls. This compatibility path
does not recreate the removed solver.

## Evidence

Compact historical summaries and figures remain under
[`docs/developer/figures/`](/docs/developer/figures/). The current experiment
inventory is in [`joint-experiments.md`](joint-experiments.md), and the active
benchmark surface is documented in [`joint-benchmark.md`](joint-benchmark.md).
The retired implementation and raw campaign wrappers remain recoverable from
Git history; they are not rebuilt or rerun by current validation.
