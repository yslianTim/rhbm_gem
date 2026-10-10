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

## Historical decisions

| Topic | Conclusion |
| --- | --- |
| OperatorPcg and Schwarz | Investigated for large connected components, but the required promotion campaign did not complete all large-case gates. They were not adopted. The current estimator is FixedNeighbor with local LegacyCompact and EIGEN. |
| Factor ownership and residency | Factor-sharing and persistent-workspace variants did not justify a current route or factor cache. Some 2048 Schwarz diagnostics exceeded the 4 GiB process limit while constructing the SPQR factor, before search PCG began. |
| Bounded rank budget | Small cases matched dense controls, but representative large SPQR prototypes returned Unavailable before reconstruction completed. No large-case rank or scaling guarantee was established. |
| Projected width and tail | Completed cases showed numerical parity and some faster assessment phases, but the projected-width candidate exceeded the 4 GiB limit on both 512-atom topologies; projected-tail symbolic work exceeded the limit on cube-512. Neither was promoted. |
| Operator and search scaling | Available Schwarz/PCG points showed no clear iteration growth, but incomplete large-case coverage and factor-memory limits did not establish a viable production route. |
| Cube ordering | Ordering choices reduced fill and RSS on completed 512/1024 cases, but no ordering completed cube-2048 within the 4 GiB gate. Production ordering was unchanged. |

## SPQR versus EIGEN

Both backends passed bounded numerical parity. SPQR was faster in the
completed 128/256 resource probes, while EIGEN used less peak RSS in every
completed case:

| Case | EIGEN seconds | SPQR seconds | EIGEN RSS MiB | SPQR RSS MiB |
| --- | ---: | ---: | ---: | ---: |
| chain-128 | 24.20 | 10.71 | 68.48 | 82.61 |
| cube-128 | 42.81 | 22.60 | 65.23 | 81.98 |
| chain-256 | 52.37 | 24.75 | 114.33 | 142.75 |
| cube-256 | 99.76 | 55.71 | 107.64 | 142.89 |

Both bounded 512 probes timed out. The original wall-time-first qualification
selected SPQR based on the completed cases. Production later prioritized
memory footprint, dependency and installation simplicity, and maintenance;
EIGEN became the sole backend. SPQR did not fail numerical parity, and its
historical speed advantage is not erased by the later decision.

Closed experiment implementations and machine-readable receipts are
recoverable from Git history and are intentionally not retained in the current
tree. This document keeps only the decisions needed to understand the current
estimator and backend.
