# Joint local A/C KKT validation

## Captured 6Z6U local problem

The failure-only fixture at `tests/fixtures/joint_local_ac_kkt/block-82.json`
replays the initial fixed-width A/C profile from the 6Z6U run. It contains the
local support geometry, effective response, widths, row and atom identities,
and numerical policy; it does not contain the CIF or map. The recaptured block
is component 1, sweep index 2, block 82, with 12 core atoms. The prior summary
called this Sweep 1; the recaptured one-based field says 2, so the sweep-label
discrepancy remains to be resolved against the full run log.

The affected observation count is 881. The local design is **881 rows × 24
columns**, with 3,368 nonzeros and numerical rank 24. The initial constrained
solve has free rank 22 because A amplitudes for local atoms 8 and 11 are active
(amplitude columns 16 and 22). Both A constraints are feasible, their dual
violations are zero, and the solver released no active variables. The solver
was valid with reason `solved`, but the profile failed `kkt-failed` after three
linear solves.

The local response norm is 12.6016184156; the parent response norm and
observation scale are 346.783950262, a ratio of 27.519. The full and free
column-normalized designs are both full rank, with condition indicators
5.317884380e9 and 5.119199361e9. Their singular-value ranges are 1.94836 to
3.66379e-10 and 1.87597 to 3.66457e-10, respectively. Thus this is a full-rank
but ill-conditioned local problem, not structural or numerical rank loss.

## Solver and active-face comparisons

| Method | Objective | Scaled projected KKT | Result |
| --- | ---: | ---: | --- |
| Production SparseQR | 1.56083675341943 | 4.65661e-10 | Fails `1e-10` |
| Dense rank-revealing QR | 1.56083675373191 | 4.07454e-10 | Fails |
| Dense SVD | 1.56083675371222 | 4.08145e-10 | Fails |
| Best objective from 4,096-face enumeration | 1.56083675385495 | 4.07454e-10 | Same active atoms 8, 11; fails |
| Direct dense QR on that free face | 1.56083675425541 | 5.82077e-10 | Fails |
| Alternate `long double` free-face QR | 1.56083675330617 | 5.82077e-11 | Passes |

The `long double` result is not higher precision on this host:
`std::numeric_limits<long double>::digits` is 53. It is an alternate
factorization and column-scaling path, and is treated as a sensitivity check,
not as a high-precision oracle.

Enumeration found 4,096 full-rank faces and 1,251 primal-feasible candidates.
The lowest-objective face has active atoms 8 and 11. One other face passes the
current scaled KKT check (`1.35128e-13`) but has objective 1.56492424545, about
0.00409 above the minimum feasible candidate. This shows why this ill-conditioned
case must compare KKT and objective together: the scaled KKT threshold alone
does not distinguish that higher-objective face from a constrained optimum.
The minimum-objective face and the active set selected by production agree.

## Release-tolerance experiment

The release tolerance is
`release_factor * epsilon * max(1, ||Z|| (response_norm + ||Z|| ||b||))`.
The local context copies the parent context, including its response norm, while
recomputing the local rank threshold. This appears to be inherited policy
behavior; the implementation does not document a separate local-norm contract.
On this fixture, parent-norm, local-effective-norm, and consistently scaled
local-RHS cases all retain the same active set, make zero releases, and return
the same coefficients. The raw release tolerance changes from
1.68299378910e-4 to 1.68299332380e-4, which does not explain the KKT failure.
The active A dual violations are zero, so an incorrect release decision is not
the cause.

## Bounded correction experiment

Starting from the production coefficients, a dense QR solve restricted to the
enumerated minimum-objective free face was used only in this offline, 881 × 22
experiment. Two residual-correction solves reduced the scaled KKT residual to
the existing threshold while leaving the same active A set. The first
correction passed KKT but increased the computed objective by 3.89e-10; after
the second correction KKT was 5.82e-11 and the objective was 6.29e-10 lower
than the original production result. Further corrections made the computed
objective fluctuate around this value, so this experiment does not justify
unconditionally applying an arbitrary number of corrections in production.
The production fix must select only a KKT-passing candidate whose measured
objective does not exceed the current candidate, and then be checked through
the unchanged `CertifyLinear` path.

This experiment rules out changing the release norm or replacing SparseQR with
plain dense QR/SVD as sufficient fixes. It supports a bounded free-face
residual-refinement step, but the production implementation and regression
evidence are tracked separately from this assessment.
