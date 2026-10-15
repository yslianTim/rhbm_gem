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
| Production bounded free-face refinement | 1.56083675296648 | 7.56219e-12 | Passes; active atoms 8, 11 |

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

## Production correction and regression result

The production solve now triggers a column-normalized dense QR refinement only
for a sparse, unit-weight local A/C problem whose initial full-rank face fails
the existing KKT certificate. It materializes the free-face matrix only when
`rows × free_columns <= 1,000,000` (at most 8,000,000 bytes for its double
entries), uses the existing rank threshold, and performs at most eight
residual-correction solves. Each candidate is checked through the unchanged
`CertifyLinear` path; it is adopted only when feasibility and the original
`1e-10` KKT gate pass and its objective does not exceed the initial SparseQR
objective. Weighted inputs, reference-SVD solves, rank-deficient faces and
larger local matrices skip refinement.

On the captured 881 × 24 problem, production returns the same active A atoms
8 and 11, preserves feasibility, and reduces scaled projected KKT to
`7.56219e-12`. Its objective is below the original SparseQR candidate. The
independent 4,096-face enumeration selected the same minimum-objective active
face; the corrected production objective is slightly lower than the enumerated
face's computed objective, consistent with rounding sensitivity in this
ill-conditioned problem. Interior, active-boundary, one-release, and
rank-deficient regressions pass; the rank-deficient case remains
`rank-deficient`.

These results identify a full-rank but ill-conditioned free face as the cause
of the original local KKT failure. The active-set release decision and release
scale were not the cause: the active dual violations were zero, all release
scale variants retained the same face, and none released a variable. Plain
dense QR and SVD did not meet the gate either. The bounded refinement improves
the certified local solve without changing the objective, rank threshold, or
KKT threshold.

## Later 6Z6U local KKT captures

The replay correction moved the full-data run past the old sweep 2 / block 82
replay disagreement. The block 82 trial fixture now passes its local replay
comparison with prediction difference 0, width-gradient difference
`7.5871e-15`, and KKT difference 0. The primary profile's scaled projected
KKT is `5.7795e-12`; the original failed decision is retained in the fixture as
historical capture evidence.

The bounded double-only refinement also recovers the captured sweep 2 / block
46 problem: its 2,889 × 24 design has rank 24 and scaled KKT `1.7344e-12`.
The older sweep 4 / block 19 fixture remains a valid expected rejection at
`2.3283e-10`; it is not the final full-data trajectory.

The final production 6Z6U run instead stops at component 1, one-based sweep 4,
block 21. Its failure-only fixture is
`tests/fixtures/joint_local_ac_kkt/block-21-sweep-4.json`: 930 rows × 24
columns, 2,894 nonzeros, design rank 24, free rank 23, condition indicator
`6.6235e10`, feasible A/C coefficients and no active-set releases. The solver
is valid with reason `solved`, but the scaled projected KKT is
`4.6566e-10`, above the unchanged `1e-10` gate. The regression requires this
state to remain rejected. This later KKT failure is unresolved; endpoint
verification and numerical qualification were not reached. See the
[full-data validation](joint-halo-pruning-validation.md#full-data-validation-6z6u).

## FullABC width and precision follow-up

The Block 21 local design is numerically full rank (24/24), with a 23/23 free
face and condition indicator `6.6235e10`. Its production constrained A/C
candidate remains rejected at scaled KKT `4.656612873e-10`; the unchanged gate
is `1e-10`. Dense QR, SVD, active-face enumeration and free-face QR also miss
that gate. A bounded 50-digit, column-normalized Householder QR reference keeps
the active constraint on atom index 5 and reaches KKT `2.9989e-66` on the
rounded double design. Evaluating the production coefficients at 50 digits
gives KKT `4.8461e-10`. This identifies a binary64 solve/certificate error in
a full-rank but ill-conditioned problem; it does not establish structural rank
deficiency.

The constrained reference objective differs from the production objective by
only about `2.34e-19` after objective scaling, while the maximum
column-norm-scaled coefficient difference is `3.63e3`. The atom 11721 G/K
terms have magnitudes about `+9.23e8` and `-9.23e8`; atom 11650's captured
Gaussian column norm is `2.10e-38`. A separate conditional width scan shows a
flat large-width profile and small-width binary64 information loss. Those
findings are detailed in the [FullABC width-degeneracy study](joint-fullabc-width-degeneracy.md).

Rounding the 50-digit fixed-design solution to double yields a binary64 KKT
of `5.47e-11` (pass), while its exact rounded coefficients and design give
50-digit KKT `1.1120e-10` (fail). Binary64 agreement therefore does not prove
precision-independent stationarity. The permanent Block 21 regression keeps
the production result rejected and compares it to the bounded constrained QR
reference without normal equations.
