# Block 21 A/C solve decision

## Decision

Do not add a production multiprecision fallback or change the A/C
parameterization for the current 6Z6U Block 21 failure. Keep the existing
SparseQR solve, bounded free-face QR correction, active-set constraints, and
`1e-10` KKT gate. A candidate that cannot pass the unchanged formal checks
remains rejected.

## Evidence reused

`LinearSolve.cpp` already normalizes free columns and, after a failed KKT
certificate, can materialize a bounded free-face matrix when
`rows × free_columns <= 1,000,000` (8 MB for its double entries). It applies
at most eight QR residual corrections. This is local to one block; it does not
form a component-wide dense A/C system.

The bounded correction passes the captured Block 82 case at scaled KKT
`7.56219e-12` and Block 46 at `1.7344e-12`, while Block 21 still has scaled
KKT `4.65661e-10`. Block 21's 930 × 24 design is full rank, but its free-face
condition indicator is about `6.62e10`; SparseQR, dense QR, SVD, active-face
enumeration, and direct free-face QR all miss the existing gate. These
comparisons and the correction's regression bounds are recorded in
[`joint-local-ac-kkt-validation.md`](joint-local-ac-kkt-validation.md).

The 50-digit constrained QR reference reaches KKT about `3.0e-66` on the
rounded double design. However, its rounded double coefficients pass the
binary64 certificate at `5.47e-11` while the same rounded candidate checked
with 50-digit functional kernels has KKT `1.112e-10`. The production
candidate's 50-digit fixed-design KKT is `4.846e-10`; its scaled objective is
within roughly `2.34e-19` of the reference while its column-scaled
coefficients differ materially. A binary64 pass is therefore insufficient
evidence for a formal fallback result. The full comparison is in
[`joint-fullabc-width-degeneracy.md`](joint-fullabc-width-degeneracy.md).

## Scope and follow-up

The offline audit already uses `boost::multiprecision::cpp_dec_float_50` with
column-normalized Householder QR and an active-set loop. Boost headers are
available through the project's configured dependencies. This establishes a
bounded research reference, but not a production fallback that returns a
double candidate passing both primary and independent numerical
certification. No such candidate is established for Block 21, so integrating
one would risk publishing an unqualified state.

Column normalization remains in the existing solve. An additional A/C basis
rotation was not added: an unconstrained orthogonalization would also have to
preserve the atom-wise `A >= 0` feasible set exactly, and no tested transform
does that while improving Block 21's formal certificate. Parameter values
remain in the original FullABC model, with signed C and nonnegative A.

This is a research-only disposition. It does not relax the objective, active
constraints, KKT threshold, global width-gradient threshold, or independent
endpoint verification. The unresolved Block 21 solve remains an explicit
numerical limitation for bounded recovery and final 6Z6U validation.
