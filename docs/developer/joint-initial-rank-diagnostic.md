# 6Z6U initial-rank failure: retained findings

The investigation at base commit `5f06b027` identified a **support/geometry
limitation in full-parameter estimation**, rather than an initialization problem.
Only this textual summary is retained; the experiment-specific tools, tests,
captures, builds and data bundles have been removed. Production initialization
and solver policy are unchanged.

## Observed failure

The selected domain contained 2,192 contributors and 262,801 observation rows
in two components. The largest component contained 2,167 atoms, 262,282 rows
and 4,334 A/C columns, matching the historical census. Historical per-atom B₀
mapping was unavailable, so this was a new baseline from the same inputs,
not a claim of bitwise reproduction.

Both EIGEN and SPQR rejected the **first all-free A/C face**, reporting rank
**4314 / 4334**, before outer search. Four fixed SPQR starts—production B₀,
0.9 × B₀, 1.1 × B₀ and per-atom simulation truth—gave the same deficient rank
and no admitted initial state. Independent no-tolerance QR reduction followed
by SVD also reported rank 4314 for all four starts; weak-direction actions
were checked against the original sparse matrix.

## Why it fails

**Twenty halo atoms each contribute to exactly one selected observation voxel.**
For each such atom i, its A and C columns are proportional:

\[
X_{A_i}=g_i(B_i)e_{r_i},\qquad X_{C_i}=h_i(B_i)e_{r_i}.
\]

The parameter change `(delta A_i, delta C_i) = (h_i, -g_i)` therefore leaves
the predicted observations unchanged. Each atom uses its own pair of parameter
coordinates, giving 20 independent null directions and the structural bound
`rank <= 4334 - 20 = 4314` for any legal positive B with this fixed support.
Nonnegative A does not restore uniqueness: an increase of A can be compensated
by unconstrained C. The measured weak subspace was concentrated in these halo
A/C coordinates.

The optional Jacobi numerical cross-check timed out, and the conservative
rank-boundary diagnostic remained unconfirmed. This limits the numerical
certificate, but does not invalidate the separate structural argument above.

## Implications

Changing B₀, accelerating factorization or introducing operator LM cannot remove
this structural deficiency while retaining the same observation domain and
full-parameter contract. Resolve that issue before increasing data scale.

One candidate is to retain all contributors but represent each one-voxel halo's
observable contribution as `lambda_i = A_i g_i(B_i) + C_i h_i(B_i)`, explicitly
leaving its individual A/C/B unidentified. This needs separate validation of
target identifiability, objective equivalence, derivatives and result semantics.
Another candidate is to expand the observation domain and recompute contributor
closure; this changes the statistical problem and can introduce further halos,
so it does not guarantee full rank.

The finding does not justify deleting 20 atoms, relaxing rank thresholds or
claiming identifiability from a minimum-norm solution. Neither target-width
identifiability nor full 6Z6U convergence was established.
