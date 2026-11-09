# Joint component numerics

This guide specifies the local block solve. The [runtime guide](joint-component-runtime.md)
defines the model, parent scale, block acceptance and global qualification.

## Constrained A/C least squares

At fixed log-widths `eta = log(tau)`, the Gaussian and charge-kernel columns
of `X(eta)` define an exactly linear model. For the block effective response
`y_eff`, solve `min_beta ||X beta - y_eff||² / 2` with nonnegative amplitude
entries and signed, unconstrained charge entries in the interleaved A/C vector
`beta`. The solver scales each design column by its Euclidean norm and maps
coefficients back to physical units.

The active set controls which amplitudes can move. Each face is solved by
Eigen `SparseQR` with `COLAMDOrdering` after a 1,024-row tiled orthogonal
reduction. Discarded orthogonal response energy, including zero-design rows,
still counts in the original residual and objective. All C columns stay free.
An infeasible step moves to the first blocking amplitude face, and the dual
gradient releases a blocked amplitude when moving into the feasible region
reduces loss. A released zero amplitude is immediately eligible for the QR
solve. The certificate separately records amplitudes at geometric equality.
Each solve starts with a fresh working set and factors.

QR must have full free-column rank. Its absolute pivot threshold is the
configured relative threshold times the maximum normalized free-column norm
used by the solver. The block relative threshold is double epsilon times
`max(block rows, 2 * block atoms)`. Zero columns, deficient rank, failed
solves and active-set exhaustion return explicit reasons. A deficient face
cannot supply a qualified linear minimizer for the profile derivative.

Let `d_j = ||X_j||₂`, `u_j = d_j beta_j / s`, and
`g_j = X_jᵀ(X beta - y_eff) / (d_j s)`, with the fixed parent scale `s`.
Local feasibility and `||u - P(u - g)||∞ ≤ 1e-10` are required, where `P`
clips amplitude coordinates at zero and leaves charge coordinates alone.
When the certificate is available but fails, the solver may refine on the
same normalized face with full-rank pivoted dense QR and up to eight
corrections. It adopts a refined point only if it passes the certificate
and does not increase the original loss. A remaining failure rejects the
profile evaluation even under the practical global convergence policy.

## Variable projection and the profile derivative

Re-solving constrained A/C at each feasible width leaves only the log-widths
in the nonlinear search. For `r = X(eta) beta*(eta) - y_eff`, the profile
objective is `phi(eta) = ||r||₂² / (2 s²)`. On a fixed, full-rank endpoint
face, let `Z` contain the column-normalized free design (all C and positive
A), `P_Z` be its orthogonal projector, `W_j = (∂X/∂eta_j) beta*`, and
`T_j = (∂X_free/∂eta_j)ᵀ r` with each entry divided by the corresponding
current free-column norm. Then the residual Jacobian column is

```text
J_j = ((I - P_Z) W_j - Z (Zᵀ Z)^(-1) T_j) / s
```

The second term accounts for the coefficient response at nonzero residual.
The endpoint face is rebuilt from the evaluated coefficients, independently
of the A/C solver's mutable working set. The kernel derivative uses both
Gaussian and charge terms. `TiledDerivative` uses 8,192-row Householder
reductions and triangular solves for the projection and correction. It retries
in reference accumulation order when projection nearly cancels a raw column.
It rejects unavailable or rank-deficient full derivatives and reduces valid
`J` and `r/s` to a compact model, keeping the full residual norm for trial
comparisons. The projected-only diagnostic can use a rank-truncated projection
to report deficient-rank evidence.

## Bounded log-width LM

[ProfileLmSearch.cpp](../../src/core/detail/joint_component/ProfileLmSearch.cpp)
owns the JCE problem adapter. It supplies the qualified A/C evaluations, full
residual norms, reduced derivatives and cached accepted stationarity. It
requires a finite `eta_max = log(tau_max)`, validates initial and trial
widths, and passes `eta_max - eta` as the step's upper bounds. A proposed
trial must pass local replay of the profile and its feasibility and agreement
checks. A step exceeding one log-width decade is rejected when at least 95%
of its norm lies in weak right-singular directions of the trial Jacobian and
its nonnegative improvement fits within `1e-12 + 2e-12 |phi_before|`.

[BoundedLmDriver.hpp](../../src/core/detail/joint_component/BoundedLmDriver.hpp)
owns the outer LM state and retry loop. It scales by Jacobian column norms,
substituting one for a zero initial norm and retaining columnwise maxima
across accepted iterations. Its initial radius is `0.1 ||D||₂` (the norm of
the scaling diagonal), independent of the log-width origin. It compares the
full residual norms for actual reduction and uses the reduced model and the
returned damping for predicted reduction, both for the actual feasible step.
It proposes acceptance at ratio `≥ 1e-4` and updates the accepted width and
residual norm only after the adapter trusts the proposal. An untrusted trial
shrinks the radius by four and multiplies the damping seed by four while retry
budget remains. A trusted low-ratio trial follows the ratio-based trust update
instead. The adapter publishes acceptance separately from its trial cache.
The block caller still performs global candidate acceptance.

The driver reports reduction, step, stationarity, budget, precision or problem
stops. The adapter records JCE reasons such as invalid inner solve, replay
failure or `bounded-subproblem:<reason>`. Production uses `ftol=1e-14`,
`xtol=1e-12`, `gtol=1e-12`, a 200-profile-evaluation budget and at most one
accepted width update per block visit, after the separately checked initial
profile. A tiny physical log-width step requires accepted bound-aware
stationarity. Local stopping criteria control the iteration. Component
convergence requires the global checks in the runtime guide.

[ActiveSetBoundedLM.cpp](../../src/core/detail/joint_component/ActiveSetBoundedLM.cpp)
owns the inner subproblem. With reduced model `J_r, r_r`, diagonal `D`,
radius `Delta` and damping `lambda`, it solves damped bound-constrained
least squares and selects damping to meet the trust radius.

```text
min_p  (||J_r p + r_r||² + lambda ||D p||²) / 2
       p_j <= eta_max - eta_j,  ||D p||₂ <= Delta
```

At fixed damping, the working set holds upper bounds while scaled, pivoted QR
solves the augmented free problem. Steps stop at blocking faces and release
a bound when its damped gradient is positive beyond tolerance. At a valid
upper face the gradient is nonpositive and the multiplier is nonnegative.
A separate stationary-null mask freezes only directions with negligible
Jacobian norm and gradient. Working bounds, frozen directions and geometric
equality therefore describe different sets. The solver verifies feasibility and
damped KKT in the original step coordinates. It tries zero damping first,
then brackets and bisects a feasible damping when needed, including recovery
from a deficient undamped free subproblem. A failed bounded solve reports its
reason without replacing accepted state.

## Compact spectra and Local Conditioning

Eigen is the production numerical backend. `CompactSvd` uses JacobiSVD when
the smaller matrix dimension is below 16 and BDCSVD otherwise. A failed,
nonfinite or near-threshold BDC result gets one Jacobi retry. The caller supplies
the observation and decision dimensions for the relative threshold or a
nonnegative absolute override, and selects the native SVD or strict
`sigma > threshold` boundary. Failed decompositions remain unavailable.

Local Conditioning uses block-sized compact QR/SVD on the endpoint's
conditional free A/C design, projected width directions and full profile
Jacobian where available. It reports rank, sensitivity near the rank threshold,
and normalized spectra. Full-rank spectra warn when threshold-sensitive or
`condition * epsilon ≥ 1e-8`. Active amplitude faces retain their separate
classification. These read-only diagnostics assess local numerical geometry.
They supply evidence to [Target Point Quality](joint-stage-integration.md)
and carry no global rank, covariance or statistical precision certificate.

See [ConstrainedAcLeastSquares.cpp](../../src/core/detail/joint_component/ConstrainedAcLeastSquares.cpp),
[TiledDerivative.cpp](../../src/core/detail/joint_component/TiledDerivative.cpp),
[CompactSvd.cpp](../../src/core/detail/joint_component/CompactSvd.cpp) and
[LocalConditioning.cpp](../../src/core/detail/joint_component/LocalConditioning.cpp).
The focused [numerical tests](../../tests/core/joint_component) cover active
faces, replay, derivative and rank boundaries, and feasible bounded steps.
