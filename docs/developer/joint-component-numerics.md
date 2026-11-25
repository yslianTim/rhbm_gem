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
corrections. The dense face is capped at one million elements. The solver
adopts a refined point only if it passes the certificate and does not
increase the original loss. A remaining failure rejects the profile
evaluation even under the practical global convergence policy.

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

[ProfileLmSearch.cpp](../../src/core/estimator/joint_component/ProfileLmSearch.cpp)
owns the JCE problem adapter. It supplies the qualified A/C evaluations, full
residual norms, reduced derivatives and cached accepted stationarity. It
requires a finite `eta_max = log(tau_max)` and accepts an optional
`eta_min = log(tau_min)`. It validates initial and trial widths and translates
the width interval into step bounds. A proposed
trial must pass local replay of the profile and its feasibility and agreement
checks. A step exceeding one log-width decade is rejected when at least 95%
of its norm lies in weak right-singular directions of the trial Jacobian and
its nonnegative improvement fits within `1e-12 + 2e-12 |phi_before|`.

[BoundedLmDriver.hpp](../../src/core/estimator/joint_component/BoundedLmDriver.hpp)
owns the outer LM state and retry loop. It scales by Jacobian column norms,
substituting one for a zero initial norm and retaining columnwise maxima
across accepted iterations. Its initial radius is `0.1 ||D||₂` (the norm of
the scaling diagonal), independent of the log-width origin. It compares the
full residual norms for actual reduction and uses the reduced model and the
returned damping for predicted reduction, both for the actual feasible step.
It proposes acceptance at ratio `≥ 1e-4` and updates the accepted width and
residual norm only after the adapter trusts the proposal. An untrusted trial
shrinks the radius by four and multiplies the damping seed by four while retry
budget remains. For trusted trials, a ratio at or below 0.25 shrinks the
radius and increases damping. A ratio at or above 0.75, or zero damping with
a ratio above 0.25, sets the radius to twice the scaled step norm and halves
damping. The adapter publishes acceptance separately from its trial cache.
The block caller still performs global candidate acceptance.

The driver reports reduction, step, stationarity, budget, precision or problem
stops. The adapter records JCE failures such as an invalid inner solve or
untrusted replay. The driver can also stop with a detailed cause, including
`bounded-subproblem:<reason>` from the inner solver. Production uses
`ftol=1e-14`, `xtol=1e-12`, `gtol=1e-12`, a 200-profile-evaluation budget
and at most one accepted width update per block visit, after the separately
checked initial profile. A tiny physical log-width step requires accepted bound-aware
stationarity. Local stopping criteria control the iteration. Component
convergence requires the global checks in the runtime guide.

[ActiveSetBoundedLM.cpp](../../src/core/estimator/joint_component/ActiveSetBoundedLM.cpp)
owns the inner subproblem. With reduced model `J_r, r_r`, diagonal `D`,
radius `Delta` and damping `lambda`, it solves damped bound-constrained
least squares and selects damping to meet the trust radius.

```text
min_p  (||J_r p + r_r||² + lambda ||D p||²) / 2
       eta_min - eta_j <= p_j <= eta_max - eta_j,  ||D p||₂ <= Delta
```

An omitted `eta_min` leaves the lower side unbounded. At fixed damping,
the working set holds lower or upper faces while scaled, pivoted QR solves
the augmented free problem. Steps stop at blocking faces. A lower face
requires a nonnegative damped gradient, and an upper face requires a
nonpositive one. The solver releases a working face when its gradient
violates that sign beyond tolerance.
A separate stationary-null mask freezes only directions with negligible
Jacobian norm and gradient. Working bounds, frozen directions and geometric
equality therefore describe different sets. The solver verifies feasibility and
damped KKT in the original step coordinates. It tries zero damping first,
then brackets and bisects a feasible damping when needed, including recovery
from a deficient undamped free subproblem. The trust feasibility allowance is
relative to the radius (`64 * epsilon * Delta`). A failed bounded solve
reports its reason without replacing accepted state.

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
They supply evidence to [Target Point Quality](joint-stage-integration.md#target-point-quality)
and carry no global rank, covariance or statistical precision certificate.

See [ConstrainedAcLeastSquares.cpp](../../src/core/estimator/joint_component/ConstrainedAcLeastSquares.cpp),
[TiledDerivative.cpp](../../src/core/estimator/joint_component/TiledDerivative.cpp),
[CompactSvd.cpp](../../src/core/estimator/joint_component/CompactSvd.cpp) and
[LocalConditioning.cpp](../../src/core/estimator/joint_component/LocalConditioning.cpp).
The focused [numerical tests](../../tests/core/joint_component) cover active
faces, replay, derivative and rank boundaries, and feasible bounded steps.

## A/C compensation at a qualified endpoint

The explicit research entry point `AssessAcCompensation` reads a block's
accepted parameters and affected Fitting Domain rows. It uses the production
kernel and support incidences. Other blocks' predictions are held fixed in
the effective response. The diagnostic owns its scratch state. Production
fitting, Point Quality and result persistence have no dependency on it.

When the sampled radii occupy a narrow band, a Gaussian and a charge column
can be almost proportional. For small `r/tau`, their center-normalized
shapes start with `1 - (r/tau)^2/2` and `1 - (r/tau)^2/6`. A large width
reduces shape differences within support, and the fitting mask can further
narrow the sampled radii. The global width bound restricts this flattening.
Even an interior width can have substantial compensation.

| Evidence | Definition and scope |
| --- | --- |
| Raw correlation | Cosine between the target's Gaussian and charge columns at its accepted width, with both Euclidean norms in original design units. A zero or subnormal norm prevents a cosine. |
| Conditional correlation and sensitivity | Residualize each free target column against the other free A/C columns in its block. Report the cosine and both residualized norms in original design units. The residualized-to-raw norm ratios measure how much independent sensitivity remains. A column exhausted to projection roundoff has a norm but no meaningful cosine. |
| Cancellation | `(||A g|| + ||C k||) / ||A g + C k||`. A single nonzero contribution gives 1. The routine reports zero contributions and a net signal within eight double epsilons of the summed contribution norms as distinct unavailable cases. |
| Weak directions | Singular values in normalized free-coefficient coordinates `u = D beta`, where `D` contains the free-column norms. For an isolated weakest direction `v`, apply `D^(-1) v` to obtain the physical coefficient change. For clustered or numerically null modes, report the subspace dimension and per-target squared projection mass divided by that dimension. Report the projection length of the unit opposing A/C contribution vector onto the weak subspace when applicable. |

The diagnostic reduces normalized sparse columns with tiled QR. The retained
`R` preserves the column norms, inner products and row-space geometry.
The right vectors of the transposed compact nuisance design give its
orthogonal span after the standard strict rank cutoff
`epsilon * max(observation rows, nuisance columns) * sigma_max`.
The same policy applies to the compact free spectrum. Numerically discarded
directions and singular values indistinguishable from the minimum within
compact roundoff enter the weak subspace. Threshold proximity is reported
separately. These cutoffs concern finite-precision arithmetic.

At `A = 0` the Gaussian column is included in the raw shape comparison and
excluded from the free-face projection and perturbation. The charge's
residualized norm remains meaningful on that face. A width at the upper
bound is tagged with the production bound tolerance. All A/C directions hold
widths fixed. Read the projected Width and profile spectra alongside them
when investigating coupling to width. A joint A/C/Width tangent-cone or
cross-block analysis would require additional directions.

Endpoint certification checks numerical convergence. The spectra describe
roundoff amplification and the effect of small local perturbations on
prediction. A large cancellation ratio reflects the actual opposing
contributions at the endpoint. A statistical uncertainty statement
requires additional assumptions or repeat data. The parent objective scale
normalizes residual energy and carries no noise variance interpretation.

The optional research executable uses the normal InitialStage and bounded
Joint fit, then analyzes selected atom serial IDs after qualification. Its
output is for offline analysis.

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/release --parallel 4 --target joint_ac_compensation
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1
build/release/bin/joint_ac_compensation \
  /path/to/6z6u.cif /path/to/emd_11103_unsharpened.map A 1 175 build/6z6u-ac.json 735 1032 441 500
build/release/bin/joint_ac_compensation \
  /path/to/8rqb.cif /path/to/emd_19436_unsharpened.map A 50 59 build/8rqb-ac.json 441 465
```

The executable normalizes each map by its map-value standard deviation.
For a half-map comparison, repeat the same selection and compare predictions
at identical fitting row IDs on the same grid. Multiply A, C and predictions
by the recorded map divisor before comparing values in input-map units.
The optional fixed-face perturbation uses `||D delta_beta|| = 0.01`,
reduced if needed to keep every positive amplitude interior. With fixed
widths the prediction change `X delta_beta` is exact and the objective change
includes both the residual cross term and the squared prediction change.

Each requested target uses one tiled reduction and up to two compact SVDs.
Only sparse basis data, a row tile and small compact factors are needed.
The existing Local Conditioning memory budget is checked before forming
the basis. Use selected targets for offline work, especially when many
targets share a block and the per-target calls would repeat its reduction.
