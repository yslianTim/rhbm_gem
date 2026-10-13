# Joint target identifiability

This guide defines the current target evidence, target runtime status, and its
separate covariance behavior.

## Derived Target Identifiability

Production no longer builds a target-specific `N × 3P` Jacobian or computes its
QR, SVD, right null space, or selector leakage. It derives target evidence from
the already-qualified component rank records. For

`y_hat = X(eta) beta`, `beta = (A1, C1, ..., AP, CP)` and `eta_i = log(B_i)`,
let `D = d y_hat / d eta` and `D_perp = (I - P_X) D`. When every A/C parameter
is free and interior, the exact rank identity is

`rank([X, D]) = rank(X) + rank(D_perp)`.

Therefore, if `rank(X) = 2P` and `rank(D_perp) = P`, then the complete
`(A, C, log B)` Jacobian has rank `3P` and a zero right null space. Every
FullABC selector, including a FullABC halo atom, is locally identifiable under
these conditions. The result is a first-order local statement, not a guarantee
of global uniqueness or practical precision.

The derived certificate requires a trusted component state; passing numerical
qualification; complete Design Rank and Projected-width Rank evidence; matching
state, parameter layout, observation rows, selection domain, scope and numerical
rank context; and no unresolved rank-threshold sensitivity. Rank records carry a
context identity so evidence for a different state or row domain cannot be
reused. Every amplitude participating in the implication must be strictly
interior (`A_i > 0`). An active `A_i = 0` boundary remains `Unavailable`.

Missing evidence and rank deficiency do not prove that every individual target
is unidentifiable. They produce `Unavailable` target evidence, while the
component's existing `Failed` rank qualification remains intact. No selection
metadata keeps `TargetRuntimeConvergence` at `NotRun`. Observable-contribution-
only nuisance groups have no independent A/C/B coordinates and remain
`Unavailable` with reason `observable-contribution-only`; they are distinct
from FullABC halo atoms.

Per-atom evidence reports `Passed` with source
`derived-from-component-identifiability`. Since no null-space calculation is
performed, `null_space_leakage` and target `column_scales` are absent. The
certificate uses the component Design Rank and Projected-width Rank
normalizations and thresholds. That numerical-rank contract is not equivalent
to the retired full-Jacobian, column-scaled SVD contract: near-threshold results
may differ. The current evidence contract is `target-estimability-v3`; tests
retain a small direct-Jacobian oracle to make this difference visible.

## Target convergence and trusted state

Target convergence remains separate from component runtime convergence,
assembly consistency, and trusted-state availability. Component target status
requires a trusted state and passing conditional A/C optimality, independent
replay, width stationarity, Design Rank, Projected-width Rank, numerical
qualification, and derived target-identifiability evidence. The assembled
result additionally requires a complete assembled state and prediction, all
observation rows, every component target result, and passing assembly mapping
and prediction/objective consistency. Search completion is reported
separately. Target evidence assesses the saved state; it does not reprofile or
replace it.

For an exact single-component match, the already-qualified component evidence
is reused only when row identity, parameter layout, saved state, and rank context
match the assembled result. Multiple components retain their local and global
evidence scopes; one component's ranks cannot certify another component.

## Covariance and point estimates

Identifiability does not provide `(JᵀJ)⁻¹` or a parameter covariance.
Covariance is computed independently by the existing tiled full-component
linearization path, with its own `TargetCovariance` memory preflight and
full-rank, interior-amplitude, residual-variance, and degrees-of-freedom gates.
Small problems can retain covariance. A preflight rejection or any numerical
gate returns `Unavailable`; rank evidence is never used to synthesize
covariance. This covariance calculation is not part of target
identifiability. An unavailable covariance leaves an available point estimate
intact, while covariance-dependent group evidence is `Unavailable`.

Observable-contribution-only nuisance atoms have no complete point estimate or
covariance. A FullABC halo point is available when its derived per-atom
identifiability evidence passes; its covariance still follows the independent
covariance gates.

## Persistence and limits

Joint JSON schema 9 stores `target-estimability-v3` and the rank context
identity. It stores no fabricated null-space leakage, full-Jacobian singular
values, or column scales. Schemas 3 through 8 remain readable. Historical
`target-estimability-v1` and `target-estimability-v2` evidence keeps its saved
meaning and statuses; loading never converts an old failed geometry result into
a derived pass or recomputes runtime or covariance status.

The rank derivation core is `O(P + T)` for `P` FullABC atoms and `T` targets.
The persisted per-atom evidence also records each observable-contribution-only
nuisance atom as `Unavailable`; if there are `U` such atoms, total additional
evidence storage is `O(P + T + U)`. Resident input and mapping storage still
contribute to its estimate. DesignSpectrum and Projected-width remain the
large-component rank bottlenecks and can still be rejected by their own
preflights. A target preflight pass does not imply that either rank stage ran or
that target evidence passed.

Observable-profile and partial-selection tests cover target selectors,
FullABC halo mapping, nuisance-only groups, rank-context mismatches, active
boundaries, missing and deficient rank evidence, threshold sensitivity,
covariance, and historical schema compatibility. Synthetic tests do not
establish qualification on a physical cryo-EM model or map.
