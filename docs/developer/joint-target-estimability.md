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
reused. Production obtains Design Rank from the same normalized free-design
QR/SVD used by projected-width preparation when every amplitude is interior, so
the free design has all `2P` columns and equals the full A/C design. If any
`A_i = 0`, the factor has fewer columns and remains evidence about that free
face only; it cannot certify full Design Rank. That full-design evidence is
`Unavailable`, and target derivation remains `Unavailable` at the active
boundary.

Missing evidence and rank deficiency do not prove that every individual target
is unidentifiable. They produce `Unavailable` target evidence, while the
component's existing `Failed` rank qualification remains intact. No selection
metadata keeps `TargetRuntimeConvergence` at `NotRun`. Current production
problems contain selected targets and retained FullABC halos only. Excluded
halos have no contributor index, per-atom evidence, point estimate, or
covariance. The eligibility gate does not imply target accuracy under the
changed additive model.

Per-atom evidence reports `Passed` with source
`derived-from-component-identifiability`. Since no null-space calculation is
performed, `null_space_leakage` and target `column_scales` are absent. The
certificate uses the component Design Rank and Projected-width Rank
normalizations and thresholds. Both preserve the original observation-row
count and use a strict `singular_value > threshold` decision; values at the
threshold do not pass, and a value within one representable step of the
threshold is marked sensitive and cannot qualify. The contract is
`target-estimability-v4`, which records the new source and active-boundary
availability semantics. Tests retain a small direct-Jacobian oracle to make
the difference from the retired full-Jacobian SVD visible.

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

Assembled target evidence is derived from the assembled state and its global
rank context. Component-local evidence cannot certify a different assembled
state or row domain.

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

A retained FullABC halo point is available when its derived per-atom
identifiability evidence passes; its covariance still follows the independent
covariance gates.

## Persistence and limits

Joint JSON schema 11 stores `target-estimability-v4` evidence, rank context
identity, and halo-pruning provenance. It stores no fabricated null-space
leakage, full-Jacobian singular values, or column scales. Schemas 3 through 10 remain readable;
historical v1, v2 and v3 evidence keeps its saved meaning and statuses. Loading
never converts missing or failed historical rank evidence into a pass or
recomputes runtime or covariance status.

The rank derivation core is `O(P + T)` for `P` retained FullABC atoms and `T`
targets. Exclusion provenance adds `O(E)` records for `E` excluded halos and is
included in problem memory preflight. Resident input and mapping storage still
contribute to its estimate. Projected-width preparation is the remaining
large-component rank bottleneck: it includes the free-design QR/SVD, projection
coefficients and projected-width QR/SVD, and can be rejected by its 3.2 GB
planning-budget preflight. The 4 GB hard limit is a resource ceiling, not proof
that full qualification fits. A target preflight pass does not imply that the
rank stage ran or that target evidence passed.

Partial-selection and target-evidence tests cover target selectors, retained
FullABC halo mapping, rank-context mismatches, active boundaries, missing and
deficient rank evidence, threshold sensitivity, covariance, and historical
schema compatibility. Synthetic tests do not
establish qualification on a physical cryo-EM model or map.
