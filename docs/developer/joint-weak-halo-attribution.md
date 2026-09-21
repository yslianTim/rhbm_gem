# Weak-halo attribution within a fixed problem

Status: **complete with numerical and inference limitations**. No solver defect
was established and no solver, kernel, objective, budget or tolerance was changed.
The observations, memberships, identities and scale are the existing `weak`
partial-selection fixture, not the historical `weak-1e-4` dataset.

## Scope and evidence

The [snapshot](figures/joint-validation/weak-snapshot.json) has 2,469 rows, one
target and one halo. Its canonical JSON SHA-256 is
`aa606e91fc7530baa0637fc1f9963f9c232992ff54bd6a9f10e7c8ef2129a3ea`.
Truth is target (A, B, C) = (2, 0.5, 0.2), halo = (2.3, 0.3, 0.15).
All three starts used exactly this input. Formal public outcomes are retained
separately from the offline assessments and restarts in the
[diagnostic records](figures/joint-validation/weak-halo.json).
These final records use the corrected precise JSON wrapper; the original attempt
and its one-ULP representation issue remain indexed in the
[recovery comparison](figures/joint-validation/serialization-recovery.json).
Internal coefficients, input snapshots and runtime statuses were unchanged.

| Initial B, target / halo (Å) | Formal stop | State | Local correction infinity norm |
| --- | --- | --- | --- |
| Production 0.5887477919 / 0.3589682426 | native-lm-stop | available; convergence failed | 0.03493886484 |
| 0.5 / 0.4 | untrusted-trial | unavailable | unavailable |
| 0.6 / 0.5 | no-trustworthy-descent-step | available; convergence failed | 0.03859433435 |

An unavailable state has no endpoint certificate or endpoint scan. This is an
explicit guarded-search limitation, not a failed attempt to manufacture a state.
It is distinct from the already documented narrower-width zero-column case.
The original runner receipt has `audit_complete=false` for this start because
there is no endpoint to audit; the compact report records collection completeness
separately and preserves that original receipt and the unavailable runtime status.
The precise cause of the rejected initial trial at (0.5, 0.4) remains
**not distinguished by this endpoint-only audit**. No diagnosis from the other
two endpoints is transferred to that missing state.

## Attribution

**Search termination is not full convergence.** At the production endpoint,
target A/B/C is accurate to approximately machine precision, but halo
(A, B, C) is approximately (0.0132040, 0.448622, 0.14999999954).
Reprofiling every A/C at a 0.1 step in the normalized correction direction
lowers the normalized objective by about `1.37157e-23`. Both double and the
independent 50/100-digit profiles resolve the decrease. One bounded restart
(1,000 profile evaluations / 500 updates) lowers the objective, but the correction
remains about 0.0343188. The second available start also remains nonconverged
after its one restart (correction about 0.0380240). Increasing the budget alone
has not established convergence.

**Resolution limits also matter, but do not explain away the correction.**
The 50/100-digit endpoint calculations agree to a maximum scaled difference below
`6e-48`. Their correction norms remain approximately 0.0349388653 and 0.0385943,
so the runtime failure is not simply an inaccurate double correction.
Dense/tiled same-state comparisons pass. However, the weak-direction derivatives
have double-versus-high-precision relative differences of approximately `9.76e-6`
and `1.46e-6`; both fail the existing offline `1e-6` derivative agreement check.
These failures remain visible and no offline regular certificate is claimed.

**The halo is extremely insensitive in this observation domain.** The production
endpoint's weakest projected-width singular value is about `1.50e-10`, versus
`0.548` for the strong direction; numerical rank nevertheless passes. The same
0.1 log-width step changes the reprofiled prediction by only about `4.94e-11` in
Euclidean norm. This is evidence of practical sensitivity limitations, not proof
of exact rank deficiency, global non-identifiability, or a reliable target-only
inference rule. A low objective is compatible with a large halo parameter error.

The three explanations are therefore not exclusive: early termination of the
search, limited double derivative resolution, and very weak halo sensitivity
coexist. The verified decrease rules out the claim that *all* remaining descent
is unresolvable in double precision. The high-precision correction rules out
relabeling this endpoint as converged.

## Reproduction and boundaries

Build the opt-in offline `joint_validation` executable and run stage A of
`tests/integration/joint_validation.py` as described in the
[delivery guide](joint-capabilities-limitations.md). Each available state uses
both signs of log-width steps `1e-6` through `1e-1`, along the correction and weakest
width direction. Every point reprofiles all A/C; steps `±1e-3` and `±1e-1` also
have independent 50/100-digit reprofiled controls. High precision promotes the
same stored double observations and distances rather than regenerating data.

This is a bounded diagnostic of one geometry. It does not justify dropping halo,
fixing its coefficients, changing runtime gates, or asserting global optimality.
