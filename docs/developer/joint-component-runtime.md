# Joint Component Estimator runtime

Joint Components fits selected non-hydrogen Targets. Read
[Numerics](joint-component-numerics.md) for the block solver and
[Stage Integration](joint-stage-integration.md) for result publication,
postprocessing and persistence.

## Model and objective

For distance `r` in Å and physical width `tau`, a Target contributes

```text
f_i(r) = A_i G(r; tau_i) + C_i K(r; tau_i)
G(r; tau) = (2 pi tau²)^(-3/2) exp(-r² / (2 tau²))
K(r; tau) = erf(r / (sqrt(2) tau)) / r
```

A is nonnegative Gaussian amplitude. Signed C multiplies the charge kernel.
Both terms share `0 < tau ≤ tau_max` and the fixed 2.5 Å support. For
`r < 1e-5 Å` the charge basis uses its center value `sqrt(2/pi)/tau`.
The solver uses `eta = log(tau)`. Inside support,
`dG/deta = G (r²/tau² - 3)` and
`dK/deta = -sqrt(2/pi)/tau * exp(-r²/(2 tau²))`, with the center branch
using `-sqrt(2/pi)/tau`. Thus `dF/deta = tau * dF/dtau`.

For observed voxel value `y_v`, the Target-only, equal-weight objective is

```text
F = sum_{v in V_fit} (sum_i f_i(r_vi) - y_v)² / (2 s²)
s = max(1, ||y_original||₂)
```

The Original observations own the fixed parent scale, including rows later
excluded from fitting. All components and block comparisons use that scale.
Map intensity normalization is a separate operation. A and C retain their
kernel coefficient units. The objective has no prior or regularization.

## Original and Fitting Domains

The Original Domain is the union of map voxels within the inclusive 2.5 Å
support of the selected Targets, deduplicated by global grid index. It records
raw voxel values, atom identities and squared-distance support incidences.
An unselected non-hydrogen atom is a Halo if it has support on these rows.

Pure Halo-Support Exclusion removes every Original row with recorded Halo
support, including the `r² = 6.25 Å²` boundary. The retained observations form
the Fitting Domain with compact Target and row indices mapped back to the
Original Domain. Fitting uses the retained raw `y_v` directly. Halo supplies
the geometric mask and has no Joint coefficients or gradients. The builder
rejects an empty Fitting Domain or a Target that loses all support. The
Original Domain remains available for postprocessing.

## Structural components

Every fitting row connects all its recorded contributors, independently of
parameter values or numerical sensitivity. The connected components preserve
stable atom identities and parent row mappings. Compact component views
share immutable parent input. Rows without contributors are constant terms
`y_v² / (2 s²)`, included once in assembled loss. The partition represents
unsupported atoms explicitly and supplies no invented coefficient or profile
for them. Invalid or duplicate identities and malformed mappings are rejected.

Each component fits independently. Assembly maps qualified component states
back by identity, preserves their coefficients and replays the assembled state
on the parent domain. Failed components stay isolated and their predictions
remain unavailable. There is no assembly-wide A/C reprofile.

## Initialization and global width bound

The workflow obtains Target width seeds through
[Joint FirstStage](joint-stage-integration.md#firststage-initialization).
Calibration uses finite positive native-qualified or recovered-qualified
Target MDPDE widths, excluding median fallbacks and initialization exceptions.
It requires at least 20 qualified donors and at least 95% Target coverage.
The sole upper bound is `tau_max = 4 × Q0.95(donor widths)`, using
Hyndman–Fan Type 7 quantiles. All components receive the same
`eta_max = log(tau_max)`. Seeds above the bound are projected for the working
initialization, with original and used widths retained. An invalid, overflowing
or unavailable bound prevents Joint search.

## Gauss–Seidel BCD

Each component visits disjoint structural cores of at most 12 atoms in stable
forward order. A core's affected rows retain every contributor across its
boundary. Holding the other Targets at their latest accepted values gives
`y_eff = y - prediction_other` on those rows. A visit profiles A/C at fixed
widths and allows at most one trusted LM width update. The trusted initial
profile can itself propose new A/C, provided the visit has no width failure.

Candidates use the accepted state's unchanged rows. On affected rows the
candidate recomputes prediction and residual using all contributors. With
`R_b` the affected rows, the increment in the cached objective is
`sum_{v in R_b}(r_new,v² - r_old,v²) / (2 s²)`. A preflight-approved
full candidate replay is the fallback if the incremental workspace exceeds
budget. Candidate acceptance requires a qualified local profile, finite
feasible state, nonincreasing local objective and agreement of the global
delta with the local delta. The roundoff enclosure is
`1e-12 + 2e-12 |reference objective|`. A positive global delta within
that enclosure is allowed only when the local objective did not increase.
Only an accepted candidate changes the shared coefficients, widths and
affected prediction/residual cache before the next visit.

A failed width solve or missing trusted profile leaves accepted state intact
and schedules one retry on the next sweep. Recovery clears pending failure
only after a successful visit. A second failure ends search with
`block-recovery-exhausted`. A safely rejected candidate or a recovered
earlier failure does not by itself disqualify a later verified endpoint.

## Global convergence and qualification

Every complete sweep performs Full Replay from accepted parameters and
compares the resulting prediction and objective with the incremental cache.
Prediction uses the enclosure `2e-12 + 2e-13 max(1, ||prediction||∞)`;
objective uses the block replay enclosure. A failed audit stops search.
Only after the audit and a valid global state evaluation does the replay
objective reconcile the cache.

Global KKT is evaluated on the component's accepted state using all its
fitting rows. A/C uses the projected-gradient formula in
[Numerics](joint-component-numerics.md#constrained-ac-least-squares), with
component-wide design norms and the parent scale. Width uses
`||eta - min(eta - dF/deta, eta_max)||∞` with feasibility and bound evidence.
At an active upper face a nonpositive gradient satisfies the one-sided
condition. The bound check uses `64 * epsilon * max(1, |eta_max|)` for
feasibility and face classification, and reports multiplier and complementarity
evidence. The raw gradient remains a diagnostic. Production `practical-v2`
requires A/C KKT ≤ `1e-7` and projected log-width KKT ≤ `1e-9`.
The `strict-v2` thresholds are `1e-10` and `1e-12`, respectively.
Confirmation needs at least two completed sweeps and a confirming sweep
without block failures.

The budget is 200 sweeps with a hard cap of 300. At the nominal limit,
the last 20 values of `q = max(AC_KKT/1e-7, width_KKT/1e-9)` must all be
finite and nonnegative, with no pending recovery. Extend only when the best
of the latest ten is at most 10 and at most 95% of the earlier ten's best,
which must exceed 1. The extension continues the same state and policy.
Budget exhaustion can report `block-sweep-budget`,
`block-hard-sweep-budget`, or `block-recovery-sweep-budget`.

Endpoint Verification evaluates the accepted coefficients and widths.
Independent Replay separately traverses recorded supports to verify
prediction, residual, objective, feasibility and both KKT conditions.
It checks rowwise prediction and residual within
`2e-12 + 2e-13 |primary prediction|`, the unnormalized half-RSS within
`1e-12 + 2e-12 |primary half-RSS|`, A/C KKT agreement within `1e-13`,
and each log-width gradient within `1e-13 + 2e-9 |primary gradient|`.
Agreement and stationarity are separate requirements. When stationarity
remains confirmed, repeat attempts are spaced ten sweeps apart. A terminal
sweep can receive a final check if recovery is clear and the sweep has not
already been certified. The assembled endpoint is verified again on the
parent Fitting Domain. Local LM termination, a budget stop or conditioning
diagnostics alone cannot qualify an endpoint.

## Resources and result contract

`ResourceWork` uses checked arithmetic to estimate resident and peak
workspace for problem and structural construction, BCD preparation, global
streaming evaluation, endpoint verification and Local Conditioning. It
accounts for support incidence, affected rows, mapping storage and largest
block workspace. Production preflights enforce the 3,200,000,000-byte
planning budget. The header also declares a 4,000,000,000-byte hard estimate
constant, currently unused by the planner. Observed RSS is a separate
measurement. Overflow or over-budget stages can fail
preflight before allocation. Global evaluation streams support entries
without global design factors or global covariance. Local factors and
derivative tiles are bounded by the current block and planned workspace.

A qualified result requires valid initialization, completed search, accepted
state integrity, qualified component and assembled evidence, consistent
mappings, valid Target parameters and available fitting rows. The captured
snapshot maps identities and availability back to the Original Domain.
An incomplete or unqualified endpoint can retain candidate diagnostics and
explicit reasons. Publication and
[Point Quality](joint-stage-integration.md#target-point-quality) use separate rules.

## Implementation references

- [JointComponentEstimator.cpp](../../src/core/detail/joint_component/JointComponentEstimator.cpp),
  [ComponentPartition.cpp](../../src/core/detail/joint_component/ComponentPartition.cpp),
  [StructuralPartition.cpp](../../src/core/detail/joint_component/StructuralPartition.cpp)
  and [WidthCalibration.cpp](../../src/core/detail/joint_component/WidthCalibration.cpp) build
  the domains, components, block topology and bound.
- [ComponentOrchestration.cpp](../../src/core/detail/joint_component/ComponentOrchestration.cpp),
  [ComponentAssembly.cpp](../../src/core/detail/joint_component/ComponentAssembly.cpp),
  [GaussSeidelBcd.cpp](../../src/core/detail/joint_component/GaussSeidelBcd.cpp),
  [BcdPolicy.hpp](../../src/core/detail/joint_component/BcdPolicy.hpp) and
  [ComponentSweepBudget.hpp](../../src/core/detail/joint_component/ComponentSweepBudget.hpp)
  own execution, assembly, production settings and budget decisions.
- [Numerics.cpp](../../src/core/detail/joint_component/Numerics.cpp),
  [BoundedLmDriver.hpp](../../src/core/detail/joint_component/BoundedLmDriver.hpp)
  and [ResourceWork.cpp](../../src/core/detail/joint_component/ResourceWork.cpp)
  implement evaluation and replay, the local LM loop and memory preflight.
- [Component tests](../../tests/core/joint_component) cover domain exclusion,
  calibration, connectivity, acceptance, recovery, replay and qualification.
