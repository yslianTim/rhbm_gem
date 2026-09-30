# Joint observable-halo parameterization

This parameterization preserves the observation domain, support geometry,
contributor closure, and component identities.

## Structural motivation

The historical full-ABC 6Z6U formulation had 4,334 free A/C columns. Twenty
singleton halo atoms each touched one selected voxel, so each atom's A and C
columns were proportional. The resulting structural rank bound was 4,314;
EIGEN, SPQR, multiple legal width starts, and an independent QR/SVD control
reproduced it. This was a support/parameterization obstruction, not an
initialization failure. The full evidence and scope are in the
[canonical historical index](joint-component-evidence.md#41-historical-full-abc-6z6u-obstruction).

Singleton-halo profiling removes that diagnosed A/C obstruction while retaining
the same observations. It does not establish target-width identifiability or a
complete 6Z6U solve, and it is not evidence that the full target estimation is
converged.

## Current parameterization contract

`JointProblem::ParameterLayout()` is fixed before initialization. Only a recorded
halo with exactly one original support row becomes contribution-only. Selected
atoms, multi-row halos and inputs without selection metadata retain FullABC.
Halos sharing a row share one unrestricted amplitude; there is no per-atom
division of that amplitude. Original atom and row indices survive in every
layout, including component layouts.
Lambda has the same fitted-map units as the observations; multiply it by the
saved normalization divisor when converting back to input-map units.

For nuisance rows U, the internal solve omits their memberships and uses
`lambda[r] = y[r] - f[r]`. The objective on remaining rows is exactly the
profile objective. The immutable input and squared distances are unchanged.
The parent observation scale and release normalization survive; rank thresholds
use the original component (or assembled problem) row count and remaining
parameter count. Components are not split again. FullABC atoms with no remaining
support, deficient problems and failed searches remain failures. Nuisance-only
components are evaluated analytically; one failed component prevents complete
assembly without erasing other available component states.

The A/C/B arrays now follow `layout.full_atoms`; nuisance amplitudes follow
`layout.groups`. With no reduced halo, the existing numerical path is used.
`FitJointComponents(problem, initial_b)` still accepts one entry per original
atom and does not infer donors for caller-provided widths. Reduced positions
need not contain a valid width.

The shared first-stage workflow skips reduced halos. Its fixed donor pool
contains only FullABC records whose original reason is `valid-width` and whose
width is finite and positive. Invalid seeds receive the median (middle-pair
arithmetic mean for even counts), when donors exist. Original reason and width
are preserved alongside used width, source and donor count; fallback does not
turn the original first-stage fit into a success or create new donors.

## Result-consumer semantics

* Second-stage points exist only for FullABC. Reduced halos report
  `observable-contribution-only` with no point or covariance.
* Peeling consumes the same joint result, adds each lambda once and retains the
  interpolation stencil's domain-coverage check. It does not extrapolate a halo
  shape outside the observed row.
* Uncertainty uses the FullABC Jacobian on informative rows. At an interior,
  full-rank endpoint its residual degrees of freedom are N-q-3p. Existing
  boundary, rank, variance and availability checks remain in force.
* New Joint JSON is schema 5, with `singleton-halo-profile-v1` or
  `full-abc-v1` as its parameterization contract. Schemas 3 and 4 retain their
  original saved semantics; they do not acquire inferred reduction, target
  evidence, or stronger qualification. SQLite schema v19 stores the Joint JSON
  payload.
* Atom CSV retains all contributors. Reduced A/B/C cells are empty with a group
  ID. The sibling `.contributions.csv` has exactly one row per nuisance group;
  unreduced results produce a header-only companion.

`JointObservableProfileTest` covers grouping, selection guards, nonzero-residual
finite differences, explicit nuisance least squares, full reconstruction,
unobserved targets, seed provenance, partial component failure, augmented
Jacobian covariance, v3/v4/v5 decoding, SQLite, CSV, peeling and initialization
skip.
Existing frozen, partial-selection and CLI regressions remain enabled.
Target selector, target convergence, and target covariance semantics are owned
by the [target-estimability guide](joint-target-estimability.md).

## Numerical boundaries

Endpoint assessment and uncertainty request only the compact SVD outputs they
need. Tiled QR reuses temporary matrices in place. Rank policies, Jacobi retry,
independent reference, and cancellation behavior remain unchanged. Operator LM,
normal equations, spectrum truncation, and relaxed solver budgets are outside
this contract.

The [canonical evidence index](joint-component-evidence.md) records the
historical campaign measurements, incomplete 6Z6U run, and artifacts that are
absent. Permanent small regression coverage is in
`ObservableProfile_test.cpp`.
