# Joint observable-halo parameterization

This work starts from `d634fa0d`. It addresses the structural obstruction in the
[6Z6U initial-rank diagnosis](joint-initial-rank-diagnostic.md), without changing
the observation domain, support geometry, contributor closure or component IDs.

## Work packages 1–2: layout and profiled solve

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

## Work package 3: result consumers

* Second-stage points exist only for FullABC. Reduced halos report
  `observable-contribution-only` with no point or covariance.
* Peeling consumes the same joint result, adds each lambda once and retains the
  interpolation stencil's domain-coverage check. It does not extrapolate a halo
  shape outside the observed row.
* Uncertainty uses the FullABC Jacobian on informative rows. At an interior,
  full-rank endpoint its residual degrees of freedom are N-q-3p. Existing
  boundary, rank, variance and availability checks remain in force.
* New JSON is v4, with `singleton-halo-profile-v1` or `full-abc-v1` as its
  parameterization contract. JSON v3 is read using its original full-ABC
  semantics; it does not acquire inferred reduction or stronger evidence.
  SQLite remains schema v19 and stores the JSON payload.
* Atom CSV retains all contributors. Reduced A/B/C cells are empty with a group
  ID. The sibling `.contributions.csv` has exactly one row per nuisance group;
  unreduced results produce a header-only companion.

`JointObservableProfileTest` covers grouping, selection guards, nonzero-residual
finite differences, explicit nuisance least squares, full reconstruction,
unobserved targets, seed provenance, partial component failure, augmented
Jacobian covariance, v3/v4 decoding, SQLite, CSV, peeling and initialization skip.
Existing frozen, partial-selection and CLI regressions remain enabled.

## Numerical operations and experiment conclusions

Endpoint assessment and uncertainty request only the compact SVD outputs they
need. Tiled QR uses its temporary matrices in place. Rank policies, Jacobi
failure/rank-boundary retry, independent reference and cancellation behavior are
preserved. Operator LM, normal equations, spectrum truncation and relaxed solver
budgets are outside this change.

The [experiment conclusions](joint-observable-acceptance.md) retain the completed
measurements and remaining 6Z6U limitations as text. Experiment-only runners,
diagnostic command extensions, extra timing fields, raw evidence and temporary
builds have been removed. Production functionality and permanent small regression
tests remain.
