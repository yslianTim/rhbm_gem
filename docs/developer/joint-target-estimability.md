# Joint target estimability

This guide owns current target-specific convergence, identifiability, and
covariance behavior. The historical full-case evidence and unresolved
large-scale boundary are indexed in the
[canonical evidence document](joint-component-evidence.md).

## Target convergence

Target qualification is separate from full-parameter runtime convergence.
Inputs without recorded selection metadata have target convergence
NotRun; roles are not inferred. Component target convergence requires a saved,
fully qualified state and passing conditional A/C optimality, independent
replay, width stationarity, and target-identifiability checks. Endpoint
verification alone does not prove identifiable parameters. Target local
correction is not part of the current qualification contract.
The assembled result additionally requires a complete assembled state and
prediction, all observation rows, every component target result, and assembled
mapping plus independent prediction/objective consistency. Search completion
remains a separate status.

The target certificate assesses the saved state. It does not reprofile to
replace that state. Any active A=0 in the FullABC parameterization makes target
geometry unavailable. Inconsistent saved widths, invalid evaluation, failed
KKT, or compact factorization failure also produce unavailable evidence.

## Observation-Jacobian identifiability

For each FullABC atom, form the complete observation Jacobian in
(A, C, log B) order after observable-contribution-only rows have been
analytically profiled. Scale each column by its saved Euclidean norm; an exact
zero column remains zero and uses scale 1. Tiled QR reduces the observation
matrix, then the compact right singular vectors are checked using the existing
rank policy and original parent row context. The compact-SVD owner documents
the decomposition dispatch and retry behavior.

For atom i, target identifiability asks whether its three selector coordinates
have any component in the right null space of the full coupled Jacobian. The
Euclidean norm of those three null-space rows is the selector leakage. A target
passes at leakage <= 1e-10. Leakage above that limit fails as
nonunique-parameter-directions. Singular values within [0.5 tau, 2 tau] make
the rank boundary sensitive and the target result unavailable. This is a local
interior identifiability statement, not a global uniqueness guarantee.

Observable-contribution-only halo groups are recorded as unavailable
representatives. A FullABC halo point is withheld from the stage adapter and
CSV when its parameter coordinates are not identified. Target point states
remain available with their separate convergence and identifiability statuses.

## Numerical evidence and persistence

The target certificate is target-estimability-v2. It adds no parameter-step
correction solve: target-specific qualification is determined by target
geometry and null-space leakage, together with the saved optimization and
endpoint evidence. Missing or threshold-sensitive geometry remains
`Unavailable`; failed target identifiability remains `Failed`. Neither state is
promoted to `Passed` because another qualification check succeeds.

## Target covariance and exports

Target covariance uses the identifiable subspace of the full coupled
(A, C, log B) observation Jacobian. It is available only when the saved state
exists, every FullABC A is interior and positive, target runtime convergence
passes, rank is not threshold-sensitive, informative rows exceed effective
rank, and residual variance is positive and finite. Degrees of freedom are
informative rows minus effective rank. Only identified requested target blocks
receive covariance; boundary, missing state, rank, variance, and degrees-of-
freedom failures retain explicit unavailable reasons. Halo covariance is not
inferred by selecting an arbitrary representative.

Joint JSON schema 5 introduced target evidence and target-runtime status; schema
6 adds optional candidate diagnostics, and schema 7 records reference-free
endpoint and qualification statuses separately. Schema 8 versions the
correction-free solver contract. Readers accept schemas 3 through 8; schemas
3–7 retain their recorded checks and convergence statuses. Historical failed
target-local-correction evidence remains failed when loaded and is not
recomputed. Schema 3/4 records retain their original evidence and do not gain
target certification or inferred candidates. CSV records target/halo/not-recorded
selection role, target convergence, and parameter identifiability. It preserves
target point states independently of convergence and leaves unidentified halo
parameters blank. The command guide owns the full JSON/CSV field list and unit
conventions.

## Limits and permanent coverage

The target calculation retains dense parameter-space compact factors and SVDs.
Complete large-scale 6Z6U target-estimability scalability has not yet been
established; preparation results and small-case correctness do not establish
that resource bound. Historical runs and their exact scope are recorded in the
[canonical evidence index](joint-component-evidence.md).

Permanent observable-profile and partial-selection tests cover target selectors,
halo confounding, threshold-sensitive evidence, active boundaries, covariance,
schema compatibility, CSV, and unavailable outcomes. These checks do not stand
in for a completed large-case run.
