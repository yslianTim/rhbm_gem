# Exact-component equivalence experiment

This testing-only experiment uses the frozen `d9f22026` observations and the [version 1 contract](joint_abc_components_contract.md). It separates structural decomposition, numerical agreement at a shared state, and the outcome of independent Guarded searches. No simulation map or first-stage estimator was regenerated.

## 1. Does the partition describe the same problem?

The contributor-CSR census has one connected component in each of the nine original datasets. The 168-atom dataset retains 139,551 rows and 407,237 memberships. Connectivity uses every saved structural membership, including zero-coefficient and underflow controls. Component views preserve squared distances and both directions of the parent/local mappings. Constant rows and unobserved atoms remain explicit.

The seven extra fixtures are disjoint compositions of saved snapshots, with separate identity namespaces and row/atom mappings. They preserve source observations and all four initial-width vectors, including `mixed`. They are algebraic fixtures; their positions are not used to invent cross-source memberships or describe a regenerated physical map.

The first gate reran 216 monolithic fits and compared 432 fit/audit records with the historical results. Scientific fields matched exactly, and regular certificates remained Legacy 38, Guarded 40, Guarded-log 40. The historical reports were retained unchanged.

The shared context owns immutable parent observations, global identity and the full observation scale. Double, 50-digit and 100-digit arithmetic compute their own scales from those observations. Same-state active-set decisions share one driver; only free-face factorization changes. SparseQR pivot evidence and SVD rank evidence remain distinct. SVD thresholds use current spectra and global matrix-family dimensions.

## 2. Do the valid numerical checks agree?

The frozen-state matrix has 256 records: sixteen datasets, two precisions, four starts, and both the initial state and saved Guarded endpoint. The validated stage-two gate passed every applicable comparison; 196 records have complete same-state equivalence evidence. Others retain active-face or rank limitations. This count describes implementation equivalence, not regular endpoint qualification.

Raw A/C/B predictions and objectives are compared before profile validity is required. Profile coefficients, gradient and objective are then compared with the shared active-set driver. Fixed-face differentiation uses one canonical A/C/B state, includes the nonzero-residual term, and compares full projected-width/Jacobian matrices. Block singular values are merged and judged at a common global threshold; repeated spectra do not require arbitrary SVD vectors to coincide.

The exploratory implementation compared derivatives at separately rounded profile coefficients and treated every face as applicable. It exposed nine difficult endpoint disagreements. Those records are retained. Fresh primary/reference face probes at the original 1e-4 and 5e-5 steps identify cross-face states as unavailable for fixed-face claims. The coefficient, prediction, objective, gradient, Jacobian and spectrum tolerances were not widened.

Global high-precision calculations retain the parent normalization. The block boundary reference includes other blocks' fixed-face/profile/path losses and constant-row loss; its KKT is global. Dense-vs-block controls cover nonzero residuals, constrained release and constants. Exact-input precision-reference reuse is confined to a single audit case, and cleared before every case or isolated rerun.

## 3. What does independent search preserve?

The source-registered regular set contains 64 paired branches: 40 originals and
24 composite branches (`regular-two`, `regular-three`, `regular-near`). Their
endpoint numerical controls all pass, with the following worst differences:

| Quantity | Observed maximum | Gate |
|---|---:|---:|
| Scaled A/C difference | 2.25891e-11 | 1e-8 |
| Scaled B difference | 2.65180e-12 | 1e-8 |
| Maximum log-B difference | 7.95542e-12 | 1e-8 |
| Prediction infinity difference / parent scale | 2.85396e-14 | 1e-8 |
| Objective difference / parent scale squared | 2.08923e-25 | 1e-12 |

Both independent runs obtain the same final global certificates:

| Dataset group | Paired branches | Monolithic regular | Assembled regular |
|---|---:|---:|---:|
| Nine original datasets | 72 | 40 | 40 |
| regular-two | 8 | 8 | 8 |
| regular-three | 8 | 8 | 8 |
| regular-near | 8 | 8 | 8 |
| regular-weak | 8 | 0 | 0 |
| regular-active | 8 | 0 | 0 |
| regular-zero | 8 | 0 | 0 |
| regular-duplicate | 8 | 0 | 0 |

All 64 required regular pairs pass; none of the 72 original Guarded scientific
records regresses. The 32 difficult composite branches preserve weak-signal,
active-face, unidentifiable or rank-failure limitations. Search completion alone
does not establish a regular certificate. The [failure matrix](figures/joint-abc-components/failure-matrix.csv)
records each limitation and each component's usable-state status.

The runner additionally checks both actual search endpoints for all 128 pairs:
256 endpoint same-state records pass every applicable comparison, with 160
complete equivalence records. These are separate from the frozen-state matrix
and from final endpoint certification.

A review of the first complete audit found that composite final certification
inherited the search diagnostics' initial weak direction. For
`regular-weak/first-stage-double`, its absolute cosine with the endpoint weak
direction was only 0.03524. This could miss unresolved derivative evidence. The
raw searches were retained. Composite legacy checks were refreshed, and the
final global/component audits and boundary scans were recomputed independently
in both runs. Original single-component legacy checks already used endpoint
directions and retain their independently computed evidence. Final certification uses each
global endpoint's weakest direction; supplemental component audits restrict the
assembled endpoint's directions without renormalizing them. Scoped fit records
refresh the legacy derivative flag before certification. No optimizer setting,
endpoint, numerical kernel or acceptance tolerance changed.

The correction removes six monolithic and seven assembled composite weak-signal
qualifications from the initial audit counts. All 896 scoped state records retain
their primary parameters, objective, KKT and gradient fields. The scoped fit
preserves raw search metadata; its adjacent `context.json` is the authoritative
final audit scope. The [audit revision](figures/joint-abc-components/audit-revision.json)
records these checks. The initial-direction audit records are archived as
superseded evidence, not used in the final qualification counts. A runner test rejects initial-direction
substitution and rescaling of restricted directions.

Each component receives its paired source start and its own 200-evaluation/100-update budget. Local search rank dimensions are recorded alongside the parent scale and parent active-set response norm. Search completion, component certificate and assembled global certificate are separate records.

The assembled estimator retains the actual component A/C/B values. Global KKT, derivatives, rank, Richardson/precision evidence and local correction assess that state. A fresh global profile is only a consistency control. Failed components retain their last trusted state; missing states produce an unavailable global prediction/objective and an explicit row mask. A failed block is never filled with zeros.

Every healthy baseline block completes all eight paired starts. In
`regular-duplicate`, six duplicate-block searches have no trusted state, while
the baseline block still completes. The other two duplicate starts retain usable
states without establishing regular global qualification. Missing blocks leave
an explicit row mask and unavailable global prediction/objective.

Four standalone reruns reproduce the raw search fit, scoped fit, audit and
context exactly after timing/RSS exclusions: `baseline/first-stage-double`, the
failed and healthy blocks of `regular-duplicate/first-stage-double` (failed block
run first), and `regular-near/narrower-double`. Unit controls also reverse the
execution order within one process, and cover invalid starts and exhausted
budgets. An unobserved-atom control preserves the observed block's successful
search and marks the assembled state unavailable.

## Reproducibility and cost

The [full repeat comparison](figures/joint-abc-components/repeat-comparison.json)
finds all 3,520 scientific JSON records identical, including failures and freshly
recomputed audits. Elapsed time, process peak RSS and execution paths/logs are
excluded. Each run performs 128 monolithic and 192 component searches. This is a
same-host reproducibility result, not a cross-platform bitwise guarantee.

Each fresh run contains immutable input, implementation source snapshots, source/executable/library hashes, full trials, contexts, spectra, face probes, audits and boundary scans. Numeric kernels use one thread. Searches use three independent dataset
processes per run; refreshed audits use three independent case processes per
run. Each case starts with an empty precision cache. No cache crosses runs. Reference computations and exact-input reuses are recorded separately.

The profile-evaluation counts below sum eight paired starts per dataset. They
exclude independent reference solves and endpoint audits, which are listed
separately in the per-scope cost tables.

| Composite | Monolithic profile evaluations | Component total |
|---|---:|---:|
| regular-two | 72 | 144 |
| regular-three | 78 | 190 |
| regular-near | 130 | 175 |
| regular-weak | 78 | 136 |
| regular-active | 324 | 399 |
| regular-zero | 129 | 184 |
| regular-duplicate | 37 | 101 |

Per-component budgets are not equal aggregate budgets to one monolithic search. Costs distinguish search, endpoint assessment, audit, nested reference/precision work, and process wall time. Global aggregate counters are not added to child counters. Process peak RSS includes earlier work in the same process and is not a component allocation peak.

Across all sixteen datasets, each run uses 1,856 monolithic and 2,337 component
search profile evaluations. Context construction is separate: 56 primary and
656 reference profile calls. Search reference calls (4,203), endpoint reference
calls (5,380), post-search trust references (180), refreshed assessment references
(2,782), and final audit references (26,524) are listed separately. High-precision
reference computations and case-local reuses also have distinct counters.

| Measured quantity | Run A | Run B |
|---|---:|---:|
| Context setup time (s) | 15.32 | 15.31 |
| Search time, summed monolithic and child scopes (s) | 888.29 | 891.34 |
| Search endpoint assessment time (s) | 907.46 | 905.36 |
| Refreshed endpoint assessment time (s) | 35.76 | 35.79 |
| Final audit time, summed scopes (s) | 4,008.89 | 4,008.66 |
| Boundary scan time, summed scopes (s) | 3,788.26 | 3,788.00 |
| Final audit elapsed wall time (s) | 2,628.05 | 2,626.99 |
| Maximum per-process peak RSS (decimal GB) | 1.256 | 1.266 |

These measurements include overlapping independent processes and runs. They are
cost accounting, not a speed benchmark. The [cost summary](figures/joint-abc-components/cost-summary.json)
and per-scope tables retain all raw values, computation/reuse counts and timing
definitions. Native process totals also include same-state checks, loading,
output and other overhead; nested timers are not added again.

The original 168-atom problem remains one structural component. This experiment establishes exact decomposition and qualification behavior; its census does not imply a speedup or smaller solve for that original problem.

## Validation and records

The [scientific validation record](figures/joint-abc-components/scientific-validation.json)
collects all three gates, both full runs and the four isolated reruns. The
[record inventory](figures/joint-abc-components/README.md) links the compact tables,
complete archives, source/input/build provenance and artifact hashes. Failed
exploratory checks and superseded audits are preserved separately from final
evidence. The [contract](joint_abc_components_contract.md#runner-operations)
documents `run`, `audit`, `summarize`, `compare` and `rerun-component`.

Production build and repository lint pass. The production shared library contains no joint-component or multiprecision symbols. The new estimator, audit support and runners remain confined to testing targets.
The 72 affected C++ tests include finite but rank-deficient raw states, a nonzero
residual Jacobian, constrained blocking/release, global/local rank controls,
constant/unobserved rows, and reversed component execution with failure/budget
controls. Nine Python runner suites pass, including hash corruption, frozen mixed
starts, exact local mappings, parent-scale replay and final audit scope.
