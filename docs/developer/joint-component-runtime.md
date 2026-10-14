# Joint component estimator and regression workflow

Include `<rhbm_gem/core/JointComponentEstimator.hpp>` and link
`RHBM_GEM::rhbm_gem`. The installed API needs no testing or offline audit target.

```cpp
model.SelectAllAtoms();
auto result = rhbm_gem::core::EstimateJointComponents(map, model);
auto problem = rhbm_gem::core::BuildJointProblem(map, model);
auto repeated = rhbm_gem::core::FitJointComponents(problem, initial_b);
```

## Runtime contract

`JointProblem` owns an immutable snapshot of observations, identities and
squared-distance memberships. V1 uses equal-weight least squares over its fixed
observation domain, with A nonnegative, C signed, and log-B widths. The public
objective is `0.5 * ||prediction - observations||² / ObservationScale()²`.

Production fitting uses the FixedNeighbor route described below. Each accepted
local Profile LM update is replayed against the global state and passes the
roundoff-safe monotone acceptance rule before it is kept. At the end of each
complete sweep, global stationarity checks identify a candidate endpoint.
Certification evaluates that accepted A/C/B state with a primary streaming pass
and an independent scalar streaming replay. It does not solve a second global
A/C problem. The replay compares prediction, residual, objective, projected
KKT and width-gradient evidence at the same parameters. Constrained KKT checks
fixed-width A/C optimality; replay agreement checks numerical consistency; full
design rank and projected-width rank qualify numerical identifiability.
Endpoint verification is distinct from full numerical qualification: a
verified endpoint with failed but potentially recoverable numerical checks
keeps its candidate and may continue through the same block-coordinate
updates. While confirmed stationarity remains intact, certification retries
remain bounded to one per 10 sweeps; losing stationarity resets eligibility.
A rank or other assessment resource limit stops further qualification retries
for that unchanged problem shape and preserves the candidate as uncertified.
Invalid primary numerical state also stops the search. At budget exhaustion,
if the final sweep has not been assessed, one bounded assessment evaluates the
final accepted state. It does not change parameters or add a sweep. That
assessment only completes the search when the final sweep is confirmed
stationary and fully qualified; budget exhaustion by itself is not convergence.
Search completion, endpoint verification, full numerical qualification, state
availability and runtime convergence remain separate outcomes.

`AssessProfileSearchForDiagnostics()` is a test-support diagnostic for Profile
search history and trusted historical states. Production components use
`CertifyFixedNeighborEndpoint()` for independent endpoint certification. Both
paths use the shared evaluation, replay and assessment primitives.

### Production FixedNeighbor route

The sole production route for `JOINT_COMPONENTS` is `FixedNeighbor`. The
estimator, workflow API, CLI and Python API do not expose a Joint route
selector; `--estimator joint-components` is sufficient.

The production FixedNeighbor policy is intentionally small and fixed: disjoint
mutable cores of 12 atoms, at most 200 forward Gauss-Seidel sweeps, and one
trusted accepted local Profile LM update per block visit. Outside-core
neighbors remain fixed while a block is solved. Each accepted block is replayed
globally, checked for monotone roundoff-safe acceptance, and contributes to the
existing global A/C KKT (`1e-10`), width-gradient (`1e-12`) and eta-confirmation
(`1e-10`) criteria. Meeting these checks makes the endpoint a candidate.
Successful production stopping additionally requires endpoint verification and
complete numerical qualification, including full design rank and projected-
width rank. A verified candidate that is not qualified is never exposed as a
scientific component or assembled state. The full parameter Jacobian is
`[X_F, D]`, where `X_F` is the free A/C design on the active face and `D` is
the log-width derivative. Its rank is
`rank([X_F, D]) = rank(X_F) + rank(P_perp D)`, where `P_perp` projects away
the column space of `X_F`. The projected-width check assesses the second term;
the separate full-design check assesses `X`. Both remain required, and a
threshold-sensitive or unavailable rank never counts as full rank. Numerical
identifiability describes parameter directions; it does not replace the A/C
KKT, width-stationarity, eta-confirmation or independent-replay checks that
establish optimization and endpoint quality.
Recoverable numerical qualification failures continue through the same bounded
block updates; resource-limited assessment stops retries without reporting
convergence. An unqualified endpoint at the sweep limit reports budget
exhaustion and remains available only as a candidate diagnostic. The
at-most-one accepted-update contract is intrinsic; diagnostic controls such as
reverse order and trajectory capture are not production configuration.

The full-design rank contract remains `tau = epsilon * max(N, P) * sigma_max`,
with `N` from the original rank policy and `P` equal to the full A/C design
column count. `sigma_max` is measured after the existing column normalization,
and singular values count as rank only when they are strictly greater than
`tau`. Production obtains this evidence from the free-design QR/SVD already
needed by projected-width preparation. When every amplitude is interior, the
free design has all `P` columns and is the full design, so the same normalized
factor supplies full Design Rank and projected-width coefficients. A
threshold-sensitive result is unavailable for qualification. The structural
obstruction runs first; it does not change the threshold or remove unsupported
observation rows from `N`. If `U` is the number of rows with an actual nonzero
design entry, then `rank(X) <= min(U, P)`. A zero design column or `U < P`
therefore proves that full-design rank cannot pass. Structural matching,
streaming estimates and SparseQR pivot rank do not certify numerical full
rank. With any active `A_i = 0`, the free-face factor has fewer than `P`
columns and cannot certify full Design Rank: production retains the
free-design diagnostic but marks full Design Rank unavailable. Structural
deficiency remains a valid full-design failure in either case.

The local numerical primitive is named **Profile LM** (Levenberg–Marquardt
optimization of the profiled width objective). `SearchProfile()` optimizes
`η = log(B)` with `InstrumentedLM<Profile>`. At every profile evaluation,
`EvaluateProfile()` conditionally solves the A/C coefficients with the existing
constrained `SolveLinear()` call, which uses EIGEN SparseQR:

```text
β*(η) = arg min 1/2 ||X(η)β - y_eff||²
        subject to A_i >= 0 and C_i unrestricted
```

The full profile derivative drives the LM width update. `AssessProfile()` keeps
its assessment responsibility. Each block visit accepts at most one trusted
update and replays the candidate against the global state.

For a sparse unit-weight local A/C solve, if the initial full-rank SparseQR
solution fails the existing KKT certificate, production may run a bounded
column-normalized QR refinement on that same free face. A preflight limits the
dense face matrix to at most 1,000,000 double entries (8 MB); the factorization
workspace is also local to that bounded face. The refinement uses the existing
rank threshold and performs at most eight residual-correction solves. A result
is adopted only when the unchanged `CertifyLinear()` check confirms feasibility
and the existing `1e-10` KKT gate, with objective no greater than the initial
SparseQR result. Rank-deficient faces, weighted inputs, reference-SVD solves,
oversized faces, and solutions that already pass KKT skip the refinement. It
does not build a global dense design or add another solver route.

The production search returns an internal `FixedNeighborComponentResult` with
the search state and independent endpoint certification. Shared orchestration
turns that into the certified public `JointComponentResult`, records component
cost/status and progress, then releases the larger search result. The batch
retains those public component results, aggregate costs and only the optional
single-component assessment that can be reused by exact-match FullABC assembly;
it does not hold a second vector of solver results. The route returns through
`JointFitResult` and the standard post-processing and persistence path. It uses
the defined observation domain, contributor/halo closure, parameterization, rank
semantics, parent scale, objective and `RuntimeConvergence()` semantics. Its
progress event includes sweep, block, accepted-local-update, KKT, width-gradient,
eta-change and objective fields.

The production FixedNeighbor implementation uses `MakeComponentProblemView()`
to prepare an immutable indexed view over the shared root snapshot. The view
keeps explicit parent/local atom and row mappings, identities and rank shape;
support rows and squared distances remain in the root snapshot. Structural
partitioning then creates one prepared block view per block, including its core
and affected-row mappings, indexed `Domain` and context template. A sweep still
rebuilds the state-dependent `eta`, `beta`, effective response and numerical
design; it does not rebuild those structural views. No numeric factor, active
set or coefficient solution is carried across sweeps.

`JointProblem::ParameterLayout()` records every retained contributor as a
FullABC atom. `OrchestrateFixedNeighborComponents()` projects that layout onto
each post-pruning component and owns component traversal, width-seed validation,
search, independent certification and component progress. Current production
problems have no nuisance groups or alternate assembly route. FullABC assembly
consumes the certified public component states and performs global
reconstruction and consistency assessment. Historical layouts remain in the
JSON decoder and exporter only.

`FitFixedNeighborComponentsImpl()` owns the Eigen thread scope and solver
provenance; `FitFixedNeighborComponents()` adds target evidence after
implementation and assembly finish. Public
`FitJointComponents()` and `EstimateJointComponents()` keep their existing
contracts and do not expose assembly details to callers.

Production search retains aggregate counters, final sweep state and the latest
endpoint assessment needed by the component route. Persisted component counters
record actual sweeps, first and confirmed stationarity, and the last certification
sweep. A completed final assessment at budget exhaustion records the same sweep
as `search_sweeps`; no certification sweep is recorded when no assessment ran.
Detailed block records, local trial trajectories, reverse order and local
certification remain opt-in experiment/test diagnostics. The preparation and
mapping counters are reported as implementation attribution; they are not
convergence criteria.
When a block cannot establish a trusted update, production emits one warning
with component/sweep/block identity, core atom IDs, affected rows, the halo
support/rank summary, initial profile status, local design dimensions/rank, and
the local stop reason. A rejected trust replay also includes its prediction,
gradient and KKT discrepancies. Successful blocks do not emit this detail.

The whole-problem `RunFixedNeighborDiagnostic()` facade is compiled into test
support and called only by component tests and the fixed-neighbor experiment.
The library retains lower-level telemetry and observer hooks used by that
support; they are implementation instrumentation. Production receives only the
qualified core-size, sweep-count and block-order controls through
`FixedNeighborSearchPolicy`, plus the progress observer.

Within one immutable problem, a full single component can share its assessment
with assembly only when observations, structural support, identities, state,
scale and linear/rank policy match exactly. This avoids repeating the reference
evaluation and assessment for the same returned state. Offline audit directions
do not affect runtime assessment. Constant rows or different contexts require
separate assessments. This is scoped result reuse, not a persistent cache. The
test-support Profile-search diagnostic assesses a selected historical trial at
its saved coefficients and widths. Assessment work and timing are reported
separately from search work. Sparse factorization details can change
floating-point trajectories within the fixed endpoint contract.

### Selection-domain contract

The Map/Model builder uses selected non-hydrogen atoms as targets S. Their
2.5 Angstrom spheres fix the unique observation rows V. Selected targets remain
in the target mapping even if clipping leaves one with no observed support; no
rows are invented for it. Each unselected non-hydrogen atom intersecting V is a
halo candidate. A candidate is retained only when it has at least three support
rows and at least three exact distinct squared distances. This is the minimum
necessary FullABC gate, not a full-component identifiability certificate.
Selected targets bypass this gate. Pruning never changes V, its row IDs,
observations, or observation scale; rows with no retained contributor stay in
the objective as constant residuals. Atom centers outside the map can qualify
if their support intersects V. Empty non-hydrogen selection is rejected. The
existing `sphere-fma-v1` arithmetic and negative/zero observations are retained.
Structural partition runs after pruning, so only retained halos can couple
components or act as bridges.

Pruning changes the additive model: the fitted prediction sums selected targets
and eligible retained halos. Excluded halo signal remains in the observations
and can affect residuals and selected-target estimates. The estimator adds no
background term or omitted-signal correction.

`JointProblemInput.selection_domain` is optional for hand-built frozen inputs
(null means selection not recorded). The Map/Model builder always records
`fixed-selected-voxel-closure-v1`, ordered unique `target_indices`, both radii
(2.5 Angstrom), the `fullabc-eligible-halo-pruning-v1` policy, and eligibility
provenance. Distinct squared distances use exact equality; adjacent
representable values remain distinct, with no tolerance. Provenance records the
minimum counts, candidate and retained halo totals, and each excluded atom ID,
support count, distinct-distance count, and reason. Retained halo indices are
the non-target contributors in the final problem. This establishes completeness
only for the supplied catalogue and truncated support model. The saved outcome
retains the policy and atom ordering; stage and CSV metadata report target and
retained-contributor roles. Excluded halos do not appear as failed fitted atoms.
Hand-built inputs may keep the legacy `all-non-hydrogen` selection contract;
inputs with no selection metadata retain their explicit FullABC semantics.

`EstimateJointComponents` and the production workflow share
`PrepareJointWidthInitialization()`: it samples only FullABC contributors with
the deterministic Fibonacci sampler, then trains local alpha and computes
first-stage MDPDE width seeds for all required contributors. The workflow first
initializes the selected model analysis and applies the prepared results there.
`EstimateJointComponents` runs the same preparation on a full model copy, then
copies raw samples and first-stage results back only for successfully initialized
targets. Only B initializes joint fitting. Its initialization data scope is
`contributor-local-sampling-may-read-outside-target-domain`: interpolation and
sampling can read beyond V, while the joint objective cannot. Supplied widths
use `caller-provided-widths`. Global map normalization is also separate from V.
The standalone copy/writeback contract preserves halo history, original
selection, second-stage and group results. Initialization records retain
atom-specific reasons and nonfinite diagnostics as unavailable, not replacement
coefficients.

An incorrectly sized initial-width vector remains a whole-input failure.
Otherwise invalid widths/initialization exceptions block only their structural
component, which retains mappings, stop reason and unavailable evidence; other
components continue. `initialization.valid` describes all widths, not whether
any component can run. Missing component states still prevent global assembly.
Peeling results, truth, prior fits and certificates do not initialize the
estimator.
`costs.construction_seconds` measures the Map/Model builder; initialization timing
includes the model copy and setup. Direct frozen-input fitting has zero builder
and initializer costs.

A/C order is `[A0, C0, A1, C1, ...]`; widths follow atom identities. Public
component, assembled-state and result objectives are
`0.5 * ||prediction - observations||² / ObservationScale()²`. Components use the
same parent scale `max(1, ||parent observations||₂)`; constant rows contribute
once to the global objective.
Width gradients are with respect to log-B. Internal raw objectives retain
half-RSS semantics. This public objective correction changes returned
scale, not the optimizer or its acceptance thresholds.

Component `state` is the actual last trusted state used by assembly. Search
completion, usable state and evidence status are separate. A missing component
prevents complete prediction/objective; the row mask retains available rows,
including constant rows. Assembly does not reprofile or zero-fill missing states.
Each attempted component search also retains a parameter-only candidate
diagnostic when its last accepted state is finite and dimensionally valid.
Its candidate status records endpoint verification and numerical qualification
separately. A verified candidate with unavailable rank remains a candidate and
does not become trusted scientific state.
The diagnostic may include eta, A/C, objective, feasibility and stationarity
summaries, but never prediction/residual vectors or rank matrices. The enclosing
saved result supplies the component identity, sweep and stop reason, input hashes
when available, and solver provenance. This is a diagnostic snapshot, not a
search-resume checkpoint. Candidate values never populate `component.state`,
assembly, uncertainty or the atom CSV.
Two-step derivative, Richardson, multiprecision, boundary and regular certificates
are `NotRun` in runtime results. Offline tools produce separate evidence without
mutating the runtime result. Operational runtime checks are not regular
certification.

`component.RuntimeConvergence()` and `result.RuntimeConvergence()` derive their
status from the returned evidence; they do not run numerical work or store another
success flag. Component convergence requires an actual trusted state and current
qualification checks: A/C feasibility, conditional A/C projected KKT
(`1e-10`), independent streaming replay, width stationarity (`1e-12`), design
rank, projected-width rank and numerical identifiability. Replay prediction
enclosure remains `2e-12 + 2e-13 |prediction|`; width-gradient agreement remains
`1e-13 + 2e-9 |gradient|`; normalized KKT agreement remains `1e-13`. Full
profile derivatives, including the nonzero-residual correction, remain in the
numerical core.

Full parameter identifiability is not the same as optimization stationarity: a
stationary point may still have rank-deficient or unavailable parameter
evidence. Eligibility is a model-selection rule; the existing full numerical
qualification remains required after pruning. The current solver-provenance
contract is v5 and the FixedNeighbor policy contract is v4. Target evidence
uses `target-estimability-v4`; JSON schema 11 records its active-boundary
semantics and halo-pruning provenance. Schemas 3 through 10 remain readable;
decoding retains their saved evidence and convergence statuses, including
historical singleton-halo parameterization, without rerunning those solves.

Global convergence additionally requires every component to be numerically
qualified, complete assembly, global evidence and assembled mapping plus
independent prediction/objective consistency. Assembly does not globally
reprofile A/C coefficients.
Missing required states or assembly give `Unavailable`. Otherwise, required
checks combine in this order: any `Failed`, any missing/`Unavailable`, any
`NotRun`, then `Passed`. Only evidence from the required scope participates.
Search completion is independent: a budget/untrusted-trial stop can leave a
trusted endpoint whose runtime checks pass; its stop reason is still reported.

Runtime convergence does not imply an offline regular certificate.

## Immutable storage and component views

The input snapshot owns observations, memberships and identities once. Numerical
observations are read-only Eigen maps; contexts and component support views retain
shared ownership, including after the original problem handle is destroyed.
Partition reverse mappings are shared once across components. Local support is
mapped on access without copying squared distances. Contiguous observation ranges
are borrowed; noncontiguous ranges use one scoped gather buffer for search and
assessment. Census and offline bundles materialize legacy mappings only at their
serialization boundary.

Component, profiled-component and prepared-block views each keep their local
atom selection and precomputed row-to-local mapping over the parent view. Support
iteration applies those mappings directly in constant work per view layer, and
FixedNeighbor reuses prepared block mappings across sweeps.

## Memory feasibility and unavailable large-case checks

`ResourceWork` plans problem construction, FixedNeighbor preparation, endpoint
verification, rank assessment and target identifiability
separately. Estimates use the actual support-entry count and checked byte
arithmetic. The planning budget is 3,200,000,000 bytes (3.200 GB, about
2.980 GiB), leaving 20% below the nominal 4,000,000,000-byte limit (4.000 GB,
about 3.725 GiB). It is an engineering margin, not a numerical tolerance.
Each log reports its stage, dimensions, resident bytes, temporary bytes,
estimated peak and budget. These are estimated allocations, not observed peak
RSS or a promise about the operating system's process limit. A rejected stage
returns an unavailable resource result with a `resource-limit:<stage>` reason;
it is not recorded as numerical nonconvergence.

Before numerical rank assessment, production counts the actual nonzero support of
the full A/C design. A zero design column or fewer unique supported rows than
full design columns proves that full design cannot have full column rank. This
one-way result marks design qualification `Failed`; it reports a rank upper
bound and reason without inventing an exact rank, singular values or condition
number. Structural matching, when run by research diagnostics, can also prove
deficiency; a full matching does not prove numerical full rank. The check
retains the full A/C column domain even when amplitudes are active at the
boundary.

If those low-cost checks find no obstruction, the rank preflight covers the
single projected-width preparation, including its free-design QR/SVD, retained
projection coefficients and projected-width QR/SVD. There is no separate
production DesignSpectrum allocation. Target identifiability reuses the full
design and projected-width rank evidence and plans only rank records,
target/halo identities and per-atom statuses; it does not build or factor a
full target Jacobian. The Joint estimator does not compute a global parameter
covariance. The parent rank peak is the maximum of sequential stage peaks with
each stage's resident data; it does not add all stage temporaries together.
The assessment telemetry records
whether a stage was entered and completed, its dimensions, elapsed time,
preflight status and byte breakdown. Current RSS is sampled at stage entry and
exit; process peak RSS is the process-wide high-water mark. Neither RSS value
is the stage's attributed allocation. On platforms without an RSS source the
fields are unavailable. Production projected-width rank assessment uses a
values-only compact SVD and does not allocate the full right-singular-vector
matrix. Its planner models six square-factor equivalents of Eigen BDC SVD
workspace, plus the input factor in the peak. Explicit weak-direction
diagnostics use a separate `ProjectedWidthDiagnostics` preflight and a second
SVD with thin right vectors; the diagnostic workspace model includes the
vector matrix and copied diagnostic result. This opt-in path does not change
rank or qualification decisions.

The rank and target preflight byte formulas and the synthetic 37,406-atom
census are recorded in the [component memory census](joint-component-memory-census.md).
The free-design factor alone scales quadratically with the number of free
columns. The projected-width stage also retains sparse inputs, projection
coefficients, tiled QR workspaces and projected-width factors. Both the full
interior shape and its all-active free face exceed the planning budget in the
synthetic census; neither full qualification is attempted. These are planned
allocation estimates, not measured RSS. Derived target evidence has no factor
or right-vector workspace. This rank limitation remains independent of the
downstream local MDPDE and group information workflow.

Prepared row mappings store sorted parent/local row pairs only for each block's
affected rows, plus cached local support indices. Their structural cost follows
the sum of affected rows and support memberships, not
`component_rows × block_count`. The global sweep KKT and width-gradient pass now
traverses support entries and keeps vector state; it does not build the global
`N × 2P` design and derivative matrices. Local block solves and endpoint
assessment still use their existing sparse design matrices.

Production endpoint verification is O(N+P) and does not build a global design
matrix or run a global reference solve. The old independent solver and
`DesignSpectrum()` remain offline/test oracles. Production has no
`DesignSpectrum()` call and no second full-design QR/SVD: interior full Design
Rank comes from the free-design `R` factor already produced for projected-width
preparation. With active amplitudes that factor certifies only the free face;
full Design Rank is unavailable unless the retained structural obstruction
proves deficiency. Derived target identifiability consumes only valid
full-design and projected-width records and does not build a target Jacobian,
QR/SVD factor or null space. Its memory is linear in retained FullABC and target
records. Excluded-halo provenance is O(E) in the number of exclusions and is
included in problem memory preflight. The
serial synthetic census rejects the large rank stage, while the derived-target
preflight itself passes; no trusted rank evidence exists in that case. No
physical map/model search has been measured.

After Joint fitting, `BuildPostFitPeelingSamples` subtracts every other
retained contributor from each selected target's raw samples using the Joint
grid's tricubic interpolation. `BuildPostJointLocalSamples` then subtracts the
target's own fixed Joint charge profile `C_i^J K(r; B_i^J)` with the same
support and interpolation convention. `FitPostJointLocalTarget` trains local
alpha and runs the existing OLS initialization and local MDPDE on that
charge-corrected Gaussian data. The resulting local A/B estimates and solver
qualification are stored separately in `PostJointLocalResult`; they do not
rewrite or requalify the Joint endpoint.

`RunJointGroupPotentialFitting` trains group alpha from qualified local MDPDE
points, then uses their working observation information `P_i = X_i^T
Sigma_y_i^-1 X_i` and `h_i = X_i^T Sigma_y_i^-1 y_i`. It reuses the current
`EstimateMuMDPDE` and `EstimateWEBFromInformation` equations without fitting
the same local observations again. Each eligible target retains only its
compact 2-by-2 precision and two-value information moment. Per-sample working
weights and covariances are temporary; no component-wide or cross-atom
covariance matrix is formed. Peeling and local sample preparation are bounded
per target, and group WEB is bounded per group.

The workflow builds one immutable `JointPostprocessingContext` after the Joint
snapshot is published. Its observation-row and contributor-ID lookups are
shared through peeling and all target-local preparation. It reads the trusted
snapshot and published endpoint records without reconstructing Joint stage
estimates per target; each target only builds its own charge support and sample
workspace. Context construction checks complete snapshot contributor and row
identity once for the run. Per-target preparation then checks only its selected
target, endpoint/peeling provenance, paired samples, geometry and support
stencils. The context borrows the map, model endpoints, problem and snapshot;
they remain alive and unchanged until that run's target-local preparation ends.
`AssessPostJointLocalResult` supplies the same source, sample,
endpoint, numerical and compact-information contract to fitting, atomic batch
publication, group eligibility and persistence validation.

Group posterior members retain Gaussian posterior A/B separately from their
fixed Joint `(C_i^J, B_i^J)` profile. The physical member response evaluates
the posterior Gaussian at `B_post` and its charge kernel at `B_i^J`. Group
means and priors are Gaussian-only; the group charge median is descriptive and
does not define a single FullABC physical response. See the
[stage-integration guide](joint-stage-integration.md) for the lifecycle and
invalidation details.

Response plots distinguish these signals in their legends: Joint peeling
includes the target's fixed charge profile, while a group prior or mean is
Gaussian-only with zero offset. The curves are shown as different
representations, not as measurements and a matching fit. Full posterior member
responses use each member's own persisted Joint charge width through
`EvaluateGroupMemberResponseAtDistance`.

The fixed Joint charge profile remains unchanged when the local width differs.
Log-Gaussian fitting filters nonpositive corrected responses according to the
existing dataset builder. Local measurement information is a working model;
it does not propagate uncertainty in other Joint contributors or cross-atom
errors. Rank preflight and target qualification retain their existing behavior
and can still exceed the memory budget independently.

The Map builder collects and sorts relevant voxel indices instead of allocating
an index array for the entire map. Voxel ordering and support arithmetic remain
unchanged. Current Joint search, solve, full workflow and command
measurements are run through the single workflow in the
[benchmark guide](joint-benchmark.md).

The FixedNeighbor experiment reports structural census counts, estimated
preparation, endpoint and rank peaks, stage timing attribution, per-sweep summaries
and process RSS. It supports preparation-only, one-sweep and bounded
multi-sweep diagnostics without changing the production policy. Aggregate
counters are the default; per-block and profile-trial records require explicit
detailed mode. The 37,406-atom preparation fixture is synthetic and does not
measure a physical map/model search.

## Internal numerical path

The current numerical path is:

```text
FixedNeighbor
  -> local Profile LM over η = log(B)
  -> conditional A/C profile solve / EIGEN SparseQR
  -> full profile derivative and LM width update
  -> global replay, endpoint certification, RuntimeConvergence
```

These stages preserve the original parent rank context and full residual norm.
They do not cut component edges or use a size-based route. The
[sparse-backend guide](joint-component-sparse-backend.md) documents the fixed
EIGEN implementation.

Permanent tests check the EIGEN path and numerical contract: trusted-state
availability and required runtime checks must agree, missing states remain
explicit, and the production KKT, width-gradient, eta-confirmation, replay,
endpoint-certification and `RuntimeConvergence` thresholds remain unchanged.

## Routine regression

```sh
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --target tests_all --parallel 4
ctest --test-dir build -L joint:runtime --output-on-failure --parallel 4
```

The runtime test lane protects the estimator behavior and public API. The
integration regression runner and its fixture catalog are documented in
[tests/README.md](../../tests/README.md). Offline certification runs are a
separate lane and do not change the runtime result.

Runtime search/assembly reports use schema version 2: `runtime_convergence`,
`runtime_checks` and `runtime_failure` replace the mixed `joint_qualified` and
`qualification_*` fields. Search-endpoint identity is retained in diagnostics.
Missing evidence never means passing. Legacy `joint_qualified` is produced only
by offline assessment and retains its derivative-dependent meaning.

`joint_component_runtime` provides commands for public API runs, physical-input
smokes, and fixture cases. The Python runner exposes run, summary, comparison,
regression and physical-smoke operations. Current entry points are maintained in
`tests/integration/joint_component_runtime.py`. It is a validation harness,
not a benchmark or experiment entry point. Current experiments and validation
tools are listed in [joint-experiments.md](joint-experiments.md).

## Extended and offline checks

Both options default to `OFF` and require `BUILD_TESTING=ON`:

```sh
cmake -S . -B build -DBUILD_TESTING=ON \
  -DRHBM_GEM_ENABLE_JOINT_EXTENDED_TESTS=ON \
  -DRHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON
cmake --build build --target tests_all --parallel 4
ctest --test-dir build -L 'joint:extended|joint:offline' --output-on-failure --parallel 4
```

Extended and offline lanes are opt-in and stay separate from the default
runtime regression.

Local-audit preparation and precision code are linked only into
`joint_component_audit`, `joint_offline_diagnostic` and `joint_offline_tests`,
never the library or ordinary test executable. Kernel or Jacobian changes
require the applicable derivative
audits; constraint changes also require boundary controls. Sparse-factorization
changes require rank and precision controls. Recorded campaign results and their limits
belong in the [canonical evidence index](joint-component-evidence.md).

```sh
python3 tests/integration/joint_component_audit.py \
  --executable build/bin/joint_component_audit --work-dir build/joint-audits
build/bin/joint_component_audit local-bundle INPUT_JSON OUTPUT_DIRECTORY
```

A standalone bundle carries immutable parent observations/context, component
memberships and initial widths. Local directions come from the actual component
endpoint; sibling/global endpoints are not inputs. Modified parent context or
invalid widths are rejected. Parent normalization, actual-state preservation,
local weak-direction evidence and independently reprofiled consistency remain
required; a local certificate cannot promote a missing or failed component.

Current sparse-factorization behavior and benchmark commands are documented by
the sparse-backend and benchmark guides.

## Saved production outcomes

`potential_analysis --estimator joint-components` stores a `JointAnalysisResult`
on the model and saves it through `DataRepository`. The default estimator remains
two-stage. See the [command contract](commands/potential_analysis.md#joint-component-opt-in).

Direct `EstimateJointComponents` callers keep the return-only joint contract.
To explicitly retain a result, use `CaptureJointAnalysisResult(fit, metadata)` and
`model.EditAnalysis().SetJointResult(...)`, then `DataRepository::SaveModel`.
Read it through `model.GetAnalysisView().GetJointResult()`. Include
`<rhbm_gem/data/io/JointAnalysisFileIO.hpp>` to export a saved result directly with
`WriteJointAnalysisResult(result, json_path, csv_path)`.

The problem snapshot, prediction vector and numerical workspaces are not captured.
Trusted result values, optional parameter-only candidate diagnostics and
identity/mapping information survive. Global and component runtime-convergence statuses are captured from the existing runtime methods;
storage validates the document structure but neither recomputes convergence nor
runs an audit. Candidate diagnostics do not define a resume contract. These saved
documents are not standalone numerical audit bundles.
Model copying retains the outcome; analysis `Clear()` and `ClearJointResult()`
remove it, while `ClearTransientFitStates()` does not.

### Saved provenance and units

Production outcome JSON writes schema 11 and retains schema 3 through 10 reading
(distinct from fixture and offline report schemas). The [halo eligibility contract](joint-observable-halo.md)
describes selection, model change and provenance. See the
[metadata contract](commands/potential_analysis.md#provenance-and-map-units-joint-json-schemas-3-through-11).
`JointAnalysisMetadata` contains optional `JointMapNormalization`, input SHA-256
values, `JointSoftwareProvenance` and optional `JointSolverProvenance`.
`JointSolverProvenance` records the current FixedNeighbor estimator and
preserves solver and local-work fields found in saved records.
Current writers record the fixed EIGEN sparse backend and the FixedNeighbor policy
fields; current FixedNeighbor records use `fixed-neighbor-production-v4` with
intrinsic OneAccepted behavior and the `joint-solver-provenance-v5` contract.
Schema 8 through 10 candidate diagnostics store endpoint-verification and numerical-
qualification statuses separately. Legacy schemas preserve their saved checks;
no legacy `inner=Passed` is mapped to a new KKT result and no historical
reference status is reinterpreted as streaming replay.
Saved `joint-solver-provenance-v2` records with LegacyCompact and OperatorPcg
fields remain readable.
`CaptureJointAnalysisResult`
records the current library's version/source/configuration/build identity. The
decoder and exporter preserve saved identity without recalculation. Direct
in-memory callers leave unknown input hashes and normalization null; capture
does not open paths.

`MapObject::MapValueArrayNormalization()` returns `std::optional<double>`: the
actual divisor if applied, or nullopt when zero SD skips normalization. Existing
callers may ignore the return. A caller retaining provenance records its requested
flag, whether the returned divisor exists, and `divisor.value_or(1)`; a command
that skips normalization records `applied=false, divisor=1`.

B and geometry use Angstrom, A uses fitted-map-unit*Angstrom^3, and C uses
fitted-map-unit*Angstrom. Multiply A/C by the recorded divisor to recover input-map
scale. These are kernel coefficient units, not an assertion of physical map
calibration; C is not a constant background offset. The objective's parent
observation scale is separate. Unknown normalization cannot support conversion.

The command's default estimator remains two-stage. Joint LS is opt-in and its
runtime checks do not establish an offline regular certificate.

Permanent [partial-selection tests](../../tests/core/joint_component/PartialSelection_test.cpp)
cover the fixed observation domain, pre-construction halo pruning, recomputed
bridge topology, unavailable target states, and the numerical policy.

## Joint validation tiers

Build with at most four jobs, then run the bounded gates from the build tree:

```sh
cmake --build build/debug-tests --parallel 4
ctest --test-dir build/debug-tests -L 'joint:contract' --parallel 4
ctest --test-dir build/debug-tests -L 'joint:workflow' --parallel 4
ctest --test-dir build/debug-tests -L 'joint:numerical' --parallel 4
```

Contract covers deterministic component, linear algebra, and small assessment
checks, including all 14 `JointComponentChecksTest` cases. Workflow covers
runtime, persistence, CLI, and smoke paths. Numerical covers bounded
FixedNeighbor, stationarity, order, prepared-block, KKT, and derivative/reference
qualification. Each gate has a CTest timeout and sets `OMP_NUM_THREADS=1`.

Ordinary Joint cleanup or refactoring accepts passing all three gates. A
numerical solver change also runs each relevant scalability case. Release or
performance qualification runs Tier 4 for its campaign scope:

```sh
ctest --test-dir build/debug-tests -L 'joint:scalability' --parallel 4
ctest --test-dir build/debug-tests -L 'joint:offline' --parallel 4
```

Tier 4 gives `joint_component_gram_scalability`, the two 256-atom FixedB cases,
and optional offline audits their own CTest entries with a one-hour timeout.
These entries are outside ordinary acceptance. Offline audit entries are
registered when configured with
`-DRHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON`.

Report test states as `passed`, `failed`, `timeout`, `not-run`, or
`not-completed`. A manually stopped run without a test result is
`not-completed`; it is neither a pass nor a failure.
