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

`JointProblem` owns an immutable snapshot. `JointProblemInput` accepts frozen
observations, identities and squared-distance memberships. V1 uses equal-weight
LS, A nonnegative, C signed, log-B, structural 2.5 Angstrom support and Guarded
search (200 profile evaluations / 100 accepted updates per component). Search
reference/replay and LM settings are unchanged by endpoint assessment cleanup.

Endpoint evaluations and reference solves are passed into assessment and trust
replay instead of being repeated. A component returns both its historical search
endpoint diagnostics and the assessment of its actual trusted state; the public
API consumes the latter. Assembly uses its raw evaluation for the returned state
and retains an independently reprofiled consistency control.

Within one immutable problem, a full single component can share its assessment
with assembly only when observations, structural support, identities, state,
scale and linear/rank policy match exactly. Offline audit directions do not
affect runtime assessment. Constant rows or different contexts require separate
assessments. This is scoped result reuse,
not a persistent cache. A fallback state's evidence is assessed at that state's
actual coefficients and widths. Assessment reuse does not change search decisions or counts; assessment work
and timing are separate from search work. Orthogonal backend row reduction can
change the search trajectory within the numerical parity contract below.

The Map/Model builder includes all non-hydrogen contributors and rejects partial
non-hydrogen selection. It uses the actual Map geometry and `sphere-fma-v1`,
retaining negative/zero observations and every structural membership, including
zero coefficients and numerically invisible basis values.

`EstimateJointComponents` runs deterministic Fibonacci sampling, alpha training
and first-stage MDPDE in a model copy. Only B initializes joint fitting. Raw
samples and first-stage results are copied back; selection, second-stage and
group results are preserved. Invalid widths yield explicit initialization
failure. No Peeling result, truth, historical fit or certificate initializes the
public estimator.

A/C order is `[A0, C0, A1, C1, ...]`; widths follow atom identities. Public
component, assembled-state and result objectives are
`0.5 * ||prediction - observations||² / ObservationScale()²`. Components use the
same parent scale `max(1, ||parent observations||₂)`; constant rows contribute
once to the global objective.
Width gradients are with respect to log-B. Internal and historical raw objectives
retain half-RSS semantics. This public objective correction changes returned
scale, not the optimizer or its acceptance thresholds.

Component `state` is the actual last trusted state used by assembly. Search
completion, usable state and evidence status are separate. A missing component
prevents complete prediction/objective; the row mask retains available rows,
including constant rows. Assembly does not reprofile or zero-fill missing states.
Two-step derivative, Richardson, multiprecision, boundary and regular certificates
are `NotRun` in runtime results. Offline tools produce separate evidence without
mutating the runtime result. Operational runtime checks are not regular
certification.

`component.RuntimeConvergence()` and `result.RuntimeConvergence()` derive their
status from the returned evidence; they do not run numerical work or store another
success flag. Component convergence requires an actual trusted state and all five
checks: inner/reference agreement, A/C KKT (including feasibility), log-B gradient,
local correction and numerical identifiability. Existing thresholds remain
`1e-10` for coefficient agreement, KKT and local correction, and `1e-12` for both
primary/reference width gradients. Full profile derivatives, including the
nonzero-residual correction, remain in the numerical core.

Global convergence additionally requires every component to converge, complete
assembly, global evidence and independently reprofiled assembly consistency.
Missing required states or assembly give `Unavailable`. Otherwise, required
checks combine in this order: any `Failed`, any missing/`Unavailable`, any
`NotRun`, then `Passed`. Only evidence from the required scope participates.
Search completion is independent: a budget/untrusted-trial stop can leave a
trusted endpoint whose runtime checks pass; its stop reason is still reported.

For the frozen representative cases, baseline and near-0.02 pass runtime
convergence. Weak-1e-4 also passes runtime convergence while its offline derivative
evidence remains resolution-unverified. Active-a fails local correction,
zero-signal fails width identification, and duplicate has no usable state. None
of these runtime results claims an offline regular certificate.

## Immutable storage and component views

The input snapshot owns observations, memberships and identities once. Numerical
observations are read-only Eigen maps; contexts and component support views retain
shared ownership, including after the original problem handle is destroyed.
Partition reverse mappings are shared once across components. Local support is
mapped on access without copying squared distances. Contiguous observation ranges
are borrowed; noncontiguous ranges use one scoped gather buffer for search and
assessment. Census and offline bundles materialize legacy mappings only at their
serialization boundary.

The Map builder collects and sorts relevant voxel indices instead of allocating
an index array for the entire map. Voxel ordering and support arithmetic remain
unchanged. `joint_component_benchmark DATASET CASE` measures the public API in a
fresh process and reports construction, search, assessment, assembly, total time
and process peak RSS. It uses only the installed public API and can also be built
against the baseline library.

## Tiled numerical backend

All production components use 8192-row tiles. The first derivative pass reduces
the normalized free design and raw width derivative. The second generates the
projected derivative and full residual-corrected Jacobian by tile, retaining only
compact QR factors, transformed residuals and column norms. LM pivots the compact
Jacobian factor and retains the original full residual norm for its objective,
actual reduction and stopping calculations. Rank thresholds still use the
original context rows/columns, not the reduced factor dimensions.

Derivative/LM dense workspace is O(tile * atoms + atoms²), in addition to sparse
designs, sparse factorization storage and O(rows) residual vectors. This is not a
bound on all sparse fill-in or on the parent-global assembled assessment as the
number of atoms grows. No component edges are cut and no dense runtime fallback
is used. Dense derivatives exist only in test support as a parity reference.

The fixture runner defaults to backend numerical comparison: trusted-state
availability, runtime checks and limitations remain fixed; converged identifiable
A/C/B endpoints use scaled 1e-10 and normalized objectives use 1e-12. For a failed
convergence check, objective may not worsen by more than 1e-12 and same-state
parity is checked at both historical and actual endpoints. Changed active faces
are reported explicitly; a free-face rank is compared across endpoints only
when their active faces match. Dense/tiled ranks at each identical state must
always match. Search counts, native stop codes and trial sequences are recorded
but are not required to be identical. Guarded acceptance and budgets remain
mandatory. `--strict-history` retains the older exact-trajectory comparison for
baseline/ownership-only validation.

Each same-state comparison checks projected/full derivatives at relative 1e-8,
spectra at 1e-10 relative to their largest singular value, local correction at
scaled 1e-10 and gradient at 1e-13 + 2e-9 * abs(reference). Rank/availability must
agree. Unavailable corrections and missing trusted states remain explicit
limitations. No fixture packages or numerical tolerances are regenerated.

## Routine regression

```sh
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --target tests_all -j
ctest --test-dir build -L joint:runtime --output-on-failure
```

The ordinary executable tests numerical, support, partition, same-state,
actual-state assembly, local-context and failure-isolation contracts. Eight
frozen representative cases run through the public API and numerical core;
physical double/float32 and CIF/MRC tests exercise fresh initialization. The
[fixture catalog](../../tests/fixtures/joint_component/README.md) is self-contained
and uses independent pre-extraction records plus scalar reference controls.
The Python runner sets numerical thread limits before importing NumPy.

Runtime search/assembly reports use schema version 2: `runtime_convergence`,
`runtime_checks` and `runtime_failure` replace the mixed `joint_qualified` and
`qualification_*` fields. `search_endpoint_eta` preserves the identity needed
for historical endpoint comparisons. Runtime directional audit evaluations are
zero. Immutable fixture packages/hashes remain unchanged: the ordinary runner
compares search/state and non-derivative assertions; the optional two-step runner
checks the original derivative/qualification conclusions separately. Missing
evidence must never be interpreted as passing. Legacy `joint_qualified` is
produced only by offline assessment and retains its derivative-dependent meaning.

`joint_component_runtime` provides `run MODEL MAP OUTPUT`, `physical OUTPUT`,
`physical-inputs OUTPUT`, and `fixture DATASET CASE OUTPUT_JSON`. The Python
runner exposes `run`, `physical`, `summarize`, `compare`, `regression` and
`physical-smoke`. Fresh-input commands use one first-stage start, not the old
four-start/audit matrix. Numerical equality is checked against a monolithic
solve of the actual loaded problem; MRC header geometry is never replaced by
generation coordinates.

## Extended and offline checks

Both options default to `OFF` and require `BUILD_TESTING=ON`:

```sh
cmake -S . -B build -DBUILD_TESTING=ON \
  -DRHBM_GEM_ENABLE_JOINT_EXTENDED_TESTS=ON \
  -DRHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON
cmake --build build --target tests_all -j
ctest --test-dir build -L 'joint:extended|joint:offline' --output-on-failure
```

Extended tests add heterogeneous-168/first-stage-float32 and the remaining
catalog datasets. For initialization changes, select additional frozen starts:

```sh
python3 tests/integration/joint_component_runtime.py regression \
  --executable build/bin/joint_component_runtime --work-dir build/joint-starts \
  --dataset baseline --all-starts
```

Offline two-step, local-audit preparation and precision code are linked only into
`joint_component_audit` and `joint_offline_tests`, never the library or ordinary
test executable. Tests
cover baseline, near-0.02/narrower, weak-1e-4 and active-a; qualification failures
must match their historical scopes. Kernel/Jacobian changes require derivative
audits; constraint changes also require boundary controls. Backend changes
require rank/precision controls and the extended lane.

```sh
python3 tests/integration/joint_component_audit.py \
  --executable build/bin/joint_component_audit --work-dir build/joint-audits
build/bin/joint_component_audit local-bundle INPUT_JSON OUTPUT_DIRECTORY

# The eight default historical two-step/qualification controls, without precision scans:
python3 tests/integration/joint_component_audit.py --two-step-only \
  --executable build/bin/joint_component_audit --work-dir build/joint-two-step
```

A standalone bundle carries immutable parent observations/context, component
memberships and initial widths. Local directions come from the actual component
endpoint; sibling/global endpoints are not inputs. Modified parent context or
invalid widths are rejected. Parent normalization, actual-state preservation,
local weak-direction evidence and independently reprofiled consistency remain
required; a local certificate cannot promote a missing or failed component.

## Retired workflows

The [evidence index](joint-component-evidence.md) records the retired research
workflows, limitations and retrieval commits. No ordinary regression requires
the historical 72/216/128-case chains or 5,008-report replay. Guarded is the only retained search branch. First-stage `mdpde_experiment solve`, `forward`
and `refine`, production second-stage and the separate fold-168 regression
remain supported. Formal workflow adoption and result persistence use the opt-in command described below. See [tiled backend acceptance](joint-component-tiled-backend.md) for
validation and measured costs.

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

The snapshot, prediction vector and numerical workspaces are not captured. Only
result values and identity/mapping information survive. Global and component
runtime-convergence statuses are captured from the existing runtime methods;
storage validates the document structure but neither recomputes convergence nor
runs an audit. These saved documents are not standalone numerical audit bundles.
Model copying retains the outcome; analysis `Clear()` and `ClearJointResult()`
remove it, while `ClearTransientFitStates()` does not.
