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
reference/replay, LM settings and endpoint assessment are unchanged by cleanup.

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
Richardson, multiprecision, boundary and regular certificates remain `NotRun`
until an offline tool supplies evidence. Operational runtime checks are not
regular certification.

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

Offline precision code is linked only into `joint_component_audit` and
`joint_offline_tests`, never the library or ordinary test executable. Tests
cover baseline, near-0.02/narrower, weak-1e-4 and active-a; qualification failures
must match their historical scopes. Kernel/Jacobian changes require derivative
audits; constraint changes also require boundary controls. Backend changes
require rank/precision controls and the extended lane.

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

## Retired workflows

The [evidence index](joint-component-evidence.md) records the retired research
workflows, limitations and retrieval commits. No ordinary regression requires
the historical 72/216/128-case chains or 5,008-report replay. Guarded is the only retained search branch. First-stage `mdpde_experiment solve`, `forward`
and `refine`, production second-stage and the separate fold-168 regression
remain supported. Endpoint assessment deduplication and memory/backend changes
are separate future work.
