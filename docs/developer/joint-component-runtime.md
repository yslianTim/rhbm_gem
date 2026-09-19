# Component-local certification and runtime parity

This work starts at `6f30510c`. The exact-component v1 and certification v2
contracts remain the numerical baseline. Local certificates supplement, and
never replace, monolithic/global and globally restricted certificates.

## Component-local evidence

`joint-abc-components-local-audit DATASET RUN OUTPUT [CASE]` reads frozen
component fits and their immutable search context. It does not read an assembled
fit or a final global context. Isolated reruns use a self-contained component
bundle with its immutable parent context and initial widths.

The Python component runner exposes these through `audit --local-only` and
`rerun-component --local-only`. The latter exports a self-contained bundle and
executes `joint-abc-rerun-local-bundle BUNDLE OUTPUT`: only local memberships/B₀,
immutable parent observations/context and the registered audit policy are read. Ordinary full audits additionally produce
`local-components/`, preserving the existing `components/` restricted evidence.

Each local scope records the stable component ID, parent snapshot identity and
exact eta/beta of the last trusted state used by assembly. Assessment uses those
coefficients, with a separately reprofiled consistency control. Normalization
and active-set response norm still come from all parent observations. Rank
sizes remain the existing independent-component sizes.

Inherited directions are discarded. The local endpoint supplies normalized
all-ones, alternating and weakest projected-width directions, with direction
norms and Jacobian response norms saved. Missing trusted states or width spectra
are explicit limitations, not unit-axis substitutes. Search completion, usable
state, derivative evidence and regular qualification remain separate.

The two-step derivative checks and registered Richardson/50-100 digit audits
retain their existing thresholds. Truth never selects directions, intervals or
qualification. Local certificates additionally require a usable actual state,
weak-direction evidence and agreement with the profiled control. Missing
components cannot be promoted by another component's certificate.

Independence means identical component input and immutable parent context:
sibling execution order and sibling endpoint availability cannot affect local
evidence. It does not assert invariance when parent observations/scale change.

## Implementation and validation records

The implementation is delivered in three stages: local certification, typed
runtime extraction, and a fresh Map/Model entry point. Validation results are
recorded under `figures/joint-component-runtime/`; historical archives are
read-only reference inputs. A successful local scope does not change any
historical global qualification.

## Typed runtime numerical core

`src/core/detail/joint_component/` owns the basis kernel, mixed constrained
linear solver, independent reference QR, complete variable-projection
Jacobian, Guarded LM, structural partition and actual-state assembly. Its
interfaces contain typed states, spectra, trial/trust evidence and assessments;
there is no JSON, file access, truth, dataset registry or multiprecision audit
orchestration in that module. The adapted Eigen LM retains its original license.

Testing entry points now serialize this core through `JointRuntimeJson.hpp`.
Legacy and Guarded-log are retained as regression controls in the internal
search API. The public estimator uses Guarded only. Shared dense/sparse linear
solver callers and the older fixed-B experiments use the same extracted
implementation.

Component search/assessment durations are measured separately. The component
adapter records process high-water RSS after the component returns; it is not a
per-phase allocation measurement. Numerical certificates and accepted-state
semantics are unchanged by the extraction.

## Public C++ entry points

Include `<rhbm_gem/core/JointComponentEstimator.hpp>` and link
`RHBM_GEM::rhbm_gem`. No testing target or offline audit library is needed.

```cpp
model.SelectAllAtoms();
auto result = rhbm_gem::core::EstimateJointComponents(map, model);
// Or freeze a problem and provide explicit positive finite widths:
auto problem = rhbm_gem::core::BuildJointProblem(map, model);
auto repeated = rhbm_gem::core::FitJointComponents(problem, initial_b);
```

`JointProblem` owns an immutable snapshot. `JointProblemInput` also permits a
caller to supply frozen observations, identities and structural squared-distance
memberships directly. The version-one numerical policy is fixed: equal-weight
LS, nonnegative A, signed C, log-B, 2.5 Angstrom support and Guarded search with
200 profile evaluations / 100 accepted updates per component. The public fit
uses the existing scoped Eigen thread guard to retain single-thread numerical
behavior, then restores the caller's thread setting.

The Map/Model builder rejects a partial non-hydrogen selection. It constructs a
unique voxel union from the actual Map geometry using `sphere-fma-v1`, retaining
negative and zero observations. Every structural membership is retained,
including zero coefficients and numerically invisible basis values.

`EstimateJointComponents` runs deterministic Fibonacci sampling, first-stage
alpha training and first-stage MDPDE fitting in a model copy. It copies only raw
samples and first-stage fitting results back to the supplied model; selection,
second-stage estimates and group results are preserved. Only B enters joint
initialization. No Peeling result, truth, historical fit or certificate is read.
Invalid widths produce explicit initialization failure instead of a fallback.

Result coefficient order is `[A0, C0, A1, C1, ...]`; width vectors follow the
problem atom identities. First-stage OLS/MDPDE diagnostic arrays are `[A, B, C]`.
Objectives are `0.5 * ||prediction - observations||² / ObservationScale()²`;
component objectives use the same parent scale, and width gradients are with
respect to log-B. Component `state` is the last trusted state actually
used by assembly. Search completion, state availability and evidence status are
separate. An unobserved or failed component prevents a complete prediction and
objective; `available_row_mask` retains usable rows, including constant rows.
There is no zero-filled substitute for missing component predictions.

Runtime performs search trust/reference/replay and operational local/assembled
assessments. Richardson, multiprecision, boundary and regular certificates are
explicitly `NotRun` until an offline adapter supplies that evidence. Runtime
operational checks must not be reported as offline regular certification.

## Fresh-input experimental entry points

The testing executable provides:

```
mdpde_experiment joint-component-runtime MODEL MAP OUTPUT
mdpde_experiment joint-component-physical OUTPUT
mdpde_experiment joint-component-physical-inputs OUTPUT
```

The Python `tests/integration/joint_component_runtime.py` runner exposes `run`,
`physical`, `summarize` and `compare`. Both numerical precisions in the physical
fixture use the same generated Map geometry; float32 observations are rounded
from that map. The fixture contains two copies of the heterogeneous 12-atom
baseline separated by 12 Angstrom and is generated by the production forward
model. The `physical-inputs` mode writes a minimal CIF and float32 MRC for a
separate test of the ordinary file-reading path.

Generation geometry and MRC header geometry are recorded separately. A loaded
map is compared to a fresh monolithic solve of its own observations and support;
it is never forced onto historical generation coordinates. Each fresh case
records four matched starts, first-stage diagnostics, structural census,
observations, local/global audits, actual-state reassembly checks and costs.

The frozen-input regression uses records saved before extraction, not two
wrappers of the extracted core as its oracle. Full end-to-end comparisons use
the shared core with different monolithic/component orchestration and are
supplementary to that independent frozen regression.

An isolated local rerun bundle contains `search_context`, `component_initial_b`,
`parent_observations` and `component_input` (own support and parent mappings).
The standalone executable reconstructs and verifies the parent normalization and
policy, validates the component graph, and records the bundle SHA-256. Its executable does not open historical fit files or
sibling/global endpoint records. The immutable parent observations still
provide the same normalization and active-set response norm.

`joint_component_runtime.py frozen-parity --reference OLD_RUN --reference-audits
OLD_AUDITS --run NEW_RUN --output REPORT` compares all historical search and audit
scientific JSON against the new core, while permitting additional local-scope
files. It requires both reference populations to exist. Timing/RSS and process
completion manifests are excluded; endpoint states, spectra, derivatives,
accepted trials, stopping decisions and existing certificates are compared.

The frozen monolithic regression input is reconstructed from the existing
`joint-abc-components/prerequisite-records.tar.gz` stage-one records plus the
original snapshots in `search-records-a.tar.gz`. This uses the already archived
pre-extraction records and does not require a missing certification archive.

## Completed validation

| Gate | Result |
| --- | --- |
| Frozen monolithic regression | 216 branches; all 432 fit/audit records unchanged; regular counts 38/40/40 for Legacy/Guarded/Guarded-log |
| Frozen component numerical parity | 128 paired branches; 1,248 search and 2,096 existing audit scientific records unchanged |
| Global qualification | All 64 required regular pairs pass; originals 40 and composites 24 on both monolithic and assembled endpoints |
| Component-local certification | 192 scopes, 128 regular; all 1,168 local scientific records unchanged by core extraction |
| Independent full-matrix repeat | Two fresh runs; all 5,008 scientific records agree |
| Physical disconnected fixture | 8 regular paired branches, 16 regular local scopes; 228 repeated scientific records agree |
| CIF/MRC entry | 4 regular paired branches, 8 regular local scopes; 114 repeated scientific records agree, including final API verification |
| Engineering | 72 focused C++ tests, 43 Python tests; testing-disabled production build and installed consumer pass |

Historical globally restricted component scopes have 96 regular certificates;
new component-local scopes have 128. The 32 additional local qualifications occur
in regular-active (8), regular-near (8), regular-weak (8), regular-three (4),
regular-two (2) and regular-zero (2). These are separate scopes; historical
certificates were neither overwritten nor relabeled.

The same-state experiment retains 196/256 complete frozen-state equivalence
records and 160/256 complete endpoint equivalence records. Applicable checks all
pass. The remaining rank, active-face and local-correction limitations are
preserved, rather than counted as complete regular evidence.

The physical fixture's maximum paired scaled A/C difference is
`3.299e-15`, maximum log-B difference `2.22e-16`, normalized prediction difference
`1.405e-16`, and normalized objective difference `9.566e-26`.
The actual MRC-header problem is validated separately against its own fresh
monolithic reference; no historical geometry is substituted.

Isolated validation includes eight difficult-composite baseline reruns in both
precisions with historical fit directories absent. A further single-file bundle
run reproduces the same scientific records in a directory containing only its
input bundle; changed parent observations and invalid initial widths are rejected.
The C++ tests cover before/after sibling execution, reversed fitting order,
missing states, zero inherited directions, and standalone bundle execution.

Detailed certificates, failure matrices, costs, hashes, source provenance and
compressed records are indexed in
[the evidence directory](figures/joint-component-runtime/README.md).
