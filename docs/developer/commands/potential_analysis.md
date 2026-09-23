# Potential Analysis Command

## Registration and Implementation Files

The command membership entry lives in
the internal command catalog in
[`src/core/command/detail/CommandCatalog.hpp`](/src/core/command/detail/CommandCatalog.hpp).

Public request and entrypoint:

- [`include/rhbm_gem/core/CommandSystem.hpp`](/include/rhbm_gem/core/CommandSystem.hpp)

Internal schema and wiring:

- [`src/core/command/detail/CommandCatalog.hpp`](/src/core/command/detail/CommandCatalog.hpp)
- [`src/core/command/CommandSystem.cpp`](/src/core/command/CommandSystem.cpp)
- [`src/python/CommandSystemBindings.cpp`](/src/python/CommandSystemBindings.cpp)

Concrete implementation:

- [`src/core/command/PotentialAnalysisCommand.cpp`](/src/core/command/PotentialAnalysisCommand.cpp)
- [`src/core/MapSampler.cpp`](/src/core/MapSampler.cpp)

## Request Surface

`PotentialAnalysisRequest` is a plain DTO that inherits shared fields from `CommandRequestBase`.

Shared fields:

- `job_count`
- `verbosity`
- `output_dir`

Command-specific fields:

- `estimator` (default `PotentialEstimator::TWO_STAGE`)
- `database_path`
- `model_file_path`
- `map_file_path`
- `simulation_flag`
- `map_normalization_flag`
- `exclude_hydrogen`
- `simulated_map_resolution`
- `saved_key_tag`
- `asymmetry_flag`
- `sampling_method`
- `enable_second_stage_failed_only_refinement`

`sampling_method` uses the shared `SphereSamplingMethod` enum and is exposed
through the `--sampling-method` CLI flag.
Alpha training is always enabled for potential analysis. Initial local/group
alpha values are internal command defaults set to `0.0` and are not exposed as
request fields or CLI options.

## Execution Contract

`rhbm_gem::core::RunCommand(PotentialAnalysisRequest{...})` returns `CommandResult`.

Expected result contract:

- `result.succeeded == true` when the command completes
- `result.succeeded == false` when normalization, validation, preflight, or execution stops the command
- `result.issues` reports public validation diagnostics as option/message pairs

## Command Behavior

The anonymous-namespace `NormalizeAndValidateRequest(...)` phase handles field validation:

- validates required model and map paths
- leaves map normalization as an execution-time choice: simulation requests always skip it,
  and non-simulation requests use `map_normalization_flag`
- rejects invalid scalar inputs that this command cannot safely recover from
- rejects an empty `saved_key_tag`

`ValidatePreparedRequest(...)` performs semantic checks after normalization:

- `--simulation` requires positive simulated resolution
- fit ranges must be ordered correctly

`ExecutePreparedRequest(...)`:

- loads command-owned model and map objects through `ReadModel(...)` and `ReadMap(...)`
- optionally switches the model object into simulation mode
- optionally runs map normalization, then runs model preprocessing
- delegates atom sampling to `MapSampler`, then performs alpha training and fitting
- persists the prepared model through `DataRepository::SaveModel(...)`
- writes to the repository using `request.saved_key_tag` as the persisted key
- clears sampled local-potential distance/value buffers after persistence to keep runtime state lean

`CommandRunner` creates `output_dir` during filesystem preflight when needed.

## Failed-only endpoint refinement

Second-stage shape fitting enables failed-only refinement by default. Disable it
with `--second-stage-failed-only-refinement false`, or set the request field
`enable_second_stage_failed_only_refinement` to `false` in C++ or `False` in
Python. Direct fitting callers use the matching `FitOptions` field.
First-stage fitting, alpha training, group estimation and the generic
`EstimateLocalGaussian` / `EstimateBetaMDPDE` entrypoints keep their native behavior.

A native successful solve is used unchanged, without a fresh-equation evaluation.
A failed solve may be refined from its own endpoint. Acceptance requires fresh
scaled equations at `1e-8`, a valid positive-variance Gaussian, weighted rank and
denominator checks, and agreement with a fixed-point continuation from the
original endpoint. The reference must reach `1e-10`; transformed coordinates and
weights must agree within `1e-6`, with identical weight-floor activation.

The candidate budget is 128 equation evaluations including verification. Reference
continuation is separate work, capped at 10,000 total fixed-point updates including
the native iterations. These constants are not public tuning options. Exact-fit
variance boundaries are not promoted to successful solves. Rejection preserves
the original numerical result; the reference is never used as a fallback.

`RHBMBetaEstimateResult.status` and `diagnostics` retain native solver history.
Use `Qualification()` to distinguish `NativeSuccess`, `RefinedSuccess` and
`Unqualified`. The optional `refinement` diagnostics record acceptance, rejection
reason, residuals, branch differences and candidate/reference work. Audit shape
records retain effective `status` for existing consumers and additionally expose
`native_status`, `qualification` and optional `refinement` evidence. No database
schema change is required.

Recovery continues to require qualified inner endpoints. Refinement does not
relax outer convergence, offset IRLS or the best-objective bound, and does not
itself establish convergence. See the [production validation](../failed-only-refinement.md).

## Internal fitting ranges

Fitting ranges are internal constants in `src/core/detail/gaussian_fit/FittingRanges.hpp`.
Signal fitting and alpha training use `[0, 1.0]` Å; the second-stage tail
objective uses `[1.2, 2.0]` Å. Both boundaries are inclusive. The two ranges
are independent and may overlap or leave a gap when the constants are changed.
An overlapping sample contributes to both objective terms, each with its own
scale and sample-count normalization. Samples in neither region contribute to
neither residual objective term.

Offset fitting still uses all raw samples, and joint polish keeps its existing
sample sources. The gap is therefore not excluded from every estimation step;
tail is an objective constraint region, not a strictly held-out validation set.
Neither PotentialAnalysis nor RHBMTest accepts `--fit-min` or `--fit-max`.

## Tests to Update When Behavior Changes

- [`tests/core/command/CommandScenarios_test.cpp`](/tests/core/command/CommandScenarios_test.cpp)
- [`tests/core/contract/CommandExecutionContract_test.cpp`](/tests/core/contract/CommandExecutionContract_test.cpp)
- [`tests/core/contract/CommandCatalog_test.cpp`](/tests/core/contract/CommandCatalog_test.cpp)
- [`tests/integration/CommandApiPipeline_test.cpp`](/tests/integration/CommandApiPipeline_test.cpp)
- [`tests/integration/python_bindings_runtime_smoke.py`](/tests/integration/python_bindings_runtime_smoke.py)

## Joint component opt-in

`--estimator joint-components` selects the map-aware `RunPotentialFittingWorkflow`.
It builds a fixed contributor workset, shares sampling and formal First fitting,
then passes only First B to `FitJointComponents`. The standalone
`EstimateJointComponents` convenience API remains available, but the command does
not call it. Joint points are mapped by identity into estimator-neutral Second
records, followed by the target summary, grid-consistent post-fit peeling and
parameter-evidence group inference. The model-only overload remains two-stage;
Joint requires map geometry.

The solver contracts remain `guarded-joint-ls-v1`,
`parent-normalized-half-rss-v1`, and `sphere-fma-v1`.

Joint uses selected non-hydrogen targets and the complete fixed-domain halo
closure, deterministic Fibonacci initialization and one worker. Hydrogen is
excluded regardless of `--exclude-hydrogen`. Backbone-only and asymmetry flags
select targets; unselected non-hydrogen atoms remain eligible contributors.
Non-Fibonacci initialization requests are rejected. `-j` values greater
than one produce a notice that joint uses one worker. The second-stage refinement
option applies only to the two-stage estimator. Map normalization and Q-score
preprocessing retain the existing command semantics: simulation skips map
normalization; other requests honor `--map-normalization`.

```cpp
rhbm_gem::core::PotentialAnalysisRequest request;
request.estimator = rhbm_gem::core::PotentialEstimator::JOINT_COMPONENTS;
request.model_file_path = "model.cif";
request.map_file_path = "map.mrc";
request.database_path = "joint.sqlite"; // new v18 database
request.saved_key_tag = "example";
auto completed = rhbm_gem::core::RunCommand(request);
```

```python
import rhbm_gem_module as gem
request = gem.PotentialAnalysisRequest()
request.estimator = gem.PotentialEstimator.JOINT_COMPONENTS
request.model_file_path = "model.cif"
request.map_file_path = "map.mrc"
request.database_path = "joint.sqlite"
request.saved_key_tag = "example"
assert gem.RunCommand(request).succeeded  # execution/persistence, not convergence
export = gem.ResultDumpRequest()
export.database_path = request.database_path
export.model_key_tag_list = ["example"]
export.printer_choice = gem.PrinterType.JOINT_ESTIMATES
export.output_dir = "results"
assert gem.RunCommand(export).succeeded
```

A returned initialization failure or missing/unconverged component is a saved
outcome, not an execution error. `CommandResult.succeeded` is true after successful
persistence. The log and saved record separately expose search completion, stop
reasons, available components, runtime convergence and regular certificate.
Invalid input/problem construction and persistence errors fail the command.
`NotRun` offline evidence never implies passing certification.

`result_dump --printer joint` writes `joint_result_<sanitized-key>.json` and
`joint_atoms_<sanitized-key>.csv`. JSON schema 3 includes metadata, identities,
row mappings/mask, initial values, actual states, objectives, cost counters,
checks, ranks and captured convergence. CSV has one row per fitted contributor (including unavailable targets):
`AtomID,ComponentID,A,B,C,StateAvailable,SearchCompleted,StopReason,RuntimeConvergence,RegularCertificate,SelectionRole`.
A/C follow the joint kernel convention, with signed C; B is the width, not log-B.
Available but unconverged states are retained. Missing estimates have empty CSV
fields and null JSON states. Initialization diagnostics with nonfinite numbers
use null; they are not converted to zero. Numeric serialization preserves finite
double values on reload.

Export needs only the saved model. It fails on missing joint records, colliding
sanitized keys or write errors. Gaussian/outlier export, display, comparison and UMAP consume the common stage
contract and enable only available data. Joint values are not written to OLS or
local MDPDE columns. UMAP excludes rows missing required features and needs at
least three valid targets.

SQLite v18 is created for new databases. Reading v17 does not modify it; the first
Save upgrades and writes in one transaction, rolling back on failure. Older
versions, including v16, remain unchanged on rejection. A saved key holds one joint outcome;
saving a model without a joint result over that key removes the previous outcome.

### Provenance and map units (joint JSON schema 3)

`metadata.model_sha256` and `map_sha256` fingerprint the original file bytes,
checked before and after loading. Paths remain descriptive, not content identity.
`metadata.software` records the library version and source, configuration and
build SHA-256 fingerprints from the existing build-time generator. Loading or
exporting a saved outcome preserves these values; it does not stamp the reader's
version or access either input file.

`metadata.map_normalization` is `{requested, applied, divisor}`. The operation is
`fit_map = input_map / divisor`, with no mean subtraction. An applied operation
records the actual pre-normalization standard deviation, including SD=1. Disabled
normalization, zero-SD maps and simulation requests use `applied=false, divisor=1`;
`requested` preserves the user's flag, and `simulation` identifies that bypass.
Pure in-memory API callers may leave the normalization record and input hashes
null; unknown provenance does not imply an identity transform.

The fixed `metadata.units` contract is `joint-kernel-map-units-v1`. For fitted map
unit U, the volume-normalized Gaussian has units Angstrom^-3 and the charge basis
`erf(r/(sqrt(2)*B))/r` has units Angstrom^-1. Thus A has units U*Angstrom^3, C has
units U*Angstrom, and B and geometry use Angstrom. U denotes the supplied fitting
scale, not an assumed physical calibration. To return to input-map scale,
multiply both A and C by `divisor`; B is unchanged. `observation_scale` normalizes
the objective only and must not be used for this conversion. CSV retains fitted
coefficients; its appended `SelectionRole` distinguishes target, halo and not-recorded; keep its companion JSON for units and
provenance. No conversion is performed during export.

Only production joint JSON schema 3 is accepted. Older joint JSON is rejected
with a request to regenerate the outcome, without migration or database writes.
The Joint snapshot JSON remains schema 3; SQLite v18 additionally stores neutral
stages, uncertainty, posterior evidence and sample geometry in `model_stage_result`.
Legacy sample BLOBs have unavailable geometry; old Joint snapshots do not gain
recomputed uncertainty, peeling or posterior.
Consumers must rebuild against the updated public C++ value types.

### Selected-domain and initialization metadata

`selection_domain` records `contract: fixed-selected-voxel-closure-v1`, ordered
unique `target_indices` into `atom_ids`, `observation_radius: 2.5`,
`support_radius: 2.5`, and `contributor_policy: all-non-hydrogen`. All other
contributors are halo atoms. `row_ids` are the fixed target observations. A
hand-built input without selection metadata retains null; no roles are inferred.
CSV `SelectionRole` is `target`, `halo`, or `not-recorded`. All contributor states
and diagnostics are saved, regardless of role or convergence.

`initialization.data_scope` is `caller-provided-widths` or
`contributor-local-sampling-may-read-outside-target-domain`. Each initialization
atom retains its reason; a failed component retains an unavailable state while
independent valid components can run. `initialization.valid` still means all
initial widths are valid. `costs.construction_seconds` records builder time;
initializer time includes model-copy setup. Neither storage nor export recomputes
roles, initialization or numerical evidence. Original selection and halo legacy
analysis are preserved; successful target first-stage data can be updated.
