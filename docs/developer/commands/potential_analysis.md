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
- `selected_chain_ids` (empty means all chains)
- `simulation_flag`
- `simulated_map_resolution`
- `saved_key_tag`

The request catalog generates both CLI options and Python bindings. Potential
analysis fixes map normalization, hydrogen exclusion, sampling and failed-only
refinement as command policies; those policies are not request fields. The
shared `SphereSamplingMethod` enum remains available to generic fitting APIs.
Alpha training is enabled, with initial local/group alpha values set internally
to `0.0`.

## Execution Contract

`rhbm_gem::core::RunCommand(PotentialAnalysisRequest{...})` returns `CommandResult`.

Expected result contract:

- `result.succeeded == true` when the command completes
- `result.succeeded == false` when normalization, validation, preflight, or execution stops the command
- `result.issues` reports public validation diagnostics as option/message pairs

## Command Behavior

The anonymous-namespace `NormalizeAndValidateRequest(...)` phase handles field validation:

- validates required model and map paths
- trims and deduplicates selected chain IDs and rejects empty IDs
- rejects invalid scalar inputs that this command cannot safely recover from
- rejects an empty `saved_key_tag`

`ValidatePreparedRequest(...)` performs semantic checks after normalization:

- `--simulation` requires positive simulated resolution
- loads the model once, verifies requested chain IDs and requires at least one
  non-hydrogen target before map preprocessing or fitting

`ExecutePreparedRequest(...)`:

- reuses the validated model and loads the map through `ReadMap(...)`
- optionally switches the model object into simulation mode
- skips map normalization for simulation requests and always attempts it for
  other requests
- computes Q-scores and fits the target selection without removing atoms from
  the model
- delegates atom sampling to `MapSampler`, then performs alpha training and fitting
- uses Fibonacci deterministic sampling, permanently excludes hydrogen, and
  enables failed-only second-stage refinement
- persists the prepared model through `DataRepository::SaveModel(...)`
- writes to the repository using `request.saved_key_tag` as the persisted key
- clears sampled local-potential distance/value buffers after persistence to keep runtime state lean

`CommandRunner` creates `output_dir` during filesystem preflight when needed.

## Target chains and fixed policies

With no chain selection, all non-hydrogen atoms in the complete model are
targets, including side-chain atoms and chains that share an entity with other
chains. `--select-chain` narrows the target set while preserving the full
`ModelObject`. Chain IDs are compared case-sensitively using
`AtomObject::GetChainID()`, which for mmCIF input prefers `label_asym_id` and
falls back to `auth_asym_id`. Unknown IDs, empty CSV tokens and selections with
no non-hydrogen targets fail before Q-score preprocessing and fitting. Repeated
IDs are deduplicated; token whitespace is trimmed.

C++ callers can assign `request.selected_chain_ids = {"A", "B", "C"}`;
Python callers can assign `request.selected_chain_ids = ["A", "B", "C"]`.
An empty vector or list selects every non-hydrogen atom in the model.

For Joint Components, selected chains define targets and their voxel
observations. Other non-hydrogen atoms remain eligible as halo contributors
when their support overlaps those observations. Chain selection does not crop
the contributor model or change the target/halo distinction.

Simulation maps skip normalization. Every non-simulation map is normalized by
its standard deviation when that operation succeeds; no mean is subtracted. A
zero-SD map records `requested: true`, `applied: false`, and `divisor: 1`.

```sh
# Default: all non-hydrogen atoms in all chains
rhbm_gem potential_analysis \
    --estimator joint-components \
    --model model.cif \
    --map map.mrc

# One target chain
rhbm_gem potential_analysis \
    --estimator joint-components \
    --model model.cif \
    --map map.mrc \
    --select-chain A

# Several target chains
rhbm_gem potential_analysis \
    --estimator joint-components \
    --model model.cif \
    --map map.mrc \
    --select-chain A,B,C
```

## Failed-only endpoint refinement

Potential analysis always enables failed-only refinement for second-stage shape
fitting. The PotentialAnalysis request and CLI do not expose a switch. Direct
fitting callers retain the generic
`FitOptions::enable_second_stage_failed_only_refinement` control. In the Joint
Components workflow only, a failed FirstStage fit for a fixed halo also uses
bounded recovery during initialization; target FirstStage fits, alpha training,
group estimation and generic `EstimateLocalGaussian` / `EstimateBetaMDPDE`
entrypoints keep their native behavior.

A native successful solve is used unchanged, without a fresh-equation evaluation.
A failed solve may be refined from its own endpoint. Acceptance requires fresh
scaled equations at `1e-8`, a valid positive-variance Gaussian, weighted rank and
denominator checks, and agreement with a fixed-point continuation from the
original endpoint. The reference must reach `1e-10`; transformed coordinates and
weights must agree within `1e-6`, with identical weight-floor activation.

The endpoint-refinement candidate budget is 128 equation evaluations including
verification. Reference continuation is separate work, capped at 10,000 total
fixed-point updates including the native iterations. These constants are not
public tuning options. Exact-fit variance boundaries are not promoted to
successful solves. Rejection preserves the original numerical result; the
reference is never used as a fallback.

For failed fixed-halo FirstStage fits, Joint initialization first attempts that
same endpoint refinement. If it is rejected and the native status is
`MAX_ITERATIONS_REACHED`, initialization retries the native fixed-point solver
with a hard 300-iteration cap and the existing tolerances and weight floor. The
retry is accepted only when that solver reports its normal `SUCCESS`; its
diagnostics record retry status, iterations, beta and variance stopping changes,
equation residual and runtime. The original native status and diagnostics remain
unchanged, and the result is identified as `RefinedSuccess`. Other failures stay
unqualified. Successful FirstStage halo fits run no recovery work. Fixed halos
remain immutable at FirstStage A/B with C=0.

`RHBMBetaEstimateResult.status` and `diagnostics` retain native solver history.
Use `Qualification()` to distinguish `NativeSuccess`, `RefinedSuccess` and
`Unqualified`. The optional `refinement` diagnostics record acceptance, rejection
reason, residuals, branch differences and candidate/reference work. Failed-only
recovery diagnostics separately record endpoint-refinement evidence and any
bounded native retry. Joint snapshots persist fixed-halo qualification and
recovery provenance in schema 15; older Joint result schemas remain readable.
Audit shape records retain effective `status` for existing consumers and
additionally expose `native_status`, `qualification` and optional
refinement/recovery evidence.

Recovery continues to require qualified inner endpoints. Refinement does not
relax outer convergence, offset IRLS or the best-objective bound, and does not
itself establish convergence. See the [production fitting contract](../production-fitting.md).

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

The map-aware `RunPotentialFittingWorkflow(map, model, options)` is the canonical
complete workflow. It dispatches by estimator, builds the corresponding fitting
workset, and passes an explicit atom list to the shared sampling engine.

For TWO_STAGE, the workset contains the selected atoms. Their samples feed batch
First-stage fitting, second-stage iterations, and group fitting. The configured
worker count continues to apply to sampling and the existing parallel fitting
operations.

For JOINT_COMPONENTS, the problem fixes the observation voxel domain and the
workset contains its target and contributor atoms. Only FullABC contributors are
sampled and fitted locally; contribution-only halo atoms retain initialization
provenance without local sampling or fitting. Joint initialization consumes
those prepared samples in isolation and supplies First B plus seed provenance to
`FitJointComponents`. Joint estimates then flow through the second-stage
mapping, target summary, post-fit peeling, fixed Joint charge-profile
correction, and fresh local alpha/MDPDE fitting for each selected target.
Group alpha and group WEB consume the qualified local estimates and their
compact MDPDE-derived information. Local fitting does not update Joint
estimates or construct a cross-atom covariance.

`FixedNeighbor` is the sole qualified bounded-memory production route with the fixed
core-12, forward, one-accepted-local-update policy documented in the
[joint runtime contract](../joint-component-runtime.md). Its local profile
solver is Profile LM, the existing profiled-width Levenberg–Marquardt
computation; it is not an automatic fallback or size-based route. Use
`--estimator joint-components`; no separate
Joint search option is required.

The standalone `EstimateJointComponents` convenience API remains available and
keeps its model-copy and target writeback behavior.

The solver contracts remain `guarded-joint-ls-v1`,
`parent-normalized-half-rss-v1`, and `sphere-fma-v1`.

Both estimators use the selected non-hydrogen atoms as targets. Sampling is
always Fibonacci deterministic, and PotentialAnalysis always enables
failed-only second-stage refinement. Joint Components uses one worker; `-j`
values greater than one produce a notice. Selected targets fix the observation
rows. Unselected non-hydrogen halo candidates enter Joint fitting only with at
least three support rows and three exact distinct squared distances. Pruning
keeps selected observation rows fixed and can change fitted target values
because excluded signal remains in the observations.

```cpp
rhbm_gem::core::PotentialAnalysisRequest request;
request.estimator = rhbm_gem::core::PotentialEstimator::JOINT_COMPONENTS;
request.model_file_path = "model.cif";
request.map_file_path = "map.mrc";
request.selected_chain_ids = {"A"}; // omit this or use {} for every non-hydrogen chain
request.database_path = "joint.sqlite"; // new v19 database
request.saved_key_tag = "example";
auto completed = rhbm_gem::core::RunCommand(request);
```

```python
import rhbm_gem_module as gem
request = gem.PotentialAnalysisRequest()
request.estimator = gem.PotentialEstimator.JOINT_COMPONENTS
request.model_file_path = "model.cif"
request.map_file_path = "map.mrc"
request.selected_chain_ids = ["A"]  # [] selects all non-hydrogen atoms
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
`joint_atoms_<sanitized-key>.csv`. It writes the sibling
`joint_atoms_<sanitized-key>.contributions.csv` only for a historical result
with nuisance groups. Joint JSON schema 14 includes metadata, identities,
contributor-selection policy and provenance, target-only parameterization,
FirstStage fixed-halo snapshots, seed provenance, row mappings/mask, initial
values, target-only trusted states, optional parameter-only
candidate diagnostics, objectives, cost counters, checks, legacy global-rank
evidence, local block-conditioning and parameter-stability diagnostics, and captured
convergence. Target-estimability rank evidence remains under target evidence.
CSV has one row per fitted contributor (including unavailable targets):
`AtomID,ComponentID,A,B,C,StateAvailable,SearchCompleted,StopReason,RuntimeConvergence,RegularCertificate,SelectionRole,Parameterization,ContributionGroupID,TargetRuntimeConvergence,ParameterIdentifiability`.
A/C follow the joint kernel convention, with signed C; B is the width, not log-B.
Available but unconverged trusted states are retained. Uncertified candidate
diagnostics appear only in JSON and never fill scientific state or CSV values.
Missing estimates have empty CSV fields and null JSON states. Initialization diagnostics with nonfinite numbers
use null; they are not converted to zero. Numeric serialization preserves finite
double values on reload.

Export needs only the saved model. It fails on missing joint records, colliding
sanitized keys or write errors. Gaussian/outlier export, display, comparison and UMAP consume the common stage
contract and enable only available data. Joint values are not written to OLS or
local MDPDE columns. UMAP excludes rows missing required features and needs at
least three valid targets.

SQLite v20 is created for new databases. Reading v17/v18/v19 does not modify them; the first
Save upgrades and writes in one transaction, rolling back on failure. Older
versions, including v16, remain unchanged on rejection. A saved key holds one joint outcome;
saving a model without a joint result over that key removes the previous outcome.

### Provenance and map units (joint JSON schemas 3 through 14)

`metadata.model_sha256` and `map_sha256` fingerprint the original file bytes,
checked before and after loading. Paths remain descriptive, not content identity.
`metadata.software` records the library version and source, configuration and
build SHA-256 fingerprints from the existing build-time generator. Loading or
exporting a saved outcome preserves these values; it does not stamp the reader's
version or access either input file.

`metadata.solver` records the current `fixed-neighbor` estimator. FixedNeighbor records
`fixed_neighbor_core_atoms: 12` and the intrinsic at-most-one trusted accepted
Profile LM update per block visit. Current writers use
`fixed-neighbor-production-v4`; `fixed_neighbor_local_work` is absent or null,
because it is legacy provenance rather than a current policy input. They also
persist the `joint-solver-provenance-v5` contract, fixed EIGEN sparse
backend, and FixedNeighbor's maximum sweeps, forward order and
`fixed_neighbor_local_search: profile-lm` local search. Readers accept
`legacy-compact` as the previous v3 value for that local-search field.
Operator, rank, preconditioner and Schwarz
controls are historical provenance only and are not emitted by current
writers. Historical `joint-solver-provenance-v2` records, including
top-level `search_method` values `legacy-compact` and `operator-pcg`, remain
readable through the compatibility decoder.

Older `fixed-neighbor-production-v1` records remain readable and preserve their
explicit `full`, `one-accepted`, or `two-accepted` local-work provenance. A v2
record may omit that field or carry transitional `one-accepted`, but explicit
v2 `full` and `two-accepted` claims are rejected; the decoder never silently
converts historical provenance into the current contract.

`metadata.map_normalization` is `{requested, applied, divisor}`. The operation is
`fit_map = input_map / divisor`, with no mean subtraction. An applied operation
records the actual pre-normalization standard deviation, including SD=1.
`requested` records the fixed command policy: false for simulation requests and
true otherwise. `applied` records whether normalization succeeded; a zero-SD
non-simulation map therefore records `requested=true, applied=false, divisor=1`.
Simulation requests record `requested=false, applied=false, divisor=1`.
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

Production joint JSON schemas 3 through 14 are accepted. Earlier joint JSON is rejected
with a request to regenerate the outcome, without migration or database writes.
New Joint result JSON uses schema 14. Schemas 3 through 13 preserve their saved
evidence and convergence statuses. Schema 10 stores `target-estimability-v4`
and rank-context identities without computed null-space leakage or full-Jacobian
column scales. Schema 11 adds the retired FullABC halo pruning policy and
exclusion records. Schema 13 identifies the FirstStage fixed-halo contributor
policy; schema 14 stores validated fixed `(A, B, 0)` halo snapshots separately
from target component states using `target-only-fixed-halo-v1`. In particular,
failed historical local or
target-local correction checks remain failed on load; the reader does not
rerun those solves. Schemas 8 through 11 record the correction-free numerical
contract and store endpoint-verification and numerical-qualification
diagnostics separately. Legacy reference agreement is not mapped to conditional KKT, and
no reference status is mapped to independent replay.
See the [target contract](../joint-target-estimability.md)
for target convergence and retained halo availability. SQLite v20 stores
estimator-neutral stage format 4, including `PostJointLocalResult`, paired
peeling, compact local information, posterior evidence and sample geometry in
`model_stage_result`. Legacy sample BLOBs have unavailable geometry; migrated
Joint results whose group evidence used the previous measurement source are
marked unavailable until the new local workflow runs. Neutral stage format 3
is read and converted to format 4 using each available local result's fixed
Joint charge profile.
Consumers must rebuild against the updated public C++ value types.

### Selected-domain and initialization metadata

Current `selection_domain` records `contract: fixed-selected-voxel-closure-v1`,
ordered unique `target_indices` into `atom_ids`, `observation_radius: 2.5`,
`support_radius: 2.5`, and `contributor_policy:
first-stage-fixed-halo-v1`. Every non-hydrogen halo with at least one support
row on the unchanged selected-target observation domain is retained. Its
validated FirstStage amplitude and width, fixed zero charge, sample count,
status, and source contract are stored in the result's `fixed_halos` array.
Observation support counts remain separate from FirstStage sampling counts.
Older payloads carrying `fullabc-eligible-halo-pruning-v1` and exclusion
records retain their original free-halo meaning. A hand-built input with the legacy `all-non-hydrogen` selection
contract remains readable; an input without selection metadata retains null and
no roles are inferred. CSV `SelectionRole` is `target`, `halo`, or
`not-recorded`; halos in current fixed-halo results are FirstStage-sourced and
are not Joint-optimized estimates.

`initialization.data_scope` is `caller-provided-widths` or
`contributor-local-sampling-may-read-outside-target-domain`. Each initialization
atom retains its reason; a failed component retains an unavailable state while
independent valid components can run. Current fixed-halo results require valid
FirstStage halo amplitudes and widths; an invalid halo is reported as
`fixed-halo-initialization-unavailable`, never replaced by zero or a median
width. Target initial-width fallback diagnostics retain their existing
provenance. See the
[halo eligibility contract](../joint-observable-halo.md). `costs.construction_seconds` records builder time;
initializer time includes model-copy setup. Neither storage nor export recomputes
roles, initialization or numerical evidence. Original selection and halo legacy
analysis are preserved; successful target first-stage data can be updated.
