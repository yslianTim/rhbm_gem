# Experiment retirement baseline

This document records the responsibility boundary before any experiment, validation runner, benchmark, or audit is retired. It is an inventory, not a removal proposal. The snapshot was taken from `develop` at `875f7aeaa3177ba0bde8227ddec27372a497639c` (2026-09-29); `origin/develop` was already at that commit and the pre-change worktree was clean.

No production source, numerical tolerance, CMake test registration, fixture, runner, or historical result was changed to make this baseline. The machine-readable build and test record is in [`figures/experiment-retirement-baseline/baseline.json`](figures/experiment-retirement-baseline/baseline.json); the CTest and generated target inventories are in the same directory.

## Classification and reading rules

- `KEEP`: permanent correctness coverage, an independent numerical oracle, an active scientific experiment, or an important property with no safe substitute.
- `MIGRATION_REQUIRED`: the campaign can eventually leave the default branch, but named checks, helpers, or provenance must first move to a permanent test or shared support. The row lists those moves.
- `READY_FOR_RETIREMENT`: no unique correctness responsibility remains; permanent owners exist; remaining use is historical campaign replay or a bounded timing harness whose results are already retained. This only marks a later candidate. It does not authorize removal in this baseline change.

The matrix distinguishes numerical/scientific correctness, workflow and persistence, resource measurement, and historical reproducibility. A Python `*_test.py` registered with CTest usually tests the runner's input/provenance/report rules; it does not mean CTest runs the full scientific campaign. External campaigns that need untracked CIF/map inputs are recorded as not run, not passed.

## Experiment responsibility matrix

| ID | Tool / entry point | Primary purpose | Numerical / scientific property protected | Workflow / performance / historical property protected | Fixture / input | Independent oracle | Current automated test | Unique vs historical-only responsibility | Candidate successor | Retirement readiness | Blocking reason / required migration | Artifact / external-data dependency | Notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S2-01 | Fold-168: `tests/integration/fold_168_regression.py`; gated CTest `fold_168_simulation_regression` | Fixed end-to-end external regression | Simulation output, estimator quality, and stored baseline comparison | Full CLI/database workflow and reproducibility of the fixed 168-atom case | `tests/benchmarks/fold_168_simulation_baseline.json`; external CIF, map, and sidecar manifest | Frozen simulation truth/scorer; not an independent production implementation | `fold_168_regression_runner_test` checks runner contracts; full regression is optional | External scientific regression remains unique; exact historical numbers are campaign-specific | Keep the external regression gate and frozen input provenance | KEEP | External input is not present in this checkout; do not infer a pass | Baseline JSON present; CIF/map/manifest external | `RHBM_GEM_ENABLE_FOLD_168_REGRESSION=OFF` in this baseline |
| S2-02 | MDPDE / forward experiment: `tests/integration/mdpde_experiment.py`, `tests/experiments/mdpde_experiment.cpp`, `tests/support/MDPDEExperiment.*`, `ForwardModelExperiment.*` | Compare forward-model and robust-fit behaviors; capture/replay solver failures | Forward discrepancy, root/branch behavior, MDPDE convergence and captured failure classification | Reproducible campaign capture, provenance, and failure replay | `tests/fixtures/mdpde/`; optional external model/map/manifest | Forward replay and stored branch evidence; not a general independent scientific truth set | `mdpde_experiment_test.py`; C++ `MDPDEExperiment_test` | Research comparisons and real-input results are historical; failure replay helpers are reused by refinement | Retain as research tool; keep failure replay in permanent tests | KEEP | None for synthetic tests; full external-input experiments remain input-dependent | `figures/mdpde-experiment/` present; external inputs vary by campaign | Shared by endpoint and failed-only refinement runners; do not remove helpers with the orchestration |
| S2-03 | Endpoint refinement: `tests/integration/endpoint_refinement.py`; `EndpointRefinementExperiment.*`; `MDPDEEndpointRefinement.*` | Compare bounded branch-preserving endpoint refinement against saved fold captures | Candidate root qualification, wrong-root rejection, fresh residual/weight consistency, covariance and budget outcomes | Historical policy comparison and exact solve-count/cost reports | `tests/fixtures/mdpde/` and saved solver captures; larger captures are external | Branch/reference replay is useful but not a broad independent oracle | `endpoint_refinement_runner_test`; `ProductionFitting_test`; `MDPDEExperiment_test` | All named production policy checks have permanent C++ owners; exact call counts and fold-wide outcomes are campaign-only | Keep `ProductionFitting_test` and `MDPDEExperiment_test`; retain shared fixture support | READY_FOR_RETIREMENT | No numerical migration blocker found. Keep the C++ tests and fixture records; exact campaign call counts are not a correctness contract | Endpoint campaign outputs removed from the worktree; hashes and source commits are recorded in the artifact manifest; captures may be external | Runner imports neutral `mdpde_experiment_support.py`; hash calls use `experiment_io.py` |
| S2-04 | Failed-only refinement: `tests/integration/failed_only_refinement.py` | Compare production failed-only policy against a preserved fold run | Successful native solves remain unchanged; failed solve is refined only with branch/residual qualification | Full-run records, saved-state equality, and historical counts | Preserved fold-168 records and `tests/fixtures/mdpde/` | Paired saved-run comparison; permanent C++ tests cover policy behavior | `ProductionFitting_test`; `MDPDEExperiment_test`; Python full-capture runner is not a CTest entry | Full historical population and exact counts are campaign-only; native-success, failed-only, branch, fresh-weight, and budget contracts have permanent tests | `ProductionFitting_test` / `MDPDEExperiment_test` | READY_FOR_RETIREMENT | No numerical migration blocker found. Keep fixture records and permanent policy tests; exact campaign counts are not a correctness contract | Failed-only campaign output removed from the worktree; its hash and source commit are recorded in the artifact manifest; preserved captures may be external | Imports neutral `fold_168_support.py` and `mdpde_experiment_support.py` |
| S2-05 | Second-stage audit: `resources/tools/developer/second_stage_audit.py` | Parse and summarize decision-evidence records | Does not decide fit correctness or alter estimator decisions | Log completeness, explicit unavailable/skipped outcomes, and diagnostic presentation | Captured audit logs or user-provided logs | None; this is an observer/parser | `second_stage_audit_parser_test` CTest | Parser contract is useful; historical traces are reproducibility evidence | Keep the parser test while the audit format is supported | KEEP | None | Historical traces in `docs/developer/audit-history/` | The audit is not an independent numerical oracle |
| S2-06 | Audit neutrality: `tests/integration/second_stage_neutrality_test.py`; `NumericalAuditNeutrality_test.cpp` | Compare second-stage numerical behavior with audit observation enabled/disabled | Audit ON/OFF numerical neutrality | Full-run persisted records and observer side effects | Same saved baseline inputs in two builds | Paired execution, not an independent solution oracle | C++ `SecondStageNumericalProbe` tests; Python runner is manual and not CTest-registered | C++ tests cover local observer neutrality; the full two-build comparison is not currently an automatic test | Permanent ON/OFF workflow test | MIGRATION_REQUIRED | Promote the full persisted-output ON/OFF comparison before retiring the Python runner | Saved captures/builds may be external | Python script currently lives beside integration tests but is not registered in `tests/CMakeLists.txt` |
| JC-01 | Runtime: `joint_component_runtime` executable and `tests/integration/joint_component_runtime.py` | Exercise the public joint runtime on synthetic, catalog, and physical cases | Runtime convergence, state/evidence parity, backend behavior | Database and complete regression workflow | `tests/fixtures/joint_component/catalog.json` and its archives; optional physical data | Fixture replay/cross-build parity; the separate offline audit is the independent oracle | `joint_component_regression`, `joint_component_physical_smoke`, `joint_component_runtime_runner_test`, joint C++ tests | Permanent runtime/workflow regression | Keep runtime regression as permanent owner | KEEP | Unique live runtime contract | Local compressed fixtures; physical case input is optional | Cost/timing fields are excluded from scientific equality; benchmark timing is not a correctness assertion |
| JC-02 | Offline audit: `joint_component_audit` and `tests/integration/joint_component_audit.py` | Verify joint derivatives, fixed-state evidence, and certification | Independent high-precision/dense reference comparisons and derivative checks | Reproducible local/two-step audit path | Standalone bundle or `joint_component` fixture catalog | Yes: multiprecision/dense global reference checks in `joint_offline_tests` | `joint_component_offline_tests`, `joint_component_offline_regression`, `joint_component_two_step_regression` when enabled | Independent oracle is unique and must survive any runner change | Keep `joint_offline_tests` and audit support | KEEP | None | Local fixture bundles; output directories are generated | Available only with `RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON` |
| JC-03 | Workflow CLI smoke: `tests/integration/joint_workflow_cli_smoke.py` | Exercise the user-visible joint CLI through save/reload/export | No independent solve oracle; checks returned/persisted values and errors | CLI, SQLite/JSON/CSV, reload, export, and target-only summary workflow | Small generated database and local synthetic fixtures | No; workflow round-trip is the contract | `joint_workflow_cli_smoke_test` | Permanent workflow owner | Keep the CTest | KEEP | None | Generated temporary databases only | Separate from full-scale resource benchmarks |
| JC-04 | Partial selection measurement: `joint_partial_selection` and `tests/integration/joint_partial_selection.py` | Record full/partial/bridge/weak synthetic selected-domain behavior and runtime/RSS | Selection, contributor, halo, and target behavior is owned by permanent C++ tests | Small synthetic timing/RSS report; explicitly not a size/resource guarantee | `tests/support/JointPartialSelection.hpp` | No separate oracle in the measurement executable; dense/reference checks live in permanent tests | No CTest invokes this executable; `PartialSelection_test.cpp`, `ObservableProfile_test.cpp`, and CLI smoke own the contract | Numerical/workflow checks have permanent owners; only bounded campaign reproduction/timing remains | `PartialSelection_test.cpp`, `ObservableProfile_test.cpp`; keep shared fixture support | READY_FOR_RETIREMENT | Before a later removal, retain the acceptance evidence and keep the shared fixture header; do not remove fixture support used by unit tests | `figures/joint-component-partial-selection-acceptance/` present | The target is built by `tests_all`, but the Python measurement runner is not a CTest entry |
| JC-05 | Postprocessing benchmark: `joint_postprocessing_benchmark` and `tests/integration/joint_postprocessing_benchmark.py` | Measure workflow/save and fixed-endpoint postprocessing/save phases | Endpoint, covariance, peeling, and target diagnostics are checked for exact equality in paired runs | Timing/RSS only; each case uses independent processes and three repetitions | Synthetic full, halo, and multi-component maps | No; equality is a paired workflow check, not an independent solve | No CTest executes the benchmark; C++ partial-selection/observable-profile tests and CLI smoke cover behavior | Numerical/workflow properties have permanent owners; only the bounded timing campaign is historical | Preserve phase definitions and report, then consolidate or remove | READY_FOR_RETIREMENT | Retain `joint-analysis-consolidation-benchmark.json` and its timing/RSS scope if retiring; do not present it as a resource guarantee | Historical report present; benchmark creates temporary SQLite files | Timing is excluded from correctness status; RSS includes fixture construction |
| JC-06 | Public API benchmark: `joint_component_benchmark` | Measure public joint API construction/search/assembly and peak RSS on supplied datasets | Output availability/objective are reported; no reference comparison is performed | Performance/resource measurement across named dataset/case inputs | Dataset bundle with `dataset.json`, `cases.json`, `voxels.csv`, `contributors.csv` | No | No CTest; runtime C++/Python regressions cover correctness separately | Performance measurement remains the primary purpose | Retain as research benchmark or replace with a defined successor | KEEP | No permanent performance owner; do not mistake successful timing execution for numerical validation | Input bundle is supplied externally by campaign | Public API only; executable is intentionally usable against pre-change libraries |
| JC-07 | Sparse benchmark driver: `joint_sparse_benchmark` | Run timed/fixed/search/rank-oracle actions used by validation campaigns | Cross-backend parity and rank evidence are checked by campaign comparators; the driver itself is not an independent oracle | Instrumentation, build fingerprinting, timing, work counts, watchdog and RSS receipts | Joint fixture or archived campaign input bundles | Dense rank oracle for bounded cases; otherwise paired backend/build comparisons | Joint runtime and validation-runner contract tests; no direct CTest of all actions | Driver and shared comparison behavior are dependencies of several historical campaigns | Extract reusable parity/provenance checks to permanent support before consolidation | MIGRATION_REQUIRED | Migrate parity, work-count, rank-status, and fingerprint helpers; keep fixtures and establish benchmark timing semantics | Sparse/fixed/operator manifests and archives in figures | It is an experiment executable, not a correctness test by virtue of being built |
| JH-A | Joint validation Stage A: `tests/integration/joint_validation.py --stage a`; `joint_validation` offline executable | Diagnose weak-halo attribution and identifiability | Weak spectrum, target-vs-halo identifiability, and local numerical diagnosis | Reproducible diagnosis over named cases | Joint fixture catalog; optional campaign input bundle | Offline audit can provide independent local references; the statistical campaign itself is not an oracle | `joint_validation_runner_test`; offline tests when enabled | Numerical diagnosis has research value independent of Stages B/C | Keep Stage A separate and retain its numerical checks | KEEP | Do not retire because the same Python runner also launches Stages B/C | Local fixtures plus optional external campaign inputs | Stage classification is by scientific purpose, not by shared runner |
| JH-B | Joint validation Stage B: `tests/integration/joint_validation.py --stage b` | Measure noise and position-mismatch behavior | Statistical/scientific outcomes under perturbation | Repeated experiment summaries and uncertainty | Generated noise/mismatch conditions and named joint cases | No deterministic numerical oracle; compare distributions and predeclared summaries | Runner contract CTest only; campaign itself is manual | Statistical experiment is a distinct research responsibility | Keep Stage B separate from A/C | KEEP | Statistical design and conclusions need an explicit successor before retirement | `figures/joint-validation/noise-runs.json` and summary artifacts present; case inputs vary | Do not infer numerical correctness from a statistical pass flag |
| JH-C | Joint validation Stage C: `tests/integration/joint_validation.py --stage c` and `joint_validation_profile.py` | Measure complete-command resource envelope | Workflow completion and persisted result correctness are separate from timing | Wall time, process-tree RSS, timeout/killed/incomplete semantics | Full command model/map/manifest inputs, often external | No; C++/CLI workflow tests own correctness, the harness owns measurement semantics | `joint_validation_runner_test` checks status semantics with mocked RSS; `joint_workflow_cli_smoke_test` owns workflow behavior | Complete-command resource measurement is not duplicated by unit tests | Consolidated resource benchmark with fixed timing/RSS semantics | MIGRATION_REQUIRED | Preserve preparation/solve/postprocessing/command timing definitions, RSS/timeout semantics, and never count incomplete/killed/not-run as pass | `figures/joint-validation/resources.json` and campaign receipts; full external inputs not present | Do not retire together with Stage A/B merely because one runner dispatches them |
| JH-01 | Reference campaign: `joint_reference_validation.py` | Compare endpoint/reference behavior against frozen and incremental runs | Reference-solve parity, complete endpoint and zero-reference-work assertions | Historical campaign and source/input provenance | Archived/input bundles and exported SQLite/JSON/CSV | No separate implementation oracle; baseline-vs-candidate replay | `joint_validation_runner_test` checks missing/censored/incomplete cases | Commit/source-pinned campaign comparison; not a permanent correctness suite | C++ `Components_test`, `Numerics_test`, offline reference tests | MIGRATION_REQUIRED | Preserve reference-specific acceptance and provenance checks in permanent tests before retiring | Six reference archives and manifests present; see artifact manifest | Pins baseline `f1ac45c…` and candidate `62d178d…`; uses neutral comparison and compact-summary support |
| JH-02 | Compact SVD campaign: `joint_compact_validation.py` | Compare compact-SVD modes and measure reference/free-design SVD work | Rank/status, spectrum, threshold and solution parity | Fixed-case timing and campaign provenance | `single-128`, `heterogeneous-168`, `single-512` bundles | Comparisons use saved/reference mode results; no standalone independent oracle in the runner | `joint_validation_runner_test` exercises compact/reference report contracts | Historical source fingerprint and exact timing comparisons | Permanent compact-SVD/rank tests plus offline dense reference | MIGRATION_REQUIRED | Preserve compact-SVD acceptance and provenance checks in permanent tests before retiring | Compact/reference acceptance artifacts present | Pins source fingerprint `97185310…`; uses neutral comparison and compact-summary support |
| JH-03 | Sparse backend acceptance: `joint_sparse_validation.py` | Compare EIGEN/SPQR sparse backend runs and costs | Solution/status parity, derivative and repeated-primary behavior | Backend build provenance, timing, RSS and source/input fingerprints | Fixture/catalog or campaign bundles | Paired backend results; not a separate high-precision oracle | No direct CTest runner entry; shared comparison primitives have focused tests | Shared parity/status comparison is now in neutral support and is also used by compact/operator validation; not tied to a literal historical commit constant | Permanent backend/parity tests and shared test support | MIGRATION_REQUIRED | Preserve sparse-specific acceptance, backend coverage, and input/build provenance in permanent tests; keep fixtures | Sparse acceptance artifacts present; SPQR backend not run in current baseline | No fixed commit constant in this file; provenance calculations are neutral support |
| JH-04 | Operator campaign: `joint_operator_validation.py` | Validate profile operator and preparation/resource behavior | Jv, Jᵀw, normal action, fixed-state parity, derivative and action preservation | Source fingerprints, watchdog/RSS, campaign time and archived receipts | Named synthetic/archived dataset cases | Paired baseline/candidate behavior; permanent operator tests provide independent algebraic checks | `joint_operator_validation_runner_test`; `ProfileOperator_test.cpp` | Fingerprint calculation and shared parity are in neutral support; the runner keeps its expected baseline policy | `ProfileOperator_test.cpp`, `Numerics_test.cpp`, shared provenance helper | MIGRATION_REQUIRED | Preserve any missing work-count and campaign provenance assertions in permanent tests | Operator receipts archive present and manifest hash matches | Hard-pins baseline commit `15e71e3…`; no other runner imports this runner |
| JH-05 | Search/Schwarz campaign: `joint_search_validation.py` | Compare bounded search/preconditioner candidates | Fixed-state scientific parity, search status, accepted state and forbidden extra work | Campaign watchdog, resource limits, source and input fingerprints | Named chain/single/heterogeneous/large cases | Paired runtime variants; C++ search tests provide algorithmic owner | `joint_search_validation_runner_test`; `Search_test.cpp`, `OperatorSearch_test.cpp` | Finite/state/scientific parity and build freshness are in neutral support; the runner retains its campaign acceptance policy | `Search_test.cpp`, `OperatorSearch_test.cpp`, permanent resource benchmark | MIGRATION_REQUIRED | Preserve search-promotion and forbidden-work assertions in permanent tests before retiring | Operator/search manifests and logs present | Hard-pins baseline commit `c6c869c…`; no other runner imports this runner |
| JH-06 | Fixed-action campaign: `joint_fixed_validation.py` | Measure frozen-state normal/composed actions and work | Normal-action parity and no unexpected factor/reference work | Fixed-state hashes, historical A/B/C timing and RSS | Frozen states generated from named campaign inputs | Paired action comparison; permanent `ProfileOperator_test` checks algebraic parity | `joint_fixed_validation_runner_test`; `ProfileOperator_test.cpp` | Shared fixed-action parity is in neutral support; fixed campaign retains A/B/C policy and exact artifacts | `ProfileOperator_test.cpp`, `Numerics_test.cpp`, consolidated fixed-work benchmark | MIGRATION_REQUIRED | Preserve action/work invariants in permanent tests; keep factor-reuse checks permanent and historical exact counts non-contractual | Fixed-actions receipt archive and JSON manifest present | Hard-pins baseline commit `ec24f30…`; no other runner imports this runner |
| JH-07 | Bounded search/rank campaign: `joint_bounded_validation.py` | Compare bounded rank prototype, dense rank oracle, and search behavior | Rank interval/status, threshold, inputs and fixed-step parity | Bounded runtime, watchdog, source/input/state hashes and resource ceilings | Synthetic topologies and named cases; dense oracle is run only on bounded sizes | Dense rank-oracle action for eligible sizes; no dense oracle for largest cases | `joint_bounded_validation_runner_test`; `FreeDesignRank_test.cpp`, `Search_test.cpp` | Rank comparison uses neutral primitives; bounded oracle and large-case policies remain campaign-specific | `FreeDesignRank_test.cpp`, `Search_test.cpp`, independent dense reference tests | MIGRATION_REQUIRED | Preserve rank-oracle status and eligibility semantics in permanent tests before retiring | Bounded/search/fixed artifacts and receipts present where manifested | Records `base_commit=5f61bb6e`; no other runner imports operator, search, or fixed runners |

### PR 3 build-policy update

This table records the build and CTest owner after PR 3. “Default” means
`BUILD_TESTING=ON` with the new benchmark and research options OFF. The earlier
target inventory remains a historical snapshot; optionalizing a tool does not
change its retirement status.

| ID / entry point | Previous default status | New category and build owner | Default built / registered? | Enable with |
|---|---|---|---|---|
| S2-01 Fold-168 | Runner contract test registered; full external regression gated | CORE runner contract; EXTERNAL full regression | Contract test yes; campaign no | `RHBM_GEM_ENABLE_FOLD_168_REGRESSION=ON` plus the existing model/map paths |
| S2-02 MDPDE / forward experiment | `mdpde_experiment` was created in every testing build | RESEARCH executable; permanent MDPDE C++ tests remain CORE | Executable no; C++ tests yes | `RHBM_GEM_BUILD_RESEARCH_TOOLS=ON` |
| S2-03 Endpoint refinement | Permanent C++ tests and runner contract registered; helper sources were compiled into `rhbm_gem` in testing builds | CORE policy tests and runner contract; test helpers now live in non-installed `rhbm_gem_test` | Yes | `BUILD_TESTING=ON` |
| S2-04 Failed-only refinement | Manual Python campaign; permanent C++ policy tests | CORE policy owners; RESEARCH/HISTORICAL manual runner | C++ tests yes; campaign is not registered | Run the retained Python runner with its existing inputs |
| S2-05 Second-stage audit parser | Parser CTest registered by default | CORE parser contract | Yes | `BUILD_TESTING=ON` |
| S2-06 Audit neutrality | C++ probes and helpers were compiled into `rhbm_gem` in testing builds; paired Python campaign was manual | CORE C++ probe tests in `rhbm_gem_test`; paired campaign remains manual | C++ tests yes; campaign no | `BUILD_TESTING=ON` |
| JC-01 Joint runtime | Runtime executable and regressions were in the default build | CORE runtime executable and regressions | Yes | `BUILD_TESTING=ON` |
| JC-02 Joint offline audit | Existing offline option | OFFLINE audit executable, support, and CTests | No | `RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON` |
| JC-03 Joint workflow CLI | CLI smoke registered by default | CORE CLI workflow contract | Yes | `BUILD_TESTING=ON` |
| JC-04 Partial-selection measurement | Measurement executable was a `tests_all` dependency | BENCHMARK measurement; permanent selection tests remain CORE | Executable no; tests yes | `RHBM_GEM_BUILD_BENCHMARKS=ON` |
| JC-05 Postprocessing benchmark | Benchmark executable was a `tests_all` dependency | BENCHMARK executable | No | `RHBM_GEM_BUILD_BENCHMARKS=ON` |
| JC-06 Public API benchmark | Executable was created in every testing build | BENCHMARK executable | No | `RHBM_GEM_BUILD_BENCHMARKS=ON` |
| JC-07 Sparse validation driver | Benchmark driver was a `tests_all` dependency; validation-runner CTests were default | BENCHMARK + RESEARCH + HISTORICAL driver; runner CTests are category-gated | No | `RHBM_GEM_BUILD_BENCHMARKS=ON` or `RHBM_GEM_BUILD_RESEARCH_TOOLS=ON` |
| Stage A/B/C `joint_validation` | Available through the offline-audit option | Shared OFFLINE + BENCHMARK + RESEARCH campaign executable; remains a manual campaign | No | Offline audits, benchmarks, or research tools |

No row's KEEP / MIGRATION_REQUIRED / READY_FOR_RETIREMENT value changes because
of build optionality. The permanent protection matrix below retains its owners;
the independent offline oracle stays behind its existing option. No runner,
benchmark source, fixture, artifact, or numerical behavior was removed.

### Current candidate summary

- `READY_FOR_RETIREMENT`: endpoint refinement (`S2-03`), failed-only full-capture orchestration (`S2-04`), partial-selection measurement (`JC-04`), and postprocessing timing harness (`JC-05`). Their remaining use is historical campaign replay or bounded timing; permanent correctness/workflow owners and retained evidence are listed in the rows. Keep their fixtures and permanent tests.
- `MIGRATION_REQUIRED`: audit ON/OFF full-workflow comparison (`S2-06`), instrumented sparse benchmark dependencies (`JC-07`), complete-command resource measurement (`JH-C`), and campaign-specific acceptance protections for reference/compact/sparse/operator/search/fixed/bounded (`JH-01` through `JH-07`).
- `KEEP`: external fold regression, MDPDE research/replay support, second-stage audit parser, Joint runtime and independent audit, CLI workflow smoke, public API benchmark, and Joint Stages A/B. Stage A and B stay separate even though one Python runner dispatches all three stages.

The PR 2 dependency disentanglement changes the owner of reusable I/O, hashing, provenance, watchdog, comparison, numerical replay, Fold-168, and MDPDE support. It removes runner-to-runner imports but does not migrate the campaign-specific scientific acceptance matrix into permanent tests. Therefore no `MIGRATION_REQUIRED` classification is upgraded here.

### Historical campaign pinning and runner dependency graph

The validation runners use hard-coded historical pins as follows: reference pins baseline `f1ac45c…` and candidate `62d178d…`; compact pins source fingerprint `97185310…`; operator pins commit `15e71e3…`; search pins `c6c869c…`; fixed-actions pins `ec24f30…`; bounded search writes `5f61bb6e` as its base-commit label and also requires matching source fingerprints. Sparse validation has no literal historical commit constant in its source, but records build/input provenance. These are campaign reproduction anchors, not current production-version requirements.

Arrows below mean “the left-hand Python module imports the right-hand module.” This is why deleting one script can break another campaign even if no CTest directly invokes it.

```text
endpoint_refinement ─┬─> mdpde_experiment ─> fold_168_regression
failed_only_refinement ─────────────────────> fold_168_regression
                         └─────────────────> mdpde_experiment

joint_reference_validation -> joint_compact_validation -> joint_sparse_validation
joint_operator_validation --------------------------------> joint_sparse_validation
joint_search_validation -> joint_operator_validation
joint_fixed_validation -> joint_operator_validation, joint_search_validation
joint_bounded_validation -> joint_operator_validation, joint_search_validation,
                            joint_fixed_validation

joint_component_runtime -> joint_runtime_support, joint_fixture_records,
                           joint_offline_support
joint_component_audit -> joint_runtime_support, joint_component_runtime,
                         joint_fixture_records, joint_offline_support
```

The reference → compact → sparse and bounded → fixed → search → operator chains are real Python import dependencies. `joint_runtime_support.py`, `joint_validation.py`, `joint_fixture_records.py`, and `joint_offline_support.py` are shared support modules; do not classify them as disposable just because they are not top-level campaigns.

## Permanent protection matrix

“Duplicate owner” means a second automated owner exists; it does not mean either owner is redundant. A property owner must remain after any historical runner is removed.

| Property | Current owner(s) | Duplicate owner? | Future owner if a historical runner leaves | Owner remains? |
|---|---|---:|---|---:|
| A/C solution correctness | `Numerics_test.cpp`, `Components_test.cpp`, `joint_component_offline_tests` | Yes | Same permanent tests and independent reference | Yes |
| Active-face correctness | `ProfileOperator_test.cpp`, `Components_test.cpp` | Yes | `ProfileOperator_test.cpp` | Yes |
| KKT / gradient | `Components_test.cpp`, offline precision/derivative tests | Yes | Permanent component and derivative tests | Yes |
| Rank determination | `FreeDesignRank_test.cpp`, `Components_test.cpp`, offline rank checks | Yes | `FreeDesignRank_test.cpp` plus dense small-case oracle | Yes |
| Weak spectrum | `ObservableProfile_test.cpp`, offline precision tests | Yes | Permanent spectrum/rank tests | Yes |
| SVD/reference parity | `Components_test.cpp`, `OfflineComponents_test.cpp`, compact/reference runners | Yes | Permanent dense/multiprecision reference tests | Yes |
| Derivative correctness | `joint_offline_tests`, `Precision_test.cpp`, operator tests | Yes | Offline oracle tests | Yes |
| Profile Jacobian | `ProfileOperator_test.cpp`, `ObservableProfile_test.cpp` | Yes | `ProfileOperator_test.cpp` | Yes |
| Jv | `ProfileOperator_test.cpp` | No | `ProfileOperator_test.cpp` | Yes |
| Jᵀw | `ProfileOperator_test.cpp` | No | `ProfileOperator_test.cpp` | Yes |
| Normal action | `ProfileOperator_test.cpp`, fixed-action runner | Yes | `ProfileOperator_test.cpp` | Yes |
| Factor reuse | `Numerics_test.cpp`, `ProfileOperator_test.cpp` instrumentation cases | Yes | Permanent operator tests | Yes |
| Bounded search status | `Search_test.cpp`, `FreeDesignRank_test.cpp`, bounded runner tests | Yes | Permanent search/rank tests | Yes |
| Selected voxels define observation domain | `PartialSelection_test.cpp` | No | `PartialSelection_test.cpp` | Yes |
| Contributor closure | `PartialSelection_test.cpp` | No | `PartialSelection_test.cpp` | Yes |
| Non-recursive halo inclusion | `PartialSelection_test.cpp` | No | `PartialSelection_test.cpp` | Yes |
| Target-vs-halo role distinction | `PartialSelection_test.cpp`, `ObservableProfile_test.cpp` | Yes | Those permanent tests | Yes |
| Singleton halo | `PartialSelection_test.cpp`, `ObservableProfile_test.cpp` | Yes | Those permanent tests | Yes |
| Bridge contributor | `PartialSelection_test.cpp`, `Components_test.cpp` | Yes | Those permanent tests | Yes |
| Unobserved target | `ObservableProfile_test.cpp` | No | `ObservableProfile_test.cpp` | Yes |
| Zero-amplitude target | `PartialSelection_test.cpp`, `ObservableProfile_test.cpp` | Yes | Those permanent tests | Yes |
| Target identifiability | `ObservableProfile_test.cpp`, offline target evidence tests | Yes | Permanent observable-profile tests | Yes |
| Halo non-uniqueness | `ObservableProfile_test.cpp` | No | `ObservableProfile_test.cpp` | Yes |
| Target convergence vs full-parameter convergence | `ObservableProfile_test.cpp`, `Components_test.cpp` | Yes | Permanent observable-profile tests | Yes |
| `FittingStage::Second` writeback | `PartialSelection_test.cpp`, `joint_workflow_cli_smoke_test` | Yes | Those permanent workflow tests | Yes; one C++ case failed in this run and is reported below |
| Persistence | `DataObjectPersistence_test.cpp`, workflow CLI smoke | Yes | Permanent data/workflow tests | Yes |
| JSON schema | `DataObjectSchemaValidation_test.cpp`, workflow CLI smoke | Yes | Permanent schema and workflow tests | Yes |
| SQLite | `DataObjectFileIO_test.cpp`, workflow CLI smoke | Yes | Permanent data/workflow tests | Yes |
| CSV | Workflow CLI smoke and observable-profile persistence tests | Yes | Permanent workflow tests | Yes |
| Save/reload | `DataObjectPersistence_test.cpp`, workflow CLI smoke | Yes | Permanent workflow tests | Yes |
| Export | Workflow CLI smoke, observable-profile persistence test | Yes | Permanent workflow tests | Yes |
| Target-only summary | `PartialSelection_test.cpp` | No | `PartialSelection_test.cpp` | Yes |
| Peeling coverage | `PartialSelection_test.cpp` | No | `PartialSelection_test.cpp` | Yes |
| Covariance availability | `ObservableProfile_test.cpp`, `PartialSelection_test.cpp` | Yes | Permanent uncertainty tests | Yes |
| Native-success unchanged | `ProductionFitting_test.cpp`, `MDPDEExperiment_test.cpp` (`FailedOnlyLeavesNativeSuccessUntouchedDespiteFreshResidual`) | Yes | Those permanent tests | Yes |
| Failed-only refinement policy | `MDPDEExperiment_test.cpp` (`PoliciesDoNotConfuseNativeSuccessWithFreshQualification`, `FailedOnlyLeavesNativeSuccessUntouchedDespiteFreshResidual`, `ProductionRefinesCapturedFailuresWithoutRewritingNativeStatus`, `FailedOnlyDoesNotPromoteAnIneligibleRankDeficientEndpoint`) | Yes | Those permanent policy tests | Yes |
| Wrong-root rejection | `MDPDEExperiment_test.cpp` (`KnownOtherRootFailsBranchComparison`) | No | Same permanent branch test | Yes |
| Fresh-weight behavior | `MDPDEExperiment_test.cpp` (`AcceptedResultHasFreshWeightsAndExistingCovarianceFormula`) | Yes | Same permanent fresh-weight test | Yes |
| Covariance consistency after refinement | `MDPDEExperiment_test.cpp` (`AcceptedResultHasFreshWeightsAndExistingCovarianceFormula`), `ProductionFitting_test.cpp` | Yes | Same permanent tests | Yes |
| Budget exhaustion | `MDPDEExperiment_test.cpp` (`BudgetRejectionPreservesNativeHistoryButIsUnqualified`), `ProductionFitting_test.cpp` recovery exhaustion tests | Yes | Same permanent tests | Yes |
| Invalid variance/state | `MDPDEExperiment_test.cpp` (`InvalidEndpointsAreNeverPromoted`, `NegativeAndNonfiniteVarianceAreNeverPromoted`) | Yes | Same permanent endpoint tests | Yes |
| Rank-deficient refinement | `MDPDEExperiment_test.cpp` (`InvalidEndpointsAreNeverPromoted`, `FailedOnlyDoesNotPromoteAnIneligibleRankDeficientEndpoint`) | Yes | Same permanent endpoint tests | Yes |
| Weak-halo diagnosis and returned-state identity | Offline `joint_offline_diagnostic_smoke`; historical `weak-halo.json` equivalence | No | Offline diagnostic, separate from numerical pass/fail | Yes |
| Stage B RNG stream and generated input | `joint_statistical_experiment_contract_test` fixed-seed hashes | No | Same research contract test | Yes |
| Stage B attempted/completed/qualified/failed/unavailable denominators | `joint_statistical_experiment_contract_test` | No | Same research contract test | Yes |
| Stage C command save/reload/export and resource status | `joint_benchmark_smoke`, `joint_benchmark_contract_test` | Yes | Unified command profile and workflow correctness tests | Yes |
| Audit ON/OFF numerical neutrality | `NumericalAuditNeutrality_test.cpp`; manual two-build runner | Partial | Permanent full-workflow paired test | Partial; migrate first |
| EIGEN backend | Joint C++ tests and runtime CTest (this snapshot) | Yes | Backend-specific permanent CTest | Yes |
| SPQR backend | Backend-conditional C++ tests and archived SPQR CTest receipts | Yes historically | Run SPQR permanent CTest in an SPQR build | Not run in this baseline |
| `BUILD_TESTING=OFF` | CMake option guards and historical testing-disabled build receipt | Yes historically | Keep a small configure/build smoke in release validation | Not run in this baseline |
| Installed consumer | `tests/cmake/consumer_smoke`, `lint_install_smoke`, historical receipts | Yes | Installed consumer smoke | Yes; `lint_all` installed, configured, built, and ran the consumer |
| Offline-audit configuration | `joint_offline_tests`, offline runner CTests | Yes | Same conditional CTest entries | Yes; enabled and run here |
| Extended tests | Conditional `joint_component_extended_regression` | No in this build | Conditional extended CTest | Not run; option OFF |
| External regression | Conditional fold-168 CTest plus runner contract test | Partial | External fold-168 CTest with inputs | Not run; inputs unavailable |

## Artifact state

The artifact manifest records path, bytes, current-tree existence, SHA-256, and whether an existing manifest hash was checked. PR 4 rechecked three paths that the PR 1 manifest had recorded as present: `joint-sparse-acceptance/measurements.json.gz`, `joint-sparse-acceptance/command-exports.tar.gz`, and `joint-operator-search/receipts.tar.gz`. All three are absent from this checkout; their expected sizes and hashes remain in the original verification manifests but could not be rechecked.

The older certification report names `scientific-records.tar.gz` and `frozen-diagnostics/frozen-endpoint-records.tar.gz` without giving a containing directory. Neither name appears in the current figures tree or the checked historical tree. All five archive references are recorded as missing; no replacement was generated. `joint-abc-components/search-records-a.tar.gz` is absent at HEAD but recoverable from commit `ba2449f`; its historical size and SHA-256 match the retired-artifact inventory.

The full fold-168 CIF/map/sidecar and some large campaign inputs are external data dependencies, not missing archives. They were not supplied to this run. Other campaign result files in `docs/developer/figures/` were not re-generated or hash-audited; an unlisted file is **not verified**, not implicitly present or missing.

See [`artifact-manifest.json`](figures/experiment-retirement-baseline/artifact-manifest.json) for the checked archives and explicit missing paths.

## Build and test inventory

The captured configuration is a fresh Debug build using SYSTEM dependencies, `BUILD_TESTING=ON`, EIGEN backend, offline audits ON, extended tests OFF, external fold regression OFF, UMAP OFF, Python bindings OFF, and experimental features OFF. Environment: AppleClang 21.0.0, CMake 4.4.3, Ninja 1.13.2, Boost 1.92.0, GoogleTest 1.18.0, SQLite 3.54.0; OpenMP 5.1 was found and ROOT support was enabled by system detection.

`ctest-list.txt` is the actual generated CTest registration list; `build-targets.txt` is the generated named target list for this configuration. CTest count is a snapshot descriptor only. Future retirement must preserve property coverage, not an identical test count.

The generated list contains 14 C++ CTest entries (13 `rhbm_tests` groups plus `joint_offline_tests`), 16 Python entries, and one shell smoke test. The named project build inventory contains two C++ test executables and eight experiment/audit/benchmark executables for this feature configuration.

The initial full CTest run produced one existing `rhbm_tests_joint_component` failure and one runner-contract error because the sandbox denied the runner's `ps` process-table query. The offline C++ audit executable passed. The complete summary and any permitted focused rerun are recorded in `baseline.json`; neither failure was changed or marked expected.

The external fold-168 regression was not configured because CIF/map inputs are not present. Joint extended tests, SPQR backend, Python bindings, and a separate `BUILD_TESTING=OFF` project build were not exercised by this configuration. `lint_all` did exercise the installed-consumer smoke against the configured build.

Repository lint and install-consumer smoke both passed. The permanent CTest documentation-layout group was rerun after adding this document and passed.

## Future PR preconditions

### Historical Joint campaign retirement

- Reusable parity, fingerprint, I/O, and watchdog mechanisms have neutral support owners with focused tests outside campaign modules.
- Watchdog/process execution and incomplete/killed/not-run status semantics do not depend on a validation runner.
- Sparse/compact/reference have no runner-level import dependency.
- Operator/search/fixed/bounded have no runner-level import dependency.
- Each fixture has an explicit permanent owner; shared fixtures are not removed with a runner.
- Historical artifact paths and provenance hashes are fixed and retrievable as documented.

### Endpoint and failed-only refinement retirement

- Native-success unchanged, failed-only selection, wrong-root rejection, fresh weights, covariance consistency, budget exhaustion, invalid variance, and rank-deficient behavior each have permanent tests.
- Historical exact call counts are no longer used as correctness contracts.

### Benchmark consolidation

- Preparation, solve, postprocessing, save, and complete-command timing boundaries are fixed.
- RSS, timeout, killed, incomplete, and not-run semantics are explicit.
- Incomplete, not-run, timed-out, or killed measurements cannot be counted as pass.
- Timing-only harnesses are identified as performance tools, never correctness tests.

## Scope result

This baseline changes no production numerical behavior. It removes no tool, target, test, fixture, or historical artifact. Future simplification should begin only with rows marked `READY_FOR_RETIREMENT`, and still preserve the permanent owner and artifact conditions stated above.

## PR 2 dependency disentanglement update

The original inventory snapshot above remains the PR 1 record from commit `875f7aeaa3177ba0bde8227ddec27372a497639c`. PR 2 started from a clean `develop` worktree at `4cf9daa183996efa7fdd1be7e49f99c7622cb51e`. The changes are currently uncommitted, so the ending `HEAD` is still `4cf9daa183996efa7fdd1be7e49f99c7622cb51e`.

The pre-change AST graph contained these runner-to-runner edges:

```text
endpoint_refinement -> fold_168_regression, mdpde_experiment
failed_only_refinement -> fold_168_regression, mdpde_experiment
mdpde_experiment -> fold_168_regression
simulation_contract -> fold_168_regression
joint_component_audit -> joint_component_runtime
joint_sparse_validation -> joint_validation
joint_compact_validation -> joint_sparse_validation, joint_validation
joint_reference_validation -> joint_compact_validation, joint_validation
joint_operator_validation -> joint_sparse_validation, joint_validation
joint_search_validation -> joint_operator_validation, joint_validation
joint_fixed_validation -> joint_operator_validation, joint_search_validation, joint_validation
joint_bounded_validation -> joint_fixed_validation, joint_operator_validation, joint_search_validation, joint_validation
joint_validation_profile -> joint_validation
joint_validation_report -> joint_validation
```

Those edges are now cut. The affected runners import neutral `experiment_*`, `joint_*_support`, `joint_validation_checks`, `joint_fixture_records`, or `joint_numerical_reference` modules. AST guards check runner imports, literal dynamic imports/launches, support-to-runner imports, and integration-test-to-test imports. Tests may still import a runner as the subject of that runner's contract test.

No `MIGRATION_REQUIRED` row was promoted: the extraction gives shared mechanisms neutral owners, but it does not move every campaign-specific scientific acceptance check into permanent tests. No historical runner, existing CTest registration, fixture, artifact, C++ support file, production estimator, or numerical threshold was removed or changed. Two focused Python CTest entries were added.

The PR 2 Debug/EIGEN build completed. Full CTest ran 33 entries: 32 passed, including the offline regression; the sole failure was the same `JointComponentPartialSelectionTest.StageAdapterUsesIdentityAndClearsMissingStates` failure already recorded in the PR 1 baseline. All 12 focused Python runner/support CTests passed. `joint_validation_runner_test` now tests timeout and RSS status classification with a mocked process-tree measurement, so it no longer needs sandbox permission to inspect the host process table.

## PR 4 Joint campaign retirement update

PR 4 began on `develop` at `e5a9e7922f063743df71b6e2561ee3c43c781ff5` (2026-09-29), after a fast-forward pull from `origin/develop`; the starting worktree was clean. The PR 1 inventory, PR 2 neutral support extraction, and PR 3 benchmark/research build separation were present. The original classification and inventory above remain historical snapshots. This section records the final PR 4 outcome without changing those earlier classifications.

| Retired runner | Historical purpose | Permanent property owner | Replacement measurement profile | Historical-only logic | Outcome |
|---|---|---|---|---|---|
| `joint_sparse_validation.py` | EIGEN/SPQR solution, status, rank, derivative, and resource comparison | `Numerics_test.cpp`, `FreeDesignRank_test.cpp`, backend-specific Joint CTests against dense/reference cases, and `experiment_support_test.py` parity/status primitives | `fixed`, `solve`, `rank` | Paired campaign acceptance, commit/build fingerprints, and exact historical resource comparisons | Retired |
| `joint_compact_validation.py` | Compact-SVD spectrum/rank/solution parity and SVD work measurements | `Numerics_test.cpp`, `FreeDesignRank_test.cpp`, `experiment_support_test.py` | `fixed`, `rank` | Legacy/compact campaign gates, source fingerprints, and timing thresholds | Retired |
| `joint_reference_validation.py` | Frozen baseline versus incremental endpoint/reference acceptance | `joint_numerical_reference.py`, `Precision_test.cpp`, `OfflineComponents_test.cpp`, offline CTests, and retained Stage A/B/C runner contracts | `solve`, `workflow`, `command` | Commit-pinned replay gates, campaign exports, and source fingerprints | Retired |
| `joint_operator_validation.py` | Profile-operator preparation, action, and fixed-state measurements | `ProfileOperator_test.cpp`, `Numerics_test.cpp`, `experiment_support_test.py` | `prepare`, `fixed` | PR0 paired resource campaign and archived receipts | Retired |
| `joint_search_validation.py` | Legacy/operator/Schwarz search comparison and promotion checks | `Search_test.cpp`, `OperatorSearch_test.cpp`, `experiment_support_test.py` | `solve` | Baseline promotion gates, campaign schedule, and large-local resource sweep | Retired |
| `joint_fixed_validation.py` | Frozen-state composed/normal actions, factor reuse, and work counts | `ProfileOperator_test.cpp`, `Numerics_test.cpp`, `experiment_support_test.py` | `fixed` | Historical A/B/C policy, three-run comparisons, and campaign receipts | Retired |
| `joint_bounded_validation.py` | Bounded rank/search outcomes, watchdog, RSS, and not-run semantics | `FreeDesignRank_test.cpp`, `Search_test.cpp`, `OperatorSearch_test.cpp`, `experiment_support_test.py`, and `experiment_process.py` | `rank`, `solve` | PR5 schedule, case eligibility campaign, and fixed resource thresholds | Retired |

The separate Python `joint_postprocessing_benchmark.py` wrapper was also retired because the retained C++ `joint_postprocessing_benchmark` is now called by the unified `workflow` and `postprocess` profiles. The seven candidate runners, their four runner-only self-tests, `joint_compact_support.py`, and `joint_validation_test_data.py` are removed. `joint_validation.py` and its Stage A/B/C workflow remain in place.

Permanent numerical protection includes sparse solve/dense-reference parity; compact spectrum, rank, and threshold boundaries; `Jv`, `Jᵀw`, adjoint and normal-action parity; factor reuse and forbidden-work counters; search endpoint/status; bounded rank and unavailable semantics; and independent reference/precision/derivative certification. The original dense and offline reference implementations remain independent. Shared comparison checks now live in `experiment_support_test.py`, and timeout/RSS/incomplete-process semantics remain in `experiment_process.py` and its focused contract tests.

Current measurements use [`joint-benchmark.md`](joint-benchmark.md), with explicit `prepare`, `fixed`, `solve`, `rank`, `workflow`, `postprocess`, and `command` profiles. The C++ drivers remain separate: `joint_sparse_benchmark`, `joint_postprocessing_benchmark`, and `joint_component_benchmark`. The first two serve the unified profiles; no benchmark algorithm was moved into Python. `BUILD_TESTING=ON` without benchmark/research options no longer registers or depends on the retired campaign entries.

Historical source remains retrievable from Git history. The seven validation source files were last changed at `7f84f931215f7b5a52948f7ecdbfc15fad79d7c5`; the postprocessing wrapper was last changed at `47514581f13e6544d2800cac59ff445948d9043c`. Historical timing gates and source pins remain evidence only. At PR 4 completion, seven hash-checked archives were present and had no active-tool dependency. Three additional archive paths previously recorded as present were absent from that checkout; two older references whose paths were never located also remained recorded as missing. PR 6 records the later removal of the seven present archives below; no archive was recreated.

PR 4 changes no production source, estimator default, numerical tolerance, rank threshold, convergence policy, oracle, or Stage A/B/C methodology. Build and smoke outcomes for this checkout are reported with the PR 4 change review; the PR 1–3 test records above are unchanged snapshots.

## PR 5 refinement and Joint responsibility update

PR 5 began on `develop` at `465168565d5bf3c046c12e471e22cd3098c59bad`
(`2026-09-29T21:47:43+08:00`) after `git pull --ff-only`; the worktree was
clean and `origin/develop` matched. PR 1–4 records above remain historical
snapshots. The current owners below supersede their endpoint/failed-only and
Stage A/B/C rows.

### Current owner mapping

| Retired historical entry | Permanent/current owner | Result and retained evidence |
|---|---|---|
| `endpoint_refinement.py`, `endpoint_refinement_test.py` | `MDPDEExperiment_test.cpp`, `ProductionFitting_test.cpp`; shared refinement test support | Campaign/report orchestration retired; fixtures retained; report and figures removed from the worktree with provenance in the artifact manifest |
| `failed_only_refinement.py` | `MDPDEExperiment_test.cpp`, `ProductionFitting_test.cpp` | Full fold comparison retired; production failed-only policy remains unchanged; figures removed from the worktree with provenance in the artifact manifest |
| Stage A weak-halo campaign | `joint_offline_diagnostic.py` and `joint_offline_diagnostic` | Numerical/identifiability diagnosis retained as OFFLINE; old `weak-halo.json` and `weak-snapshot.json` retained |
| Stage B noise/mismatch campaign | `joint_statistical_experiment.py` and `joint_statistical_experiment` | KEEP as RESEARCH with the same noise, seeds, case order, mismatch, roles and metrics; historical outputs retained |
| Stage C command/resource campaign | `joint_benchmark.py --profile command` | Old watchdog/profile orchestration retired; CLI smoke and unified benchmark own command/save/reload/export and resource status; old results retained |

The corresponding protection matrix rows now name permanent refinement tests,
Stage B deterministic RNG/denominator contracts, the offline diagnostic smoke,
and the unified command benchmark. Exact historical call totals and work-count
fingerprints remain campaign evidence, not correctness assertions.

### Joint validation function inventory

This classifies the top-level functions as they existed before the PR 5 split.
The mixed modules are removed after their active responsibilities moved to the
single-purpose owners above.

| Source function | Responsibility | Current owner / disposition |
|---|---|---|
| `joint_validation.py:run_a` | Stage A process loop and receipts | `joint_offline_diagnostic.py`; shared budget ledger and watchdog retired |
| `clean_map`, `noise_field`, `conditions`, `metrics`, `run_b` | Stage B scientific inputs, per-fit metrics and execution | `joint_statistical_experiment.py`; methods and iteration order preserved |
| `run_c` | Stage C generated cases, census, process runs and repetitions | Retired; command profile in `joint_benchmark.py` owns current complete-command measurement |
| `main` | Stage dispatch, shared campaign budget/provenance ledger | Removed with mixed runner |
| `joint_validation_report.py:wilson`, `parameter_stats`, `statistical_summary` | Stage B reporting | `joint_statistical_experiment.py`, with explicit denominator fields |
| `compact_a` | Stage A diagnostic compaction and state checks | `joint_offline_diagnostic.py`; retained-state and snapshot checks remain |
| report `main` | Merged A/B/C reports and artifact manifest | Removed; Stage A and B write their own output; Stage C uses benchmark schema |
| `joint_validation_profile.py:main` | Stage C stack sampler and campaign accounting | Retired; no current owner |
| `joint_validation_support.py:resource_run` | Stage C analysis/export watchdog and resource report | Retired; shared process support and command benchmark remain |
| C++ `Weak` | Stage A numerical diagnosis | `joint_offline_diagnostic.cpp` |
| C++ `TwoAtoms`, `Statistical` | Stage B deterministic generated-input check and numerical fits | `joint_statistical_experiment.cpp` |
| C++ `Inspect`, `Generate` | Stage C historical resource case preparation | Retired with Stage C campaign matrix |
| C++ `Read`, `Write`, `Outcome`, `Census` | Small serialization/census helpers | Kept only in the owner that needs them; no all-stage dispatcher |

The new integration architecture guard checks that Stage A, Stage B and the
benchmark runner do not import or launch one another. `experiment_process.py`
remains shared process infrastructure; the retired Stage C profile/watchdog
framework is not duplicated.

### Permanent refinement properties

| Property | Permanent owner | Covered after PR 5 |
|---|---|---:|
| Native success leaves parameters, status, covariance and iterations unchanged | `FailedOnlyLeavesNativeSuccessUntouchedDespiteFreshResidual`; production test | Yes |
| Native failure enters failed-only policy; rank-ineligible endpoint is rejected before root continuation | `ProductionRefinesCapturedFailuresWithoutRewritingNativeStatus`; `FailedOnlyDoesNotPromoteAnIneligibleRankDeficientEndpoint` | Yes |
| Wrong mathematical root fails branch comparison | `KnownOtherRootFailsBranchComparison` | Yes |
| Refined state uses fresh weights | `AcceptedResultHasFreshWeightsAndExistingCovarianceFormula` | Yes |
| Covariance follows the existing formula and rejected candidates preserve it | `AcceptedResultHasFreshWeightsAndExistingCovarianceFormula`; `ProductionFitting_test` | Yes |
| Budget exhaustion stays unqualified and preserves the native endpoint | `BudgetRejectionPreservesNativeHistoryButIsUnqualified`; `ProductionFitting_test` | Yes |
| Invalid/negative/NaN/Inf variance cannot be promoted | `InvalidEndpointsAreNeverPromoted`; `NegativeAndNonfiniteVarianceAreNeverPromoted` | Yes |
| Rank deficiency and invalid denominator stay explicit | `InvalidEndpointsAreNeverPromoted`; `FailedOnlyDoesNotPromoteAnIneligibleRankDeficientEndpoint` | Yes |

The rank-ineligible property reflects the current API exactly: native `SUCCESS`
returns before the failed-only entry point; other native statuses enter it, and
the refinement routine rejects unsuitable data before attempting root/reference
continuation. No new status allowlist or production policy change was introduced.

### Historical artifact ownership

| Artifact group | Current meaning | Active input? |
|---|---|---:|
| Endpoint campaign outputs (removed; retrieval provenance in the artifact manifest) | Endpoint campaign historical evidence | No |
| Failed-only campaign output (removed; retrieval provenance in the artifact manifest) | Failed-only fold comparison historical evidence | No |
| `figures/joint-validation/weak-snapshot.json`, `weak-halo.json` | Stage A diagnosis evidence | No |
| `figures/joint-validation/noise-runs.json`, `noise-summary.json`, `noise-table.md` | Stage B research results | No; full-run output is newly generated |
| `figures/joint-validation/resources.json`, `resource-preflight.json`, `diagnostic-profiles.json`, persistence records and logs | Stage C historical benchmark evidence | No |
| `figures/joint-validation/campaign.json`, `manifest.json`, `source-index.json`, build metadata | Shared historical provenance | No |

No historical figure, archive, result, fixture or source history was physically
removed. The machine-readable [artifact manifest](figures/experiment-retirement-baseline/artifact-manifest.json)
records the PR 5 owner mapping and no active dependency.

### Files removed and build categories

Removed active campaign files are `endpoint_refinement.py`,
`endpoint_refinement_test.py`, `failed_only_refinement.py`,
`joint_validation.py`, `joint_validation_test.py`,
`joint_validation_profile.py`, `joint_validation_report.py`,
`joint_validation_support.py`, and the multiplexed C++
`tests/experiments/joint_validation.cpp`. Their owners are listed above. No
fixture, independent MDPDE reference, production source or historical report was
removed.

| Responsibility | Build category | Registered smoke/check |
|---|---|---|
| Permanent MDPDE/refinement correctness | CORE | `rhbm_tests_utils_hrl`, `rhbm_tests_core_contract` |
| Stage A numerical diagnosis | OFFLINE | `joint_offline_diagnostic_smoke`, plus offline numerical tests |
| Stage B research experiment | RESEARCH | `joint_statistical_experiment_contract_test`, `joint_statistical_experiment_smoke` |
| Stage C complete command measurement | BENCHMARK | `joint_benchmark_smoke` and `joint_benchmark_contract_test` |

No per-stage CMake option was added. SPQR and external Fold/6Z6U cases are
reported as not run or unavailable when their backend or inputs are absent.

### PR 5 equivalence and verification record

| Check | Before / after result |
|---|---|
| Stage A weak-halo start 0 | Rebuilt the pre-PR5 `joint_validation.cpp` runner from the starting revision and ran the same fixture/start as the new offline diagnostic. Start B and returned state match; design/projected-width/normalized-width/Jacobian ranks are `4/2/2/2`; target runtime convergence is `passed`, full runtime convergence is `failed`; the weak direction is the halo-width direction (up to eigenvector sign). Independent 50/100-digit values agree in both runs; derivative and correction audits both remain `failed`. Numerical diagnosis fields match; elapsed-time/RSS and expected source/build fingerprints differ. |
| Stage B fixed-seed case | Historical and smoke records for iid seed `20260921` match exactly: input SHA-256 `508655679930ee0d406fde4ea1ecb68b2eaff5d2133bc7071ef176638d79857a`, observation SHA-256 `c5147ebda1fe72d1012b49d91838743ce4e8ca00abf4b3f8d3dcfa233fba02e1`, completed/qualified status, parameters, errors, prediction RMSE `0.0003535161322428431`, and residual RMSE `0.0031937884385234948`. |
| Stage C command smoke | Existing PR 4 `command` profile is the unchanged owner. Its small command smoke completed analysis and reload/export, wrote SQLite and JSON/CSV exports, and passed the separate process-status contract (`completed`, `process_error`, `timeout`, `rss_limit`, `not_run`). Timeout/RSS-limit classifications are contract-tested, not induced in this smoke. |
| Core and full CTest | Full run: 32 tests, 31 passed after rerunning the two process-inspection smokes with permission; one existing failure remains in `JointComponentPartialSelectionTest.StageAdapterUsesIdentityAndClearsMissingStates` (summary-column formatting assertions). The MDPDE/refinement group passed again after adding explicit native-success status and iteration-preservation assertions. |
| Optional environment | Debug build uses `EIGEN`; SPQR was not run. External Fold and 6Z6U inputs were unavailable in this checkout. |

The archived Stage A JSON carries a different source fingerprint from the PR 5
starting checkout, so it was retained as historical evidence but not used as
the before-run for the equivalence check above. Historical artifacts remained
physically present at the PR 5 starting checkout.

## PR 6 historical artifact and documentation update

PR 6 began on `develop` at `a8cd8d1465705b1d7d5aa28d68b47018ff0c0f71`
(`2026-09-29T22:53:05+08:00`) after `git pull --ff-only`; the starting worktree
was clean and `origin/develop` matched. PR 1–5 sections above remain historical
snapshots. This section records the current ownership and worktree state.

### Active dependency scan and artifact decisions

The source scan covered `src/`, `include/`, `tests/`, `resources/`, `cmake/`,
`CMakeLists.txt`, and developer Markdown. It checked the candidate filename,
relative path and distinctive stem. No production, permanent test, current
benchmark/research tool, or current documentation command reads the seven files.
The direct archive links and retired-runner reproduction instructions have been
removed from current documentation. The exact pre-removal hashes and compressed
file sizes were rechecked against the starting worktree before removal.

| Artifact | Classification | Active callers | Decision | Reason |
|---|---|---:|---|---|
| `joint-reference-acceptance/historical-inputs.tar.gz` | HISTORICAL_EVIDENCE | 0 | REMOVE_FROM_WORKTREE | Retired campaign input bundle; provenance and retained summary/verification are in the manifest/Git history. |
| `joint-reference-acceptance/historical-audits.tar.gz` | HISTORICAL_EVIDENCE | 0 | REMOVE_FROM_WORKTREE | Retired audit matrices; exact compressed hash and source commit are recorded. |
| `joint-reference-acceptance/historical-results.tar.gz` | DUPLICATE_HISTORICAL_EVIDENCE | 0 | REMOVE_FROM_WORKTREE | Summary, verification, and compact decisions remain; raw campaign outputs are in Git history. |
| `joint-reference-acceptance/incremental-inputs.tar.gz` | HISTORICAL_EVIDENCE | 0 | REMOVE_FROM_WORKTREE | Retired experiment input bundle; not used by current tests or benchmark. |
| `joint-reference-acceptance/incremental-audits.tar.gz` | HISTORICAL_EVIDENCE | 0 | REMOVE_FROM_WORKTREE | Retired fixed-state audits; compact conclusions and provenance remain. |
| `joint-reference-acceptance/incremental-results.tar.gz` | DUPLICATE_HISTORICAL_EVIDENCE | 0 | REMOVE_FROM_WORKTREE | Historical command results are summarized and retrievable from Git history. |
| `joint-fixed-actions/receipts.tar.xz` | OBSOLETE_RAW_LOG | 0 | REMOVE_FROM_WORKTREE | Raw per-process receipts/logs; summary, comparison, and experiment metadata remain. |

### Retained artifact decisions

| Artifact group | Classification | Active owner/caller | Decision | Reason |
|---|---|---|---|---|
| `tests/fixtures/joint_component/` | ACTIVE_FIXTURE / ACTIVE_REFERENCE | Joint permanent C++ and Python tests | KEEP_ACTIVE | Deterministic regression and oracle inputs; hashes/sizes remain in the fixture inventory. |
| `tests/fixtures/mdpde/` | ACTIVE_FIXTURE / ACTIVE_REFERENCE | MDPDE and production-fitting permanent tests | KEEP_ACTIVE | Includes captured offset-IRLS and shape maximum-iteration cases used by permanent protection. |
| `tests/benchmarks/fold_168_simulation_baseline.json` | ACTIVE_BENCHMARK_INPUT | Optional Fold-168 regression | KEEP_ACTIVE | A compact reference baseline, not an output archive. |
| `figures/joint-validation/noise-runs.json` and summaries | ACTIVE_RESEARCH_OUTPUT | `joint_statistical_experiment.py` produces the current result schema; no runtime reader | KEEP_HISTORICAL | Preserves the full fixed-seed Stage B outcomes and denominator alongside the current research method. |
| `figures/joint-validation/weak-halo.json`, `weak-snapshot.json` | HISTORICAL_EVIDENCE | Offline diagnostic owner; no historical-result reader | KEEP_HISTORICAL | Records the weak-halo conclusion and original input; its old source fingerprint differs from current code. |
| `figures/joint-validation/resources.json`, preflight and diagnostic records | HISTORICAL_EVIDENCE | Joint benchmark owns current resource measurement; no historical-result reader | KEEP_HISTORICAL | Retains historical timeout/resource outcomes and their limitations. |
| `figures/joint-reference-acceptance/summary.json`, `verification.json`, `testing-disabled-result.json` | HISTORICAL_EVIDENCE | No active reader | KEEP_HISTORICAL | Compact acceptance summary plus the distinct testing-disabled full-command outcome. |
| `figures/joint-fixed-actions/summary.json`, `comparison.json`, instrumentation patch | HISTORICAL_EVIDENCE | Permanent operator tests and unified benchmark own current behavior | KEEP_HISTORICAL | Preserves the decision record, detailed comparisons and measurement-only patch. |
| `figures/joint-sparse-acceptance/late-profile/stack-*.txt` | HISTORICAL_EVIDENCE | No active reader | KEEP_HISTORICAL | Small timeout stack samples have distinct diagnosis value and remain hash-linked by their receipt. |

The 10.3 MB `heterogeneous-168.tar.gz` Joint fixture remains active in
`tests/fixtures/`; the 5 MiB figure threshold does not apply to fixture inputs.
No permanent fixture, benchmark input, Stage B output, or independent oracle
was removed.

The machine-readable [artifact manifest](figures/experiment-retirement-baseline/artifact-manifest.json)
is the sole source for each removed path, status, compressed bytes, SHA-256,
source commit, last-present commit, purpose, and replacement. `removed_from_worktree`
means the current branch no longer tracks that file; Git history is unchanged.
The earlier missing paths (`scientific-records.tar.gz`,
`frozen-diagnostics/frozen-endpoint-records.tar.gz`, sparse measurements and
exports, and operator receipts) remain recorded as missing; none was recreated.

### Retained evidence and ownership

- Joint component and MDPDE fixtures, the independent numerical references, and
  the fixed-seed Stage B research outcomes remain in the tree. The Stage B tool
  can still run from its fixed-seed design without the removed archives.
- Compact Joint reference summaries, verification records, fixed-action
  summary/comparison, and small decision evidence remain historical evidence.
- The remaining `docs/developer/figures/` tree contains curated results,
  diagnostics, baseline records, and manifests; generated full benchmark and
  research outputs belong under `build/` or an explicit work directory.

### Current and historical documentation

| Document | PR 6 role | Status |
|---|---|---|
| `joint-benchmark.md` | Sole current Joint benchmark instructions, including command profile | Current |
| `joint-noise-mismatch-validation.md` | Stage B research design, metrics, denominator, smoke and full run | Current research |
| `joint-fixed-actions.md` | Current fixed-action benchmark entry point and numerical contract | Current |
| `joint-component-reference-acceptance.md` | Retained acceptance conclusions and Git-history retrieval | Historical |
| `joint-fixed-actions-validation.md` | Retained campaign results and limitations | Historical |
| `figures/joint-reference-acceptance/README.md` | Small evidence index and removed-archive provenance | Historical |

Developer guidance now assigns permanent inputs to `tests/fixtures/`, generated
outputs to build/work directories, and `figures/` to curated compact evidence.
The old archive checksums were removed from campaign-local metadata; the central
artifact manifest owns current path/status/hash/provenance records.

### Guards added or retained

| Guard | Owner | Coverage |
|---|---|---|
| Current Markdown local-link validation | `DocumentationLayoutTest.MarkdownLinksResolveWithinRepository` | Current links must resolve; historical removed files are no longer linked. |
| Manifest and fixture validation | `historical_artifact_manifest_test` | Unique artifact IDs/paths, allowed statuses, removed-file absence, missing-file records, source/hash metadata, and fixture sizes/hashes. |
| Active dependency and large artifact ownership | `historical_artifact_manifest_test` | Removed archive names/paths/stems cannot appear in current code or developer commands; figure files at least 5 MiB must be listed in the manifest. |

The 5 MiB ownership threshold is above all retained figure files and below each
of the seven removed archives; it makes future large files require an owner
without banning large active fixtures elsewhere in the repository.

### PR 6 behavior boundary

No production source, estimator behavior, numerical tolerance, rank or
convergence threshold, permanent regression fixture, or independent oracle was
changed or removed. Stage B noise generation, seeds, mismatch definition and
metrics are unchanged. Joint benchmark and Stage B research inputs remain
available without historical archives. Git history was not rewritten; the
external archive service remains unused.

### Verification

| Check | Result |
|---|---|
| Production-only Release configure/build (`BUILD_TESTING=OFF`) | Passed from a fresh build directory. |
| Default Release CTest | 25/26 passed. The sole failure is the PR 5 baseline issue `JointComponentPartialSelectionTest.StageAdapterUsesIdentityAndClearsMissingStates`; the other core, integration, documentation-link and manifest checks passed. |
| Artifact manifest/ownership test | Passed: seven unique removed entries, retained fixture hashes/sizes, no active references, and no unowned figure at or above 5 MiB. |
| Offline diagnostic smoke | Passed with a fresh work directory: 3/3 starts completed; weak-halo rank/precision/restart output was produced. |
| Unified Joint benchmark smoke | Passed for all seven profiles, including command persistence/export. Timeout and RSS classifications remain covered by the benchmark contract test. |
| Stage B statistical smoke | Passed. Seed `20260921` input hash `508655679930ee0d406fde4ea1ecb68b2eaff5d2133bc7071ef176638d79857a`; observation hash `c5147ebda1fe72d1012b49d91838743ce4e8ca00abf4b3f8d3dcfa233fba02e1`; 3/3 summarized attempts completed and qualified. |
| Optional backend/external cases | Four SPQR-specific tests were skipped in this EIGEN build. Fold/6Z6U external inputs were unavailable and not run. |

The full tracked-byte accounting is stored under `pr6_update.size_accounting_*`
in the artifact manifest. It distinguishes current-tree removal from Git
history size; the seven removed archives total 252,805,334 compressed bytes.

### Current-tree size

| Measure | Before PR 6 | After PR 6 |
|---|---:|---:|
| Tracked files | 679 | 673 |
| Current tracked tree bytes | 277,957,447 | 25,171,196 |
| Tracked `figures/` files | 115 | 108 |
| Tracked `figures/` bytes | 258,695,460 | 5,898,770 |

The seven removed compressed files account for 252,805,334 bytes (241.094 MiB).
The current tracked tree is smaller by 252,786,251 bytes (241.076 MiB), after
including PR 6's manifest, guard, and documentation additions. The recorded
active fixture/reference inventory remains 20 files totaling 14,181,744 bytes.
These are current tracked-tree figures; Git history and clone history size are
unchanged.
