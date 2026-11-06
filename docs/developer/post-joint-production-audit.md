# Post-Joint production consolidation audit

## Scope and baseline

The starting local and fetched `origin/develop` HEAD is
`98b101cd6670451e1f4f7f37abbd51ef4eb96d56`. The working tree was clean.
This cleanup keeps Original-grid collection at 2 Å, non-self Peeling support
at 2.5 Å, direct own-charge correction, and Local fitting at 1 Å. Joint,
FirstStage, Two-Stage, Gaussian MDPDE and Group WEB algorithms are outside the
change boundary. Stage JSON 5, SQLite 21 and Joint JSON 24 retain their shapes.

`PostJointLocalFittingRadius` lives beside the persisted Local result type,
so core fitting and data validation share one contract value. Strict equality
is intentional: 1.0 is exactly representable in binary64 and survives current
JSON/SQLite round-trips exactly. All recorded Local statuses require it.

The previous delivery report was available in the ignored
`build/release/grid-domain-validation/` directory. Both existing build caches
had `BUILD_TESTING=OFF`; baseline builds explicitly enable testing with the
existing Release/Debug presets, Apple Clang, Eigen, ROOT and OpenMP. No
dependency is added. Builds use four jobs; CTest is serial with
`OMP_NUM_THREADS=2`, `OMP_THREAD_LIMIT=4`, and single-thread BLAS/Accelerate.
Explicit per-test OpenMP settings remain bounded by that limit.

Baseline evidence and logs are retained locally under ignored
`build/release/production-consolidation/`. Release passes all 21 CTest entries
in 40.38 s; Debug passes all 21 entries in 458.22 s, including the unchanged
numerical/replay group in 240.71 s. Registration covers all 1,071 enabled
GoogleTests with zero missing registrations.

Three independent Release evidence runs compare all saved Joint and Stage
values, including grid identities/responses, Local alpha/OLS/MDPDE,
information, qualification, fixed charge and Group posterior/covariance.
Comparison preserves binary64 values exactly and excludes only generated run
UUIDs and elapsed-time fields; source/run relationships remain checked by the
normal regression tests. All three baseline captures agree. A separate
four-eligible-member Group fixture supplies a deterministic WEB reference
without random Group alpha splits. The varied 20-Target workflow also agrees
across the three baseline captures and exercises JSON/SQLite round-trips.
Temporary capture code is removed after collecting evidence.

Representative baseline Target 1 has 1,237 candidate/corrected/positive rows,
171 retained rows, alpha 0, OLS `(A,tau)` =
`(1.9999999999922033, 0.6499999999986469)`, MDPDE `(A,tau)` =
`(1.9999999999921967, 0.6499999999986452)`, and information condition
`64.36153474356654`. Joint assembled objective is
`2.5239616261232942e-23`. Across the workflow, all 24,740 Peeling observations
are available. Median postprocessing is 0.036063750 s; estimated sparse
workspace is 2,372,456 B and peak payload is 4,748,808 B. Process RSS measured
with full evidence serialization includes capture overhead; clean workload
measurements are reported separately in the delivery evidence.

## Repository-wide reference audit

Searches covered the entire tracked source, public headers, tests, CMake and
documentation. Definitions, callers and negative tests were read before
classification. “Production” below includes required validation, persistence
and unavailable semantics, not only the numerical solver.

| Symbol/file | Current usage | Production dependency | Test dependency | Decision / proposed change | Risk and validation |
| --- | --- | --- | --- | --- | --- |
| `FitPostJointLocalTarget` | Workflow and test adapter in `PotentialFittingWorkflow` | Local dataset, fit, rank and information | Range, nonpositive, rank and width-sensitivity cases | Refactor: expose the existing internal production entry and remove the runtime radius argument | Exact Local parity and dataset reference |
| `FitPostJointLocalTargetForTesting` | Header/wrapper; Pure Exclusion and partial-selection tests | None beyond forwarding | Direct Local qualification cases | Remove wrapper; call production function | Retain all caller coverage |
| `kPostJointFittingDistanceMax` | Core fitting range | 1 Å Local policy | Dataset consistency | Share one value with the persisted result type and data validation | Exact 1.0 equality; invalid-radius tests |
| `PostJointCandidateRadius` | Original support collector in `PostFitPeeling` | Inclusive 2 Å candidates | Grid boundary/identity tests | Keep, separate from fitting policy | Original rows/counts unchanged |
| `PreparedLocalGaussianDesign` and dataset functions | FirstStage, Two-Stage, Post-Joint | Shared Gaussian preparation/estimation | Solver and empty-dataset regression | Keep generic range support; remove only Post-Joint experimental selection | FirstStage/Two-Stage sentinel and numerical tests |
| `PostJointLocalResult::fitting_radius` | Stage JSON; SQLite stores Stage payload | Persisted estimation evidence | Editor, assessment, persistence and Group tests | Keep field; require exact production radius in all statuses | Reject 1.5/2/zero/nonfinite values without migration |
| `AssessPostJointLocalResult` | Workflow, editor, Stage validation and Group fitting | Source/count/charge/diagnostic/information qualification | Partial-selection and persistence negative tests | Tighten range only | Available/Unavailable/Ineligible and Group eligibility |
| `StageResultJson` / `PersistenceSchema` | JSON and repository I/O | Current formats and strict rejection | Round-trip, malformed/legacy formats | Keep versions and field encoding | JSON 5 / SQLite 21 / Joint 24 regression |
| `AtomLocalPotentialView` / `ModelAnalysisEditor` | UI response views, atomic updates, invalidation | Independent grid lifetime and source validation | Map/geometry/run/snapshot invalidation | Keep shared validation paths | Workflow and persistence tests |
| Five research CSV/JSON files | Only the grid validation report links them | None; no CMake/CTest/runtime readers | None | Remove after retaining representative conclusions | Repository reference scan and documentation guard |
| `post-joint-grid-validation.md` | Stage integration links the report | Decision context | Documentation links | Replace phase history with production decision record | Verify ranges, equations, results and current versions |
| Five mapping/identity metrics counters | Increments in context/workflow; fixed-count assertions in `PartialSelection_test` | No logger/UI/persistence/benchmark consumer | Assertions only | Remove instrumentation and assertions, keep actual checks | Malformed mapping/source/geometry regressions remain |
| Preparation call/row metrics (`target_preparation_count`, `total_prepared_sample_count`, `valid_local_sample_count`, `invalid_local_sample_count`) | Increments and a conservation assertion in partial-selection tests | No logger/UI/persistence consumer; Local results already retain actual dataset counts | Shared-versus-standalone output comparison | Remove redundant counters and their counting loop | Keep per-sample identity/response/reason comparisons and Local result counts |
| `estimated_reconstruction_workspace_bytes` | Sparse Original row payload estimate, debug logger, workflow test property | Memory diagnosis | Workflow memory evidence | Rename to `estimated_peeling_workspace_bytes`; formula unchanged | Exact payload comparison |
| Remaining `JointPostprocessingMetrics` fields | Availability, phase timing, RSS, Halo/tail diagnostics, payload | Active diagnostic responsibilities | Workflow numerical and resource checks | Keep | Paired phase/RSS measurements |
| `JointPeelingContributor` | Sorted ID/index lookup; captured role, eligibility and finite geometry validation | Contributor completeness and provenance | Malformed universe, role, index, foreign and nonfinite cases | Keep: reducing captured invariants would require a different validation design | Preserve existing negative tests; negligible payload compared with sparse rows |
| `PostJointLocalSample` | Grid identity, distance, correction/failure state; Gaussian sample conversion | Solver dataset and failure identity | Independent peeling/charge and unavailable references | Keep: cached intermediates avoid recomputation and retain reference coverage | Exact selected rows and direct Gaussian/erf checks |
| `AssignGridObservationsForTesting` | Abstract editor/Group/UI fixtures adapt a raw sample buffer | No runtime use | Model analysis, UMAP and partial-selection fixtures | Replace with direct grid-observation construction | Keep real-map collectors and Original/global identities in integration tests |
| Post-Joint interpolation/node cache/raw pairing | No remaining runtime branch found | None | Old-format rejection is still meaningful | No numerical cleanup needed; retain malformed/legacy rejection | Final symbol scan |
| `MapInterpolation.hpp` and sampler tests | FirstStage interpolation remains active | FirstStage | Closed-form tricubic weight/clamp references | Keep | Sampler CTest coverage |
| `RHBM_GEM_TEST_INSTRUMENTATION` | FirstStage hooks, editor invalidation, MDPDE and Joint diagnostics | Separate test-library build | Multiple unrelated suites | Keep project-wide switch | Build and registration guard |
| `tests/CMakeLists.txt` | Group filters, strict registration coverage, Debug numerical timeout | Test execution | All 21 CTest groups | Keep registrations and existing 900 s Debug numerical timeout | Zero uncovered enabled GoogleTests |

No branch or worktree is created, no Notion document is accessed or changed,
and no temporary benchmark executable is added to the repository. Existing
unrelated worktrees are left intact.
