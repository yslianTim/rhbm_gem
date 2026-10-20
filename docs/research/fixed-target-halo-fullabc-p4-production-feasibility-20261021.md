# Fixed-target Observation-driven Halo FullABC — P4 Production gate

Date: 2026-10-21
Branch: `develop`
Decision: **do not promote**. The Production Halo reconstruction remains `geometry-local-target-median-v1`.

## Decision basis

The research estimator is implemented behind an explicit optional postprocessing result. Callers that do not provide it retain the existing geometry-local-target-median-v1 path. The P2 change does not alter the Pure Exclusion Target Joint solver, its endpoint contract, or Target A/B/C; no SQLite schema or persisted Halo snapshot semantics changed.

The P3 physical synthetic validation establishes a useful signal/Peeling result, not a Production-wide reliability claim. Ten of ten first Target Joint endpoints qualified and downstream Target A/B/C drift was exactly zero. Eight of ten FullABC endpoints qualified; the exact-three and exact-four observation free Halo stress cases failed A/C feasibility, KKT/replay, endpoint and numerical qualification. Among the eight valid comparisons, Halo-only signal RMSE and Peeling RMSE improved in all eight. Local/Group estimates improved inconsistently, and FullABC has a separated high-noise A-parameter outlier (A RMSE 1.413, maximum absolute A error 3.150). The P3 tables and raw cases are in [the physical validation report](fixed-target-halo-fullabc-p3-physical-validation-20261021.md).

The complete requested 5,938,087-row / 37,406-contributor Halo FullABC workflow has not been measured. The previously recorded 37,406-target synthetic FixedNeighbor memory test is a zero-sweep preparation and zero-state endpoint preflight: it measured 1,471,479,808 bytes peak RSS, but it did not execute FirstStage, free-Halo FullABC Active-set LM, or downstream work. See [the existing large-component resource report](joint-bounded-giant-component-resource-validation-20261020.md). That result cannot certify this Halo estimator's 4 GB ceiling.

## Promotion gates

| Gate | Evidence | Status |
| --- | --- | --- |
| First-round Target A/B/C unchanged | The second estimator accepts only the qualified Pure Exclusion snapshot; P3 independently measured zero postprocessing drift in all paired cases. | Pass for implemented path |
| Sole preclassification rule is distinct `n_h < 3` | P1 deterministic tests and real P0 census; n=3/4 synthetic cases entered the free set. | Pass |
| Reliable FullABC on normal geometries | Eight valid physical cases qualified; both exact n=3/4 free-Halo stress cases failed numerical qualification. Rank-deficiency recovery is not implemented. | **Insufficient** |
| Halo signal and Peeling improve | Both improved in all eight paired synthetic cases; the result is not a production-data comparison. | Partial |
| Downstream Local/Group improve | Local A/B and Group posterior RMSE improved in only four of eight paired cases; Group uncertainty changes are small. | **Not established** |
| Failures have a production-safe handling policy | Unqualified free-Halo results are unavailable and do not silently fall back. A recovery strategy for n≥3 non-identifiable/ill-conditioned fits remains separate work. | **Incomplete** |
| Reconstruction-only contributors retained | P2 tests and P3 map-edge n=0 fixed Halo cover this behavior. | Pass for in-memory research path |
| Chemical-label invariance | P1/P2 permutation and relabeling tests passed. | Pass for tested paths |
| Full-scale 4 GB limit | P3 RSS is 50,151,424 B on a small synthetic map; no large Halo FullABC process RSS or complete workflow run exists. | **Not measured** |
| Deterministic results | A repeated 4 Å experiment returned identical Halo profiles and accuracy values; P1/P2 also test deterministic ordering. | Partial |
| Persistence, UMAP, Two Stage regression | No new persisted Halo contract was promoted, so those pipelines were not exercised against FullABC snapshots. Existing regression suite is run separately. | **Not established for new contract** |
| Numerical acceptance criteria and complete regression | P3 preserved all failure records; all 31 CTest test names passed across the full run and corrected cleanup-test rerun. Production-scale qualification remains absent. | Pass for repository regression; insufficient for promotion |

## Large-component resource scope

| Stage | Evidence available | Gate result |
| --- | --- | --- |
| Real map/model construction and Halo count at 5.94M rows | Not available for the requested large component; the actual Halo count must be counted from the selected model. | Not measured |
| Target prediction subtraction and fixed-median Halo prediction | Correctness tested on small OriginalInput domains; not profiled at 5.94M rows. | Not measured at scale |
| Sparse free-Halo support and FixedNeighbor workspace | P3 uses 211–300 free support entries and a maximum estimated workspace of 4,163,594 B; no large Halo component. | Not measured at scale |
| Constrained A/C and Active-set LM | Passed on eight small physical cases; unqualified on n=3/4 free-Halo stress cases. | Not certified at scale |
| Endpoint replay, Halo snapshot, Peeling, Local, Group | Synthetic P3 run covered these with 4,000/4,000 Peeling samples and 20/20 Local and Group members for qualified scenarios. | Small-case only |
| SQLite save/reload, UMAP and Two Stage | Not run with a FullABC Halo snapshot; the Production path and schema remain unchanged. | Not measured |
| Process RSS | P3 small-case peak 50,151,424 B; previous Target-only synthetic preflight peak 1,471,479,808 B. | Neither is a full-scale Halo FullABC RSS result |

P3 estimates memory and records process RSS separately. It retains all OriginalInput rows, including residual-constant zero-sensitivity rows, and uses sparse Halo support incidence. This is a construction-level design benefit; it does not remove the need to run the complete large workflow under a process-level measurement.

## Production action

No promotion commit is made. The standard Peeling reconstruction contract, Target-only Pure Exclusion Joint fitting, Target qualification and downstream behavior remain the existing Production behavior. `EstimateHaloFullABCJoint` remains an internal research estimator with explicit unavailability on solver failure. There is no FullABC SQLite schema version and no implicit interpretation of existing kNN-4 snapshots as FullABC results.

The next promotion attempt requires, at minimum: a complete 5.94M-row workload with the actual Halo count and stage-by-stage RSS; successful qualification or an explicit, tested handling policy for free-Halo failures; real-data paired Peeling/Local/Group outcomes using the same frozen Target endpoint; persistence/reload and UMAP/Two Stage verification if introducing a saved contract; and full regression evidence. This implementation does not meet those gates, so the current kNN-4 Production default is retained.

## Validation

The complete build passed with `cmake --build build/debug --parallel 4`. The first full `ctest --test-dir build/debug --output-on-failure -j1` run completed 30/31 tests successfully in 2,566.94 seconds. Its only failure was `joint_experiment_cleanup_test`, which correctly rejected the two new experiment sources because they were absent from its whitelist. The whitelist was updated, and `ctest --test-dir build/debug -R '^joint_experiment_cleanup_test$' --output-on-failure -j1` then passed. Thus all 31 tests passed across the full run plus the targeted rerun; the full 31-test suite was not repeated after this test-only whitelist correction.

The long numerical/scalability items passed: `rhbm_tests_joint_component_numerical` in 748.64 seconds, `joint_component_gram_scalability` in 1,276.40 seconds, `joint_component_mapping_memory_scalability` in 6.15 seconds, `joint_component_regression` in 151.11 seconds, and `joint_component_physical_smoke` in 73.32 seconds. All remaining tests except the corrected whitelist check passed in the full run. The P3 runner was built with `cmake --build build/debug --target joint_halo_fullabc_validation --parallel 4`, and all ten physical scenarios wrote their JSON/CSV results. The P2 targeted tests and component contract validation are documented in [the P2 report](fixed-target-halo-fullabc-p2-implementation-20261020.md). No Notion documents were changed.
