# Joint component historical evidence

This is the canonical index of retired Joint Component experiments and their
historical results. It records what was tested, what the evidence established,
and which questions remain open. Current implementation and numerical behavior
belong in the current-owner documents below; stable mathematical definitions
belong in the two formal contracts. Historical benchmark results are not current
performance measurements, and historical success does not establish
cross-platform behavior or real-data generalization.

No historical Markdown report or evidence artifact is deleted by PR-H1 or PR-H2.
Any later deletion is a separate PR-H3 decision.

## 1. Scope and ownership

| Concern | Current owner |
| --- | --- |
| Runtime estimator and selection domain | [Runtime guide](joint-component-runtime.md) |
| Sparse backend and bounded SPQR rank | [Sparse backend](joint-component-sparse-backend.md) |
| Compact SVD | [Compact SVD](joint-component-compact-svd.md) |
| Fixed factors and normal actions | [Fixed actions](joint-fixed-actions.md) |
| Tiled/profile operator | [Profile operator](joint-profile-operator.md) |
| Operator search | [Operator search](joint-operator-search.md) |
| Observable halo | [Observable-halo parameterization](joint-observable-halo.md) |
| Target estimability | [Target estimability](joint-target-estimability.md) |
| Offline weak-halo diagnosis | [Weak-halo guide](joint-weak-halo-attribution.md) |
| Statistical research protocol | [Noise/mismatch experiment](joint-noise-mismatch-validation.md) |
| Current benchmark commands | [Benchmark guide](joint-benchmark.md) |
| Workflow integration | [Stage integration](joint-stage-integration.md) |
| Data and persistence | [Data I/O architecture](architecture/dataobject-io-architecture.md) and [Joint JSON/CSV contract](commands/potential_analysis.md#provenance-and-map-units-joint-json-schemas-3-and-4) |
| Certification mathematics | [Certification contract](joint_abc_certification_contract.md) |
| Component mathematics | [Component contract](joint_abc_components_contract.md) |

## 2. Algorithm-formation evidence

These findings describe the historical experiments that shaped the current
contracts. The four retired experiment pages are short retrieval pointers; this
section is the summary of their scientific conclusions.

| Historical question | Evidence and limitation | Contract retained today | Provenance |
| --- | --- | --- | --- |
| Do sampled observations match the frozen support and sampler? | On the truth-neighbor control, compact prediction versus sampler RMSE was **5.17001e-16**; float32 residual was **2.91318e-8**. This checks observation consistency, not convergence with unknown neighbors. | Preserve the sampler, strict support boundary, contributor membership, and post-sum quantization controls. | [Profile experiment](joint-abc-profile-experiment.md); [coverage experiment](joint-abc-coverage-experiment.md); [fixture catalog](../../tests/fixtures/joint_component/catalog.json). |
| Does an estimated-neighbor sweep improve every atom? | Across eight checkpoints, aggregate A/B/C errors improved by about **19–21 / 25–26 / 4.7–4.9 times**, while C worsened for **70–76 of 168 atoms**. Aggregate improvement is not per-atom recovery. | Keep per-atom evidence and frozen support provenance; do not infer individual recovery from aggregate error. | [Coverage experiment](joint-abc-coverage-experiment.md); historical source and artifact hashes are indexed by the retired-artifact inventory below. |
| Does matched joint A/C solve qualify from full design rank? | All 16 joint runs exhausted 2,000 evaluations without qualification despite design rank 336; frozen-neighbor controls qualified 1,335/1,344 blocks. Rank and lower C error do not establish stationarity. | Search completion, a usable state, and numerical qualification are distinct. | [Component experiment](joint-abc-components-experiment.md); [runtime contract](joint-component-runtime.md). |
| What did Experiments A/B/C establish about composite MDPDE? | A: baseline-best28 at alpha 0 had A/C RMSE **0.011471954 / 0.00063732864**; alpha 0.1 qualified at **0.0091138914 / 0.00049314988**, while alpha 0.5 and 1 exhausted their budgets. B: alpha 0 qualified at **0.0091707311 / 0.00032313699**; positive-alpha runs exhausted 100 evaluations. C: two starts both exhausted 100 evaluations; one endpoint had A/C RMSE **0.017916446 / 0.00011271024** and effective sample fraction **0.8983%**. The C comparison had no independent repeat or common-alpha control, so it neither isolates the cause nor disproves heterogeneous MDPDE. | These runs are not prerequisites for current joint LS and do not justify changing its objective. | [Coverage experiment](joint-abc-coverage-experiment.md); [retired artifact inventory](figures/joint-component-evidence/retired-artifacts.json). |
| Does fixed-B variable projection support the joint profile? | The historical profile work retained independent double/float32 recovery and full profile-Jacobian checks at nonzero residual. It is finite-fixture evidence, not a general convergence result. | Keep the constrained linear subproblem, full nonzero-residual correction, and independent reference checks in the current numerical contract. | [Profile experiment](joint-abc-profile-experiment.md); [certification contract](joint_abc_certification_contract.md). |
| Did the coverage set make weak and degenerate cases successful? | Nine immutable datasets remain in the fixture catalog. Weak, active-face, zero-signal, and duplicate cases retain their original resolution, correction, width-identification, and state-availability limitations. | Keep these cases as distinct regression outcomes; do not turn unavailable or failed evidence into passes. | [Fixture catalog](../../tests/fixtures/joint_component/catalog.json); [fixture guide](../../tests/fixtures/joint_component/README.md); [v1 acceptance](joint-component-v1-acceptance.md). |
| Did guarded certification qualify all cases? | The frozen baseline recorded 128 paired branches, 64 required global regular pairs, and 192 component-local scopes (128 regular); 96 regular globally restricted scopes use a separate certification scope. These counts apply only to the recorded fixtures and policy. | Keep global, restricted, and local certification scopes separate; local certification cannot replace global qualification. | [Certification experiment](joint-abc-certification-experiment.md); [certification contract](joint_abc_certification_contract.md); [component contract](joint_abc_components_contract.md). |
| Did exact component decomposition preserve parent evidence? | Transition receipts compare 1,248 migrated search/state/context records, 2,096 global/restricted audit records, and 3,712 baseline records including separate local scopes. Four representative local certificates were rerun; the migrated offline driver did not rerun the full audit matrix. | Preserve structural partition, shared parent scale/rank context, same-state assembly, and failure isolation. | [Transition receipt](figures/joint-component-evidence/transition-acceptance.json); [global pairs](figures/joint-component-evidence/global-pairs.csv); [scope matrices](figures/joint-component-evidence/global-and-restricted-scopes.csv) and [local scopes](figures/joint-component-evidence/local-scopes.csv). |

Public objectives are normalized by the parent observation scale squared;
archived experiment objectives use raw half-RSS. The transition receipt is
historical comparison evidence, not a new runtime measurement. Its 22 default
and 3 extended/offline test counts came from the recorded macOS/compiler
configuration and do not establish cross-platform behavior.

## 3. Backend and numerical evolution

| Milestone | Historical evidence | Retained limitation | Current owner |
| --- | --- | --- | --- |
| Tiled derivative and assessment | For heterogeneous-168 / first-stage-float32 on the recorded macOS arm64 host, peak RSS changed from **1,252.14 MiB to 359.14 MiB (71.3% lower)**; total time changed from 28.7429 s to 27.6702 s (3.7% lower), with assessment slightly slower. Fixed-state parity passed at the same states. | This establishes a storage reduction for that host and case, not general scalability or a universal speedup; sparse factor storage and compact quadratic factors remain. | [Profile operator](joint-profile-operator.md); [tiled evidence](joint-component-tiled-backend.md); [measurement JSON](joint-component-tiled-backend.json). |
| Sparse EIGEN/SPQR backend | On completed fixed states, backend comparisons retained status, rank, active face, coefficients and existing Guarded checks. The single-128 complete command measured 19.20 s versus 11.78 s (1.63x); all six runs converged. Heterogeneous-168 and single-512 fixed-state timings are preserved in the [acceptance report](joint-component-sparse-acceptance.md). | Both 6Z6U attempts stopped during reference at 600 s after primary rank 4314/4334; neither had an accepted state. Both single-512 complete-command attempts stopped at 600 s without persisted output. Missing archives are listed in §9. | [Sparse backend](joint-component-sparse-backend.md). |
| Compact SVD | The historical compact campaign passed its recorded numerical and fixed-state performance gates. The single-512 candidate completed three full commands (median 289.238 s); the ce58c897 command baseline stopped at 600 s, so no completed-command speedup is established. | The acceptance is tied to its source, backend, host, and cases. It is not a current performance guarantee. | [Compact SVD](joint-component-compact-svd.md); [compact acceptance](joint-component-compact-svd-acceptance.md). |
| Fixed factors and normal actions | The recorded campaign completed **72/72 processes** and **24/24** three-repetition A/B, B/C, and A/C comparisons. Normal actions used two Q/Q′ calls versus six for the composed control. | These fixed-state/fixed-step results do not promote a different full-search default or establish large-workflow scalability. Some original step subtotals omitted gradient setup; corrected bounds are not exact end-to-end step times. | [Fixed actions](joint-fixed-actions.md); [validation report](joint-fixed-actions-validation.md); [summary](figures/joint-fixed-actions/summary.json) and [comparison decisions](figures/joint-fixed-actions/comparison.json). |
| Endpoint reference | The 43.81-minute historical campaign completed all 18 fixed-state audits and all 128/512 analysis-export comparisons. Single-128 changed from 5.229 s to 4.912 s; single-512 from 307.805 s to 300.564 s, with peak RSS from 2.500 GiB to 1.872 GiB for single-512. Search reference evaluations changed from six to zero. | These are incremental historical results on the recorded host; search and endpoint assessment remained the dominant costs. | [Compact SVD](joint-component-compact-svd.md); [endpoint-reference report](joint-component-reference-acceptance.md); [archive inventory](figures/joint-reference-acceptance/archives.json). |
| Profile operator | All 18 fixed-state audits and 20 preparation-only workloads completed. On 6Z6U, baseline and candidate reached the sampled RSS watchdog on both backends and had no numerical verdict. The 10,000-atom runs stopped after input/layout/basis preparation; they did not solve, rank, search, assess, or estimate uncertainty. | Operator construction and action parity do not prove full search or full-workflow scaling. | [Profile operator](joint-profile-operator.md); [acceptance record](figures/joint-profile-operator/acceptance.md). |
| Operator LM and Schwarz search | The baseline campaign ended at its 80-minute deadline with 41 completed numerical runs, one deadline termination, and 16 entries not run. Single-512 Schwarz medians were 419.562 s vs legacy 229.232 s (Eigen) and 350.104 s vs 217.933 s (SPQR), exceeding the 1.10 performance ratio on both backends. | Required comparisons remained incomplete; 6Z6U and large-local campaign cases were not run. Production remains LegacyCompact; no promotion was made. | [Operator search](joint-operator-search.md); [validation report](joint-operator-search-validation.md); [manifest](figures/joint-operator-search/manifest.json). |
| Bounded SPQR rank prototype | Small controls established both full-rank and deficient decisions against the dense oracle. The representative chain-128 and cube-128 SPQR prototypes exhausted 100 million charged work units and returned Unavailable before completing the reconstruction certificate. | Unavailable is not deficient. No 128-atom definitive rank result, 512/2,000-atom rank result, or scaling guarantee was established. The campaign was not repeated after its final deadline-guard correction. | [Sparse backend](joint-component-sparse-backend.md); [bounded-rank report](joint-bounded-search-rank.md); [validation](joint-bounded-search-rank-validation.md). |

The sparse and compact-SVD acceptance reports preserve exact completed-case
tables, fingerprints, and retained test logs. Their fixed-state evidence must not
be confused with the full-command resource envelope in §6.

## 4. Structural and identifiability boundaries

### 4.1 Historical full-ABC 6Z6U obstruction

At historical base 5f06b027, the selected domain had 2,192 contributors and
262,801 rows. Its largest component had 2,167 atoms and 4,334 free A/C columns.
Twenty singleton halo atoms each contributed to exactly one selected voxel, so
their A and C columns were proportional:

    X_Ai = g_i(B_i)e_ri, X_Ci = h_i(B_i)e_ri.

Each pair supplies one independent null direction
    (delta A_i, delta C_i) = (h_i, -g_i),
giving the structural bound rank <= 4334 - 20 = 4314 for this fixed support and
any legal positive B. EIGEN, SPQR, four legal B starts, and an independent
no-tolerance QR/SVD reduction all reproduced rank 4314 before outer search.
Changing initialization or accelerating factorization could not remove this
full-ABC structural obstruction.

This was a historical full-ABC parameterization result. It does not establish
failure of the later observable-halo or target-estimability formulation. The
diagnostic's optional Jacobi cross-check timed out and its conservative
rank-boundary diagnostic remained unconfirmed; those numerical limits do not
invalidate the separate structural argument.

### 4.2 Observable-halo result

Singleton-halo profiling removed the diagnosed initial A/C obstruction without
changing observation rows, support closure, or component IDs. On the prepared
6Z6U input there were 27 singleton halos on 26 rows, leaving 2,165 FullABC atoms
and 262,775 informative rows. Two initial constrained A/C solves had KKT
certificates and ranks 4,194 and 33. Seventeen other FullABC halos had exactly
two informative rows, so a valid initial A/C solve did not establish their
width identifiability.

The complete observable-halo 6Z6U analysis did not finish. The SPQR correctness
attempt stopped at 3,600.02 s and 4.29 GiB sampled RSS inside outer search under
the 3,600 s / 8 GiB gate; no endpoint, convergence, uncertainty, peeling,
database save, or export was completed. A later optimized attempt was stopped
by the user at 180 s / 3,198 MiB and is incomplete evidence. The 600 s reviewed
baseline also stopped, at 600.03 s / 3.70 GiB. Single-128 and single-512 completed
historical comparisons, but do not turn the 6Z6U result into an acceptance.

The completed Single-128 and Single-512 analysis-export comparisons each had
three runs per version. Reviewed stage-3 versus optimized medians were 4.68 s
versus 3.34 s (1.40x) and 289.79 s versus 122.93 s (2.36x), respectively.
All completed runs passed runtime convergence and paired endpoint checks.
Maximum sampled RSS was 0.608 / 0.551 GiB for Single-128 and 2.964 / 3.020 GiB
for Single-512. These measurements predate removal of experiment instrumentation;
they are not a current performance claim.

See the [observable-halo report](joint-observable-acceptance.md) for the paired
performance measurements and the [current parameterization contract](joint-observable-halo.md)
for the retained behavior.

### 4.3 Target estimability and weak halos

Target convergence and full-parameter convergence are separate. The complete
large 6Z6U target-estimability solve has not been established; truth replay,
small-case correctness, and preparation-only census results are not a completed
large target certificate. The local interior target certificate is not a global
uniqueness proof, and its large dense parameter-space compact factors have no
demonstrated 6Z6U resource bound.

The recorded preparation-only target study reproduced 1,559 targets, 2,192
contributors, 262,801 rows, 866,216 memberships, and components of 2,167 and
25 atoms. Its independent scalar truth replay had maximum absolute error
2.382330182015835e-7, RMSE 1.8633138857787364e-8, and normalized objective
3.185058819046246e-16. Truth was used for validation only. The large solve
budget remained unset, and the manifests, truth files, receipts, and replay tool
were removed.

Weak-halo controls show why passing numerical rank is not enough for practical
recovery. In the retained weak fixture the weakest projected-width singular
value was about 1.50e-10 versus 0.548 for the strong direction, while numerical
rank passed. The production endpoint still had local correction 0.03493886484;
a bounded restart lowered objective but did not establish convergence, and
independent 50/100-digit controls confirmed the unresolved correction. This is
evidence of practical sensitivity limitation, not proof of exact rank
deficiency. See the active [weak-halo diagnostic](joint-weak-halo-attribution.md)
and its [input snapshot](figures/joint-validation/weak-snapshot.json) and
[diagnostic record](figures/joint-validation/weak-halo.json).

## 5. Workflow, selection, and persistence evidence

### Partial selection and halo closure

Selected non-hydrogen atoms fixed the historical observation rows. Every
eligible contributor intersecting those rows was included once; halo did not
expand rows or recursively add its neighbors. A bridge halo had one parameter
vector. The exhaustive small-domain oracle covered outside-map halos, unobserved
selected targets, hydrogen exclusion, exact cutoff, clipping, anisotropic
spacing, and nonrecursive closure.

The matched/omitted-halo synthetic control used the same observation rows.
Complete closure gave normalized objective 8.139650468e-32 and target A/B/C
errors below 3e-15. Omitting the halo gave objective 0.07347280397 and target
errors delta-A=1.909360106, delta-B=0.1355602156 Angstrom,
delta-C=0.1243884602. This demonstrates parameter compensation in that
synthetic control; it is not a real-data accuracy claim. A separate weak
partial-selection case retained available states but failed local correction in
all three runs. Small synthetic success does not establish large-scale resource
bounds.

See the [partial-selection report](joint-component-partial-selection-acceptance.md)
and [partial-selection receipt](joint-component-partial-selection-acceptance.json).
The current closure rules are owned by the [runtime guide](joint-component-runtime.md).

### Stage and persistence evolution

The stage-integration work retained target/halo roles and explicit unavailable
states, kept Joint points out of OLS/MDPDE diagnostics, and made stage provenance
and downstream invalidation explicit. Post-fit peeling used the fixed Joint rows,
included halo contributions, and retained stencil coverage failures. Missing
evidence did not produce fabricated uncertainty or a substitute group posterior.

The S1-S5 consolidation paired pre/post numerical results without relaxing frozen
tolerances. The final persistence milestone used SQLite v19 and neutral analysis
document v2, moved v17/v18 adapters to load/migration, rejected conflicting
duplicate representations, and rolled back schema and records on failed
migration or writes. The SQL migration fixture was generated by the v18
implementation at df138da5. These are historical acceptance facts; current
workflow and persistence semantics are owned by the [stage integration](joint-stage-integration.md)
and [data I/O architecture](architecture/dataobject-io-architecture.md).

## 6. Historical resource envelope

The full-command resource study used an Apple M1 iMac (8 cores, 16 GiB,
macOS 26.6.2 arm64), Release builds, one numerical worker, and a nominal 100 ms
sampled 4 GiB process-tree watchdog. Input hashes, sampled-memory gaps and
per-case phase records are retained in
[resources.json](figures/joint-validation/resources.json). The analysis/export
pair shared a 600 s deadline.

| Historical workload | Result | Interpretation |
| --- | --- | --- |
| Connected 128 | 3/3 completed in 17.693-19.740 s; all passed, maximum sampled RSS 368 MiB. | Demonstrates only this historical case/host. |
| Connected 512 | 0/1 completed; 600.030 s time limit, 759 MiB sampled RSS; no exported result. | Numerical outcome unknown. |
| Connected 2,167 | Not run after the 512 stop. | No connected-2,167 result. |
| Multi-component, total 2,167 | 0/1 completed; 600.032 s time limit, 1,570 MiB sampled RSS; no exported result. | Does not establish many-small-component scalability. |
| 6Z6U, 2,192 contributors | 0/1 completed; 600.037 s time limit, 1,956 MiB sampled RSS; no exported result. | Does not establish 6Z6U convergence or failure. |

The 37,406-atom catalogue control completed only with two contributors and 681
rows; it is not evidence for running the full catalogue domain. A process
timeout/resource stop is not a numerical failure. Process completion is not
convergence or scientific qualification. Sampled RSS is not a hard OS limit;
the capped cases lack final OS peak readings. No capped case was retried with
more time, split into artificial components, or run with relaxed checks.

## 7. Statistical evidence

The active [noise/mismatch protocol](joint-noise-mismatch-validation.md) retains
the complete design and reproduction commands. Its fixed-seed pilot generated
326 inputs and 450 fitted outcomes: all retained a state, 65 passed the current
runtime checks, and 385 failed. Recovery reruns were not pooled as more
replicates. The matrix is deliberately unbalanced and tests one noise-correlation
length and one position-mismatch type; it is not a population success-rate,
coverage, or real cryo-EM validity study. Process completion and numerical
qualification are reported separately.

## 8. Explicitly unproven boundaries

- Full 6Z6U end-to-end convergence has not been demonstrated.
- Large connected-component scalability has not been established.
- Observable-halo profiling removes the diagnosed full-ABC structural
  obstruction, but does not establish complete 6Z6U acceptance.
- A complete large 6Z6U target-estimability solve has not been established.
- Numerical full rank does not guarantee practical parameter sensitivity or
  recovery for a weak halo.
- Search termination does not establish runtime convergence.
- The historical operator/Schwarz campaign is insufficient for production
  promotion; the production default remains LegacyCompact.
- The bounded SPQR rank prototype can return Unavailable on a large case and
  does not guarantee definitive rank within its configured budget.
- A timeout or resource stop is neither a numerical failure nor a convergence
  result.
- Low objective does not establish parameter recovery or a global optimum.
- Small synthetic controls do not establish real cryo-EM statistical validity.
- The Stage-B pilot is not a population-level precision or coverage study.
- Historical benchmark values are not current performance guarantees.
- Process completion alone is not numerical or scientific qualification.

## 9. Historical retrieval

The canonical index and adjacent inventories identify the surviving evidence:

- The final historical runtime baseline was
  ba2449faf1e1b163f1b4209c3ef386ff6516de43; the pre-extraction component
  reference was 6f30510c06def630226a3dab25e9bcc579d72687.
- The accepted v1 implementation baseline was
  7c0b7a4543fdc43390ca59b3aa3b5d9267fbaac7; its report records runtime source
  fingerprint e601fd9b727d59ec9668e08a1a3949ddef2365eca859c688ae4737767435c450.
- The 6Z6U initial-rank investigation starts at 5f06b027; observable-halo
  profiling starts at d634fa0d.
- Sparse-backend acceptance records baseline
  ceb6155992b9c891ab3fa878e875c6b10aa0590f and candidate production
  fingerprint 9718531067e72cd3ecf180d12ae8033b8e39360b7fe1d49d24ae847de4240fa0.
- Compact-SVD acceptance records baseline ce58c89747d4091f91967e7f4bd9c7890d51c203
  and candidate e32919f3; endpoint-reference comparison is e32919f3 against
  c9e0f8c9.
- Operator search uses pristine baseline
  c6c869cd4f4939104dfde96984ae17f736a8ce4a; bounded rank work starts at
  5f61bb6e.
- Analysis consolidation starts at c5619e0c; its migration fixture was
  generated by v18 at df138da5.
- Partial selection was tested from base
  bb72080c4878ee2333ad3e084f5e7a9d687abf33, with production source SHA-256
  e346827a2d7ffe05f4d83bcf06c8b92f416d85825f9d4cf3d8f22a20dee73e34.

For a historical source file or artifact, retrieve the recorded commit/path in
an isolated checkout, for example:

    git show ba2449faf1e1b163f1b4209c3ef386ff6516de43:docs/developer/figures/joint-abc-components/search-records-a.tar.gz > /tmp/search-records-a.tar.gz

The [retired-artifact inventory](figures/joint-component-evidence/retired-artifacts.json)
records former paths, byte sizes, source commit, and SHA-256 values for the
retired Joint tree. The [endpoint-reference archive inventory](figures/joint-reference-acceptance/archives.json)
records its former paths, retrieval commits, sizes, and hashes. The tiled,
transition, fixed-action, operator, sparse, compact, resource, and statistical
receipts remain under figures/; verify retrieved bytes against their
inventories before use.

Known unavailable material is explicitly **absent**:

- The older certification report's scientific-records.tar.gz and
  frozen-diagnostics/frozen-endpoint-records.tar.gz are absent from the final
  baseline tree and are not promised recoverable artifacts.
- Sparse acceptance's measurements.json.gz and command-exports.tar.gz are
  absent; their previously recorded hashes remain in
  [sparse verification metadata](figures/joint-sparse-acceptance/verification.json).
- Operator-search receipts.tar.gz is absent; its recorded hash remains in
  [the operator manifest](figures/joint-operator-search/manifest.json).
- Initial-rank experiment-specific tools, tests, captures, builds, and data
  bundles were removed; only the text finding remains.
- Observable-halo experiment runners, diagnostic extensions, raw measurements,
  source archives, and temporary builds were removed; only the report's text
  summaries remain.
- Target-estimability manifests, truth files, receipts, replay entry point, and
  reproduction tools were removed; only the summarized results remain.
- Bounded-rank reproduction instructions and experiment artifacts were removed;
  the report explicitly records that only textual results remain.
- Older certification scientific-records archives are not required inputs; the
  fixture catalog identifies the actual retained component archive.

The retired reports and archives remain in Git history where the inventories
identify them. Their former paths may not exist in the current checkout.
