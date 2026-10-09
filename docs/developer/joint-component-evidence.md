# Joint component historical evidence

This is the canonical index of retired Joint Component experiments and their
historical results. It records what was tested, what the evidence established,
and which questions remain open. Current implementation and numerical behavior
belong in the current-owner documents below; stable mathematical definitions
belong in the two formal contracts. Historical benchmark results are not current
performance measurements, and historical success does not establish
cross-platform behavior or real-data generalization.

Historical entries are not uniformly reproducible. Entries whose
experiment-specific artifacts have been retired are conclusion-only historical
records. They are not promised as reproducible experiments from the current
checkout.

## 1. Scope and ownership

| Concern | Current owner |
| --- | --- |
| Runtime estimator and selection domain | [Runtime guide](joint-component-runtime.md) |
| Sparse backend and shared factorization | [Sparse backend](joint-component-sparse-backend.md) |
| Compact SVD | [Compact SVD](joint-component-compact-svd.md) |
| Fixed factors and normal actions | [Fixed actions](joint-fixed-actions.md) |
| Historical profile/operator evidence | [Retired profile-operator record](joint-profile-operator.md) |
| Historical OperatorPcg evidence | [Retired OperatorPcg record](joint-operator-search.md) |
| Observable halo | [Observable-halo parameterization](joint-observable-halo.md) |
| Target estimability | [Target estimability](joint-target-estimability.md) |
| Offline weak-halo diagnosis | [Weak-halo guide](joint-weak-halo-attribution.md) |
| Statistical research protocol | [Noise/mismatch experiment](joint-noise-mismatch-validation.md) |
| Current benchmark commands | [Benchmark guide](joint-benchmark.md) |
| Historical assessment and end-to-end evidence | [Benchmark guide](joint-benchmark.md#historical-evidence) |
| Historical operator factor ownership and lifetime | [Retired OperatorPcg record](joint-operator-search.md) and [ownership campaign](figures/joint-factor-ownership-r1/campaign-manifest.json) |
| Workflow integration | [Stage integration](joint-stage-integration.md) |
| Data and persistence | [Data I/O architecture](architecture/dataobject-io-architecture.md) and [Joint JSON/CSV contract](commands/potential_analysis.md#provenance-and-map-units-joint-json-schemas-3-4-and-5) |
| Certification mathematics | [Certification contract](joint_abc_certification_contract.md) |
| Component mathematics | [Component contract](joint_abc_components_contract.md) |

### Current scalability axes

The active benchmark separates FixedNeighbor search-only attribution, production
solve with endpoint certification and `RuntimeConvergence`, workflow/persistence
smoke, and complete-command smoke. Current measurements use the shared EIGEN or
SPQR backend and the unchanged production policy: core 12, Forward order, 30
maximum sweeps, local `LegacyCompact`, and at most one trusted local update per
block visit.

Rank certification, Krylov iteration scaling, operator throughput, Schwarz
geometry and operator-factor residency are closed historical dimensions. Their
compact results remain below and in the figures, but no current driver or
runtime route depends on them.

## 2. Algorithm-formation evidence

These findings describe the historical experiments that shaped the current
contracts. This section is the canonical summary of their scientific
conclusions.

| Historical question | Evidence and limitation | Contract retained today | Provenance |
| --- | --- | --- | --- |
| Do sampled observations match the frozen support and sampler? | On the truth-neighbor control, compact prediction versus sampler RMSE was **5.17001e-16**; float32 residual was **2.91318e-8**. This checks observation consistency, not convergence with unknown neighbors. | Preserve the sampler, strict support boundary, contributor membership, and post-sum quantization controls. | [Fixture catalog](../../tests/fixtures/joint_component/catalog.json); [retired path and hash inventory](figures/joint-component-evidence/retired-artifacts.json). |
| Does an estimated-neighbor sweep improve every atom? | Across eight checkpoints, aggregate A/B/C errors improved by about **19–21 / 25–26 / 4.7–4.9 times**, while C worsened for **70–76 of 168 atoms**. Aggregate improvement is not per-atom recovery. | Keep per-atom evidence and frozen support provenance; do not infer individual recovery from aggregate error. | Historical source and artifact hashes are indexed by the [retired-artifact inventory](figures/joint-component-evidence/retired-artifacts.json). |
| Does matched joint A/C solve qualify from full design rank? | All 16 joint runs exhausted 2,000 evaluations without qualification despite design rank 336; frozen-neighbor controls qualified 1,335/1,344 blocks. Rank and lower C error do not establish stationarity. | Search completion, a usable state, and numerical qualification are distinct. | [Runtime contract](joint-component-runtime.md); [retired path and hash inventory](figures/joint-component-evidence/retired-artifacts.json). |
| What did Experiments A/B/C establish about composite MDPDE? | A: baseline-best28 at alpha 0 had A/C RMSE **0.011471954 / 0.00063732864**; alpha 0.1 qualified at **0.0091138914 / 0.00049314988**, while alpha 0.5 and 1 exhausted their budgets. B: alpha 0 qualified at **0.0091707311 / 0.00032313699**; positive-alpha runs exhausted 100 evaluations. C: two starts both exhausted 100 evaluations; one endpoint had A/C RMSE **0.017916446 / 0.00011271024** and effective sample fraction **0.8983%**. The C comparison had no independent repeat or common-alpha control, so it neither isolates the cause nor disproves heterogeneous MDPDE. | These runs are not prerequisites for current joint LS and do not justify changing its objective. | [Retired path and hash inventory](figures/joint-component-evidence/retired-artifacts.json). |
| Does fixed-B variable projection support the joint profile? | The historical profile work retained independent double/float32 recovery and full profile-Jacobian checks at nonzero residual. It is finite-fixture evidence, not a general convergence result. | Keep the constrained linear subproblem, full nonzero-residual correction, and independent reference checks in the current numerical contract. | [Certification contract](joint_abc_certification_contract.md); [retired path and hash inventory](figures/joint-component-evidence/retired-artifacts.json). |
| Did the coverage set make weak and degenerate cases successful? | Nine immutable datasets remain in the fixture catalog. Weak, active-face, zero-signal, and duplicate cases retain their original resolution, correction, width-identification, and state-availability limitations. | Keep these cases as distinct regression outcomes; do not turn unavailable or failed evidence into passes. | [Fixture catalog](../../tests/fixtures/joint_component/catalog.json); [fixture guide](../../tests/fixtures/joint_component/README.md); [v1 freeze receipt](joint-component-v1-acceptance.json). |
| Did guarded certification qualify all cases? | The frozen baseline recorded 128 paired branches, 64 required global regular pairs, and 192 component-local scopes (128 regular); 96 regular globally restricted scopes use a separate certification scope. These counts apply only to the recorded fixtures and policy. | Keep global, restricted, and local certification scopes separate; local certification cannot replace global qualification. | [Transition receipt](figures/joint-component-evidence/transition-acceptance.json); [certification contract](joint_abc_certification_contract.md); [component contract](joint_abc_components_contract.md). |
| Did exact component decomposition preserve parent evidence? | Transition receipts compare 1,248 migrated search/state/context records, 2,096 global/restricted audit records, and 3,712 baseline records including separate local scopes. Four representative local certificates were rerun; the migrated offline driver did not rerun the full audit matrix. | Preserve structural partition, shared parent scale/rank context, same-state assembly, and failure isolation. | [Transition receipt](figures/joint-component-evidence/transition-acceptance.json); [global pairs](figures/joint-component-evidence/global-pairs.csv); [scope matrices](figures/joint-component-evidence/global-and-restricted-scopes.csv) and [local scopes](figures/joint-component-evidence/local-scopes.csv). |

Public objectives are normalized by the parent observation scale squared;
archived experiment objectives use raw half-RSS. The transition receipt is
historical comparison evidence, not a new runtime measurement. Its 22 default
and 3 extended/offline test counts came from the recorded macOS/compiler
configuration and do not establish cross-platform behavior.

The v1 freeze accepted fixed-position, all-non-hydrogen, structural 2.5 Angstrom
support, equal-weight Guarded joint LS as an opt-in estimator; it did not change
the default two-stage estimator or package version. Historical baseline,
near-0.02, weak-1e-4, active-a, zero-signal, and duplicate controls retained
their distinct convergence, derivative-resolution, width-identification, and
state-availability outcomes.

## 3. Backend and numerical evolution

| Milestone | Historical evidence | Retained limitation | Current owner |
| --- | --- | --- | --- |
| Tiled derivative and assessment | For heterogeneous-168 / first-stage-float32 on the recorded macOS arm64 host, peak RSS changed from **1,252.14 MiB to 359.14 MiB (71.3% lower)**; total time changed from 28.7429 s to 27.6702 s (3.7% lower), with assessment slightly slower. Fixed-state parity passed at the same states. | This establishes a storage reduction for that host and case, not general scalability or a universal speedup; sparse factor storage and compact quadratic factors remain. | [Profile operator](joint-profile-operator.md); [measurement JSON](joint-component-tiled-backend.json). |
| Sparse EIGEN/SPQR backend | The experiment-specific receipts have been retired from the current tree. This row is the historical record for the pre-promotion comparisons: completed fixed-state EIGEN/SPQR comparisons agreed within recorded tolerances; historical 6Z6U primary rank was 4314/4334; and a historical single-512 complete-command run timed out. The later current FixedNeighbor promotion is recorded separately in `figures/joint-fixed-neighbor-backend-qualification-r1/`. | These fixed-state comparisons do not establish the current backend choice or a completed large historical workflow. | [Sparse backend](joint-component-sparse-backend.md). |
| Compact SVD | The historical compact campaign passed its recorded numerical and fixed-state performance gates. For SPQR fixed states, the auto/legacy median ratios were 0.268861 for Single-128 reference plus derivative (maximum 1.10), 0.730425 for Heterogeneous-168 (maximum 1.10), 0.009068 for Single-512 reference SVD, 0.005008 for Single-512 free-design SVD, and 0.184551 for Single-512 reference plus derivative (the latter three maximum 0.70). Three samples per backend passed the numerical audits. The Single-512 candidate completed three full commands (median 289.238 s); the ce58c897 command baseline stopped at 600 s, so no completed-command speedup is established. Baseline `ce58c89747d4091f91967e7f4bd9c7890d51c203` had production fingerprint `9718531067e72cd3ecf180d12ae8033b8e39360b7fe1d49d24ae847de4240fa0`; candidate `e32919f3` had fingerprint `f1ac45c36945e8e89e588af5c305158a7723b268a21c5c7ab0b22a3ae676a755`. | These results are tied to their source, backend, host, and cases; they are not current performance guarantees. | [Compact SVD](joint-component-compact-svd.md); [verification receipt](figures/joint-compact-acceptance/verification.json). |
| Fixed factors and normal actions | The recorded campaign completed **72/72 processes** and **24/24** three-repetition A/B, B/C, and A/C comparisons. Normal actions used two Q/Q′ calls versus six for the composed control. | These fixed-state/fixed-step results do not promote a different full-search default or establish large-workflow scalability. Some original step subtotals omitted gradient setup; corrected bounds are not exact end-to-end step times. | [Fixed actions](joint-fixed-actions.md); [summary](figures/joint-fixed-actions/summary.json) and [comparison decisions](figures/joint-fixed-actions/comparison.json). |
| Endpoint reference | The 43.81-minute historical campaign completed all 18 fixed-state audits and all 128/512 analysis-export comparisons. Single-128 changed from 5.229 s to 4.912 s; single-512 from 307.805 s to 300.564 s, with peak RSS from 2.500 GiB to 1.872 GiB for single-512. Search reference evaluations changed from six to zero. Detailed comparisons and archive provenance remain in the [summary](figures/joint-reference-acceptance/summary.json) and [archive inventory](figures/joint-reference-acceptance/archives.json). | These are incremental historical results on the recorded host; search and endpoint assessment remained the dominant costs. | [Compact SVD](joint-component-compact-svd.md). |
| Profile operator | The experiment-specific receipts have been retired from the current tree. The summary below is the canonical historical record. All 18 fixed-state audits passed, establishing operator/action parity on finite recorded controls. The 10,000-atom workload was preparation-only and does not prove solve, search, or full-workflow scalability. Historical 6Z6U runs hit the sampled RSS watchdog and established no numerical verdict. | Fixed-state and preparation-only results do not establish full-workflow scalability or a 6Z6U numerical outcome. | [Profile operator](joint-profile-operator.md). |
| Operator LM and Schwarz search | The experiment-specific receipts have been retired from the current tree. The summary below is the canonical historical record. The full promotion campaign was incomplete and required comparisons were not all run. The historical Single-512 performance gate did not support promotion; the qualified FixedNeighbor route was subsequently adopted as the production default. No 6Z6U promotion conclusion was established. | The campaign is closed and does not define a current solver route. | [Operator search](joint-operator-search.md). |
| Bounded SPQR rank prototype | Small controls established full-rank and deficient decisions against the dense oracle. Two representative SPQR prototypes returned Unavailable before completing the reconstruction certificate. | Unavailable is not deficient; no definitive 128-atom or 512/2,000-atom rank result or scaling guarantee was established. A separate Single-512 Schwarz attempt exceeded the sampled RSS threshold on both backends and had no endpoint. The rank prototype is historical; shared sparse factorization remains current. | [Sparse backend](joint-component-sparse-backend.md). |

These fixed-state results must not be confused with the full-command resource
envelope in §6.

## 4. Structural and identifiability boundaries

### 4.1 Historical full-ABC 6Z6U obstruction

At historical base 5f06b027, the selected domain had 2,192 contributors,
262,801 rows, and two components. Its largest component had 2,167 atoms,
262,282 rows, and 4,334 free A/C columns. Historical per-atom B0 mapping was
unavailable, so its starts were not a bitwise reproduction of the earlier fit.
Twenty singleton halo atoms each contributed to exactly one selected voxel, so
their A and C columns were proportional:

    X_Ai = g_i(B_i)e_ri, X_Ci = h_i(B_i)e_ri.

Each pair supplies one independent null direction
    (delta A_i, delta C_i) = (h_i, -g_i),
giving the structural bound rank <= 4334 - 20 = 4314 for this fixed support and
any legal positive B. Both EIGEN and SPQR rejected the first all-free A/C face
before outer search. Four SPQR starts (production B0, 0.9 B0, 1.1 B0, and
per-atom simulation truth) and an independent no-tolerance QR/SVD reduction all
reproduced rank 4314; weak-direction actions were checked against the original
sparse matrix. Nonnegative A does not restore uniqueness because a change in A
can be offset by signed C. Changing initialization or accelerating
factorization could not remove this full-ABC structural obstruction. The
measured weak subspace was concentrated in these halo A/C coordinates.

This was a historical full-ABC parameterization result. It does not establish
failure of the later observable-halo or target-estimability formulation. The
diagnostic's optional Jacobi cross-check timed out and its conservative
rank-boundary diagnostic remained unconfirmed; those numerical limits do not
invalidate the separate structural argument.

Two remedies were considered at the time: represent each one-row halo by its
single observable contribution while leaving its individual A/C/B unidentified,
or expand the observation domain and recompute contributor closure. Observable
profiling became the maintained parameterization. Domain expansion changes the
statistical problem and can create further halos, so it has no guaranteed
full-rank outcome.

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
historical comparisons, but do not establish a 6Z6U numerical verdict.

The completed Single-128 and Single-512 analysis-export comparisons each had
three runs per version. Reviewed stage-3 versus optimized medians were 4.68 s
versus 3.34 s (1.40x) and 289.79 s versus 122.93 s (2.36x), respectively.
All completed runs passed runtime convergence and paired endpoint checks.
Maximum sampled RSS was 0.608 / 0.551 GiB for Single-128 and 2.964 / 3.020 GiB
for Single-512. These measurements predate removal of experiment instrumentation;
they are not a current performance claim.

The paired performance measurements and incomplete 6Z6U outcomes above are
historical; the retained parameterization behavior is owned by the
[current contract](joint-observable-halo.md).

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

The experiment-specific receipts have been retired from the current tree. The
summary below is the canonical historical record. Selected atoms defined the
historical observation rows, and halo contributor closure was required for
correctness. In the matched synthetic control, omitting the halo caused strong
target compensation; complete closure produced near-zero objective and target
errors. A weak partial-selection state could remain available while local
correction failed. These small synthetic cases do not establish large-scale
scalability.
The current closure rules are owned by the [runtime guide](joint-component-runtime.md).

### Stage and persistence evolution

The stage-integration work retained target/halo roles and explicit unavailable
states, kept Joint points out of OLS/MDPDE diagnostics, and made stage provenance
and downstream invalidation explicit. Post-fit peeling used the fixed Joint rows,
included halo contributions, and retained stencil coverage failures. Missing
evidence did not produce fabricated uncertainty or a substitute group posterior.

The S1-S5 consolidation paired pre/post numerical results without relaxing frozen
tolerances. Its postprocessing comparison retained exact endpoint A/C/B,
objective, convergence, target peeling, and covariance in all 18 paired runs;
the bounded fixtures did not establish a general speedup or size limit. The
final persistence milestone used SQLite v19 and neutral analysis document v2,
moved v17/v18 adapters to load/migration, rejected conflicting
duplicate representations, and rolled back schema and records on failed
migration or writes. The SQL migration fixture was generated by the v18
implementation at df138da5. The retained paired postprocessing measurements are
in [joint-analysis-consolidation-benchmark.json](joint-analysis-consolidation-benchmark.json).
The earlier v1 result format was schema 2; the
current writer is schema 5 and retains older v3/v4 reads without inferring
target evidence. These are historical workflow results; current
workflow and persistence semantics are owned by the [stage integration](joint-stage-integration.md)
and [data I/O architecture](architecture/dataobject-io-architecture.md).

## 6. Historical resource envelope

The experiment-specific receipts have been retired from the current tree. The
summary below is the canonical historical record. The campaign ran on an Apple
M1 iMac (8 cores, 16 GiB, macOS 26.6.2 arm64). A nominal 100 ms sampled 4 GiB
process-tree RSS watchdog was used, with a 600 s analysis/export deadline.

| Historical workload | Result | Interpretation |
| --- | --- | --- |
| Connected 128 | Completed. | Historical-host evidence only. |
| Connected 512 | 600 s timeout. | Numerical outcome unknown. |
| Connected 2,167 | Not run / unavailable. | No scalability conclusion. |
| Multi-component, total 2,167 | 600 s timeout. | No many-small-component scalability proof. |
| 6Z6U | 600 s timeout. | No convergence or failure verdict. |

Timeout is not a numerical failure; a missing result is not failed convergence;
completion is not scientific qualification. Sampled RSS was a watchdog
observation, not a hard OS limit.

## 7. Statistical evidence

The active [noise/mismatch protocol](joint-noise-mismatch-validation.md) retains
the complete design and reproduction commands. Its fixed-seed pilot generated
326 inputs and 450 fitted outcomes: all retained a state, 65 passed the current
runtime checks, and 385 failed. Recovery reruns were not pooled as more
replicates. The matrix is deliberately unbalanced and tests one noise-correlation
length and one position-mismatch type; it is not a population success-rate,
coverage, or real cryo-EM validity study. Process completion and numerical
qualification are reported separately. A serialization-recovery comparison
found identical input snapshots and numerical statuses, with at most a one-ULP
parameter representation difference caused by the historical experiment
wrapper's default JSON parsing; the production decoder was already precise.
This conclusion is retained without the experiment-specific recovery record.

## 8. Explicitly unproven boundaries

- Full 6Z6U end-to-end convergence has not been demonstrated.
- Large connected-component scalability remains unestablished. The current
  one-level pilot stopped before an eligible PCG point when bounded SPQR rank
  reached `rank-work-budget` at chain-128; see the
  [current campaign record](joint-benchmark.md#current-one-level-scaling-evidence).
- Observable-halo profiling removes the diagnosed full-ABC structural
  obstruction, but does not establish complete 6Z6U qualification.
- A complete large 6Z6U target-estimability solve has not been established.
- Numerical full rank does not guarantee practical parameter sensitivity or
  recovery for a weak halo.
- Search termination does not establish runtime convergence.
- The historical operator/Schwarz campaign is closed; the later FixedNeighbor
  qualification adopted FixedNeighbor as the production estimator.
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

### Current-tree retained historical evidence

- The transition receipt, scope matrices, compact-SVD verification, and
  fixed-action comparisons remain in their linked locations above.
- The endpoint-reference package retains its README, archive inventory,
  summary, and verification files.
- Frozen fixtures remain catalogued in `tests/fixtures/joint_component/`.
- Current weak-halo and noise/mismatch results remain under
  `figures/joint-validation/` and are described by their active research guides.

### Conclusion-only retired campaigns

The experiment-specific current-tree artifacts were intentionally removed for:

- Sparse backend campaign (H06).
- Operator/Schwarz promotion campaign (H07).
- Profile-operator audit campaign (H08).
- Partial-selection campaign (H09).
- Historical resource envelope (H12).

The summaries in §§3, 5, and 6 are the canonical historical records for these
campaigns; they are not promised as reproducible experiments from the current
checkout.

### Historical source references

The canonical index and adjacent inventories identify historical source points:

- The final historical runtime baseline was
  ba2449faf1e1b163f1b4209c3ef386ff6516de43; the pre-extraction component
  reference was 6f30510c06def630226a3dab25e9bcc579d72687.
- The v1 implementation baseline was
  7c0b7a4543fdc43390ca59b3aa3b5d9267fbaac7; its report records runtime source
  fingerprint e601fd9b727d59ec9668e08a1a3949ddef2365eca859c688ae4737767435c450.
- The 6Z6U initial-rank investigation starts at 5f06b027; observable-halo
  profiling starts at d634fa0d.
- Compact-SVD historical results use baseline ce58c89747d4091f91967e7f4bd9c7890d51c203
  and candidate e32919f3; endpoint-reference comparison is e32919f3 against
  c9e0f8c9.
- Analysis consolidation starts at c5619e0c; its migration fixture was
  generated by v18 at df138da5.

For a historical source file or artifact, retrieve the recorded commit/path in
an isolated checkout, for example:

    git show ba2449faf1e1b163f1b4209c3ef386ff6516de43:docs/developer/figures/joint-abc-components/search-records-a.tar.gz > /tmp/search-records-a.tar.gz

The [retired-artifact inventory](figures/joint-component-evidence/retired-artifacts.json)
records former paths, byte sizes, last-present commits, and SHA-256 values for
historical files removed from the current tree. The [endpoint-reference archive
inventory](figures/joint-reference-acceptance/archives.json) records former
paths, retrieval commits, sizes, and hashes for that retained evidence package.

Known unavailable material is explicitly **absent**:

- The older certification report's scientific-records.tar.gz and
  frozen-diagnostics/frozen-endpoint-records.tar.gz are absent from the final
  baseline tree and are not promised recoverable artifacts.
- The sparse-backend campaign's `measurements.json.gz` and `command-exports.tar.gz` were
  already absent before its remaining receipts were retired.
- Initial-rank experiment-specific tools, tests, captures, builds, and data
  bundles were removed; the structural finding is preserved in §4.1.
- Observable-halo experiment runners, diagnostic extensions, raw measurements,
  source archives, and temporary builds were removed; the conclusions and
  limitations are summarized in §4.2.
- Target-estimability manifests, truth files, receipts, replay entry point, and
  reproduction tools were removed; only the summarized results remain.
- Bounded-rank reproduction instructions and experiment artifacts were removed;
  the retained prototype results and limits are summarized in §3.
- Older certification scientific-records archives are not required inputs; the
  fixture catalog identifies the actual retained component archive.

The retired reports and archives remain retrievable from Git history where the
inventories identify a last-present commit. Their former paths may not exist in
the current checkout.
