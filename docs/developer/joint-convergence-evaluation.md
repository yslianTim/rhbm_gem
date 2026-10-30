# FixedNeighbor convergence evaluation — 6CVM Chain A

Date: 2026-10-30 (Asia/Taipei). Initial revision: `f4bce1ea21ea9d6878436670ff88447cbb00638f`.
Scope: R0 → R1 → R2 → R4 hotspot correction only → R6. No partition/order/scheduling experiments.

## Reproduction contract

- Existing `develop`, initially clean; no new branch/worktree; no Notion writes.
- Release preset, AppleClang, Eigen 5, existing ROOT/OpenMP/UMAP dependencies; no added dependency.
- All builds use `--parallel 4`; CTest runs serially with internal thread limit 4.
- Model: `6cvm.cif`, SHA-256 `8483d9096eadd9fb6538ff88b78daf010957f9718779ee613c97764a985d7a1d`.
- Map: `emd_7770_unsharpened.map`, SHA-256 `e4b6e7bf92d0452f14e2f2f8c2dcca0e21818de62ebc3808040f7086ab584da7`.
- Non-hydrogen Chain A targets; map normalization divisor 0.1593869396082229.
- 8,193 targets; 508,449 fitting rows; 1 connected component; 742 fixed cores of at most 12 atoms.
- Original-domain halo 1,645; freely fitted halo 0. Original rows 622,982; excluded rows 114,533.
- Original FirstStage initialization, globally calibrated Bmax=5.608794450554618 Å; 11 projected starts.
- Research runs rebuild the actual model/map fitting domain and require exact saved fitting row identities.
  They reuse the captured original FirstStage initialization, apply `MakeBoundedInitialization`,
  and start with zero A/C, exactly as production search. They never restart from candidate widths.
- Practical-v2 thresholds stay A/C ≤ 1e-7, projected log-width ≤ 1e-9.
- Research executables, full state files, SQLite and verbose logs live in `/tmp/rhbm-convergence`.
  Historical diagnostic evidence is in `/tmp/rhbm-6cvm-diagnosis`; these paths are local scratch assets.

## R0 flow audit

1. `PotentialFittingWorkflow.cpp` builds the pure-halo-support-exclusion domain; only targets are free.
   FirstStage donors calibrate one global bound; projection changes working starts, not donor estimates.
2. `Components.cpp` partitions by shared observations and maps component rows/atoms to compact local views.
   `StructuralTopology::BuildCores(12)` fixes the structural cores before search. Metadata and local row
   mappings are prepared once; no per-sweep repartition/reinitialization occurs.
3. `SearchFixedNeighborComponent` creates complete accepted eta/beta/prediction/residual/objective state.
   A/C start at zero. Every complete sweep visits every original core in forward order.
4. A block builds effective observations from accepted residual plus the old core prediction.
   `SearchProfile` reuses constrained A/C and bounded active-set LM, with at most one accepted local
   width update per visit. Width failures cannot fall back to a merely trusted initial profile.
5. Missing trusted candidates set a pending recovery. A second failed visit stops with
   `block-recovery-exhausted`; pending recovery at budget uses `block-recovery-sweep-budget`.
   These differ from normal `block-sweep-budget` exhaustion.
6. Local objective increase leaves the accepted state intact. Otherwise the code copies eta/beta,
   overwrites only the core, replays affected rows (full replay on memory fallback), checks finite
   parameters/prediction/residual/objective and local/global objective consistency, then commits.
   Rejected candidates never partially write the accepted state. There is no new solver path.
7. Each complete sweep performs full prediction/objective replay and an independent streaming global
   accepted-state evaluation. Cache consistency must pass before objective reconciliation.
8. A/C uses column-scaled projected KKT. Width stationarity uses bound-aware projected log-width KKT;
   infeasible upper bounds force nonstationarity. Raw width gradient is auxiliary and can remain
   nonzero on an active upper face. Complete-sweep confirmation is required.
9. Stationarity triggers endpoint certification; valid but unqualified endpoints can be retried
   after 10 sweeps. Terminal budget also certifies the actual final accepted state.
   Certification preflights memory, recomputes the primary endpoint, independently replays it,
   checks all original feasibility/stationarity/consistency conditions, and checks endpoint trust.
10. Only qualified certification plus accepted-state integrity creates `component.state`.
    Assembly checks mappings, component qualification and global replay. Unqualified candidates
    remain diagnostics. `PotentialFittingWorkflow` returns before halo/peeling/local/group fitting
    when the Joint result is unqualified. CLI preserves diagnostics and returns exit 3.

Per-candidate full eta/beta copies exist, but preserve the established transaction boundary.
They are measured, not refactored in this work. No evidence yet establishes ill-conditioning,
non-identifiability or a local-solver defect as the cause of the slow tail.

### Historical endpoint, independently reproduced before this task

The existing same-revision CLI reproduction exited 3 after 200 sweeps (428.123 s search).
A/C KKT=3.7578613380675596e-7; projected width KKT=1.6868359486998941e-9;
raw width gradient=3.2339581943321713e-8; objective≈0.01977142389929956.
Endpoint/replay qualification failed stationarity, while consistency and state integrity passed.
Independent KKT discrepancy=1.0842021724855044e-19; width gradient discrepancy=2.0466996070848568e-18.

| Serial | Label residue / atom | Zero-based core | Failed condition |
|---:|---|---:|---|
| 5497 | LEU 684 C | 478 | A/C |
| 5509 | PRO 685 CD | 463 | A/C, width |
| 5496 | LEU 684 CA | 478 | A/C, width |
| 5498 | LEU 684 O | 478 | A/C |
| 5508 | PRO 685 CG | 463 | A/C |
| 5495 | LEU 684 N | 439 | width |
| 5654 | ALA 704 C | 493 | width |

Label 684/685/704 correspond to author 685/686/705. These identities must not be interchanged.
LEU/PRO involve distinct original cores; ALA is a separate candidate region.

An initial offline attempt was stopped and excluded: its vectorized Eigen `log` differed from
production's per-atom `std::log`, perturbing the search trajectory. The accepted study uses the
production scalar conversion and checks the complete 200-sweep eta/A/C endpoint against the CLI record.
The first CTest invocation used OMP_NUM_THREADS=1, which conflicts with one existing HRL thread-restore
fixture's two-thread assumption. That group passes with OMP_NUM_THREADS=2; the full baseline suite
is rerun in that environment (joint groups retain their own one-thread setting).

### Verified R0 result

The new 200-sweep run matches **all 200** historical rounded A/C/raw-width/accepted-block records.
All 8,193 final eta values and all 16,386 A/C coefficients are **exactly equal** to the historical
CLI candidate. Endpoint A/C/width KKT and independent discrepancies above are reproduced exactly.
The stopping reason is `block-sweep-budget`, accepted-state integrity passes, unresolved recovery=0.
Certification and independent replay remain unqualified because stationarity fails; no official state
or assembly is available. This establishes feasibility of a continuous R1 extension.

Full-precision sweep evidence: [R0 CSV](evidence/joint-convergence-r0.csv).
The CSV distinguishes rejected updates from genuinely unchanged cores (the existing solver's
`unchanged_blocks` counter includes safe objective rejections). Across 200 sweeps: 130,433 accepted,
17,898 rejected, 69 unchanged, 0 failed. Sweep 200: 312 accepted, 428 rejected, 2 unchanged.

| Measured phase | Seconds |
|---|---:|
| Complete sweep search | 392.970465 |
| Local profile search (included above) | 299.813554 |
| Affected candidate replay (included above) | 34.298110 |
| Candidate parameter copies (included above) | 0.639624 |
| Full sweep replay (included above) | 12.336593 |
| Global streaming KKT (included above) | 24.573978 |
| Endpoint certification | 0.252105 |
| Independent replay (part of certification) | 0.129337 |
| Local conditioning assessment | 1.610276 |
| Solver wall time incl. preparation/certification/conditioning | 395.258723 |

Peak process RSS=873,742,336 bytes. This offline process also retains the decoded historical record,
map and test-only block records; it is not a solver-workspace measurement. Global search preflight
passes at 219,762,986 bytes; endpoint preflight passes at 235,957,873 bytes; candidate scratch peak
198,232 bytes. `/usr/bin/time -l` could not access `kern.clockrate` under process restrictions;
RSS comes from the study executable's successful `getrusage` call, not that wrapper.

Full baseline Release build and **19/19 CTest groups pass** (34.01 s). The only production change in
R0 exposes the already-computed projected width KKT as `width-kkt` beside raw `width-grad` in CLI
progress. Its regression verifies distinct values; no numerical decision or serialized field changes.

The final 20-sweep window's best A/C KKT improves by 32.75%, projected width by 32.58% between
its two halves. This supports testing an extension; it does not establish eventual qualification.

Shared-row incidence confirms local coupling: serial 5495 shares contributor incidences with cores
438/439/463/478/479; 5496 also with 468; 5497–5498 involve 439/463/468/478/479;
5508–5509 involve 439/463/468/478/479 (5509 also 461). Serial 5654 instead involves
252/493/494/495. Counts are contributor-row incidences, not unique rows across multiple neighbors.
Final-minus-initial eta changes for 5495/5496/5497/5498/5508/5509/5654 are respectively
0.02216217 / -0.71758790 / 0.57404160 / -0.38523613 / -0.30884112 / 0.06191469 / -0.10479786.
All seven are interior to the calibrated width upper bound at this endpoint. In particular, serial
5497 has A=0 but nonzero C, so its width still affects the charge contribution.

## R1 continuous original search — Qualified

A single run with a 600-sweep cap, identical scalar-log initialization and unchanged local solvers,
cores/order/acceptance/thresholds, stopped at **sweep 250**. Its first 200 sweeps match R0 exactly
in objective, both KKT metrics, raw gradient and block counts; the 200-sweep eta/A/C arrays again
match the historical CLI candidate exactly. No width-only restart or reconstructed continuation.

| Checkpoint budget | Actual sweeps | A/C KKT | Projected width KKT | Qualified | Search seconds | Peak RSS bytes |
|---:|---:|---:|---:|---|---:|---:|
| 200 (R0 independent run) | 200 | 3.7578613380675596e-7 | 1.6868359486998941e-9 | Failed | 392.970465 | 873742336 |
| 300 | Not reached: qualified at 250 | — | — | Not executed | — | — |
| 400 | Not reached: qualified at 250 | — | — | Not executed | — | — |
| 600 (configured cap) | 250, early stop | 4.2690769063921663e-8 | 9.699454295297727e-10 | Qualified | 487.347409 | 986546176 |

The first actual A/C pass is sweep **241**, first width pass **250**, first simultaneous pass and
complete endpoint certification **250**. Independent A/C KKT=4.2690769069993195e-8;
KKT discrepancy=6.0715321659188248e-18; gradient discrepancy=2.5283204579171602e-18.
All feasibility, replay and endpoint checks pass; integrity=true, unresolved recovery=0.
The separately recomputed endpoint also passes the existing production `AssembleComponents`
checks (mapping, numerical qualification and independent replay).

Full trajectory: [R1 CSV](evidence/joint-convergence-r1.csv). Final objective=0.019770945620759782.
Accepted=144,252; rejected=41,039; genuinely unchanged=209; failed=0. Solver wall=489.667256 s,
certification=0.252499 s, independent replay=0.128226 s. Time-to-qualified excludes model/map loading,
FirstStage initialization, assembly and downstream processing, and includes complete sweep work;
the reported solver wall additionally includes preparation and final conditioning diagnostics.
The higher research-process RSS includes retained test-only records for 50 more sweeps. Global
preflight remains exactly 219,762,986 bytes and candidate scratch remains 198,232 bytes; no new
production workspace was introduced in this experiment.

The tail is not purely geometric: after sweep 200, objective descends into a lower region while
KKT rebounds (sweep 225 A/C≈1.81e-6, width≈5.52e-9), then both settle by 250. There are short-term
and odd/even variations. This run provides no evidence of a machine-precision barrier or persistent
local failure. It does not establish statistical identifiability or global optimality.

Decision: **Outcome A**, original search qualifies within 300. R2 will evaluate a 200 nominal / 300
hard budget, using recent normalized projected-KKT progress. Running to 400/600 after qualification
would have no scientific value and was intentionally avoided.

## R2 adaptive budget

Production policy: nominal 200, hard 300. No public tuning option. A constant-size 20-sweep ring stores
q=max(A/C KKT / 1e-7, **projected** width KKT / 1e-9). At the nominal boundary, all window metrics
must be finite and nonnegative, no recovery may remain pending, the best q in the most recent
10 sweeps must be ≤10, and it must improve at least 5% over the preceding 10-sweep best (which
must exceed 1). A granted extension runs the same search up to the hard cap; transient rebounds
inside that extension do not revoke the allocation. Qualification still stops search immediately.
No extra KKT evaluation, accepted-state copy, local solve or replay is introduced by this decision.

For 6CVM the boundary compares q=5.587652802821822 against 3.75786133806756 (32.75% improvement),
so it extends. A denied nominal extension keeps `block-sweep-budget`; an exhausted granted
extension uses `block-hard-sweep-budget`. Pending recovery keeps `block-recovery-sweep-budget`.
All three paths retain the original endpoint-certification and unqualified-state safeguards.

The policy is versioned `fixed-neighbor-production-v5`; the existing maximum-sweeps provenance
field now records hard cap 300. No JSON schema or diagnostic tuning field is added. The reader
also accepts the exact existing v4/200 historical pair so the actual baseline database remains
readable; mismatched v4/300 is rejected. No numerical compatibility branch is added.

Budget regression coverage includes easy early qualification with identical state/work, progress,
rebound, odd/even oscillation, separately normalized width, distant/stalled endpoints, nonfinite/
negative metrics, unresolved/recovered failures, independent component history, fixed/zero budgets,
and exact hard-cap behavior. A coupled Gaussian pair requires an actual extension (nominal 30),
qualifies with the same parameters as the original fixed larger budget, and passes production
component assembly and replay. A 31-sweep hard cap remains unqualified and publishes no state;
a one-byte memory budget fails before search/state creation. Existing recovery/candidate rollback
and preflight tests remain unchanged.

### Verified production CLI result

The complete real CLI, original files/selection and `-j 4`, exits **0** with policy v5/300,
qualifies at **250**, and saves a qualified component and assembled global state. The final eta/A/C
arrays are exactly equal to R1. Search=500.783628 s; construction=0.693987 s;
FirstStage initialization=10.846225 s; assessment=2.867448 s; assembly=0.276636 s.
Whole command incl. loading/Q-score/sampling/downstream/SQLite=575.361841 s; peak process
RSS=2,811,363,328 bytes. These whole-workflow timing/RSS values must not be directly compared
with the offline solver-only process. The budget's storage is fixed-size and independent of atoms,
rows and extension length; endpoint preflight stays 271,418,201 bytes, exactly the original CLI value.

SQLite verification: runtime and target convergence both passed; endpoint certification sweep 250;
1,645 halos reconstructed, none unavailable; 8,193 peeling entries and 8,193 post-joint local entries;
819,300/819,300 peeling samples valid; 167 group summaries. This demonstrates execution and
persistence of downstream stages, not that every separate statistical fit is scientifically identified.
Machine-readable evidence: [R2 CLI summary](evidence/joint-convergence-r2-cli.json).

Final R2 Release build and **19/19 CTest groups pass** (35.44 s), including 13 budget regressions,
actual independent multi-component extension, historical-provenance validation, and a full map/model
workflow regression whose coupled component extends past nominal 35 and qualifies at 40 before
halo reconstruction, peeling, local MDPDE and group-stage processing. Failure workflows still skip
these stages. The full-workflow test uses the existing scoped test facade, not a production switch.

An interleaved 20-pair easy-case measurement (8 atoms, identical 3-sweep qualified results) gives
median solver wall 0.0144783 s fixed versus 0.0144969 s adaptive (about 0.13% difference).
The budget object is 208 bytes, with no heap allocation. No material easy-case overhead was observed.
Decision: **promote adaptive budget**. It repairs the actual 6CVM qualification failure without
changing its accepted numerical trajectory or weakening any qualification check.

## R4 hotspot correction study

The controlled study used a test-build-only callback after all ordinary blocks and before the normal
full sweep replay. No production CLI option, block reordering/repartition, active-block scheduling or
new numerical solver is used. The callback was removed after the promotion decision; solver-study code never entered the repository.

Detection recomputes the actual accepted endpoint and uses the same column-scaled A/C projected
residual and bound-aware projected width residual as production. No A-based width exclusion is
applied. Violating seeds sharing fitting rows form separate connected regions. Each region adds
neighbors ranked by shared seed-row incidence (stable identity tie-break), up to 24 atoms; regions
with overlapping selected patches safely skip the overlap. More than 12 seeds, >6,000 patch rows,
>16 MiB local-buffer screening estimate or failed existing memory preflight falls back.
The first real-data trigger is one correction round at sweep 200, after unchanged ordinary work.

A patch fixes every external contributor and constructs effective observations from accepted residual
plus old core prediction. It reuses `SearchProfile`, constrained A/C, bounded active-set LM and the
calibrated width bound, allowing at most 20 trusted local width updates for the joint correction.
Candidates must pass local trust/feasibility, finite checks, the existing local objective rule, full
global objective evaluation and the same local/global objective enclosure. State mutation occurs
only after acceptance. The ordinary full sweep replay, global KKT and endpoint certification then
run normally; a local improvement alone never publishes a component.

All ten controlled synthetic contracts are **Verified**:
[synthesized results](evidence/joint-hotspot-synthetic.json). Single hotspot, separated regions,
overlapping contributors, active upper bound and A=0/nonzero-charge width cases pass endpoint
certification. Near-singular local geometry safely rejects the correction (unqualified endpoint).
Injected objective rejection and width failure preserve all accepted arrays/objective; memory
fallback performs no patch solve; no-hotspot performs no correction. Independent numerical replay
consistency is checked in every case. These are patch/endpoint unit experiments, not claims that
an isolated harness published fully assembled components.

### Real-data candidate — numerically Qualified, promotion Rejected

Before correction, all 8,193 eta and 16,386 A/C values exactly match the R0 sweep-200 state.
Detection reproduces exactly the seven known serials. LEU/PRO form a 24-atom patch; the separate
ALA region forms a 20-atom patch. Largest patch=2,426 rows. Both corrections accept; none rejects,
fails recovery or falls back. Total correction round incl. detection/assembly/solving=0.628395 s.
Normal full replay and the following ordinary sweep qualify the endpoint at **201**.

| Run | Sweeps | Search seconds | Objective | A/C KKT | Projected width KKT | Endpoint / replay |
|---|---:|---:|---:|---:|---:|---|
| A: original continuous path (offline) | 250 | 487.347409 | 0.019770945620759782 | 4.26908e-8 | 9.69945e-10 | Passed |
| B: adaptive (complete production CLI) | 250 | 500.783628 | 0.019770945620759765 | 4.26908e-8 | 9.69945e-10 | Passed |
| C: same budget + correction (offline) | 201 | 398.824218 | 0.019771423894066865 | 8.07636e-8 | 4.26290e-10 | Passed |

The matched offline A/C comparison observes **18.16%** less search time in this one candidate run.
B confirms the adaptive allocation preserves A's path and parameters; its CLI timing is reported
separately because process/orchestration instrumentation differs. No matched offline B timing
repeat or production Hotspot CLI run is claimed. C wall incl. preparation/certification/conditioning
=401.124044 s, peak RSS=982,286,336 bytes versus A's 986,546,176 bytes. The local patch-buffer
formula estimates 4,110,592 bytes; it excludes global streaming evaluation and retained profile-trial
objects and therefore is **not** a proven bound on the complete correction workspace. The 16 MiB
research screen and injected fallback test must not be presented as a complete production preflight.

C passes complete endpoint certification, independent replay, accepted-state integrity and separate
production assembly replay. Recomputed existing local-conditioning/point validation gives the same
aggregate 7,766 Available / 427 Warning as B; serials 4523/4524 remain Adequate. This is a valid
stationary endpoint, not a numerical failure or a tolerance relaxation.

However C's objective is 4.78273307e-7 higher (0.002419% relative) and some point estimates differ
materially, even away from the detected hotspot. For serial 4524, original-path A/B are
12.52779165 / 0.169727533 Å versus correction A/B≈2.215e-9 / 0.659342950 Å; serial 4523 changes
from A=1.59566024, B=0.382712383 Å to A≈2.213e-9, B=0.544678905 Å. Different trajectories and
parameters are allowed; these differences do **not** invalidate the endpoint's qualification.
But the same coarse conditioning classification is not evidence that these point estimates have
equivalent statistical accuracy. No ground-truth/equivalence evidence establishes that tradeoff here.

Decision: retain the original-path adaptive solver. In the requested correctness-first scope,
one observed speed gain with a higher-objective different solution, inconclusive parameter-equivalence
and incomplete correction-workspace accounting does not justify a second production mechanism.
The study does **not** prove hotspot correction generally ineffective. Its numerical contracts are
Verified/Qualified, its production promotion is Rejected, and its broader precision/performance
reproducibility remains Inconclusive. Alternative late/persistence triggers and production downstream
Hotspot validation are Not executed after that gate decision. No update-order/partition/scheduling
experiments were executed.

Evidence: [full candidate trajectory](evidence/joint-hotspot-candidate.csv),
[candidate metrics and parameter comparison](evidence/joint-hotspot-candidate.json).

## R6 cleanup and final audit

- Removed the temporary mutable-state hotspot study hook from both declaration and call site.
  No patch solver/detector/trigger, numerical experiment switch, trace option or extra CMake target
  remains in production. Offline study programs were never installed or tracked.
- Retained the budget state machine, projected-KKT progress display, provenance and required tests.
  The sole formal numerical path is original forward FixedNeighbor + progress-qualified 200→300
  budget → unchanged certification/replay → qualified assembly → existing downstream workflow.
- `git diff` against the initial revision is empty for ConvergenceCriteria, Constrained A/C
  (`LinearSolve.cpp`), bounded LM, local Search, Numerics/independent replay, Components/assembly
  and PotentialFittingWorkflow. No tolerance/constraint/acceptance/replay/halo-domain change.
- Existing A=0, width donors, failure rollback, recovery, memory fallback, accepted/candidate integrity,
  endpoint, assembly and downstream tests all run in the complete CTest suite.
- Full Release build (all default targets) passes. Final CTest **19/19 groups pass**, 0 failed,
  37.13 s. Build parallelism=4; CTest parallelism=1, internal OpenMP limited to 4 (joint groups 1).
  No new dependency, branch, worktree or Notion write. `git diff --check` passes.
- Retained repository evidence is bounded (approximately 108 KiB before the final CLI summary),
  rather than committing large SQLite/log/state files. Detailed scratch artifacts are separate.

For performance interpretation, a conservative component-solver time-to-qualified is the measured
whole `SearchFixedNeighborComponent` call including preparation, certification and final conditioning:
R1=489.667256 s, C=401.124044 s (18.08% observed reduction). Neither includes input loading,
FirstStage, global assembly or downstream. The sweep-only times in the tables use their own explicitly
labeled narrower scope. No predicted qualification times are substituted for actual runs.

### Final cleaned production acceptance — Qualified

The final all-target Release build is followed by a fresh real 6CVM CLI run after hook removal.
CLI exit=**0**, sweep/certification=**250**, component/target/global qualification all pass, SQLite
`quick_check=ok`, SQLite schema 20 / Joint schema 22. Input SHA-256 and every fitting row identity
match the verified inputs. Final eta/A/C arrays exactly match both R1 and R2. Search=502.349581 s;
workflow time through qualified global assembly incl. construction/FirstStage/assessment=518.188925 s;
whole CLI incl. loading/Q-score/downstream/SQLite=578.104026 s; whole-process peak RSS=2,850,586,624 bytes.
Final objective=0.019770945620759765, A/C KKT=4.269076906392166e-8,
projected width KKT=9.699454295297727e-10. Endpoint preflight stays 271,418,201 bytes.
[Final machine-readable acceptance](evidence/joint-convergence-final-cli.json).

Final downstream persistence includes 1,645 reconstructed halos, 8,193 peeling entries, 8,193
post-joint local entries and 167 group summaries. Local statuses are 8,013 Available / 180 Unavailable:
176 `invalid-local-parameters`, 3 `local-mdpde-unqualified`, 1 `insufficient-positive-responses`.
These existing downstream eligibility outcomes are separate from the now-qualified Joint component;
this work does not claim universal statistical identifiability or successful inference for every group.
The 427 point-quality Warnings likewise remain visible rather than being promoted to Available.

Cleanup removed 27 owned temporary programs/source/object/duplicate-data files, including the unused
patch harness and intermediate R2 SQLite/JSON duplicates. The final SQLite and bounded raw audit
logs/states/trajectories remain under `/tmp/rhbm-convergence`; the pre-existing diagnosis directory
and the user's original inputs/databases were not modified. No temporary numerical interface remains.

Final recommendation: retain adaptive budget and the honest projected-KKT display. The original
6CVM path needs 250 sweeps under these conditions; other difficult real datasets are not guaranteed
to qualify by 300 and were not evaluated here. Nonconvex stationary solutions, statistical precision
and the downstream unavailable fits remain relevant limits. Hotspot's observed speed/solution tradeoff
is documented without replacing the lower-objective validated original path.

Git phase commits before final handoff:

- `562a05ac38ec4b8e631c2701c144873ff15657c2` — R0 baseline and projected-width progress.
- `37eced56c6ca1a85418e53088d1763030e56fc9c` — R1 extended-search evidence.
- `43f01b8c848ceb199bea33ff510aad77d129bed3` — R2 production adaptive budget and regressions.
- `3602a9b4a246e7dc58bdf7e842a1ce1da3bbd039` — R4 evidence and rejection of production promotion.

The final closeout SHA, clean-tree status and verified remote push are recorded in the handoff report
after committing this evidence. An existing detached worktree was observed and left untouched; this
work created no branch/worktree and remained on `develop` throughout.
