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
