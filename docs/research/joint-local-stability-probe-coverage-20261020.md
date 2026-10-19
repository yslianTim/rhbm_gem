# Bounded local stability probe coverage — 2026-10-20

## Decision

The representative probe policy now prioritizes Width-bound-active blocks,
blocks near the practical KKT limits, poor conditioning, and active A/C faces.
It samples ordinary well-conditioned blocks across core-atom positions. The
existing limit of eight perturbation blocks per fit and the per-block CPU,
evaluation, and memory preflights remain in force. One global slot is reserved
for a later Width-bound-active component so earlier ordinary components cannot
consume the entire fit budget.

The bounded optimizer, accepted endpoint, global qualification, and independent
replay are unchanged. The 6Z6U reruns retained the same reported objectives,
KKT values, sweep counts, and profile counts as the P1 reruns. The synthetic
fixture selected its one bound-active component. Unavailable and NotAssessed
results remain explicit; neither is counted as Stable.

## Probe coverage and cost

| Workload | Stage | Blocks | Selected | Assessed Stable / Sensitive / Unstable | Unavailable | NotAssessed | Perturbation solves / evaluations | Diagnostic CPU / wall (s) | Estimated workspace peak (B) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 6Z6U simulation | P1 | 138 | 1 | 1 / 0 / 0 | 0 | 137 | 6 / 24 | 0.046322 / not emitted | not emitted |
| 6Z6U simulation | P2 | 138 | 5 | 2 / 0 / 3 | 0 | 133 | 30 / 120 | 0.314277 / 0.315784 | 106,629,766 |
| 6Z6U unsharpened | P1 | 140 | 1 | 0 / 0 / 0 | 1 | 139 | 3 / 31 | 0.080683 / not emitted | not emitted |
| 6Z6U unsharpened | P2 | 140 | 6 | 0 / 0 / 1 | 5 | 134 | 19 / 159 | 0.365126 / 0.366709 | 110,724,617 |
| Bound-active synthetic | Before cross-component reservation | 25 | 8 | 8 / 0 / 0 | 0 | 17 | 16 / 78 | 0.011524 / 0.011766 | 3,219,315 |
| Bound-active synthetic | P2 | 25 | 8 | 8 / 0 / 0 | 0 | 17 | 15 / 71 | 0.010738 / 0.010981 | 3,219,315 |

The P1 process peak RSS was 618,905,600 B (simulation) and 1,434,451,968 B
(unsharpened). The P2 reruns measured 708,755,456 B and 1,389,789,184 B,
respectively. These whole-process RSS values include the complete fit and are
not estimates of local-diagnostic-only allocations. The bounded local
workspace preflight estimates are reported separately above.

The synthetic fixture has 25 disconnected components, and its sole bound-active
atom is in the last component. Before the global reservation, the first eight
ordinary component probes consumed the fit budget. P2 retained the eight-probe
cap and selected the later bound-active component while still assessing eight
blocks. Its bound-active block was Stable; bound feasibility, active-set KKT,
global practical KKT, RuntimeConvergence, and independent replay passed. The
truth-width outside the global cap remains a known source of parameter bias in
this fixture.

## 6Z6U numerical comparison

| Metric | Simulation P1 → P2 | Unsharpened P1 → P2 |
|---|---:|---:|
| Targets / observations | 1,418 / 180,478 | 1,418 / 189,377 |
| Bmax (Å) | 2.2759153354042887 | 3.2084356818726385 |
| Search sweeps / profile evaluations | 83 / 22,908 → 83 / 22,908 | 124 / 34,848 → 124 / 34,848 |
| Objective | 1.2330854370810358e-13 → 1.2330854370810358e-13 | 0.018417539368903803 → 0.018417539368903803 |
| Global A/C KKT | 9.434574186999167e-8 → unchanged | 8.174157334426418e-8 → unchanged |
| Projected Width KKT | 8.776510069988752e-10 → unchanged | 9.402726930396454e-10 → unchanged |
| Active-bound targets | 0 → 0 | 2 → 2 |
| Active-set subproblem failures | 0 → 0 | 0 → 0 |
| RuntimeConvergence / Independent Replay | Passed / Passed → unchanged | Passed / Passed → unchanged |
| Peak RSS (B) | 618,905,600 → 708,755,456 | 1,434,451,968 → 1,389,789,184 |
| Wall time (s) | 303.731 → 303.694 | 435.112 → 435.770 |

On the unsharpened map, selected blocks `1/block-5` and `1/block-132` each
contain one bound-active target. Their baseline A/C KKT values are
`3.5208046721135755e-10` and `7.528810352930293e-10`, respectively. These are
above the strict local `1e-10` certificate but below the unchanged practical
`1e-7` A/C criterion. Their projected Width KKT values are
`5.4610760358286825e-12` and `8.681722007963799e-12`, both below the practical
`1e-9` limit. Both baselines are eligible and independently replayed. Their
perturbation endpoints still fail the strict independent replay gate, so the
reports remain Unavailable rather than Stable. One other selected block is
Unstable based on an independently verified direction whose prediction change
exceeded the existing warning limit.

## Implementation and validation

The bounded context is now passed into block conditioning and stability
assessment, so local diagnostics can report active-bound counts and projected
Width KKT using the same upper bound as the accepted global endpoint. Probe
selection ranks those measured local values, skips free-A/C rank-deficient
blocks unless they contain an active Width bound, and samples well-conditioned
blocks over spatial positions. The global budget reservation prevents early
ordinary components from starving a later bound-active component.

Persistence adds optional local KKT and active-bound metrics plus an explicit
`probe_selected` flag. Readers accept older records where these keys are
absent. No persistence schema migration, public CLI option, warning threshold,
solver certificate, global threshold, or production default changed.

Validation used the Release non-instrumented candidate runner with internal
OpenMP/BLAS threading set to one:

- Direct local-conditioning, parameter-stability, and persistence filter: 23/23 passed.
- Focused CTest group covering numerical, contract, workflow, regression,
  bounded synthetic, persistence, and runtime runner: 9/9 passed.
- Bound-active synthetic run: 25 blocks, eight selected and assessed, including
  the bound-active block; 15 perturbation solves and 71 profile evaluations.
- 6Z6U simulation and unsharpened bounded reruns both retained
  `RuntimeConvergence=Passed`, bound feasibility, projected Width KKT, and
  independent replay.

## Remaining qualification limits

The two maps show increased probe representation under fixed resource budgets,
not exhaustive per-block stability coverage. The simulation reported three
Unstable blocks and the unsharpened map reported five Unavailable blocks; these
diagnostic outcomes do not fail the separately qualified Joint endpoint. The
whole-process RSS changed between reruns, so those differences are not
attributed solely to local diagnostics. The 37,406-target workload and full
bounded CLI/persistence/downstream workflow remain separate promotion gates.
