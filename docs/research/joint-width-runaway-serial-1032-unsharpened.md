# Serial 1032 width runaway on the unsharpened 6Z6U map

## Scope and provenance

This report is the P1 bounded replay and reproducibility audit for Chain A,
GLU 65, CD, serial 1032. It uses the incident map
`emd_11103_unsharpened.map`, not the simulation map used by the older FullABC
and Block 21 studies. The input artifacts are external to the repository and
are identified here by SHA-256:

| Input | SHA-256 |
| --- | --- |
| `6z6u.cif` | `349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa` |
| `emd_11103_unsharpened.map` | `fa9a6b15ad86d471130fc83495754cbf0a2e4084adc5327bf09f21d7e5ec60d1` |

The production CLI and both replay executables used the current `develop`
implementation, Pure Halo-Support Exclusion, 12-atom forward FixedNeighbor
blocks, a 200-sweep cap, one compute thread, and the same model/map selection.
The instrumented replay additionally used a 2,048-event bounded capture:

```sh
build/bin/joint_width_trajectory \
  "$DATA_DIR/6z6u.cif" \
  "$DATA_DIR/emd_11103_unsharpened.map" \
  /private/tmp/joint-unsharpened-width-1032-profile-final \
  1032 2048 pure
```

The tool writes only bounded JSONL/JSON diagnostics to the output directory;
the CIF/MAP and a full database are not repository artifacts. All three runs
used 1,418 targets, 241,178 original observations, 51,801 excluded
observations, and 189,377 fitting observations. The map normalization divisor
was exactly `0.005236894964862669`.

The relevant build fingerprint was:

| Field | Value |
| --- | --- |
| Source SHA-256 | `7384f21fd92765a39b0bf44ba09256004cf654892cf6343d65c4bbf6a0b7236d` |
| Configuration SHA-256 | `737f816d2f3d78888e1b72cf27da2b9b2830337852c1a82766a6207875b0d597` |
| Build SHA-256 | `e01b538803981b638b3261f43106a346513848f3622dabc1bffa00d3d209eae2` |
| Input material SHA-256 | `661fb5ad31212153112659084c18a289e03b927350491cbbf6434a08a310d659` |
| Problem material SHA-256 | `c4ee1f3c3474f0f27463f671d72e9cdaf0c7e6bf2395a3394bdde792c2082398` |
| Initialization material SHA-256 | `1e00b58163933a21462a27a4a06111c9261d401ebb77d25860286595d230b5a6` |
| Partition material SHA-256 | `3f7446b8c5e6aa217507b04dba847f5729745db0f14a38b3af03e5220a153450` |

## Production, instrumented, and uninstrumented results

The production CLI persisted the same failure state that the diagnostic replay
observed: `search_completed=false`, runtime convergence `Unavailable`, stop
reason `block-recovery-exhausted`, 70 sweeps, 19,705 profile evaluations, and
9,760 accepted updates. Its persisted candidate objective was
`0.018417514708820186`, candidate A/C KKT
`5.2758671483044486e-6`, and candidate width-gradient infinity norm
`4.4276910662124235e-5`.

| Run | Sweeps / profile evaluations / accepted updates | Max width (Å) | Max local A/C condition | Max cancellation | Wall time (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Production CLI, persisted SQLite | 70 / 19,705 / 9,760 | not captured | not captured | not captured | search cost 271.627 |
| Instrumented replay | 70 / 19,705 / 9,760 | 6099.810066010588 | 2.3142908215317008e8 | 2.8333918707490444e7 | 414.194 |
| No-instrumentation replay | 70 / 19,705 / 9,760 | not captured | not captured | not captured | 281.969 |

The instrumented replay produced the bounded-output hashes
`sweeps.jsonl = 2ae9b37d93cf72a260b6c1e319e51cbdc2991d27dbc4f1dcdda14ee08c8ff3c4`
and `width-trajectory.jsonl =
e81adf32aa06e3e780b4cc23a77c543ab6eec0521649b188132fcb8e31140f71`.
The instrumented and no-instrumentation runs had identical input, problem,
initialization, and partition material hashes. Their policy fingerprints
correctly differed only because one build role enabled test instrumentation.
The no-instrumentation executable does not capture the trajectory, so the
evidence claim is aggregate/state parity rather than a claim of an unobserved
bit-for-bit JSONL trace.

## Resolving the historical 13-sweep discrepancy

An earlier diagnostic executable was run before P1 added the production
preprocessing lifecycle. That run omitted `MapValueArrayNormalization()` and
therefore followed a different numerical problem. It reached 13 sweeps with
3,654 profile evaluations and 1,822 accepted updates; its sweep-8 width was
`42766.18467850702 Å`, condition `1.1375911952e10`, and cancellation
`1.3792245219e9`.

Adding only `InitializeFromSelection()` to that old executable left the 13
sweep trace unchanged. Adding the missing map normalization, with the exact
divisor recorded above, changed the replay to the production 70-sweep trace.
Thus the discrepancy is resolved as a replay-input/lifecycle mismatch, not as
an instrumentation-induced optimizer change. No production objective,
qualification threshold, width bound, or estimator parameterization was
changed.

| Quantity | Old unnormalized diagnostic | Corrected production-equivalent replay |
| --- | ---: | ---: |
| Completed sweeps | 13 | 70 |
| Profile evaluations | 3,654 | 19,705 |
| Accepted updates | 1,822 | 9,760 |
| Sweep 6 width (Å) | 7.647783566 | 7.647783566434263 |
| Sweep 7 width (Å) | 1365.304216 | 1365.3042156820059 |
| Sweep 8 width (Å) | 42766.18467850702 | 6099.810066010588 |
| Stop reason | `block-recovery-exhausted` | `block-recovery-exhausted` |
| Formal runtime convergence | Unavailable | Unavailable |


## Accepted sweep 6 → 7 → 8 updates

The following rows are the corrected production-equivalent capture for the
serial-1032 block. Local objective changes are normalized profile-objective
changes; global changes are the original full-component objective changes. LM
actual and predicted reductions are retained separately because they use the
LM residual convention.

| Sweep | Accepted `B` (Å) | Δlog B | LM ratio | Local ΔF | Global ΔF | A | C | Local condition | Cancellation | Local KKT |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 6 | 1.8277316809 → 7.6477835664 | 1.4313401966 | 0.8008009 | 7.6838e−8 | 1.8636e−6 | 7.1668e4 | −9.0482e1 | 3.5692e2 | 3.6226e1 | 1.56e−17 |
| 7 | 7.6477835664 → 1365.3042160 | 5.1847166754 | 0.4282655 | 2.3820e−8 | 1.1406e−6 | 1.4283e16 | −6.0975e8 | 1.1527e7 | 1.2864e6 | 1.36e−12 |
| 8 | 1365.3042160 → 6099.8100660 | 1.4968803613 | 0.6907077 | 5.3249e−8 | 7.5540e−7 | 2.7160e19 | −5.8088e10 | 2.3143e8 | 2.6515e7 | 1.46e−11 |

The accepted LM actual/predicted reductions were respectively
`2.6253160529e−4 / 3.2783629426e−4`,
`8.2052853808e−5 / 1.9159342080e−4`, and
`1.8424957192e−4 / 2.6675477734e−4`. Thus each accepted update had a
finite, positive local/global descent signal. The sweep-7 profile also tried
`B=9.2361928202e9 Å`; that proposal was rejected as `rank-deficient` with
local rank 23 rather than 24 and local KKT `3.5474e−2`.

The fitted basis correlation moved from `0.9999532350` at sweep 6 to
`0.999999999999955` at sweep 7 and `0.999999999999999` at sweep 8. The
corresponding normalized A/C terms were approximately
`(207.42,-196.45)`, `(7.4916e6,-7.4916e6)`, and
`(1.5974e8,-1.5974e8)`. Cancellation therefore began at the first accepted
multi-decade update and grew by later updates; it was not created only by the
eventual recovery failure.

## Why the current trust checks accepted the updates

`RejectWeakWidthTrial()` is deliberately narrow. It requires all of the
following before rejecting a trial: log-width step larger than `ln(10)`, at
least 95% of the step in the weak right-singular subspace, and an objective
improvement no larger than the replay enclosure
`1e-12 + 2e-12 |F_before|`. It does not reject a direction merely because its
condition indicator is large.

For the current trace, the sweep-6 jump is below `ln(10)`. The sweep-7 and
sweep-8 accepted trials recorded `weak_width_step_fraction=0.0` and
`weak_width_guard_rejected=false`; their objective reductions were also far
above the replay enclosure. The large sweep-7 proposal was stopped by the
separate rank/trust path. Consequently the existing safeguard did not fire
for the accepted runaway steps for a specific, observable reason; changing it
to an arbitrary condition-number veto would conflate weak identifiability with
unreliable objective descent.

## Conditional width profile

After sweep 8, the serial-1032 block was profiled with external neighbors
fixed. For each width, constrained A/C was re-solved. The base state was
`B=6099.810066010588 Å`, with affected-row count 3,132 and normalized local
profile objective `2.8879621233752627e−4`; the exact base value is retained in
the JSON artifact because the reported state is not an inner profile optimum.

| B (Å) | Profile objective | Δ from sweep state | Local KKT | Target width gradient | Condition | Cancellation | Prediction Δnorm |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5 | 2.8949694001e−4 | +7.0073e−7 | 2.78e−17 | 1.6746e−6 | 5.41 | 1.00 | 2.5015 |
| 7.647783566 | 2.8863307410e−4 | −1.6314e−7 | 4.16e−17 | 8.2845e−9 | 3.60e2 | 4.02e1 | 1.2142 |
| 1365.304216 | 2.8862891084e−4 | −1.6730e−7 | 4.50e−13 | 2.1348e−10 | 1.159e7 | 1.357e6 | 1.2146 |
| 6099.810066 | 2.8862891072e−4 | −1.6730e−7 | 1.70e−11 | 1.7340e−6 | 2.314e8 | 2.7080e7 | 1.2146 |
| 8192 | 2.8862891070e−4 | −1.6730e−7 | 2.91e−11 | 6.5181e−6 | 4.174e8 | 4.8842e7 | 1.2146 |

This is a conditional slice, not a global alternative solution. It does show
that a finite, much better-conditioned width (`7.6478 Å`) has an almost equal
objective to the runaway tail. The extreme tail gains only about 4.2e−9 in
this local normalized objective relative to that finite point while its
condition and cancellation grow by roughly seven orders of magnitude. The
`0.5 Å` point is well-conditioned but has a visibly larger profile loss.

External-neighbor state changes are material: 1,406 atoms changed from sweep
6 to 7 and again from sweep 7 to 8. Their maximum external Δeta values were
`0.4909124782` and `0.2934143935`; scaled coefficient changes were
`0.3173491701` and `0.2900963209`. The captured projected-width singular
values at sweep 8 were
`1.05682, 1.04691, 1.02760, 1.02218, 1.01575, 1.00280, 1.00170,
0.98803, 0.98296, 0.95709, 0.94862, 0.94193`.

## Precision separation

The exact-map profile now records two independent references at 50 and 100
decimal digits:

1. `fixed_binary64_design`: the binary64 design entries produced by the
   production kernel are promoted exactly and re-solved with constrained
   Householder QR;
2. `functional_kernel`: squared distances and kernels are recomputed directly
   at 50/100-digit precision and then solved with the same constrained method.

The references agree between 50 and 100 digits to the displayed precision and
reach KKT below `1e−66` at the reported points. Representative profile
changes from the sweep-8 base are:

| Width | Fixed-design ΔF / prediction Δ | Functional-kernel ΔF / prediction Δ |
| ---: | ---: | ---: |
| 0.5 Å | `8.6802928548e−7 / 2.2101818965` | `8.6802928929e−7 / 2.2101818931` |
| 6099.810066 Å | `−2.0463960618e−15 / 2.8413193970e−8` | `2.2826370905e−26 / 1.2029796528e−19` |

At the moderate width the two references are effectively identical. In the
runaway tail, both references still solve the constrained A/C problem, but
the fixed-design and functional-kernel prediction changes separate by orders
of magnitude because the opposing terms are huge. This is evidence for
precision sensitivity under cancellation, not evidence that SparseQR alone is
the root cause. The conditional objective remains flat even after the A/C
subproblem is solved accurately, so practical non-identifiability and the
block-coordinate search path are both implicated.

The fixed-design/functional split is specific to this exact-map replay. The
older Block 21 report remains useful as a separate historical study, but its
simulation-map fixture is not substituted for this event.

## Conclusions and limits

* The first current exact-map runaway update is sweep 6, where `B` jumps from
  `1.8277` to `7.6478 Å`; the first severe multi-decade update is sweep 7.
* The accepted updates pass existing trust checks because they have positive,
  material LM and objective reductions, while the weak-width replay-scale
  guard conditions are not met. This is a policy limitation, not a failure to
  preserve the numerical failure state.
* A/C cancellation is already visible at the accepted sweep-6 endpoint and
  becomes extreme after the sweep-7 and sweep-8 updates. Neighbor changes
  alter the conditional problem substantially, so the later cancellation
  cannot be attributed to the width proposal in isolation.
* The profile contains finite widths with essentially the same local objective
  and much better conditioning. The evidence supports weak practical
  identifiability plus search-path sensitivity; it does not establish a unique
  finite global minimizer.
* The exact-map high-precision references do not show a simple solver-precision
  explanation. They do show that binary64 kernel rounding is amplified in the
  runaway tail. SparseQR is therefore not established as the sole primary
  problem.
* No production objective, width bound, qualification threshold, or estimator
  parameterization was changed in P1. The diagnostic path was changed to
  normalize the map and initialize model analysis through the same lifecycle as
  production, then stream bounded width records without retaining every
  block/sweep record. The paired replay and fingerprints verify that these
  diagnostics do not alter the observed search result.

The replay itself remains numerically unqualified. Its final candidate had
candidate A/C KKT `5.2758671483e−6`, candidate width-gradient infinity norm
`4.4276910662e−5`, and endpoint/independent replay status `NotRun`; those
values are diagnostics, not formal target estimates. The stop reason remains
`block-recovery-exhausted`, so no formal Joint target estimate is promoted.
