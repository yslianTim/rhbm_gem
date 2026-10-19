# Serial 1032 width runaway on the unsharpened 6Z6U map

## Scope and provenance

This report is the P1 bounded replay for Chain A, GLU 65, CD, serial 1032.
It uses the corrected incident map `emd_11103_unsharpened.map`, not the
simulation map used by the older FullABC and Block 21 studies. The input
artifacts are external to the repository and are identified here by SHA-256:

| Input | SHA-256 |
| --- | --- |
| `6z6u.cif` | `349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa` |
| `emd_11103_unsharpened.map` | `fa9a6b15ad86d471130fc83495754cbf0a2e4084adc5327bf09f21d7e5ec60d1` |

The replay used the current `develop` implementation, Pure Halo-Support
Exclusion, 12-atom forward FixedNeighbor blocks, a 200-sweep cap, one
compute thread, and a 2,048-event bounded capture:

```sh
build/bin/joint_width_trajectory \
  "$DATA_DIR/6z6u.cif" \
  "$DATA_DIR/emd_11103_unsharpened.map" \
  /private/tmp/joint-unsharpened-width-1032-profile-final \
  1032 2048 pure
```

The tool writes only bounded JSONL/JSON diagnostics to the output directory;
the CIF/MAP and a full database are not repository artifacts. The run used
1,418 targets, 241,178 original observations, 51,801 excluded observations,
and 189,377 fitting observations, matching the incident domain counts.

## Historical incident versus current replay

The replay is reproducible on the corrected map, but it is not an exact
reproduction of every historical number. The original report reached 70
completed sweeps and recorded `B=6099.810066 Å` at sweep 8. The current
`develop` replay stops after 13 completed sweeps with
`block-recovery-exhausted`; its accepted sweep-8 value for serial 1032 is
`B=42766.18460821246 Å`. Therefore the two traces must not be silently
merged. The exact map does reproduce the same qualitative transition:
ordinary width → multi-decade width jump → severe A/C cancellation.

| Quantity | Historical incident | Current exact-map replay |
| --- | ---: | ---: |
| Completed sweeps | 70 | 13 |
| Sweep 6 width (Å) | 7.647783566 | 7.647783566434263 |
| Sweep 7 width (Å) | 1365.304216 | 1365.3042156820059 |
| Sweep 8 width (Å) | 6099.810066 | 42766.18460821246 |
| Stop reason | `block-recovery-exhausted` | `block-recovery-exhausted` |
| Formal runtime convergence | Untrusted | Unavailable |

The difference is an unresolved provenance/configuration discrepancy, not
evidence that the historical trace was incorrect. All conclusions below label
which measurement they use.

## Accepted sweep 6 → 7 → 8 updates

The following rows are the current exact-map capture for the serial-1032
block. Local objective changes are normalized profile-objective changes;
global changes are the original full-component objective changes. LM actual and
predicted reductions are retained separately because they use the LM residual
convention.

| Sweep | Accepted `B` (Å) | Δlog B | LM ratio | Local ΔF | Global ΔF | A | C | Local condition | Cancellation | Local KKT |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 6 | 1.8277316809 → 7.6477835664 | 1.4313401966 | 0.8008009 | 7.6838e−8 | 1.8636e−6 | 3.7532e2 | −4.7385e−1 | 3.5692e2 | 3.6226e1 | 2.78e−17 |
| 7 | 7.6477835664 → 1365.3042156820 | 5.1847166754 | 0.4282655 | 2.3820e−8 | 1.1406e−6 | 7.4799e13 | −3.1932e6 | 1.1527e7 | 1.2864e6 | 9.09e−13 |
| 8 | 1365.3042156820 → 42766.1846082125 | 3.4443704386 | 0.6906870 | 5.3248e−8 | 7.5540e−7 | 2.4095e21 | −1.0484e11 | 1.1376e10 | 1.3033e9 | 8.83e−11 |

The accepted LM actual/predicted reductions were respectively
`2.6253160529e−4 / 3.2783629426e−4`,
`8.2052852045e−5 / 1.9159342080e−4`, and
`1.8424522786e−4 / 2.6675646365e−4`. Thus each accepted update had a
finite, positive local/global descent signal. The sweep-7 profile also tried
`B=9.2361928202e9 Å`; that proposal was rejected as `rank-deficient` with
local rank 23 rather than 24 and local KKT `3.5474e−2`.

The fitted basis correlation moved from `0.9999532350` at sweep 6 to
`0.999999999999955` at sweep 7 and `1.000000000000001` at sweep 8. The
corresponding normalized A/C terms were approximately `(1.0863,-1.0288)`,
`(39232.9,-39232.9)`, and `(4.1121e7,-4.1121e7)`. Cancellation therefore
began at the first accepted multi-decade update and grew by later updates; it
was not created only by the eventual recovery failure.

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
`B=42766.18460821246 Å`, with affected-row count 3,132 and normalized local
profile objective `2.887962134625143e−4`.

| B (Å) | Profile objective | Δ from sweep state | Local KKT | Target width gradient | Condition | Cancellation | Prediction Δnorm |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5 | 2.8949694161e−4 | +7.0073e−7 | 3.47e−17 | 1.6746e−6 | 5.41 | 1.00 | 1.3100e−2 |
| 7.647783566 | 2.8863307532e−4 | −1.6314e−7 | 3.30e−17 | 8.2845e−9 | 3.60e2 | 4.02e1 | 6.3585e−3 |
| 1365.304216 | 2.8862891206e−4 | −1.6730e−7 | 3.18e−12 | 2.4836e−8 | 1.159e7 | 1.357e6 | 6.3606e−3 |
| 6099.810066 | 2.8862891193e−4 | −1.6730e−7 | 1.46e−11 | 3.0756e−6 | 2.314e8 | 2.708e7 | 6.3606e−3 |
| 42766.184608 | 2.8862891185e−4 | −1.6730e−7 | 2.98e−11 | 2.7293e−4 | 1.138e10 | 1.331e9 | 6.3606e−3 |

This is a conditional slice, not a global alternative solution. It does show
that a finite, much better-conditioned width (`7.6478 Å`) has an almost equal
objective to the runaway tail. The extreme tail gains only about 4.2e−9 in
this local normalized objective relative to that finite point while its
condition and cancellation grow by roughly seven orders of magnitude. The
`0.5 Å` point is well-conditioned but has a visibly larger profile loss.

External-neighbor state changes are material: 1,406 atoms changed from sweep
6 to 7 and again from sweep 7 to 8. Their maximum external Δeta values were
`0.4909124782` and `0.2934143935`; scaled coefficient changes were
`0.1550909840` and `0.1894672826`. The captured projected-width singular
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
| 0.5 Å | `8.68030e−7 / 1.15745e−2` | `8.68030e−7 / 1.15745e−2` |
| 6099.810066 Å | `8.33395e−14 / 5.90240e−9` | `6.44252e−15 / 1.77808e−10` |

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
  stream bounded width records without retaining every block/sweep record; a
  regression test verifies that this does not alter the search result or
  retain the large diagnostic vectors.

The replay itself remains numerically unqualified. Its final candidate had
global projected A/C KKT `5.7641e−4`, global width-gradient infinity norm
`2.1341e2`, and endpoint/independent replay status `NotRun`; those values are
diagnostics, not formal target estimates.
