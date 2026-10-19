# Joint width stabilization study and production gate

> **Superseded:** Use [the corrected scoped and scale-aligned study](joint-width-trust-corrected-study.md) and its [paired results](joint-width-trust-corrected-results.json) for current findings. The runs below applied the width-trust callback outside the intended serial-1032 core when tracked evidence was absent, and compared normalized simulation estimates against original-domain A/C truth. Their 27k rejection counts and simulation A/C RMSE are retained as historical evidence of those defects, not valid strategy results.

## Decision

No width bound, objective penalty, qualification-threshold change, or reduced
production parameterization is promoted by this study. The P2 paired runs
reject both experimental width-trust policies: B1 and B2 fail the supplied
simulation and unsharpened-map searches, and both worsen the weak synthetic
truth error despite reducing the accepted width path. The current production
Pure Halo-Support Exclusion estimator is therefore retained unchanged, with
the exact-map result recorded as an unqualified diagnostic result.

## Compared evidence

The requested strategies are not all equivalent experiments. The exact
unsharpened map has no parameter truth table, so its RMSE and scientific
accuracy are explicitly unavailable. The simulation controls have truth, but
they do not by themselves calibrate a physical upper width bound. The table
keeps those scopes separate.

| Strategy | Qualification | A/B/C accuracy | Prediction accuracy | Width/conditioning evidence | Runtime / peak RSS | Gate result |
| --- | --- | --- | --- | --- | --- | --- |
| A. Current production baseline on exact unsharpened map | Unavailable; `block-recovery-exhausted` after 70 sweeps | N/A — no truth | Measured fitting-domain prediction RMSE `0.9260385` | Accepted max `B=6099.8101 Å`; condition `2.3143e8`; cancellation `2.8317e7` | 276.401 s / 1,425,555,456 B; test-instrumented paired runner | Failure state is preserved; no promotion claim |
| B1. Conditioning-aware width trust | Unavailable; `block-sweep-budget` after 200 sweeps; replay failed | N/A — no truth | Prediction RMSE `1.0253483` | Accepted max `B=52.4268 Å`; 27,994 policy rejections; serial-1032 scope | 436.369 s / 1,439,776,768 B; diagnostic runner | Rejected: new sweep-budget failure and prediction regression |
| B2. Conditioning + cancellation-aware width trust | Unavailable; `block-sweep-budget` after 200 sweeps; replay failed | N/A — no truth | Prediction RMSE `1.0253483` | Same accepted path and 27,994 rejections as B1; no additional cancellation benefit | 439.450 s / 1,440,301,056 B; diagnostic runner | Rejected: same failure and no improvement over B1 |
| C. Physically calibrated width bounds | Not run as a width-bound estimator | Truth-controlled Pure Exclusion control has 10/10 usable seeds, but no width-bound comparison | Existing noise control reports candidate RMSE A/B/C `0.01663 / 5.520e−4 Å / 0.001590`; this is not a bound study | Width-runaway incidence under a calibrated bound not measured | Not promotable; physical calibration is missing |
| D. Regularized/reduced or precision-assisted model | Full estimator qualification not measured | Not measured on real map | Local conditional prediction changes only | Exact-map 50/100-digit constrained references reach KKT below `1e−66`; fixed-design and functional-kernel references separate in the runaway tail | High-precision work is included in the diagnostic wall time; isolated production RSS not measured | Research reference only; does not establish a deployable model |

The exact-map baseline's measured process peak RSS is the diagnostic process
high-water mark, including detailed capture and high-precision profile work.
It must not be presented as a normal production runtime or memory benchmark.

## Strategy A: current production baseline

The exact-map replay is documented in
[`joint-width-runaway-serial-1032-unsharpened.md`](joint-width-runaway-serial-1032-unsharpened.md).
It retains the current objective, `A >= 0`, `B > 0`, signed `C`, Pure
Halo-Support Exclusion, 12-atom blocks, existing trust checks, and unchanged
practical qualification. It does not produce formal target estimates. This is
the correct baseline for the event: an unqualified numerical result is not a
failed command and must not be converted into a candidate estimate.

The bounded replay is one incident run, not a population estimate of runaway
frequency. Its observed width trajectory is nevertheless sufficient to reject
the claim that the event is caused only by a missing solver exception: the
accepted updates had material LM and global objective decreases, while the
conditional profile exposed an almost-flat, increasingly ill-conditioned tail.

## Strategy B: condition-aware trust

The source already contains `RejectWeakWidthTrial()`. It rejects only a
large-log-width trial whose proposed step is overwhelmingly in a weak singular
subspace and whose objective reduction is within the existing replay
enclosure. It does not apply an arbitrary condition-number cutoff.

On the exact-map trace, the guard rejected no accepted serial-1032 update:
the sweep-7 and sweep-8 accepted events had weak-step fraction `0.0`, and the
apparent descent was well above the replay enclosure. A separate rank-deficient
proposal was rejected by the rank/trust path. This is a meaningful negative
result: the current guard is not a fix for this event, but it is also not
evidence that every high-condition trial should be rejected.

The conditional profile supplies the calibration evidence still missing from
an end-to-end policy. `B=7.6478 Å` has local condition about `3.57e2` and
cancellation about `36`, while its profile objective is only about `4.2e−9`
higher than the `B=1365–42766 Å` tail. A useful future rule would need to
combine objective-reduction reliability, projected width information,
cancellation, and neighbor-state changes on multiple truth-controlled cases.
Before P2, no such predeclared end-to-end rule had been run. The P2 runner
below defines an evidence policy without a condition-number cutoff and tests
it on paired truth-controlled and supplied real-data fixtures. It remains a
test-only policy; it is not promoted.

## Strategy C: physically calibrated width bounds

The repository's truth-controlled Pure Exclusion noise study is relevant as a
model-validity control, not as a width-bound calibration. For the matching
support, sigma-0.02, ten-seed study, Pure Exclusion had 10/10 usable seeds and
candidate RMSE `0.01663 / 5.520e−4 Å / 0.001590` for A/B/C; its maximum
reported process RSS was 722,452,480 B. The study also records point-quality
and stability coverage separately, so the usable-seed result is not a claim
that every target is formally stable.

No distribution of physically valid widths was used to choose an upper bound,
and no paired run compared a calibrated bound against the same truth-controlled
maps. The support radius `2.5 Å` is an observation-support cutoff, not a width
upper bound. Promoting a bound without that paired study would change the
estimator policy and could remove valid broad components, so Strategy C fails
the promotion gate for now.

## Strategy D: reduced or precision-assisted model

The exact-map conditional reference separates two numerical questions. A
50/100-digit constrained solve on the promoted binary64 design is internally
consistent; a separate 50/100-digit solve recomputing the functional kernels
is also internally consistent. At `B=6099.810066 Å`, their profile objective
changes are `8.33395e−14` and `6.44252e−15`, while their prediction changes
are `5.90240e−9` and `1.77808e−10`. This shows that A/C cancellation amplifies
kernel-rounding differences, but it does not produce a formally qualified
full-component endpoint.

The historical Block 21 precision study reaches the same disposition on a
different simulation-map fixture: high-precision local references are useful
for diagnosis, but no safe production fallback candidate was established.
The current task therefore does not add a multiprecision fallback, an
orthogonal reparameterization, a width penalty, or a fixed-width substitution.

Related dynamic-halo experiments provide an additional model-validity warning:
their qualified research group improved some target truth RMSEs but worsened
field prediction and observation objective relative to the production control.
That tradeoff is not a width-stabilization result and is not used to justify a
production change.

## P2 paired end-to-end width-trust experiments

### Predeclared policies

The experiment executable is `build/bin/joint_width_trust_experiment`. It
uses the same 12-target forward FixedNeighbor cores, 200-sweep cap, map
normalization, initialization, objective, parameterization, local/global
qualification contracts, and independent replay for all three groups. The
policy hook is compiled only into test-instrumented targets; production code
has no new width rule.

* **A — production baseline:** no width-trust hook.
* **B1 — conditioning-aware:** for an accepted profile candidate, require
  tracked atoms, positive projected-width `sigma_min`, finite positive raw
  Gaussian and charge column norms, and matching baseline atom evidence. Set
  `amplification = max(kappa_normalized)`; projected `sigma_min` is a required
  positive evidence gate and is recorded, but is not turned into an
  uncalibrated multiplier.
* **B2 — conditioning + cancellation-aware:** use the B1 evidence and set
  `amplification = max(kappa_normalized * max(1, R_cancel))`.

For B1/B2, reject only when evidence is unavailable or when the measured
objective reduction is no larger than

```text
epsilon * amplification * max(1, abs(F_before))
```

where `epsilon` is binary64 machine epsilon. This is a numerical-evidence
policy, not an arbitrary condition or cancellation threshold. The diagnostic
record retains raw column norms, normalized condition, projected singular
values, cancellation ratio, objective reduction, and the rejection reason.
No width bound, regularizer, objective change, or qualification threshold is
introduced.

### Reproducibility and scope

The paired synthetic cases are generated by the committed runner with seeds
1001, 1002, and 1003. The well-conditioned case uses widths near `0.72 Å`
and noise sigma `5e-4`; the weak-identifiability case uses widths near `8 Å`
and noise sigma `1e-5`. Both use the same generated truth and support for A,
B1, and B2 within each seed.

The supplied simulation map and `.simulation.json` sidecar are one fixed
truth/noise realization, so its runner `seed` is `null`; no additional
simulation-noise seed artifact was invented. The unsharpened map has no A/B/C
truth. For both real fixtures, the policy evidence and intervention are
limited to serial 1032's block (`diagnostic_scope=serial-1032-target-only`)
while the search still contains all 1,418 selected targets. Consequently,
the real-data rows are evidence about the observed failure block, not a claim
that B1/B2 were applied as a full-target policy.

The external inputs are identified by the P1 hashes: `6z6u.cif`
`349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa`,
`emd_11103_unsharpened.map`
`fa9a6b15ad86d471130fc83495754cbf0a2e4084adc5327bf09f21d7e5ec60d1`, and
the supplied simulation map
`6d7771a9fd584e562550da473d624307b22d20396f11a2a79a064bff2551682a` plus
sidecar
`33a3c221ec3f2fff18dcd2a0ea6484f68334971fc688f7cd100e220ad5848931`.
All experiment builds were made with `cmake --build build --parallel 4`.
Runtime and RSS below are the
instrumented runner's process measurements, not a production benchmark.

Prediction RMSE is the candidate residual RMSE on the fitting-domain
observations. A/B/C RMSE is reported only where the runner has a truth table.
No cancellation quantiles were measured; only maxima are reported.

### Truth-controlled paired results

The rows preserve failures rather than filtering to successful runs. Values
are `A RMSE / B RMSE (Å) / C RMSE`, followed by fitting-domain prediction
RMSE, maximum trusted width, policy rejections, sweeps, and wall/RSS.

| Dataset / seed | Strategy | Qualification / replay | A/B/C RMSE | Prediction RMSE | Max trusted B (Å) | Rej. | Sweeps | Wall / RSS |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| well / 1001 | A | Passed / Passed | 0.00415027 / 0.000547510 / 0.000114488 | 0.000494384 | 0.830075 | 0 | 5 | 0.135 s / 27.9 MB |
| well / 1001 | B1 | Passed / Passed | 0.00415021 / 0.000547524 / 0.000114493 | 0.000494384 | 0.830075 | 1 | 5 | 0.064 s / 26.9 MB |
| well / 1001 | B2 | Passed / Passed | 0.00415021 / 0.000547524 / 0.000114493 | 0.000494384 | 1 | 5 | 0.065 s / 27.0 MB |
| well / 1002 | A | Passed / Passed | 0.00412459 / 0.000506650 / 0.0000726005 | 0.000489607 | 0.829668 | 0 | 5 | 0.362 s / 29.5 MB |
| well / 1002 | B1 | Passed / Passed | 0.00412474 / 0.000506661 / 0.0000726045 | 0.000489607 | 0.829668 | 1 | 5 | 0.254 s / 28.9 MB |
| well / 1002 | B2 | Passed / Passed | 0.00412474 / 0.000506661 / 0.0000726045 | 0.000489607 | 1 | 5 | 0.268 s / 29.4 MB |
| well / 1003 | A | Passed / Passed | 0.00418979 / 0.000532387 / 0.000180084 | 0.000507110 | 0.830606 | 0 | 5 | 0.404 s / 29.5 MB |
| well / 1003 | B1 | Passed / Passed | 0.00418981 / 0.000532383 / 0.000180085 | 0.000507110 | 0.830606 | 1 | 5 | 0.102 s / 29.1 MB |
| well / 1003 | B2 | Passed / Passed | 0.00418981 / 0.000532383 / 0.000180085 | 0.000507110 | 1 | 5 | 0.107 s / 29.1 MB |
| weak / 1001 | A | Unavailable / Failed | 31.7882 / 1.30925 / 0.0284405 | 0.00000991269 | 9.13313 | 0 | 200 | 4.658 s / 32.5 MB |
| weak / 1001 | B1 | Passed / Passed | 36.8420 / 3.12708 / 0.101431 | 0.0000101584 | 7.96090 | 1 | 97 | 1.051 s / 28.7 MB |
| weak / 1001 | B2 | Passed / Passed | 36.8688 / 3.12812 / 0.101478 | 0.0000101589 | 7.86434 | 1 | 96 | 1.015 s / 28.7 MB |
| weak / 1002 | A | Unavailable / Failed | 25.4488 / 1.08688 / 0.0255257 | 0.00000981255 | 9.09858 | 0 | 200 | 3.064 s / 31.4 MB |
| weak / 1002 | B1 | Passed / Passed | 41.6210 / 2.14831 / 0.0564232 | 0.00000984895 | 8.55727 | 1 | 101 | 1.453 s / 31.2 MB |
| weak / 1002 | B2 | Passed / Passed | 41.6210 / 2.14831 / 0.0564232 | 0.00000984895 | 8.55727 | 1 | 101 | 1.457 s / 31.1 MB |
| weak / 1003 | A | Passed / Passed | 37.8176 / 1.86699 / 0.0450866 | 0.0000101364 | 9.13306 | 0 | 150 | 2.037 s / 31.3 MB |
| weak / 1003 | B1 | Passed / Passed | 39.8809 / 2.58453 / 0.0737142 | 0.0000101664 | 8.88774 | 1 | 98 | 1.218 s / 31.0 MB |
| weak / 1003 | B2 | Passed / Passed | 41.4116 / 2.61291 / 0.0745124 | 0.0000101687 | 8.28098 | 1 | 95 | 1.179 s / 30.9 MB |

The well-conditioned cases are preserved to the displayed precision, but the
policy changes a numerically accepted update in every B1/B2 run. In the weak
cases, B1/B2 reduce the accepted maximum width and often reach a stationary
stop sooner, but they worsen all three parameter RMSEs relative to A and
slightly worsen prediction RMSE. Thus numerical convergence and parameter
identifiability are not interchangeable here.

### Supplied 6Z6U simulation and unsharpened map

| Dataset / strategy | Qualification / replay | A/B/C RMSE | Prediction RMSE | Max trusted B / condition / cancellation | Rej. / accepted updates | Sweeps / evaluations | Wall / RSS |
| --- | --- | ---: | ---: | --- | ---: | ---: | ---: |
| simulation / A | Passed / Passed | 8.48009 / 1.46484e−7 / 0.626977 | 9.20736e−7 | 0.500000 / 3.32726 / 1.13168 | 0 / 11,454 | 83 / 22,908 | 306.151 s / 581,844,992 B |
| simulation / B1 | Unavailable / Failed | 8.85286 / 0.0388439 / 0.608313 | 0.129203 | 0.501900 / 3.32876 / 1.13198 | 27,559 / 41 | 200 / 55,200 | 432.762 s / 602,308,608 B |
| simulation / B2 | Unavailable / Failed | 8.85286 / 0.0388439 / 0.608313 | 0.129203 | 0.501900 / 3.32876 / 1.13198 | 27,559 / 41 | 200 / 55,200 | 434.178 s / 556,728,320 B |
| unsharpened / A | Unavailable / Failed | N/A — no truth | 0.926038 | 6099.810 / 2.31429e8 / 2.83166e7 | 0 / 9,760 | 70 / 19,705 | 276.401 s / 1,425,555,456 B |
| unsharpened / B1 | Unavailable / Failed | N/A — no truth | 1.02535 | 52.4268 / 1.71087e4 / 1.71530e3 | 27,994 / 6 | 200 / 56,776 | 436.369 s / 1,439,776,768 B |
| unsharpened / B2 | Unavailable / Failed | N/A — no truth | 1.02535 | 52.4268 / 1.71087e4 / 1.71530e3 | 27,994 / 6 | 200 / 56,776 | 439.450 s / 1,440,301,056 B |

On the supplied unsharpened fixture, the serial-1032 accepted width is much
smaller under B1/B2, but both policies create a new `block-sweep-budget`
failure and fail independent replay. They also regress fitting-domain
prediction RMSE from `0.926038` to `1.02535`. On the simulation fixture,
baseline is qualified with near-truth prediction while both policies fail the
search and replay, with materially worse parameter and prediction accuracy.
B2 is numerically identical to B1 on both supplied real fixtures; the
cancellation term provides no additional improvement.

### P2 answers and decision

1. **Width runaway:** B1/B2 reduce the trusted serial-1032 endpoint width in
   the unsharpened paired run, but the run is not qualified and the policy was
   not applied to every target. This is not sufficient evidence of a safe
   stabilization.
2. **Cancellation-aware benefit:** none beyond B1 in the measured fixtures.
3. **Qualification:** improved on weak synthetic seeds, but failed on the
   truth-controlled simulation and real unsharpened map; no general gain.
4. **Truth accuracy:** well cases are preserved, while weak synthetic and
   simulation parameter accuracy regress under B1/B2.
5. **Rejected valuable updates:** yes. The policies reject updates that
   produce lower parameter accuracy on weak synthetic cases and lower
   simulation prediction quality.
6. **New recovery failures:** yes — both B1 and B2 reach the sweep budget and
   fail replay on the supplied simulation and unsharpened fixtures.
7. **Cost:** real-data B1/B2 runs are roughly 1.4–1.5× the baseline wall
   time; RSS is of the same order but not lower in a way that compensates for
   the failure. Synthetic timing is not a meaningful production estimate.
8. **Best-supported strategy:** A, the unchanged production baseline. B1/B2
   remain research-only diagnostics.

## Regression and normal-case controls

The existing normal synthetic FixedNeighbor fixtures remain passing, including
the bounded width-capture test and the weak-width safeguard contract. The
near-degenerate Block 21 precision regression remains intentionally
unqualified under the unchanged KKT gate while its independent 50/100-digit
references qualify. These controls preserve the distinction between solver
diagnostics, scientific accuracy, and formal qualification.

## Promotion-gate conclusion

| Gate | Result |
| --- | --- |
| Improve or preserve truth accuracy across paired cases | Rejected: B1/B2 regress weak-synthetic and simulation parameter/prediction accuracy |
| Preserve good normal cases | Yes on three well-conditioned seeds to displayed precision |
| Keep qualification standards unchanged | Yes |
| Preserve model constraints and independent replay | Rejected: both policies fail replay on both supplied real fixtures |
| Acceptable runtime/memory tradeoff | Rejected: higher real-data cost with no qualification benefit |
| Production promotion | **Rejected** |

The next scientifically useful experiment is a broader, independently
calibrated policy study over multiple simulation noise seeds and full-target
evidence. Until that exists, the production behavior remains unchanged and
the exact-map failure remains a saved, unqualified diagnostic result.
