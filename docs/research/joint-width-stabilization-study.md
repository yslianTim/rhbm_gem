# Joint width stabilization study and production gate

## Decision

No width bound, objective penalty, qualification-threshold change, or reduced
production parameterization is promoted by this study. The current evidence
supports retaining the existing Pure Halo-Support Exclusion estimator and
recording the unqualified exact-map result. The P1 conditional profile shows
that a finite width can be nearly as good locally as the runaway tail, but it
does not define a globally valid replacement search rule.

## Compared evidence

The requested strategies are not all equivalent experiments. The exact
unsharpened map has no parameter truth table, so its RMSE and scientific
accuracy are explicitly unavailable. The simulation controls have truth, but
they do not by themselves calibrate a physical upper width bound. The table
keeps those scopes separate.

| Strategy | Qualification | A/B/C accuracy | Prediction accuracy | Width/conditioning evidence | Runtime / peak RSS | Gate result |
| --- | --- | --- | --- | --- | --- | --- |
| A. Current production baseline on exact unsharpened map | Unavailable; `block-recovery-exhausted` after 13 sweeps | Not measured; no truth | Not measured against truth | Accepted sweep-8 `B=42766.18 Å`; max local condition `1.1376e10`; max cancellation `1.3792e9` | 180.063 s / 1,600,421,888 B; diagnostic run, not production benchmark | Failure state is preserved; no promotion claim |
| B. Existing condition-aware weak-width guard plus conditional profile | Same full-run status; separate alternate end-to-end policy not run | Not measured on real map | Conditional profile only | Guard rejections `0`; finite `B=7.6478 Å` has condition `3.57e2`, cancellation `36.2`, and local objective within `4.2e−9` of the runaway tail | Same run as A; alternate runtime/RSS not measured | Evidence for further research, not a calibrated threshold |
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
No such predeclared rule was run here; therefore no condition threshold was
promoted.

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
| Improve or preserve truth accuracy across paired cases | Not established for a width policy |
| Preserve good normal cases | Existing tests pass; new policy not run |
| Keep qualification standards unchanged | Yes |
| Preserve model constraints and independent replay | Yes; no production numerical changes |
| Acceptable runtime/memory tradeoff | Not measured for an alternate policy |
| Production promotion | No |

The next scientifically useful experiment is a predeclared condition-aware
policy sweep over multiple truth-controlled maps, with a paired unbounded
baseline and explicit width/cancellation/stability outcomes. Until that exists,
the production behavior should remain unchanged and the exact-map failure
should remain a saved, unqualified diagnostic result.
