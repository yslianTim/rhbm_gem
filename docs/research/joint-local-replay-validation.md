# Joint local replay validation

## Scope and reproduction

This study isolates the rejected Local Profile LM trial captured from the 6Z6U
workflow. It does not rerun the global endpoint verifier or alter the Joint
objective. The fixture contains the 881-row local response, 12-atom support,
trial widths and A/C coefficients, and both evaluation paths. It can be replayed
without the original CIF, map, or SQLite database:

```sh
cmake --build build/joint-validation --target joint_local_replay_analysis --parallel 4
./build/joint-validation/bin/joint_local_replay_analysis \
  tests/fixtures/joint_local_replay/block-82-trial.json
```

The pre-correction captured trial is component 1, one-based sweep 2, block 82,
profile evaluation 1. The local design has 881 rows, 24 columns, 1,684 stored support
entries, and A/C scaled KKT `5.820766091346741e-11`. The primary profile
evaluation is valid and feasible. `CheckReplay` rejects it as
`replay-disagreement`.

The numerical tables and error decomposition below are the baseline recorded
before the production corrections. They are retained to explain the original
failure; they are not measurements of the final arithmetic path.

The analysis used AppleClang 21.0.0 on arm64 with IEEE-754 binary64, fast math
disabled, and `boost::multiprecision::cpp_dec_float_50` (50 decimal digits) for
the reference. Fixture doubles are read directly as `double` and converted to
the multiprecision type before any kernel, residual, or gradient calculation;
the reference `erf` uses `boost::math::erf` at that precision. On this platform,
`long double` has 53 binary digits, the same precision as `double`, so it was not
used as an oracle.

## Pre-correction results

| Quantity | Primary | Local replay | 50-digit reference | Primary/replay difference |
| --- | ---: | ---: | ---: | ---: |
| Gaussian kernel, maximum absolute error | `1.96e-16` | `1.84e-16` | — | `1.11e-16` |
| Charge kernel, maximum absolute error | `4.97e-16` | `4.97e-16` | — | `0` |
| Gaussian log-width derivative, maximum absolute error | `4.93e-16` | `4.56e-16` | — | `4.44e-16` |
| Charge log-width derivative, maximum absolute error | `3.15e-16` | `3.15e-16` | — | `2.22e-16` |
| Prediction, maximum absolute error | `2.13e-8` | `1.96e-8` | — | `1.284e-8` |
| Residual, maximum absolute error | `2.13e-8` | `1.96e-8` | — | `1.284e-8` |
| Width gradient infinity norm | `9.23e-5` | `1.506e-4` | `3.163e-4` | `6.291e-5` |
| A/C scaled KKT | `5.821e-11` | `1.746e-10` | `3.143e-10` | `1.164e-10` |

The primary and replay prediction errors are both small in absolute terms, but
the high-precision comparison shows that neither is an adequate certificate for
this trial's width gradient or A/C stationarity. The replay prediction is
slightly closer to the reference on this fixture; that is an observed result,
not a general guarantee about either implementation.

## Root cause

The local rows are severely cancellation-sensitive. In the row with the largest
primary/replay difference, row ID `6132229`, ten Gaussian/charge terms have a
sum of absolute values of about `5.736e7`, while the final prediction is about
`-5.77e-2`. The cancellation indicator is `9.95e8`. The largest local
prediction discrepancy is `1.284e-8`.

The per-row error decomposition shows that two effects matter:

1. **Prediction accumulation:** the sparse product accumulates the rounded
   terms differently from the scalar compensated replay. On the largest
   differing row, the primary accumulation error is `1.147e-8`, while replay's
   is `-1.368e-9`.
2. **Kernel and product rounding under large coefficients:** individual kernel
   differences are at most about `5e-16`, but A/C coefficients magnify them.
   For the same row, the net kernel-evaluation contribution is `6.03e-9` on
   both paths, and product rounding contributes about `-1.23e-9` on both.

Across rows, the maximum absolute errors against the multiprecision prediction
are `2.13e-8` for the sparse primary path and `1.96e-8` for the compensated
replay. Kahan summation on the primary terms reduces its worst prediction error
to `1.96e-8`; Neumaier and pairwise reductions reach `1.86e-8`. These are modest
improvements because the input terms themselves already carry kernel and
product rounding at the observed coefficient scale.

Residual formation adds essentially no error: each path's stored residual
matches its own `prediction - y` reconstruction within `2.78e-17`. The
primary/replay residual difference tracks the prediction difference. This is why
the width-gradient difference is amplified: swapping the primary residual for
the replay residual changes the width gradient by `6.291e-5`; swapping primary
and replay kernel derivatives while holding the residual fixed changes it by
only `1.36e-20`. The two production reduction arrangements differ by about
`3e-14`, also negligible next to the residual effect.

The A/C column norms agree with the high-precision reference within
`8.27e-15`. The largest raw `X^T r` discrepancy between primary and replay is
`7.93e-8`, which becomes a scaled-gradient discrepancy of order `1e-10`. The
high-precision fixed-state KKT is `3.143e-10`, above the unchanged `1e-10`
qualification threshold. The primary's double-precision KKT value therefore
understates the stationarity error for this state.

## Kernel and floating-point checks

The fixture has no support entries outside the `2.5` cutoff, so a cutoff error
does not explain this 6Z6U trial. The kernel checks found:

- `std::numbers::pi` and `M_PI` are equal in this build.
- The maximum measured `sqrt(square) * sqrt(square) - square` magnitude is
  `1.08e-15`.
- Center and near-center charge values and width derivatives agree with the
  50-digit reference to below `6.4e-16` absolute error.
- The local `CheckReplay` scalar kernel does not apply the production support
  cutoff. At the next representable squared distance above `6.25`, production
  returns zero while local replay returns charge `0.399857984`; the endpoint
  independent replay applies the cutoff and returns zero. This inconsistency
  did not affect the captured fixture, but it is a real path mismatch and needs
  a regression and correction.
- The tested FMA contribution had no material difference from the ordinary
  expression on these inputs. No FMA setting explains the observed discrepancy.

## Error bounds and trust meaning

For each prediction row, the offline tool combines measured kernel/product
rounding with a standard `gamma_n = n*epsilon/(1-n*epsilon)` accumulation term.
The resulting maximum row bounds are `5.46e-7` for primary and `7.24e-7` for
replay; both cover the observed prediction errors, with no measured bound
violation. These bounds are 25–37 times larger than the worst prediction
errors, so they do not distinguish a trustworthy trial from a bad one. The
analogous gradient bounds are `4.57e-4` and `3.34e-4`, with no observed
violation, but they are too broad to certify the `1e-12` width-stationarity
requirement. These are fixture-specific a posteriori checks, not a proven
production enclosure.

The local trust gate is right to reject this trial: the prediction, width
gradient, and KKT comparisons all exceed their existing consistency limits,
and the multiprecision KKT independently exceeds `1e-10`. Explaining a
discrepancy as floating-point error does not establish local optimality. No
trust threshold change is supported by this study. The independent endpoint
verification remains a separate calculation and remains required.

The numerical correction should improve residual/prediction accuracy and make
the unchanged A/C KKT threshold reflect the true computed residual. It must
also keep an explicit independent KKT check at `1e-10`; the replay agreement
check alone is not a substitute for that qualification. The cutoff mismatch
should be corrected as a separate formula-consistency fix. These changes are
evaluated in the subsequent production-fix and regression commits.

## Pre-correction 6Z6U outcome

The first full workflow run reproduced the failure at component 1, sweep 2,
block 82. That result is superseded by the post-correction run below.

## Post-correction behavior

The captured block 82 trial now replays as `trusted` with prediction difference
0, width-gradient difference `7.587070620200531e-15`, and KKT difference 0.
Its primary evaluation is valid, feasible and has scaled projected KKT
`5.779504479164199e-12`. The test also confirms that deliberately changed
coefficients, width gradient and KKT are still rejected. The capture preserves
the original failed decision (`replay-disagreement`) as historical data; the
test recomputes the current decision from its exact trial parameters.

The supported production changes were the evidence-backed numerical and path
corrections: replay derives widths with the same vectorized `exp` evaluation as
the profile, applies the same 2.5 support cutoff, and keeps scalar kernel
evaluation and accumulation independent from sparse matrix construction.
Prediction, transpose products, derivative reductions and replay now use
compensated binary64 summation with an FMA product-error term. Local free-face
QR refinement uses the same compensated double residual for at most eight
correction solves. Accepted FixedNeighbor caches are rebased to the full replay
prediction, residual and objective already computed for that candidate.
Independent endpoint replay remains a separate scalar calculation.

The current offline decomposition still measures the alternate matrix-term
accumulation against the scalar replay: its maximum prediction difference is
`4.1796e-9`. In that analyzer, the matrix-term primary accumulator is a Kahan
diagnostic; it is not the exact `StablePrediction()` path used by the current
trust decision. The trust fixture's production `CheckReplay()` comparison is
0. The highest-difference diagnostic row is row `5931763`: the absolute
contributions sum to `5.7363e7`, the 50-digit prediction is `-0.49413534`, and
the cancellation indicator is `1.161e8`. For that matrix-term primary path,
the maximum prediction error against the 50-digit reference is `2.118e-8`;
the scalar replay's is `1.742e-8`. Kernel errors are around `1e-16` before
coefficient scaling, while their maximum net contribution error reaches
`1.742e-8`. Accumulation alone is not the only source.

The analyzer's residual cross-check compares the profile residual with that
alternate matrix-term prediction, so its `4.1796e-9` primary formation delta
includes the prediction-path difference. The scalar replay residual is formed
from its own prediction with `2.78e-17` maximum round-trip error. This does not
indicate a separate residual formula defect. In the pre-correction capture,
the `6.291e-5` primary/replay width-gradient difference was driven almost
entirely by the residual difference; swapping derivative kernels alone
changed the result by about `1.36e-20`.

On the corrected fixture, primary and replay width gradients differ by only
`7.5871e-15`, but each differs from the 50-digit functional-kernel gradient by
about `1.5367e-4` absolute. The 50-digit gradient infinity norm is
`1.5857e-4`; the width gradient is not stationary. The corrected binary64
primary and replay KKT values both equal `5.7795e-12`, while the 50-digit
functional-kernel KKT is `1.5760e-10`. Thus the unchanged cross-path trust
gate confirms agreement between its two binary64 implementations; it is not a
rigorous enclosure of the exact-function arithmetic error. The fixture is a
nonstationary LM trial, and this result does not qualify a component endpoint.

No prediction, gradient, KKT, width-stationarity, rank or endpoint threshold was
relaxed. No error enclosure was adopted: the measured fixture bounds are too
broad to establish the `1e-12` width-stationarity gate or distinguish correct
from incorrect trials. The 50-digit Boost calculations remain offline
research/test support; production does not invoke a multiprecision solver.

The high-precision reference exposes an important limit: for the original
fixed-state fixture, the exact-kernel A/C KKT estimate is `1.576e-10`, above
the production `1e-10` threshold, even though the corrected binary64 primary
and replay certificates agree below that threshold. Cross-path agreement is
therefore not a proof of exact mathematical stationarity. It remains paired
with the unchanged production feasibility/KKT checks and does not establish
component completion or endpoint qualification.

## Replay path independence and stationarity scope

The checks provide different amounts of independent evidence. `CheckReplay()`
recomputes the scalar-kernel prediction, residual and width derivative with
compensated accumulation and an FMA product-error term. For its A/C KKT check,
it reuses the primary sparse design's column norms and transpose product, then
compares that KKT with the primary certificate. `CheckIndependentReplay()`
performs a separate scalar-kernel traversal and independently rebuilds
prediction, residual, objective, column norms, A/C gradients/KKT and width
gradient; it does not solve another constrained A/C problem. Both paths use
the same fixed support, rounded parameters, objective, binary64 math library
and deterministic kernel formulas. They are independent evaluation paths, not
independent high-precision references. The Block 21 QR audit and the offline
50-digit functional-kernel calculation supply the separate accuracy evidence.

The captured Block 82 trial passes `CheckReplay()` while its width-gradient
magnitude is greater than the endpoint stationarity threshold `1e-12`, but
`CheckIndependentReplay()` returns `independent-kkt-failed` for its recomputed
A/C KKT (`1.1641532e-10`). The former reports KKT `5.7795e-12` using
primary-design columns; the latter independently rebuilds the columns from
scalar kernels. A 50-digit functional-kernel audit reports KKT `1.576e-10`.
This fixture therefore preserves evidence that replay agreement and primary
KKT are not enough to prove precision-independent local optimality. A separate
well-conditioned synthetic trial now verifies that both replay paths can pass
while the width gradient remains nonstationary, as expected for an LM trial.
Deliberately changed width gradient, KKT certificate and coefficient fixtures
continue to be rejected.

## Final 6Z6U run

The current validation run used source HEAD `b53c1d8bf90bf73b2b3920d4d6a4a94beafae5fc`,
a fresh read-only-source SQLite backup at
`/private/tmp/rhbm-6z6u-width-stability-measure2-b53c1d8b.sqlite`, and unique
key `6z6u_width_measure2_b53c1d8b`. The old block 82 `CheckReplay()`
disagreement did not recur. Search stopped at component 1, one-based sweep 4,
block 21 with `inner-kkt-failed`. Its 930 × 24 design is full rank (rank 24,
free rank 23), but the feasible local A/C candidate has scaled KKT
`4.656612873e-10`. No endpoint verification, global rank qualification,
target estimability or trusted assembly was reached.

The CLI exited 0, but search completion was false. It recorded 5,056
evaluations and 620 accepted block updates over three completed sweeps. The
wrapper measured 33.645 s wall time and macOS `ru_maxrss` reported child peak
RSS of 1,340,866,560 bytes. The detailed full-data counts, width identities
and unavailable qualification states are in the
[halo-pruning validation](joint-halo-pruning-validation.md#full-data-validation-6z6u).
