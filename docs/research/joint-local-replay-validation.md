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

The captured trial is component 1, one-based sweep 2, block 82, profile
evaluation 1. The local design has 881 rows, 24 columns, 1,684 stored support
entries, and A/C scaled KKT `5.820766091346741e-11`. The primary profile
evaluation is valid and feasible. `CheckReplay` rejects it as
`replay-disagreement`.

The analysis used AppleClang 21.0.0 on arm64 with IEEE-754 binary64, fast math
disabled, and `boost::multiprecision::cpp_dec_float_50` (50 decimal digits) for
the reference. Fixture doubles are read directly as `double` and converted to
the multiprecision type before any kernel, residual, or gradient calculation;
the reference `erf` uses `boost::math::erf` at that precision. On this platform,
`long double` has 53 binary digits, the same precision as `double`, so it was not
used as an oracle.

## Results

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

## 6Z6U baseline outcome

The full workflow was run with the provided CIF and map and a backup of the
SQLite database. It reproduced the failure at component 1, sweep 2, block 82.
The initial local profile was `qualified-inner`; local trial trust failed as
`replay-disagreement`. Search did not complete, stationarity and endpoint
qualification were not confirmed, and there was no assembled scientific state
or trusted target A/C/B estimate. The CLI's final process status is not treated
as numerical success.
