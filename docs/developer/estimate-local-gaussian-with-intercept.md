# EstimateLocalGaussianWithIntercept Algorithm Notes

This note documents the current implementation of
`rhbm_gem::core::EstimateLocalGaussianWithIntercept` in
`src/core/GaussianEstimator.cpp`. It is intended as a discussion aid for future
algorithm changes.

## Purpose

`EstimateLocalGaussianWithIntercept` estimates a local three-parameter
`GaussianModel3D`:

- amplitude
- width
- intercept

The underlying beta fit is performed on a zero-intercept response after
subtracting a fixed intercept offset. The function wraps that fixed-intercept
fit in an outer fixed-point loop that repeatedly re-estimates the intercept from
the residuals.

The public inputs are the local samples, `alpha_r`, `FitOptions`, and an optional
`intercept_initial` whose default is `0.0`. The function requires finite,
non-negative, ordered fit-range endpoints, a finite non-negative `alpha_r`, and
a finite initial intercept.

## Fixed-Intercept Fit

The inner fixed-intercept path is implemented by
`EstimateLocalGaussianWithOffsetModel`.

For each sample, it builds a zero-intercept response by subtracting the offset
part of the supplied model:

```text
model_offset(distance) = offset_model.ResponseAtDistance(distance)
                       - offset_model.SignalAtDistance(distance)

updated_response = sample.response - model_offset(distance)
```

When `offset_model.width != 0`, the offset is:

```text
model_offset(distance) = offset_model.intercept
                       * offset_model.InterceptBasisAtDistance(distance)
```

The intercept basis depends on the supplied model width. When
`offset_model.width == 0`, `GaussianModel3D::ResponseAtDistance` returns the
intercept directly, so `model_offset(distance) = offset_model.intercept`.

The updated samples are passed to `rhbm_helper::BuildMemberDataset`. The current
linearization keeps samples in the inclusive requested fit range only when
`updated_response > 0`, then builds the log-quadratic regression:

```text
X(distance) = [1, -0.5 * distance^2]
y(distance) = log(updated_response)
```

If no sample survives those filters, `BuildDatasetSeries` supplies one
all-zero basis/response row rather than leaving the dataset empty.

`rhbm_helper::EstimateBetaMDPDE` estimates the two zero-intercept beta
parameters. It receives the `RHBMExecutionOptions` produced by
`MakeExecutionOptions`, which copies `FitOptions::thread_size`, forces
`quiet_mode = true`, and leaves the other execution settings at their defaults.
`DecodeLocalGaussianResult` decodes both the OLS and MDPDE beta vectors, attaches
`offset_model.intercept`, leaves their standard-deviation models at zero, and
preserves the `RHBMBetaEstimateResult` in `LocalGaussianResult::fit_result`.

## Outer Fixed-Point Loop

The high-level flow is:

1. Run an initial zero-offset fit.
2. Build `current_model` from the initial MDPDE model and the clamped
   `intercept_initial`.
3. For each outer iteration:
   - fit the signal with `current_model` as the fixed intercept offset;
   - estimate a raw intercept from the residuals of that fixed-intercept fit;
   - clamp the raw intercept into the allowed intercept range;
   - measure fixed-point error as the absolute difference between the raw
     intercept and the intercept used for the current fit;
   - record the fitted result when it has the lowest fixed-point error so far;
   - damp and clamp the raw intercept update before building the next
     `current_model`.

The outer loop uses the `max_iterations` and `tolerance` values from the same
`RHBMExecutionOptions` shape produced by `MakeExecutionOptions` for the inner
MDPDE fits. At the time of this note, `FitOptions` does not expose either
setting, so both the inner MDPDE iterations and the outer fixed-point loop use
the `RHBMExecutionOptions` defaults: `100` iterations and `1.0e-5` tolerance.
The outer loop's maximum-iteration diagnostic is controlled separately by
`FitOptions::quiet_mode`; `MakeExecutionOptions` forces only the inner estimator
calls into quiet mode.

The raw intercept is clamped before error calculation. The outer loop converges
when the fixed-point error is below
`RHBMExecutionOptions::tolerance`:

```text
fixed_point_error = abs(raw_intercept - current_intercept)
```

The damped update is:

```text
next_intercept = current_intercept
               + 0.5 * (raw_intercept - current_intercept)
```

## MDPDE Intercept Estimation

The raw intercept is estimated by `EstimateResidualInterceptParameter`.

It decodes the current inner MDPDE beta result into a signal model, thereby
fixing the amplitude and width for the intercept fit. It then builds a
one-parameter linear regression dataset from samples whose distances are in
`[1.0, 2.0]`:

```text
X(distance) = [signal_model.InterceptBasisAtDistance(distance)]
y(distance) = sample.response - signal_model.SignalAtDistance(distance)
```

Samples with a non-finite residual or a non-finite or effectively zero basis are
omitted. The resulting model is linear in its only beta parameter:

```text
y = X * intercept
```

`rhbm_helper::EstimateBetaMDPDE` estimates this one-element beta vector using
the same `alpha_r` value and `RHBMExecutionOptions` as the inner amplitude/width
fit. The MDPDE beta value is used as the raw intercept; the corresponding OLS
beta is not attached to the returned OLS model.

Residual intercept estimation returns `current_intercept` unchanged when:

- the decoded MDPDE width is non-finite or non-positive;
- fewer than two usable regression rows remain;
- the estimated MDPDE beta does not contain exactly one finite value;
- a candidate would make any zero-intercept sample response non-finite or too
  large for the stored `float` response type.

An `EstimateBetaMDPDE` status such as `MAX_ITERATIONS_REACHED` or
`NUMERICAL_FALLBACK` does not by itself reject the intercept. The candidate is
accepted when its beta vector has the required finite one-parameter shape and
passes the zero-intercept sample safety check.

## Max-Iteration Fallback

If the outer loop reaches the maximum iteration count before convergence, the
fallback result is the recorded fixed-intercept fit with the lowest
fixed-point error:

```text
best_error = min(abs(raw_intercept - current_intercept))
```

The function normally returns the already-recorded `best_result`; it does not
refit that candidate despite the wording of the current debug message. Only if
no best candidate was recorded does it call `EstimateLocalGaussian` with
`best_intercept` as a last-resort fallback. In normal finite-data paths, the
first iteration records a candidate because residual intercept estimation falls
back to the current intercept when it cannot estimate a better one.

The raw intercept update is only used for convergence, fallback ranking, and
building the next fixed-intercept candidate. The result returned from a
converged iteration is the fixed-intercept fit for that iteration, with the
current intercept attached to the decoded OLS and MDPDE models.

## Discussion Points

Items worth revisiting before changing the algorithm further:

- The intercept clamp `[-1.0, 1.0]` is hard-coded and independent of data scale.
- Residual intercept estimation uses `[1.0, 2.0]`, independent of the configured
  fit range.
- The intercept basis changes with the fitted width, so changing the width also
  changes the response offset subtracted on the next iteration.
- Samples whose offset-adjusted response is not positive are omitted from the
  log-quadratic beta fit.
- The same `alpha_r` controls robustness for both the inner amplitude/width fit
  and the outer intercept fit.
- The fixed-point convergence and fallback ranking use only intercept update
  error, not MDPDE objective value or amplitude/width parameter movement.
