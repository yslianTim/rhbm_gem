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

The updated samples are converted into an `RHBMMemberDataset` over the requested
fit range, then `rhbm_helper::EstimateBetaMDPDE` estimates zero-intercept beta
parameters. The decoded OLS and MDPDE Gaussian models get the fixed intercept
attached before returning.

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

The loop uses the `max_iterations` and `tolerance` values from
`RHBMExecutionOptions` produced by `MakeExecutionOptions`. At the time of this
note, those values are the execution defaults: `100` iterations and `1.0e-5`
tolerance.

The outer loop converges when the fixed-point error is below
`RHBMExecutionOptions::tolerance`:

```text
fixed_point_error = abs(raw_intercept - current_intercept)
```

The damped update is:

```text
next_intercept = current_intercept
               + 0.5 * (raw_intercept - current_intercept)
```

## Residual Intercept Estimation

The raw intercept is estimated by `EstimateResidualInterceptParameter`.

It decodes the current MDPDE beta result into a signal model and builds residual
samples only for distances in `[1.0, 2.0]`:

```text
basis(distance) = signal_model.InterceptBasisAtDistance(distance)
residual        = sample.response - signal_model.SignalAtDistance(distance)
```

The intercept candidate is then estimated as a robust slope-through-origin
problem:

```text
residual ~= intercept * basis
```

The Huber slope solver starts from ordinary least squares, computes residual
scale from median absolute deviation, applies Huber weights, and iterates until
the slope change is below tolerance or the Huber iteration limit is reached.

## Max-Iteration Fallback

If the outer loop reaches the maximum iteration count before convergence, the
fallback result is the recorded fixed-intercept fit with the lowest
fixed-point error:

```text
best_error = min(abs(raw_intercept - current_intercept))
```

If no best candidate is recorded, the function calls `EstimateLocalGaussian`
with the best known intercept as a last-resort fallback. In normal finite-data
paths, the first iteration records a candidate because residual intercept
estimation falls back to the current intercept when it cannot estimate a better
one.

The raw intercept update is only used for convergence, fallback ranking, and
building the next fixed-intercept candidate. The result returned from a
successful iteration is the fixed-intercept fit for that iteration, with the
current intercept attached to the decoded OLS and MDPDE models.

## Discussion Points

Items worth revisiting before changing the algorithm further:

- The intercept clamp `[-1.0, 1.0]` is hard-coded and independent of data scale.
- Residual intercept estimation uses `[1.0, 2.0]`, independent of the configured
  fit range.
- The fixed-point convergence and fallback ranking use only intercept update
  error, not MDPDE objective value or amplitude/width parameter movement.
