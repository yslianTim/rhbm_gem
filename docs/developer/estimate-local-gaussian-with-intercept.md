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
   - read the fitted MDPDE robust objective from `RHBMBetaEstimateResult`;
   - check convergence with MDPDE parameter relative change or robust objective
     relative change;
   - estimate a raw intercept from residuals only when another iteration is
     needed;
   - damp and clamp that intercept update for the next iteration.

The first fixed-intercept iteration seeds the previous convergence state. The
combined convergence check is applied from the second fixed-intercept iteration
onward.

The outer loop converges when either value is below
`RHBMExecutionOptions::tolerance`:

```text
max_relative_change(amplitude, width, intercept)
OR relative_change(MDPDE robust objective)
```

Relative scalar change is measured as:

```text
abs(current - previous) / max(1.0, abs(previous), abs(current))
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
fallback result is the recorded iteration with the lowest finite MDPDE robust
objective.

The objective is computed in `RHBMMemberDataset` regression space, matching the
beta fit. For residual

```text
r = y - X * beta_mdpde
```

and final variance `sigma_square`, the weight is:

```text
w = exp(-0.5 * alpha * r^2 / sigma_square)
```

The minimized bounded objective is:

```text
alpha > 0: mean((1 - w) / alpha)
alpha = 0: mean(0.5 * r^2 / sigma_square)
```

If no finite-objective candidate is ever recorded, the function returns the
latest fitted result as a last-resort fallback.

The intercept update is only used to build the next fixed-intercept candidate;
it is not used for convergence or fallback ranking.

## Discussion Points

Items worth revisiting before changing the algorithm further:

- The intercept clamp `[-1.0, 1.0]` is hard-coded and independent of data scale.
- Residual intercept estimation uses `[1.0, 2.0]`, independent of the configured
  fit range.
- The robust objective is measured in log-quadratic beta-fit space, not original
  electric-potential response space.
