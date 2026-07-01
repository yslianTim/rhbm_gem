# Local Gaussian Estimation with Offset

`rhbm_gem::core::EstimateLocalGaussian` estimates the amplitude and width of a
local `GaussianModel3D` while holding the supplied offset model fixed. The
implementation is in [`src/core/GaussianEstimator.cpp`](/src/core/GaussianEstimator.cpp);
the public interface and `FitOptions` are in
[`include/rhbm_gem/core/GaussianEstimator.hpp`](/include/rhbm_gem/core/GaussianEstimator.hpp).

## Model

For distance `r`, the fitted response is:

```text
response(r) = signal(r; amplitude, width)
            + offset * offset_basis(r; width)
```

The model's offset parameter is not an additive vertical constant. It is
multiplied by a basis that depends on width:

```text
offset_basis(r; width) =
    sqrt(2 / pi) / width                              when r < 1.0e-5
    erf(r / (width * sqrt(2))) / r                    otherwise
```

The estimator therefore accepts a complete `GaussianModel3D` as the fixed
offset model. Its amplitude and width define the offset basis contribution that
is subtracted from each sample before fitting the next amplitude and width.

## Algorithm

1. Validate the fit range, `alpha_r`, and `offset_model.GetOffset()`.
2. Remove the fixed offset-model contribution from each sample:

   ```text
   adjusted_response(r) =
       sample_response(r)
       - (offset_model.ResponseAtDistance(r)
          - offset_model.SignalAtDistance(r))
   ```

3. Build the member dataset from adjusted samples inside the inclusive
   `FitOptions` distance range and with positive adjusted response.
4. Fit the two-parameter log-quadratic regression with
   `rhbm_helper::EstimateBetaMDPDE`.
5. Decode OLS and MDPDE amplitude/width parameters and attach the fixed
   `offset_model.GetOffset()` to both returned models.

If no sample survives the range and positivity filters, the dataset builder
uses one zero-valued fallback row.

## Regression Shape

The retained samples are transformed into:

```text
X(r) = [1, -0.5 * r^2]
y(r) = log(adjusted_response(r))
```

`alpha_r` controls the MDPDE fit. `FitOptions` controls the fit range and thread
count used by the execution options.

## Internal Fixed-Point Offset Estimator

`GaussianEstimator.cpp` also keeps an internal
`EstimateLocalGaussianWithOffset` helper for future workflows. It is not declared
in the public header and is not used by the current first-stage, second-stage, or
group fitting paths.

This helper alternates between two steps for a single atom:

1. hold the current offset model fixed and reuse `EstimateLocalGaussian` to fit
   amplitude and width with MDPDE; and
2. hold the fitted amplitude and width fixed, build residual offset-basis
   samples for distances in `[1.0, 2.0]`, and estimate the offset with the shared
   Huber M-estimator through origin.

The offset M-estimator uses an amplitude-scaled regularization prior scale
`max(abs(amplitude) * 0.1, 1.0e-12)`. The offset update is damped before the next
fixed-point iteration. If the loop reaches its iteration limit, the helper
returns the best fixed-point candidate tracked during the loop.

## Important Implementation Constraints

- `EstimateLocalGaussian` is the only public local Gaussian fitting entry point.
- The default offset model is `GaussianModel3D{ 0.0, 1.0, 0.0 }`.
- Second-stage local fitting passes the joint-offset snapshot model directly so
  the fixed offset basis uses that model's width.
- Non-positive adjusted responses do not participate in the logarithmic signal
  fit.
- Gaussian standard-deviation models remain zero because this path does not
  decode parameter uncertainty.
