# Local Gaussian Estimation with Offset

`rhbm_gem::core::EstimateLocalGaussianWithOffset` estimates the amplitude,
width, and offset of a local `GaussianModel3D`. The implementation is in
[`src/core/GaussianEstimator.cpp`](/src/core/GaussianEstimator.cpp); the public
interface and `FitOptions` are in
[`include/rhbm_gem/core/GaussianEstimator.hpp`](/include/rhbm_gem/core/GaussianEstimator.hpp).

## Model

For distance `r`, the fitted response is:

```text
response(r) = signal(r; amplitude, width)
            + offset * offset_basis(r; width)
```

The model's offset parameter is not an additive vertical constant. It is
multiplied by a basis that depends on the current width:

```text
offset_basis(r; width) =
    sqrt(2 / pi) / width                              when r < 1.0e-5
    erf(r / (width * sqrt(2))) / r                    otherwise
```

The estimator therefore alternates between fitting the signal parameters and
fitting the offset.

## Algorithm

1. Validate the fit range, `alpha_r`, and `offset_initial`.
2. Fit a zero-offset Gaussian to initialize amplitude and width.
3. Clamp `offset_initial` to ±10% of the initial MDPDE amplitude.
4. Repeat:
   1. keep the current width and offset fixed, subtract their offset
      contribution from every sample, and fit amplitude and width;
   2. keep the new amplitude and width fixed and estimate an offset from
      residual samples at distances in `[1.0, 2.0]`;
   3. clamp the estimated offset to ±10% of the latest MDPDE amplitude;
   4. stop when the offset change is below `1.0e-5`;
   5. otherwise update the offset with damping factor `0.5` and continue.
5. Return the converged fixed-offset fit. If the loop reaches 100
   iterations, return the iteration with the smallest offset change.

The loop uses the default `RHBMExecutionOptions` iteration limit and tolerance.
`FitOptions` currently controls only the fit range, thread count, and outer-loop
diagnostic logging.

## Signal Fit with a Fixed Offset

`EstimateLocalGaussianWithOffsetModel` removes the offset contribution using
the supplied model:

```text
adjusted_response(r) =
    sample_response(r)
    - (offset_model.ResponseAtDistance(r)
       - offset_model.SignalAtDistance(r))
```

`rhbm_helper::BuildMemberDataset` then keeps samples that:

- are inside the inclusive `FitOptions` distance range; and
- have a positive adjusted response.

The retained samples are transformed into a two-parameter log-quadratic
regression:

```text
X(r) = [1, -0.5 * r^2]
y(r) = log(adjusted_response(r))
```

`rhbm_helper::EstimateBetaMDPDE` fits this dataset with `alpha_r`. The decoded
MDPDE parameters provide the next amplitude and width; the fixed offset is
attached to both the decoded OLS and MDPDE models.

If no sample survives the range and positivity filters, the dataset builder
uses one zero-valued fallback row.

## Offset Fit from Residuals

`EstimateResidualOffsetParameter` decodes the current MDPDE signal model and
uses only samples in the fixed distance range `[1.0, 2.0]`:

```text
X(r) = [signal_model.OffsetBasisAtDistance(r)]
y(r) = sample_response(r) - signal_model.SignalAtDistance(r)
```

This is a one-parameter regression through the origin:

```text
y = X * offset
```

The same `alpha_r` and execution options used by the signal fit are passed to
`rhbm_helper::EstimateBetaMDPDE`. Its one-element MDPDE beta vector is the
offset candidate.

The current offset is retained when:

- the fitted width is non-finite or non-positive;
- fewer than two finite, non-zero-basis residual samples are available;
- the estimated beta is not one finite value; or
- subtracting the candidate offset would produce a non-finite value or
  exceed the `float` response range for any sample.

When the current offset is retained, the offset change is zero and the
outer loop converges immediately.

## Iteration and Result Selection

For current offset `c`, latest MDPDE amplitude `a`, and clamped candidate
`c_raw`:

```text
lower  = -0.1 * a
upper  =  0.1 * a
c_raw  = clamp(offset_candidate, lower, upper)
error  = abs(c_raw - c)
c_next = clamp(c + 0.5 * (c_raw - c), lower, upper)
```

Each iteration records the fixed-offset signal fit and its `error`.
Convergence returns the current iteration, whose models still contain `c`.
Maximum-iteration fallback returns the recorded fit with the lowest `error`;
it does not perform another fit in the normal fallback path.

The returned `LocalGaussianResult` preserves the signal fit's
`RHBMBetaEstimateResult`. Gaussian standard-deviation models remain zero because
this path does not decode parameter uncertainty.

## Important Implementation Constraints

- Signal fitting uses `FitOptions::distance_min` and `distance_max`; offset
  fitting always uses `[1.0, 2.0]`.
- The offset range is recalculated from the latest MDPDE amplitude on every
  iteration. A zero amplitude constrains the offset to zero.
- Non-positive adjusted responses do not participate in the logarithmic signal
  fit.
- Width changes also change the offset basis used in the next iteration.
- `alpha_r` controls robustness in both the signal and offset MDPDE fits.
- Convergence and fallback ranking consider only offset change, not
  amplitude/width movement or the MDPDE objective.
