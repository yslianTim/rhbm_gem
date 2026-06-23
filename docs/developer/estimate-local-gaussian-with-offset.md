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
   4. combine the new amplitude and width with the clamped offset candidate,
      then stop when its three-parameter L2 distance from the current model is
      below `1.0e-5`;
   5. otherwise update the offset with damping factor `0.5` and continue.
5. Return the converged fixed-offset fit. If the loop reaches 100
   iterations, return the iteration with the smallest three-parameter L2
   distance.

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

The offset is estimated with a Huber M-estimator. The solver starts from the
ordinary least-squares slope, computes a robust residual scale from median
absolute deviation, and iteratively applies Huber weights:

```text
scale  = max(1.4826 * MAD, 1.0e-12)
cutoff = 1.345 * scale
```

Iteration stops when the slope change is below `1.0e-8` or after 50
iterations. This outer offset fit does not use `alpha_r`.

The current offset is retained when:

- the fitted width is non-finite or non-positive;
- no finite Huber slope can be estimated from the usable residual samples; or
- subtracting the candidate offset would produce a non-finite value or
  exceed the `float` response range for any sample.

When the current offset is retained, its contribution to the L2 distance is
zero. The outer loop still continues if the fitted amplitude or width changed
by enough to exceed the convergence tolerance.

## Iteration and Result Selection

For current model `m_current`, latest MDPDE signal model `m_fit`, latest
amplitude `a`, and clamped candidate offset `c_raw`:

```text
lower     = -0.1 * a
upper     =  0.1 * a
c_raw     = clamp(offset_candidate, lower, upper)
m_raw     = m_fit.WithOffset(c_raw)
error     = norm(m_raw.ToVector() - m_current.ToVector())
c_next    = clamp(c + 0.5 * (c_raw - c), lower, upper)
```

The norm is the unnormalized Euclidean distance across amplitude, width, and
offset. It uses the raw clamped offset candidate before damping. Each iteration
records the fixed-offset signal fit and its `error`.
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
- `alpha_r` controls only the inner amplitude/width MDPDE fit; the outer offset
  fit uses the fixed Huber constants above.
- Convergence and fallback ranking consider the unnormalized L2 movement of
  amplitude, width, and offset, but not the MDPDE objective.
