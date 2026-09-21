# Noise and position-mismatch validation

Status: **complete with applicability limitations**. The 326 generated inputs
produced 450 fitted outcomes: all retained a state, 65 passed current runtime
convergence, and 385 failed. These totals describe this deliberately unbalanced
pilot matrix; they are not a population success-rate estimate.

These are the final 450 outcomes after fixing an experiment-only JSON parsing
issue and repeating the same matrix. The original run is retained separately;
recovery repetitions are not pooled as additional Monte Carlo samples. Seeds,
observations, geometry, selection and runtime success criteria were unchanged.

## Controlled design

Two carbon atoms have target/halo separations 1.2 and 0.6 Å. Their (A, B, C) values
are (2, 0.5, 0.2) and (2.3, 0.55, 0.15). The 25³ map has 0.5 Å spacing and origin
(-6, -6, -6); the target fixes 515 rows. Halo never expands these rows. Python's
independent scalar forward generator is checked against the existing independent
C++ scalar control for every input; a discrepancy above `1e-12` stops the runner.
There is no map normalization: A is in input-map-unit·Å³, B in Å and C in
input-map-unit·Å. These coefficients are not independently calibrated charges.

Matched noise conditions use σ = 1%, 5%, 10% of the **unshifted noiseless** target
RMS. Independent and correlated fields have the same theoretical marginal
variance. Correlated noise uses a separable seven-tap Gaussian (σ = one voxel),
L2-normalized weights, padded generation and valid cropping; individual samples
are neither demeaned nor variance-normalized. All noisy conditions use PCG64
seeds 20260921–20260940. The exact inputs and output hashes are recorded.

Mismatch moves only the generating halo by +0.05 or +0.15 Å along the atom axis,
including its generating support; the fitting coordinates and target domain stay
fixed. Each shift has a noiseless control and 5% independent-noise replicates.
This isolates a position-model mismatch, including its cutoff effect, not an
optimizer error relative to a correctly specified model.

Fixed initialization uses B₀ = (0.6, 0.6) for every condition. Production
initialization is paired on identical observations for noiseless matched cases,
5% independent/correlated noise, and 0.15 Å mismatch with/without noise. It reads
the full local map as the production initializer normally does, possibly beyond
the target domain. Initializer output does not supply known A/C background.

## Results and interpretation

Both noiseless matched geometries recover A/B/C to approximately `1e-15` with both
initializers and pass runtime checks. Noise increases parameter error, and equal
marginal variance does not make independent and correlated noise equivalent.
For the 1.2 Å geometry, fixed-init target A RMSE at 5% noise is about 0.0417 for
independent noise and 0.1194 for correlated noise. At 0.6 Å the corresponding
values are 0.0617 and 0.1415. The
[complete condition table](figures/joint-validation/noise-table.md) reports both
roles, both initialization modes and actual convergence counts.

Mismatch produces systematic parameter compensation even without noise. At
0.15 Å shift and 0.6 Å separation, fixed-init absolute A deviations are about
0.7168 for the target and 0.9014 for halo. These are deviations from generating
parameters, not proof of solver bias relative to a correctly specified LS model.
The noiseless mismatched result is also preserved as a procedural reference for
paired noisy comparisons; it is not called a proven global optimum.

Many available estimates fail `width-stationarity` and/or `local-correction`.
The report preserves those failures even when two initialization modes give very
similar parameter values. There is no evidence here for changing the runtime
criteria or presenting these cases as converged. Precision, search and statistical
uncertainty are distinct; this experiment does not certify all failed endpoints.

## Statistical reporting

[Per-fit records](figures/joint-validation/noise-runs.json) retain all 450 outcomes,
including status, stop reasons, initialization validity, parameter errors and
residual diagnostics. [Grouped statistics](figures/joint-validation/noise-summary.json)
separate target/halo, fixed/production initialization, all available estimates,
converged-only estimates, and available-but-nonconverged estimates. They include bias, RMSE, sample variance, Monte
Carlo standard errors for bias and mean squared error, prediction RMSE, residual
RMSE, neighboring-voxel residual correlation, and failure reasons.

Unavailable-rate denominators are all recorded fits; nonconvergence-rate
denominators are fits with a state. Process failures/unexecuted inputs are counted
separately. Both rates have 95% Wilson intervals. With 20 replicates, zero failures
still permits a roughly 16% upper bound: no rare-failure or precision-coverage
claim is supported. A single noiseless run has no estimated variance or MCSE.
Missing quantities stay null; they are never zero-filled. All-available parameter
statistics are conditional on state availability, and converged-only statistics
are additionally selected by the runtime checks.

Residual correlation is a descriptive nearest-neighbor statistic over the fixed
domain. It is not an independently calibrated hypothesis test. This pilot tests
one noise correlation length and one type of model mismatch, not a general
cryo-EM noise distribution or real-data scientific validity.
