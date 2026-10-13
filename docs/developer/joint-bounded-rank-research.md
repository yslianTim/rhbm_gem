# Bounded global rank research

This document records the rank contract and reproducible controls for the
memory-bounded global design rank study. The dense SVD path is the small-case
oracle; experimental methods must match its normalized matrix and strict rank
boundary before they can be considered for production.

## Existing numerical rank contract

`DesignSpectrum()` receives the full A/C design with `P = 2 * atom_count`
columns. It computes each weighted column norm, divides nonzero columns by
that norm, and leaves a zero-norm column zero by substituting scale 1. Its SVD
decision uses the original rank-policy row count `N` and full design column
count `P`, even when some observation rows have no structural support:

```text
tau = machine_epsilon * max(N, P) * sigma_max(Z)
rank = count(sigma_i > tau)
```

The comparison is strict. `sigma_max` is from the column-normalized design
`Z`; `N` is not replaced by the number of supported rows. Design rank uses the
full A/C domain. The active-face A/C domain is smaller and has its own column
count; its rank cannot stand in for full-design rank.

`ReferenceQR()` is the current dense rank-factorization primitive used before
the SVD. Its name does not mean the removed global independent reference
solve is used in production. The production estimator and endpoint checks are
unchanged by this research.

## Dense reference controls

`BoundedRank_test.cpp` fixes deterministic small matrices for well-conditioned
and ill-conditioned full rank, duplicate and zero columns, structural and
numerical deficiency, strict near-threshold values, an active amplitude face,
target/halo columns, coupled A/C columns, connected full-rank and deficient
components, and weak-signal cancellation sensitivity. The oracle compares
`ComputeSpectrum()` with a direct dense Jacobi SVD of the same normalized
matrix, and checks the original row count, column norms, threshold and strict
rank count. A connected full-rank control is included independently of the
large sparse preparation fixture.

The 37,406-atom preparation fixture is intentionally not a full-rank positive
control. It has 5,938,087 observation rows and 74,811 support memberships;
its eventual structural assessment must count unique supported rows rather
than treating affected-row entries as unique rows.

## Evidence states

Research reports distinguish an exact numerical SVD result, a structural rank
upper bound, an experimental floating-point estimate, and unavailable work.
A structural matching equal to the number of columns means only that the
support graph is structurally full rank. It is not a numerical full-rank
certificate. A strict structural upper bound below the column count proves
deficiency without identifying the exact numerical rank or singular spectrum.

Streaming bounds are sufficient only when a reliable lower bound on
`lambda_min(Z^T Z)` and upper bound on `lambda_max(Z^T Z)` satisfy
`lambda_min_lower > (tau / sigma_max)^2 * lambda_max_upper`. Floating-point
accumulation must be enclosed conservatively before this may be called a
certificate. Compensated summation alone does not establish that enclosure.

SparseQR's pivot-diagonal rank is retained as an experimental estimate. Its
threshold is not the SVD singular-value threshold, so comparisons record
disagreements rather than adjusting the production policy.
