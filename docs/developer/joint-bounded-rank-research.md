# Bounded global rank research

This document records the historical rank contract and controls for the
memory-bounded global design rank study. `DesignSpectrum()` remains its
small-case oracle. Current production obtains interior Design Rank from the
free-design factor in projected-width preparation; see the
[runtime contract](joint-component-runtime.md). Experimental methods in this
study must match the oracle's normalized matrix and strict rank boundary before
they can be considered for production.

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

## Structural rank evidence

`AssessDesignStructuralRank()` traverses the kernel support at the endpoint
widths and applies a preflight before retaining compact row/column pairs. It
uses the actual nonzero Gaussian and charge entries, then counts unique rows
and supported columns. Its rank upper bound is at most the smaller of those
counts. A zero column or an upper bound below the full-design column count
proves deficiency. Structural support patterns that are identical are counted
as diagnostics only; they do not establish equal coefficient values.

`AssessStructuralRank()` also constructs the compact column/row bipartite graph
for offline controls. Maximum matching is an upper bound on numerical rank. A
matching smaller than the column count proves deficiency; a full matching
means only structurally full rank. The production endpoint runs the basic
obstruction check and does not spend time on matching when that check is
inconclusive. Matching is preflighted separately before its arrays are
allocated and can return `Unavailable`.

For the giant synthetic fixture, the full design has 74,812 columns, 149,622
actual nonzeros, and 37,406 unique supported rows. Therefore its structural
rank upper bound is 37,406 and full-design rank qualification is `Failed`.
The active state has 37,406 free A/C columns. Those columns have the same
37,406 supported rows, so the row-count check does not prove free-face
deficiency. The C-only graph has a perfect 37,406 by 37,406 matching; this is
still structural evidence, not a numerical full-rank certificate.

## Streaming scalar Gershgorin study

The experimental streaming pass retains column norm sums, row absolute sums
and one radius per column. For normalized `Z`, it estimates
`rho_i = sum_r |Z_ri| * (sum_j |Z_rj| - |Z_ri|)`, then
`lambda_min >= min_i(1-rho_i)` and
`lambda_max <= max_i(1+rho_i)`. It uses compensated sums and O(N+P) working
arrays; `ResourceWork` preflights those arrays before allocation.

The bound calculation uses ordinary floating-point normalization and
accumulation. Compensated summation reduces error but does not enclose it, so a
passing estimate for a general matrix remains `Unavailable` and is never a
full-rank certificate. The current exact certificate is deliberately narrow:
each column has one exact `+1` or `-1`, and those entries occupy distinct rows.
The normalized columns are then a signed permutation of identity, so every
singular value is exactly 1. It qualifies as full rank only when the unchanged
relative threshold is below 1. Zero columns and near-threshold cases remain
unavailable to this experimental numerical path; structural evidence handles
zero columns separately.

On the 5,938,087 × 74,812 synthetic full design, the streaming scan completed
in 0.00718 s. It reported `lambda_min_lower = -2.41421`,
`lambda_max_upper = 4.41421` and a 6.82843 gap, so the estimate did not prove
full rank. The scan's estimated peak was 678,017,655 bytes; sampled RSS was
674,004,992 bytes before and 675,790,848 bytes after, with a 675,790,848-byte
process peak. This is a completed estimate, not a timeout or a resource
rejection, and its numerical status remains experimental because the
accumulation has no certified enclosure.

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

## SparseQR feasibility experiment

The offline helper accepts an already normalized sparse design, applies the
unchanged `tau = epsilon * max(N, P) * sigma_max` value as Eigen's pivot
threshold, and records the resulting pivot-diagonal rank. It times
`analyzePattern()` and `factorize()` separately and records R storage, fill
ratio, the conservative resource estimate, RSS and status. Natural ordering
and Eigen COLAMD run in separate serial test processes so one ordering's
allocator high-water mark does not contaminate the other. The estimate reserves
`8 * N * P` bytes for worst-case row-column workspace, additional P-square and
input-sparse terms, then a 5% allocation margin. It is a preflight bound; it
does not prove that a passing SparseQR call stays below an operating-system
limit.

### Small SVD comparisons

Eight small controls were compared under both orderings: well-conditioned full
rank, ill-conditioned full rank, exact duplicate columns, strong A/C coupling,
mixed target/halo, and three near-threshold cases. The dense SVD oracle used the
same normalized sparse matrix and original row count. The near-threshold cases
had `sigma_min / tau` of 0.5, 1.0 and 2.0. Dense SVD ranks were 1, 1 and 2;
SparseQR ranks were 1, 2 and 2 under both orderings. Thus both orderings
disagreed with the strict SVD boundary at equality (2 disagreements among 16
ordering comparisons). The other five controls matched the oracle under both
orderings. No threshold was adjusted to hide the mismatch.

### Chain and cube topology measurements

Topology rank values below are SparseQR estimates. Dense SVD was deliberately
not run on these tall matrices, so these rows do not establish exact numerical
rank. `N/C` means Natural/COLAMD. Memory values are decimal MB; the planned peak
uses the 5% margin and the observed value is process peak RSS from that isolated
test process.

| Case | Matrix rows × columns | Input nonzeros | QR rank estimate N/C | R nonzeros N/C | Fill ratio N/C |
| --- | ---: | ---: | ---: | ---: | ---: |
| chain-128 | 58,300 × 256 | 131,840 | 256 / 256 | 18,029 / 18,029 | 0.137 / 0.137 |
| cube-128 | 47,792 × 256 | 131,840 | 256 / 256 | 31,062 / 23,186 | 0.236 / 0.176 |
| chain-256 | 116,540 × 512 | 263,680 | 512 / 512 | 71,554 / 71,554 | 0.271 / 0.271 |
| cube-256 | 93,431 × 512 | 263,680 | 512 / 512 | 123,752 / 97,399 | 0.469 / 0.369 |
| chain-512 | 233,020 × 1,024 | 527,360 | 1,024 / 1,024 | 285,101 / 285,101 | 0.541 / 0.541 |
| cube-512 | 183,040 × 1,024 | 527,360 | 1,024 / 1,024 | 493,084 / 367,413 | 0.935 / 0.697 |
| chain-1024 | 465,980 × 2,048 | 1,054,720 | Unavailable / Unavailable | NotRun / NotRun | NotRun / NotRun |
| cube-1024 | 361,867 × 2,048 | 1,054,720 | Unavailable / Unavailable | NotRun / NotRun | NotRun / NotRun |

| Case | R factor bytes N/C | Symbolic seconds N/C | Numeric seconds N/C | Planned peak MB N/C | Observed process peak MB N/C | Result |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| chain-128 | 217,376 / 217,376 | 0.00023 / 0.00098 | 0.663 / 0.683 | 177.97 / 177.92 | 152.04 / 151.99 | Completed |
| cube-128 | 373,772 / 279,260 | 0.00021 / 0.00141 | 0.713 / 0.623 | 154.99 / 154.93 | 143.28 / 143.28 | Completed |
| chain-256 | 860,700 / 860,700 | 0.00040 / 0.00209 | 4.917 / 5.119 | 590.66 / 590.68 | 491.24 / 435.99 | Completed |
| cube-256 | 1,487,076 / 1,170,840 | 0.00041 / 0.00329 | 5.127 / 4.975 | 490.42 / 490.45 | 470.01 / 444.06 | Completed |
| chain-512 | 3,425,312 / 3,425,312 | 0.00082 / 0.00440 | 39.867 / 40.325 | 2,199.01 / 2,199.05 | 1,392.15 / 1,308.25 | Completed |
| cube-512 | 5,921,108 / 4,413,056 | 0.00084 / 0.01050 | 42.272 / 35.780 | 1,767.16 / 1,767.12 | 1,166.62 / 1,164.25 | Completed |
| chain-1024 | Not allocated | NotRun / NotRun | NotRun / NotRun | 8,539.39 / 8,539.46 | 99.14 / 99.21 | ResourceLimit; factorization NotRun |
| cube-1024 | Not allocated | NotRun / NotRun | NotRun / NotRun | 6,744.24 / 6,744.27 | 95.06 / 95.09 | ResourceLimit; factorization NotRun |

For completed cases, both ordering estimates were full rank, but no exact
topology SVD comparison was run. COLAMD reduced cube R nonzeros by about 25%,
21% and 25% at 128, 256 and 512 atoms; chain R nonzeros were unchanged. The
completed 512 cases stayed below the 3.2 GB planning budget and observed peaks
were 1.16–1.39 GB. The 1,024-atom chain and cube estimates exceeded budget at
about 8.54 GB and 6.74 GB, so all four ordering runs returned `ResourceLimit`
before factorization. These are expected safe resource outcomes, not rank
failures or timeouts.

The giant 5,938,087 × 74,812 full-design shape has 149,622 input nonzeros. Its
SparseQR preflight estimates 4,154,773,326,379 temporary bytes and a
4,155,449,133,611-byte peak including resident data, then rejects the call
without factorization. The observed process peak was 675,823,616 bytes while
the input and preflight were resident. This estimate is intentionally
conservative; it does not establish what unbounded Eigen SparseQR would
allocate. At 512 atoms, observed SparseQR RSS was higher than the current dense
tiled `DesignSpectrum` memory estimate for the same shape; the 1,024 cases
were not factorized. SparseQR therefore showed no production memory advantage
over the existing tiled dense rank planner in this study.

### Qualification decision

The production one-way structural obstruction is retained: a zero design
column or a proven upper bound below the full-design column count marks design
rank qualification `Failed`, without an exact rank or fabricated spectrum.
Structural matching remains research evidence; full matching is not a
numerical-rank result. The scalar streaming pass remains experimental because
compensated accumulation does not provide a rigorous floating-point enclosure;
only the exact disjoint signed-unit-column case is certified. SparseQR remains
experimental because its pivot rank disagrees at the strict SVD boundary,
large topology rows were not checked against SVD, and its conservative memory
plan rejects 1,024-atom cases. No method in this study qualifies the global
design as a full production replacement for `DesignSpectrum()`.

Projected-width assessment, local correction and target identifiability keep
their existing resource bottlenecks and were not redesigned in this study.

## Projected-width matrix-free operator controls

The small research operator uses the same normalized free-design matrix `Z`
and observation-scaled width derivative `D` as its inputs. An exact-oracle
build computes a dense thin-U SVD of normalized `Z`, applies the original
observation-row policy and strict singular-value threshold, and retains the
resulting orthonormal basis `Q`. Its actions are

```text
Apply(v)        = D v - Q (Q^T D v)
ApplyAdjoint(u) = D^T (u - Q (Q^T u))
```

This applies the same orthogonal projector in both directions. It stores a
sparse normalized `Z`, sparse `D`, and a dense `N × rank(Z)` basis; each call
uses observation and parameter vectors. It does not materialize projection
coefficients, an observation-by-width dense matrix, or a projected-width
factor. The dense SVD oracle basis is `O(N p_F)` in the full-rank case, so this
mode is limited to small correctness controls and is not memory-bounded.
Its conservative preflight includes dense SVD/basis workspace, sparse copies,
vectors, caller resident bytes and the planning budget before allocation.

The controls compare forward and adjoint actions with the dense oracle, check
the adjoint inner-product identity, projection idempotence, design
orthogonality, and the projected-width spectrum from the existing tiled-QR
reduction. They report absolute and relative action errors, projection
residual, the original width threshold, estimated workspace, operator calls,
runtime and process peak RSS. Cancellation-sensitive examples are included;
reported floating-point parity is not a rigorous operator-error enclosure.

The bounded projection preflight currently reserves only sparse input and
vector-sized action work. Even when that storage fits the 3.2 GB planning
budget, it returns numerical status `Unavailable`: this repository has no
bounded global projector with validated rank behavior. No iterative residual
or approximate Ritz value is treated as rank evidence. The 37,406-atom shape
has an action-work estimate only; no projection call or full rank assessment
was run for it.
