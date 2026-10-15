# Bounded global rank research

This document records the retired Global Rank Certification contract and
research controls. The production rank pipeline, `DesignSpectrum()` oracle,
and its experiment drivers are no longer in the current estimator. The
equations below explain historical results only; they are not current runtime
checks, gates, or resource estimates. See the
[runtime contract](joint-component-runtime.md) for the current convergence,
conditioning, and stability architecture.

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
P_Z = Z Z^+
J_perp = (I - P_Z) D
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

## Structured sparse projection feasibility route

A separate test-only route normalizes the free-design sparse matrix, forms the
global sparse Gram matrix `Z^T Z`, and uses Eigen `SimplicialLDLT` with AMD
ordering to solve `(Z^T Z) c = Z^T u`; it returns `Z c`. This is a controlled
projection experiment, not a production solver. It keeps all global column
couplings in the Gram graph and factorization. It does not drop separator edges
or combine local rank decisions into a global claim.

Before copying or multiplying matrices, its preflight reserves a worst-case
quadratic sparse Gram/factor/workspace bound, input copies, graph adjacency,
and observation/parameter vectors against the existing 3.2 GB planning
budget. For admitted small cases it records exact graph components,
articulation-column separator size, largest interior region after removing a
single articulation, interface dimension, Gram and factor nonzeros, maximum
factor column size, fill ratio, pivot-diagonal condition estimate, projection
residual, idempotence error, runtime and whole-process RSS. A separator is a
measurement only; the full global Gram matrix is still factored.

The normal equations square the condition number of `Z`. A successful sparse
factor and a small solve residual therefore establish neither the original
SVD rank contract nor a singular-value bound. Numerical breakdown or a bad
projection is reported as unavailable. No rank certificate is produced, even
when local blocks or the sparse factor appear full rank. The giant synthetic
shape is handled by the checked preflight only; its structural full-design
rank obstruction remains unchanged.

### Projected-width measurements and qualification result

The projected-width policy uses the original observation row count `N` and
width column count `P`:

```text
tau = epsilon * max(N, P) * sigma_max(J_perp)
```

Free-design columns are normalized before projection; the width derivative
uses the existing observation scaling and production column construction. The
free-design rank path uses its existing strict `>` comparison. The projected
width SVD keeps Eigen's current rank-boundary behavior, and a singular value
within one representable step of `tau` is threshold-sensitive and cannot
qualify `identified`. These controls do not substitute supported rows for
`N`, alter the threshold, or treat a structural matching as numerical rank.

The deterministic dense oracle contains 13 cases, each with eight matrix rows:
well- and ill-conditioned full rank, exact projected-width deficiency,
free-design deficiency, a near-threshold width direction, a connected
full-rank support graph, local and strong coupling, a weak separator,
target/halo selection, an active amplitude boundary, multiple components,
and cancellation-sensitive columns. It reports the normalized design column
norms, both spectra and thresholds, ranks under both the Eigen and strict
boundaries, support counts, projection residual and elapsed time. In the
serial Release run, the width oracle classified 11 of 13 controls as full rank
and two as deficient. Ten controls had both full design and full projected
width rank; three failed at least one rank check. The near-threshold case had
`N=1,000`, `P=3`, `tau=2.220446049250313e-13`, and a minimum singular value
equal to `tau`: Eigen reported rank 3 while the strict count was 2, and the
existing threshold-sensitivity guard made it unsuitable for qualification.
The ill-conditioned full-rank control's minimum singular value was
`5.0e-15`, about 2.8 times its `1.77636e-15` threshold. The exactly deficient
width control had singular values approximately
`[sqrt(3), 1, 2.05e-16]`; the weak-separator control had
`[1, 0.87287156, 2.81e-16]`. These are small dense references, not evidence
for a giant component. The isolated oracle test's process peak was 24,772,608
bytes; this is whole-process RSS, not an allocation-level measurement.

The values-only and thin-right-vector CompactSvd modes agreed on singular
values, rank, threshold and threshold-sensitivity for all 13 controls. The
below/equal/above strict-boundary controls returned ranks 2/2/3 in both modes.
The right-vector output for three width columns occupies 72 bytes in these
fixtures. Tiny-matrix timings ranged from sub-microsecond values to tens of
microseconds on first use and were not stable enough to compare; no
function-level RSS saving was measured. Production still requests the
thin right vectors because runtime JSON and offline audit tools consume the
up-to-three weak directions. The directions do not participate in the rank
decision, but removing them from the current result would change a persisted
diagnostic contract. This study therefore did not promote a values-only
production path. For `P=37,406`, the explicit `P × P` output matrix alone is
11,193,670,688 bytes. Omitting it would remove that allocation in principle,
but the same all-active production preflight is still 91,426,019,120 bytes;
this study measured no giant RSS reduction.

The exact-oracle-backed matrix-free operator uses the same orthogonal
projection for both actions:

```text
Apply(v)        = D v - Q (Q^T D v)
ApplyAdjoint(u) = D^T (u - Q (Q^T u))
```

Here `Q` is the thin-U basis from a dense SVD of normalized free design `Z`.
Across the controls, forward and adjoint absolute errors were at most
`5.65e-17` and `5.01e-17`, the adjoint inner-product difference at most
`4.45e-16`, `Z^T Apply(v)` at most `1.78e-16`, and projection idempotence
error at most `3.93e-16`. Dense and tiled projected-width singular values
matched in the controls; the near-threshold case remained marked sensitive.
These algebra checks use the dense oracle basis and cost `O(N rank(Z))` storage,
so they are not a bounded projector. The giant action-vector preflight is
751,790,728 bytes including the assumed 600,000,000 resident bytes and passes
the 3.2 GB storage budget, but the bounded global projection result is
`Unavailable`; it performs zero projection calls and supplies no rank evidence.

The test-only structured path factors the full sparse Gram matrix, preserving
all global couplings. The 8-row weak-separator positive control had one exact
block, one articulation column, a one-dimensional interface, a two-column
largest interior, seven factor nonzeros for ten Gram nonzeros, and a dense
projection difference of `2.37e-16`. The strongly coupled control increased
the diagonal pivot ratio to about 98.7 and had a `5.93e-16` projection
difference. A deficient free-design control broke down and returned
`Unavailable`. Every structured projection leaves rank status `Unavailable`,
including successful small projections.

Chain/cube scalability cases are sparse free-design projection controls only.
They do not include a width derivative or a dense SVD, so their numerical
rank is explicitly unknown. Each was run in a separate serial CTest process;
RSS is the whole-process high-water mark. MB below is decimal.

| Topology | Atoms | Rows × columns | Gram / factor nonzeros | Maximum factor column | Estimated peak MB | Process peak RSS MB | Projection error / idempotence | Runtime (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Chain | 128 | 58,300 × 256 | 1,528 / 892 | 4 | 42.41 | 36.95 | 2.14e-15 / 2.89e-15 | 0.00468 |
| Cube | 128 | 47,792 × 256 | 6,768 / 9,500 | 72 | 41.64 | 36.72 | 2.03e-15 / 2.59e-15 | 0.00436 |
| Chain | 256 | 116,540 × 512 | 3,064 / 1,788 | 4 | 68.69 | 45.61 | 2.12e-15 / 2.92e-15 | 0.00589 |
| Cube | 256 | 93,431 × 512 | 14,728 / 28,816 | 124 | 66.84 | 45.43 | 2.06e-15 / 2.76e-15 | 0.00752 |
| Chain | 512 | 233,020 × 1,024 | 6,136 / 3,580 | 4 | 146.97 | 68.76 | 2.12e-15 / 2.95e-15 | 0.00945 |
| Cube | 512 | 183,040 × 1,024 | 31,616 / 105,600 | 282 | 142.97 | 68.04 | 2.03e-15 / 2.60e-15 | 0.01624 |
| Chain | 1,024 | 465,980 × 2,048 | 12,280 / 7,164 | 4 | 403.37 | 111.87 | 2.12e-15 / 2.92e-15 | 0.01968 |
| Cube | 1,024 | 361,867 × 2,048 | 65,864 / 327,260 | 442 | 395.04 | 110.17 | 2.06e-15 / 2.68e-15 | 0.04822 |

All eight projection runs completed below the 3.2 GB planning budget. The
connected chain/cube patterns had one exact block; no separator was used to
drop edges. Cube factor fill grew from 1.40× Gram nonzeros at 128 atoms to
4.97× at 1,024 atoms, while chain factor nonzeros stayed below its Gram count.
These runs record sparse solve completion, projection residual and
idempotence; they did not compare the large projections with an independent
dense oracle. Since normal equations square `cond(Z)`, the fill and
pivot-ratio measurements do not establish stable global rank evidence.

For the giant shape, all estimates below use 5,938,087 rows, 37,406 atoms,
74,811 support memberships, 149,622 design nonzeros, and an assumed
600,000,000 resident bytes. No giant factor, Gram matrix, dense design, or
projection was allocated. The production stage values are from the current
`PlanMemory` path; the structured estimates are its conservative checked
four-times dense Gram/fill workspace bound. Peak RSS is only the preflight test
process, not a giant-data measurement.

| Free columns | Production free factor | Free SVD workspace | Projection coefficients | Width factor | Width SVD workspace | Production estimated peak | Structured estimated peak | Status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 74,812 | 44,774,682,752 | 134,324,048,256 | 22,387,341,376 | 11,193,670,688 | 33,581,012,064 | 266,691,682,368 | 359,052,341,856 | ResourceLimit |
| 56,109 | 25,185,759,048 | 75,557,277,144 | 16,790,506,032 | 11,193,670,688 | 33,581,012,064 | 167,865,180,120 | 202,337,211,624 | ResourceLimit |
| 37,406 | 11,193,670,688 | 33,581,012,064 | 11,193,670,688 | 11,193,670,688 | 33,581,012,064 | 91,426,019,120 | 90,396,764,144 | ResourceLimit |

The matrix-free action estimate is 751,790,728 bytes and passes its storage
preflight, but projection and rank are `Unavailable`. The structured preflight
exceeds budget for all three faces before allocation. The structural
row-count bound proves the giant full design deficient because 37,406
supported rows are fewer than 74,812 columns. It does not prove active-free
design rank or projected-width rank. The giant shape is not a full-rank
positive control and no real 37,406-atom cryo-EM component was tested.

### Rank-certificate contract and promotion gate

`ClassifyProjectedRankEnvelope` is a test-support contract checker, not a
spectral estimator. It accepts an envelope only when the caller asserts that
projection, matrix-vector, spectral, floating-point, normalization and
threshold-rounding errors are rigorously enclosed. It certifies full rank only
when the lower singular-value bound exceeds the upper threshold, and
deficiency only when the upper minimum-singular-value bound is at or below the
lower threshold. Equality is deficient under the existing strict comparison.
The helper rounds interval arithmetic outward; actual certification still
requires an independently justified enclosure for every error term.
Rule-level tests cover full, ill-conditioned, exact deficient, below/equal/
above threshold and overlapping intervals. The exact threshold and its
one-ULP neighbors return `Unavailable` when outward-rounded intervals overlap;
the strict production boundary itself remains checked separately by the
below/equal/above SVD control. The 13 actual fixture loops supply
no such enclosures: all 13 returned `Unavailable`, with zero full or deficient
certificates and zero false certifications. This is not a 0% empirical error
rate for a completed certification method; no projected-width rank method in
this study produced certificates.

`ResourceLimit` means a checked memory preflight rejected work before the
factor or projector allocation. `Unavailable` means the computation lacks a
valid numerical result or sufficient evidence. Neither status passes rank
qualification.

| Method | Full-rank evidence | Deficiency evidence | Memory evidence | Decision |
| --- | --- | --- | --- | --- |
| Dense SVD oracle | Dense small-matrix reference | Dense small-matrix reference | Small fixtures only; 24.8 MB process peak | Reference |
| Matrix-free operator | None from the bounded route; exact-oracle algebra only | None from the bounded route | Giant action estimate 751.8 MB; no projector ran | Experimental only |
| Structured projection | No global SVD bounds | Sparse factor breakdown on a deficient control only | Topologies observed below 112 MB RSS; giant preflight 90.4–359.1 GB | Not worthwhile as current normal-equation route |
| Projected-width spectral bounds | No rigorous enclosure implemented | No rigorous enclosure implemented | Not measured for a projected-width rank solver | Requires further research |

The dense-oracle operator algebra is correct on the small controls, but
promoting it would require a bounded, numerically reliable global projector.
The current structured normal-equation preflight rejects the giant case and
its conditioning squares. Stronger sparse factorization and/or rigorous
singular-value enclosures would be a separate research effort. No experimental
route changes production `AssessEvaluatedImpl()`, Trusted State availability,
Design Rank, Projected-width Rank, or Derived Target Identifiability.
