# Joint Component Assessment Memory Census

This census records the current production rank preflight for one deterministic
synthetic component. It is not a physical cryo-EM support graph. The planner
uses a 3,200,000,000-byte planning budget and a 4,000,000,000-byte hard limit.
Estimates model allocation peaks; rejected factors are not measured RSS.

## Synthetic shape

| Measure | Value |
| --- | ---: |
| Observations | 5,938,087 |
| Contributors / FullABC atoms | 37,406 |
| Targets | 37,406 |
| Halo atoms | 0 |
| Connected components | 1 |
| Support memberships | 74,811 |
| Prepared blocks | 3,118 |
| Full A/C design columns | 74,812 |
| Projected-width columns | 37,406 |
| Full `(A, C, log B)` column count | 112,218 |

Each atom has one or two support memberships. The shape is intentionally
sparse, and repository tests use it to exercise checked preflight arithmetic.
No physical model or map is represented.

## Current rank assessment estimates

The production rank stage is `ProjectedWidthAssessment`. It contains the
free-design tiled QR and compact SVD, projection coefficients, projected-width
tiled QR and compact SVD, sparse inputs, and concurrent QR workspace. The
planning budget is 3.2 GB, so each case below is rejected before a large factor
is allocated.

| Amplitude face | Free columns | Free-design factor | Free-design SVD workspace | Projection coefficients | Projected-width factor | Projected-width SVD workspace | Estimated peak | Status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| All amplitudes interior | 74,812 | 44,774,682,752 | 134,324,048,256 | 22,387,341,376 | 11,193,670,688 | 33,581,012,064 | 266,666,536,155 | ResourceLimit; factorization NotRun |
| Half of amplitudes at the bound (`A_i = 0`) | 56,109 | 25,185,759,048 | 75,557,277,144 | 16,790,506,032 | 11,193,670,688 | 33,581,012,064 | 167,840,033,907 | ResourceLimit; factorization NotRun |
| All amplitudes at the bound (`A_i = 0`) | 37,406 | 11,193,670,688 | 33,581,012,064 | 11,193,670,688 | 11,193,670,688 | 33,581,012,064 | 91,400,872,907 | ResourceLimit; factorization NotRun |

The stage totals include these additional sparse and tile allocations:

| Amplitude face | Sparse design workset | Tiled QR workspace peak | Temporary bytes |
| --- | ---: | ---: | ---: |
| All amplitudes interior | 108,176,144 | 265,872,936,096 | 266,091,682,368 |
| Half of amplitudes at the bound | 106,979,216 | 167,049,426,200 | 167,265,180,120 |
| All amplitudes at the bound | 105,782,224 | 90,613,257,680 | 90,826,019,120 |

The sparse workset includes the free-design and raw-width sparse matrices,
normalization/index storage, triplet buffers and their construction copy. The
planner varies the free-design dimension; it does not reuse a fixed estimate
for one active face. For the all-active row, the component's projected-width
factor is already 11,193,670,688 bytes, above the full planning budget. Thus no
free-design column count can make this 37,406-width shape pass the complete
rank-stage preflight. There is no universal safe column limit: the planner also
depends on observation rows, width columns, support entries, maximum tile size,
resident inputs and projection dimensions.

The old planner estimated the independent Dense DesignSpectrum stage at
189,484,131,163 bytes and the projected-width stage at 39,065,926,939 bytes.
Those are historical planner estimates, not observed RSS. The current stage
removes the separately resident Dense DesignSpectrum allocation but accounts
for the free-design factorization that projected-width preparation already
performed. Its estimate is larger than the old projected-width estimate
because that estimate did not fully account for free-design QR/SVD work and
concurrent retained data.

## Preflight observation and remaining factors

The serial `ProductionScaleSingleComponentPreparationStaysBelowBudget` test
constructs the synthetic input, runs endpoint verification and evaluates the
rank preflight. Its current process peak RSS was 702,464,000 bytes; it did not
allocate the rejected rank factors. This RSS is for the whole test process and
is not a per-stage allocation measurement or a claim of memory saved by
rejection. RSS data is unavailable on platforms without a supported source.

For this shape, the structural row-count check sees 37,406 supported rows and
74,812 full-design columns, proving `rank(X) <= 37,406 < 74,812`. It returns a
valid full-design deficiency result before factorization. The free-face
matching check has size 37,406, but that structural matching does not prove
numerical free-design rank. Projected-width qualification is also `NotRun` after
the preflight rejection. No exact free-design or projected-width spectrum is
claimed for this fixture.

The production Dense DesignSpectrum factorization has been removed. Remaining
quadratic allocations are the free-design `p_F × p_F` QR/SVD work, the
projection-coefficient `p_F × q` matrix, and the projected-width `q × q`
QR/SVD/right-vector work, where `p_F` is the number of free A/C columns and `q`
is the number of width columns. Tiled sparse inputs and row tiles add further
workspace. Derived Target Identifiability remains linear in retained FullABC
and target records; it adds no dense factor. Excluded-halo provenance adds
O(E) identity and reason records for E exclusions, included in problem
construction memory preflight.
Covariance is independent: for 112,218 parameters its preflight estimates
418,399,065,443 bytes and rejects before factorization.

## Bounded local A/C KKT refinement

The local constrained solver has a separate, KKT-triggered dense QR refinement
for sparse unit-weight problems. It allocates only the current free-face
matrix and runs only after a full-rank SparseQR solve fails the existing KKT
certificate. The matrix-entry preflight is `rows × free_columns <= 1,000,000`,
or at most 8,000,000 bytes of double-valued matrix entries; QR workspace is
additional but remains limited to that local face. At most eight residual
correction solves are attempted, and a candidate is accepted only through the
unchanged feasibility, KKT and objective checks.

On the captured 6Z6U block, the face is 881 × 22: 19,382 doubles, or 155,056
bytes for the dense matrix. This is a local solve allocation, separate from the
projected-width rank-stage estimates above. Weighted, reference-SVD,
rank-deficient and over-limit cases do not allocate this matrix.

## Projected-width allocation lifetime

The rank stage has two sequential factor phases, but data created during
preparation remains live as the later phase runs:

| Phase | Retained inputs | Large temporary allocations | Growth |
| --- | --- | --- | --- |
| Free-design preparation | sparse free design and raw width derivative | triplet assembly, tiled QR `R`, SVD workspace, `p_F × P` coefficients | `p_F²`, `p_F P`, support nonzeros, tile rows × columns |
| Projected-width reduction | sparse inputs and coefficients | projected derivative tiles and QR assembly workspace, dense QR `R` | `P²`, tile rows × `P` |
| Projected-width SVD | projected-width `R` and singular values | SVD workspace and current full right-vector matrix | `P²` |

`free_design` and `raw` are stored sparsely in proportion to support
nonzeros. Per-tile dense row blocks can be regenerated from those sparse
matrices. The projection coefficients are retained so each width column can
be generated as a tiled `raw - free_design * coefficients` action. The
reduction streams those row tiles into a dense `P × P` QR factor. Production
currently requests the full thin right-vector matrix before copying at most
three weak directions; those directions do not affect the rank decision, but
are emitted for offline diagnostics. Process peak RSS is the whole process
high-water mark and does not attribute memory to an individual function.

Small dense-SVD controls live in `tests/support/ProjectedWidthResearch.*` and
`tests/core/joint_component/ProjectedWidthResearch_test.cpp`. They cover
full-rank, deficient, threshold-adjacent, connected, locally/strongly coupled,
weak-separator, target/halo, active-face, multiple-component and cancellation
cases. The test reports normalized design column norms, both singular spectra,
both thresholds, support counts, projection residual, elapsed time and process
peak RSS. Its dense QR/SVD oracle is for small controls; it is not a
memory-bounded algorithm. The 37,406-atom input remains a resource preflight
and structural-deficiency shape, not a full-rank control.

The planner's 3.2 GB preflight leaves margin under the 4 GB hard limit; it does
not establish that a real 37,406-atom component can be fully qualified within
4 GB. The full interior, mixed active and all-active estimates above all
exceed the planning budget, and no large rank factorization or RSS experiment
was run.

## Projected-width bounded-research preflight

The projected-width study adds small dense-oracle fixtures, exact-oracle
matrix-free action controls, a test-only sparse normal-equation projection,
and a rank-envelope contract checker. Results and topology measurements are
in [bounded projected-width research](joint-bounded-rank-research.md).

The following repeat of the current production planner uses the same giant
dimensions, 74,811 support memberships, 149,622 structural nonzeros, and an
explicit 600,000,000 resident-byte assumption. It is a checked shape
preflight; no dense matrix or factor was allocated. Values differ slightly
from the earlier sampled-resident census because that run used the then-live
resident RSS as its input.

| Free columns | Free-design factor | Free-design SVD workspace | Projection coefficients | Projected-width factor | Projected-width SVD workspace | Sparse workset | Tiled QR workspace | Estimated peak | Status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 74,812 | 44,774,682,752 | 134,324,048,256 | 22,387,341,376 | 11,193,670,688 | 33,581,012,064 | 108,176,144 | 265,872,936,096 | 266,691,682,368 | ResourceLimit |
| 56,109 | 25,185,759,048 | 75,557,277,144 | 16,790,506,032 | 11,193,670,688 | 33,581,012,064 | 106,979,216 | 167,049,426,200 | 167,865,180,120 | ResourceLimit |
| 37,406 | 11,193,670,688 | 33,581,012,064 | 11,193,670,688 | 11,193,670,688 | 33,581,012,064 | 105,782,224 | 90,613,257,680 | 91,426,019,120 | ResourceLimit |

The test-only structured projector estimates 359,052,341,856 bytes,
202,337,211,624 bytes and 90,396,764,144 bytes for those same three free
column counts. Its checked worst-case sparse Gram/fill workspace bound rejects
all three before allocation. The matrix-free action-vector preflight estimates
751,790,728 bytes, which passes the storage budget, but has no validated
bounded projector and returns numerical status `Unavailable`. Its projection
call count is zero. The preflight test's 24,494,080-byte process peak records
only the process that evaluated these shape formulas; it is not an RSS
measurement for the giant data or an estimate of memory saved by a rejected
factorization.

The full design's structural upper bound remains 37,406 for 74,812 columns,
so that design is structurally deficient. This does not supply a numerical
rank for the free face or the projected width. The giant component remains
unqualified within both memory limits; no full-rank positive control or real
4 GB run exists.
