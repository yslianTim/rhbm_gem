# Joint Component Assessment Memory Census

This census describes the current production stages and a deterministic
single-component synthetic shape. It is not a physical cryo-EM support graph.
The memory planner uses a 3,200,000,000-byte (3.2 GB) budget and a
4,000,000,000-byte hard limit. Estimates include stage-resident data and
temporary allocations; rejected factors are not observed RSS.

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
| Design columns | 74,812 |
| Projected-width columns | 37,406 |
| Full `(A, C, log B)` column count | 112,218 |

Each atom has one or two support memberships. The shape is intentionally
sparse, and repository tests use it to exercise checked preflight arithmetic.
No physical model or map is represented.

## Current stage estimates

| Stage | Estimated peak bytes | Outcome |
| --- | ---: | --- |
| FixedNeighbor preparation | 910,595,445 | Pass |
| Global state evaluation | 972,989,100 | Pass |
| Endpoint verification | 625,051,683 | Pass |
| Structural design check | 578,519,527 | Pass; row-count rank obstruction |
| DesignSpectrum | 189,484,131,163 | ResourceLimit; NotRun |
| Projected-width assessment | 39,065,926,939 | ResourceLimit; NotRun |
| Derived Target Identifiability | about 482,240,000 | Pass; linear evidence only |
| Independent target covariance preflight | about 418,400,000,000 | ResourceLimit; factorization NotRun |

Derived Target Identifiability retains existing component rank evidence, atom
identities, active-amplitude checks, target/halo mapping and per-atom statuses.
It does not allocate a full Jacobian factor, right singular vectors or null-space
workspace. Its rank-derivation core is `O(P + T)` for `P` FullABC atoms and
`T` targets. Per-atom `Unavailable` records for `U` observable-contribution-only
nuisance atoms make total additional evidence storage `O(P + T + U)`; the
estimate also includes resident input and mapping storage. A target-memory
preflight pass does not establish that either rank stage ran or that target
evidence passed. In this synthetic case the structural check proves Design Rank
deficiency, so no trusted rank evidence is available to derive a target pass.

The covariance estimate is separate. It includes tiled Jacobian assembly,
factor and right-vector storage, and compact-SVD workspace. At 112,218
parameters, the modeled factor and right-vector allocations are each
100,743,036,192 bytes. The covariance preflight rejects this case before any
tiles or factors are allocated. No covariance factorization ran, and no RSS
savings are inferred from its estimated peak.

The retired Target Geometry planner estimated 317,422,079,259 bytes for this
shape and rejected it before allocation. That planner estimate is a historical
baseline, not an observed memory saving. The new covariance estimate has a
separate purpose and is not target-identifiability work.

The synthetic structural check counts 37,406 unique supported rows against
74,812 full-design columns, proving `rank(X) <= 37,406 < 74,812` without
forming the dense design factor. If the obstruction is absent, DesignSpectrum
still uses its dense `ReferenceQR` rank factorization. The projected-width
stage still retains its width-space factors and tiled workspace. Both remain
above the 3.2 GB budget for this shape and are the current rank-memory
bottlenecks.

## RSS and scope

The serial scalability regression records current RSS before and after each
preflight and process peak RSS where the platform supports it. RSS is a process
high-water mark, not an allocation-level stage measurement. A preflight
rejection means the corresponding dense factor was not allocated; its estimate
cannot be reported as observed RSS or actual memory saved.

The giant fixture is synthetic evidence only. A complete qualification on a
real cryo-EM model/map has not been established. Windows/MSVC RSS and build
behavior are not measured by this census.
