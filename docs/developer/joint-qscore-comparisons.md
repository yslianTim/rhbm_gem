# Joint-peeled QScore comparisons

The sequence summary compares three atom-wise Pearson correlations. Each residue
is the arithmetic mean of its available non-hydrogen atom scores within its own
chain. Blue retains the stored Standard QScore and its published endpoint
eligibility; Green and Red require the same available paired record. Halos are
contributors, never Target outputs. Painter selection changes aggregation only.

| Curve | Observations | Reference |
| --- | --- | --- |
| Blue: Standard QScore | Raw map at Standard QScore points | Stored height/offset and width 0.6 Å |
| Green: Joint-Peeled QScore (standard reference) | Joint neighbor contributions subtracted from those map observations | Same Standard reference as Blue |
| Red: Joint-Peeled QScore (Joint model reference) | Exactly the same vector as Green | Target endpoint's normalized Gaussian plus charge kernel |

Blue versus Green isolates neighbor peeling. Green versus Red isolates reference
choice. Neither comparison implies an ordering of scores. Blue can include more
atoms than the paired curves. The summary produces plots only, without an
auxiliary atom-count TSV. Missing values are omitted, never zero-filled; residues
without valid scores are omitted entirely. Residues are sorted by sequence ID.
Plot limits include finite correlations, including negative values, within
[-1,1]. Empty curves are skipped and an all-empty summary is reported.

## Sampling and arithmetic

`PrepareQScoreSamples` is shared with Standard QScore. It repeats the interpolated
center eight times, then visits shells 0.1 through 2.0 Å in 0.1 Å increments.
`GetRadialPointsForQScore` retains deterministic spiral-sphere sampling, the
existing non-hydrogen neighbor exclusion, up to 50 attempts, and **all** accepted
points from the first successful attempt. Empty shells remain empty. Map values
use the existing trilinear interpolation, including its boundary handling.

Each retained sample has equal weight, including the eight center copies.
Standard reference evaluation retains scalar, once-per-shell arithmetic and
`height + offset` at the center. All scores use the original mean-subtracted
Pearson calculation. New comparisons diagnose empty/mismatched samples,
nonfinite values, zero observation variance, zero reference variance, and
nonfinite correlation separately. The legacy Standard missing-data behavior
remains unchanged.

## Sample-centered neighbor peeling

For each sample x of Target i, the peeled observation is

```
y_peel(x) = trilinear_map(x) - sum(j != i, ||x-r_j||² <= 6.25) f_j(||x-r_j||)
f_j(r)   = A_j G(r; tau_j) + C_j K(r; tau_j)
```

A single 4.5 Å broad-phase neighbor query per Target covers the 2.0 Å sampling
radius plus 2.5 Å contribution radius; a tiny broad-phase roundoff guard cannot
expand the exact inclusive squared-distance support. Hydrogen and self are
excluded. Contributors are evaluated analytically in serial-ID order, with a
long-double neighbor sum. Joint prediction maps are not interpolated. Neither
Original-grid peeling nor own-charge correction is applied to these samples.

Other Targets must match the same run's snapshot endpoint, method, role,
component, identity, and qualification. Recorded Halos use the snapshot's
`halo_median_reconstruction`, checked against qualified geometry-local donors.
For needed Halos outside the recorded Original Domain, the existing kNN-4 Target
median policy selects donors, including its deterministic distance/ID tie-break,
and reconstructs A, tau, C. One donor index and a contributor model cache serve
the entire comparison pass. Halo representatives are never published as Joint
optimizer endpoints and never written back into the Joint snapshot. Grouping
parameters and charge medians do not enter this calculation.

Both neighbor and Target references use `GaussianModel3D::ResponseAtDistance`,
including normalization and stable origin treatment. C multiplies the charge
kernel; it is not a constant offset. Any necessary missing or invalid contributor
makes the pair unavailable. No per-point filtering can produce different masks
for Green and Red.

## Analysis, ownership, and persistence

The qualified Joint publication path computes comparisons while the original
map and problem are available, before local/group inference and SQLite save.
It verifies Original row identities/observations and map geometry. Results are
owned by `LocalPotentialEntry`, read through `AtomLocalPotentialView`, and
validated/published in one batch by `ModelAnalysisEditor`. Painter only reads
these records and does not load or interpolate a map.

Joint-peeled QScore is optional derived analysis. Only its builder and batch
publisher share a best-effort exception boundary, after qualified Joint
publication and map binding. A recoverable `std::exception` clears all QScore
records, logs one warning, and lets Original-grid Peeling, Post-Joint Local and
GroupingStage proceed under their existing policies. `std::bad_alloc` is
re-thrown. Joint solving, Halo reconstruction, Joint publication/map binding,
contributor/context construction and core downstream fitting remain outside
this boundary. In particular, unavailable formal Halo reconstruction still
prevents downstream processing and preserves `JointDownstreamUnavailable`.

The publisher validates the entire batch before updating records. QScore
metrics accumulate locally and commit only after publication succeeds, without
replacing core metrics. A batch failure leaves `post_joint_qscore = null` for
all atoms, including any previous run's results; it does not manufacture zero
scores or per-atom Unavailable records. A normally completed builder can still
publish individual `EvidenceStatus::Unavailable` records with explicit reasons.
Failure handling does not change the formal Blue/Green/Red definitions or the
shared sampling points, observations, peeling responses, order and weights.

SQLite schema **25** stores Stage Result JSON **9**; Joint Result JSON stays
**25**. `post_joint_qscore` is null when not run or when batch computation or
publication fails, otherwise contains:

- calculation contract `joint-peeled-qscore-v1` and complete Joint source;
- one optional pair `[standard, joint]`, availability status and failure reason;
- stored Standard reference height/offset, one shared sample count;
- contributor/missing counts, evaluated contribution count and Halo coverage;
- extra-domain Halo identities and qualified donor identities.

No sampling positions or peeled responses are persisted. Validation rejects
contradictory availability/reason fields, nonfinite/out-of-range scores,
impossible sample/coverage counts, endpoint/source/reference mismatches, and
extra-Halo donors that differ from the qualified geometry-local selection. Recorded Halo provenance comes from the
same immutable Joint snapshot. An unavailable record may explicitly report a
missing Halo reconstruction; it cannot carry scores. Old SQLite/JSON versions
are rejected without migration or guessed reconstruction.

Endpoint/run/snapshot replacement clears derived results through existing Joint
invalidation. Geometry edits clear analysis, and bound map revisions invalidate
analysis on access. Changing Standard reference parameters clears only these
comparisons, as do changes in Hydrogen eligibility. Selection changes preserve records. After SQLite reload, the saved
snapshot, endpoints, reference, and source identities are validated together.

## Cost and verification

Runtime metrics report Targets, available/unavailable pairs, samples, evaluated
neighbor contributions, reconstructed Halo/extra-Halo coverage, elapsed seconds,
and an estimated transient/result payload. The pass retains samples for one
Target and compact result records, not all Target sample vectors. Map/problem
storage and allocator overhead are outside the payload estimate; process RSS
also includes fitting, ROOT, and the test runtime.

`QScoreHelper_test.cpp` checks sampling order, fixed counts and representative
positions, plus six Standard score goldens captured from the independent legacy
scalar oracle (absolute tolerance 2e-14). `PostJointQScore_test.cpp` covers
independent scalar peeling/kernel oracles, shared observations, support boundaries, Hydrogen,
Target/Halo coverage, missing/nonfinite data, map boundaries, sorted averages,
negative display ranges. Persistence checks cover paired round-trip, unavailable
reasons, corrupt sources/counts,
extra-Halo donors and invalidation. Full CTest additionally covers unchanged
Joint endpoints, Original-grid peeling, Post-Joint Local and Grouping workflows.
