# Joint Peeling Label-free Policy Experiment

| Item | Value |
| --- | --- |
| Date | 2026-10-20 |
| Branch | `develop` |
| Input endpoint | Fixed synthetic qualified Target A/B/C |
| Kernel/support | Production `EvaluateKernel`, 2.5 Å |
| Scenarios | 10 deterministic geometries, 2 Halos each |
| Result rows | [joint-label-free-halo-policy.csv](joint-label-free-halo-policy.csv) |
| Runner | `tests/experiments/joint_label_free_halo_policy.cpp` |

## Experimental contract

Each scenario defines Target positions and already qualified Target A/B/C, Halo
positions, Halo truth used only to create and score the simulation, and one raw
map field. The raw field is the sum of the fixed Target predictions and the
Halo truth predictions on the same 0.5 Å node grid for every candidate policy.
All candidates use that same field, geometry, contributor universe, support
radius and kernel. No candidate reruns Joint fitting.

The ten cases cover nearby similar profiles, large Target/Halo mismatch,
misleading chemical groups, equal truth across groups, dense local donors,
sparse donors, overlapping Halos, map boundaries, mixed charge signs, and Halo
widths unlike the Target distribution. Donor counts are 12 in most cases, 24
for the dense case and 4 for the sparse case.

Compared policies are the current GroupKey median with global fallback, global
Target median, kNN medians for `k = 1, 2, 4, 8`, a 2.5 Å radius with actual
available donors, 2.5 Å with global fallback on shortage, and 5 Å with actual
donors or global fallback. Local search uses the existing KD-tree; unbounded
kNN first queries nearest neighbors and then expands the radius to include all
ties at the kth distance. Results are sorted by distance and then stable Atom
ID before truncation.
Componentwise medians keep A, B and C separate. Every output is checked for
finite A/B/C, `A >= 0` and `B > 0`. An empty local neighborhood is unavailable
unless that candidate explicitly uses global fallback.

The conditional observation-driven experiment fixes each Halo width to its
geometry-only kNN-4 donor median, subtracts the same fixed Target field from the
raw map, and fits all Halo A/C values together. It solves a linear least-squares
problem with A constrained nonnegative and C unrestricted. The fit checks full
design rank with relative threshold `1e-10`; deficient cases return unavailable
samples. It does not read the simulated Halo truth. The truth is used only to
construct `raw_map_values` and compute errors. This runner's small dense design
matrix is a research fixture; it is not evidence that a production-sized
observation-by-contributor matrix meets the project memory contract.

Each numeric estimate is recomputed after reversing Target traversal order.
Every CSV row reports `deterministic=true`. Runtime is the median of seven
measurements after one warm-up. Workspace is a conservative in-process byte
estimate for the policy lookup and temporary design, not process RSS.

## Aggregate results

Signal RMSE is computed on the shared Halo-support nodes after subtracting the
fixed Target field. Parameter columns are separate A, B and C RMSEs. The
candidate rows below exclude no signal rows except where availability is
explicitly less than complete; unavailable values are written as `NA`, never
as zero.

| Policy | Halo signal RMSE mean / worst | Median per-Halo RMSE | A/B/C RMSE | Donor count / Halo | Availability | Fallback Halos | Median CPU / case | Max workspace |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Current GroupKey + global | 0.4876 / 2.2074 | 0.4732 | 1.7539 / 0.1930 / 0.8535 | 6.2 | 20/20 | 0 | 0.019 ms | 2,912 B |
| Global Target median | 0.7101 / 2.8033 | 0.7423 | 2.5109 / 0.2577 / 1.4043 | 12.4 | 20/20 | 0 | 0.013 ms | 2,912 B |
| kNN-4, unbounded | 0.4780 / 2.2011 | 0.4260 | 1.5907 / 0.1948 / 0.8339 | 4.0 | 20/20 | 0 | 0.028 ms | 2,912 B |
| kNN-8, unbounded | 0.4665 / 2.2846 | 0.4143 | 1.5964 / 0.1941 / 0.8034 | 7.6 | 20/20 | 0 | 0.033 ms | 2,912 B |
| kNN-2, 2.5 Å + global shortage fallback | 0.4901 / 2.2011 | 0.4500 | 1.6167 / 0.1904 / 0.8744 | 2.1 | 20/20 | 1 | 0.019 ms | 2,912 B |
| Observation A/C, kNN-4 width | 0.0714 / 0.5026 | 0.0643 | 1.9685 / 0.1826 / 0.1461 | 4.0 | 18/20 | 0 | 3.774 ms | 44,112 B |

The GroupKey baseline is exact in the intentionally label-aligned nearby case
and wrong in the deliberately misleading-group case (signal RMSE 1.1559). The
geometry kNN-4 estimate lowers that latter error to 0.1381. kNN-8 lowers it to
0.0331, but raises the overlap case error to 1.0639 versus 0.0840 for the
GroupKey baseline. The global median is worse than the current baseline in 7 of
10 cases and has the largest worst-case error.

The observation-driven A/C method reduces signal error in 9 of the 10 cases,
including to 0.0042 in the misleading-group case and 0.1080 when Halo widths
differ from Target widths. It returns unavailable for the overlapping-Halo
case: the design rank is 2 of 4 and 515 shared map nodes are marked unavailable.
It does not zero-fill those nodes. Its largest available signal RMSE is 0.5026
in the large Target/Halo mismatch case, where fixing B to a Target-derived
width leaves substantial bias.

## Phase 2 decision

Global Target median is not a good production candidate for these scenarios.
Geometry-only kNN is label-free and deterministic, but its modest average
improvement over the current baseline does not remove large errors in mismatch
and overlap cases. The observation-driven A/C result is promising, but this
small noiseless two-Halo fixture does not test real sample noise, the production
Peeling interpolation path, downstream Local/Group estimates, or scaling to
large overlapping contributor sets. The dense design also cannot be carried
into production under the memory contract.

Therefore no production policy is selected in this phase. Phase 3 will compare
these candidates through actual Peeling and downstream metrics, and expand
overlap, noise and contributor coverage tests. The current median algorithm
remains the production baseline until that evidence exists.
