# Joint Label-free Peeling Production Policy

## Selected reconstruction contract

The production policy is `geometry-local-target-median-v1`. For each Halo or
eligible reconstruction-only contributor at position \(x\), take the
`min(4, N)` nearest qualified Joint Target endpoints in the same model, sorted
by Euclidean distance in Å and then Atom ID. Compute separate medians for A,
B and C. Donors must pass the existing Joint target qualification, point
validity, and local-conditioning checks. No radius limit or chemical-label
filter is used. If no qualified donors exist, reconstruction is unavailable
and reconstruction-dependent downstream stages are skipped. A shortage below
four uses all available donors.

The geometry-local snapshot records its versioned contract, fixed neighbor
count, donor eligibility/search/tie policy, Joint run ID, per-Halo donor Atom
IDs and count, source, and A/B/C. Original-Domain Halos and off-domain
contributors use the same donor index and estimator. A qualified Target keeps
its own Joint point estimate. The cache remains bounded at 4,096 reconstructed
grid nodes; donor lookup is linear in donor count plus KD-tree queries and does
not create an observations-by-contributors matrix.

## Evidence used for the choice

Phase 2 compares policies on ten fixed-endpoint synthetic cases with two Halos
per case, including misleading labels, sparse and dense donors, overlap,
boundaries, mixed charge signs and shifted widths. Signal error is scored on
the shared Halo support. The complete case-level results are in
[the Phase 2 CSV](joint-label-free-halo-policy.csv).

| Policy | Mean / worst Halo signal RMSE | Median per-Halo RMSE | Availability | Median CPU / case | Max workspace |
| --- | ---: | ---: | ---: | ---: | ---: |
| GroupKey median + global fallback | 0.4876 / 2.2074 | 0.4732 | 20/20 | 0.019 ms | 2,912 B |
| Global Target median | 0.7101 / 2.8033 | 0.7423 | 20/20 | 0.013 ms | 2,912 B |
| Geometry kNN-4 | 0.4780 / 2.2011 | 0.4260 | 20/20 | 0.028 ms | 2,912 B |
| Geometry kNN-8 | 0.4665 / 2.2846 | 0.4143 | 20/20 | 0.033 ms | 2,912 B |
| Observation A/C, kNN-4 width | 0.0714 / 0.5026 | 0.0643 | 18/20 | 3.774 ms | 44,112 B |

kNN-4 is slightly better than the old baseline on mean, median, and worst
signal RMSE in this paired synthetic sweep, with full donor coverage. Global
median is materially worse. kNN-8 has a lower mean but a worse maximum and
raises overlap-case error. Observation A/C is more accurate in available cases,
but its rank is 2/4 on overlapping Halos and it needs a dense research design;
it is not a suitable production implementation under the memory bound.

Phase 3 uses one noisy three-Target/two-Halo map and runs actual Peeling,
Post-Joint Local and Group Fitting on the same qualified endpoint. Results and
the limited scope are in [the Phase 3 report](joint-peeling-pipeline-policy-experiment-20261020.md).
The added kNN-4 row is included in the regression runner. That fixture tests
paired sample and downstream behavior but is too small to establish general
science accuracy. The observation A/C row in the original replay had lower
Halo signal RMSE (0.0199), Peeling RMSE (0.0194), and Tail Ratio MAE (0.231)
than its GroupKey baseline (0.1599, 0.0447, and 1.139), but does not remove the
overlap rank failure or dense-memory limitation.

There is no pre-existing absolute science-error threshold in the project. The
selection criterion for the current implementation is therefore paired
performance at least as close as the existing baseline in the broad sweep,
complete qualified-donor availability, deterministic results, bounded lookup,
and no change to Joint fitting. On that criterion kNN-4 is selected as the
production label-free estimator. This is evidence of baseline-level accuracy
on the tested cases, not a universal error bound. Its absolute errors remain
large in the deliberate Target/Halo mismatch and overlap stress cases; further
map-specific truth validation is still needed before claiming general
scientific accuracy.

The tests hold map values, coordinates, atom identity, roles, contributor
eligibility, endpoint, sample geometry and interpolation fixed while changing
chemical labels. They cover Halo A/B/C, original and off-domain predictions,
Peeling response, Tail Ratio, own-charge correction and local Gaussian input.
The Joint solver is not called with changed inputs and its implementation is
outside this change.

## Persistence and resource limits

New reconstruction snapshots use Joint JSON schema 18. Schema 17
`group-key-target-median-v1` snapshots remain readable and round-trip under
their saved contract. The current Peeling path rejects those legacy snapshots
if asked to recompute with them, instead of silently treating them as kNN-4.
SQLite remains schema v20; no database migration is required.

The experiment workspace estimates and cache metrics are not process RSS.
The end-to-end fixture records 148 reconstruction-only nodes, a 148-entry
cache peak, and 167,282 bytes of estimated Peeling workspace. Peak RSS for a
production-sized map has not been measured. The Phase 2 observation design is
research-only and its dense memory figure is not part of the production
workspace estimate.
