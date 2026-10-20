# Joint Peeling Policy Pipeline Experiment

## Experimental setup

The Phase 3 replay uses one deterministic synthetic map with three selected
Targets and two Halos. A single qualified Joint workflow produces the endpoint;
all policy rows then reuse that same endpoint, map, sampling geometry, and
contributor set. The map has small fixed sinusoidal noise (maximum amplitude
0.015) added before the Joint run. The target and Halo truth values are used
only to score outputs.

The experiment compares the saved GroupKey/global-fallback profiles, one global
Target median, nearest one and two Target profiles by distance, and a research
observation-driven A/C fit. The latter fixes B to the nearest-two geometry
median and fits A/C against the original-domain observations after subtracting
the fixed Joint Target field. Its 6,644 by 4 dense design uses about 208 KiB for
the design coefficients alone. This is a small-fixture research method, not a
production memory design.

For every policy, the actual Peeling path produces Target samples. The report
scores Halo profile and support-signal error, Peeling RMSE against the true
Target-only field, Tail Ratio error against the same field, Post-Joint Local
MDPDE A/B error, and the group MDPDE/posterior. It also records complete sample
coverage, effective local sample count, outlier count, reconstruction-only
node evaluations, cache peak, and the Peeling workspace estimate. The same
fixed endpoint means differences below come only from Halo reconstruction.

## Results

The full numeric table is in
[joint-label-free-peeling-pipeline.csv](joint-label-free-peeling-pipeline.csv).
Values below summarize the same rows.

| Policy | Halo A/B/C RMSE | Halo signal RMSE (median abs / max abs) | Peeling RMSE | Tail Ratio MAE | Local A/B RMSE | Group A/B RMSE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GroupKey + global fallback | 1.203 / 0.180 / 0.302 | 0.160 (0.165 / 0.306) | 0.0447 | 1.139 | 0.0622 / 0.0121 | 0.1068 |
| Global Target median | 1.203 / 0.180 / 0.302 | 0.160 (0.165 / 0.306) | 0.0447 | 1.139 | 0.0622 / 0.0121 | 0.1068 |
| Geometry kNN-1 | 1.300 / 0.250 / 0.285 | 0.154 (0.162 / 0.377) | 0.0444 | 1.130 | 0.0622 / 0.0121 | 0.1068 |
| Geometry kNN-2 | 1.237 / 0.212 / 0.251 | 0.144 (0.124 / 0.470) | 0.0336 | 0.732 | 0.0650 / 0.0132 | 0.1083 |
| Geometry kNN-4 (min(4, 3) donors) | 1.203 / 0.180 / 0.302 | 0.160 (0.165 / 0.306) | 0.0447 | 1.139 | 0.0622 / 0.0121 | 0.1068 |
| Observation A/C, kNN-2 width | 0.877 / 0.212 / 0.0047 | 0.0199 (0.00025 / 0.178) | 0.0194 | 0.231 | 0.0593 / 0.0115 | 0.1049 |

All policies had complete Peeling coverage. The fixed Joint endpoint and
100 effective samples per Target were retained. Each candidate evaluated 148
reconstruction-only nodes with a 148-entry peak cache and an estimated
167,282-byte Peeling workspace. The observation fit used a separate dense
design; its additional factorization workspace is not included in that Peeling
estimate.
The follow-up kNN-4 row is identical to the GroupKey and global rows in this
fixture because it has three qualified Targets, all in the same GroupKey, so
both kNN-4 and the historical baseline use the median of all three donors.

Observation-driven A/C had the lowest Halo signal and Peeling RMSE in this
case. It also lowered Tail Ratio error. Local A/B errors improved slightly and
the three-member group MDPDE changed little. Group alpha remained 0 because
this fixture has fewer than the ten main-chain members used for group-alpha
training; therefore the experiment does not assess group-alpha training
quality. The group posterior standard deviations were approximately 0.329 for
A and 0.078 for B under observation A/C, compared with 0.333 and 0.077 under
the baseline. No member was classified as an outlier by any policy.

This fixture deliberately has same-GroupKey Targets whose profiles differ from
the nearby Halo truth. The GroupKey median therefore has high Halo signal
error. Geometry kNN-2 improves Halo signal and Peeling error, but its largest
absolute support error is higher than the baseline. Observation A/C performs
best here, though its Halo amplitude RMSE remains 0.877; the signal metric is
more favorable because signal depends jointly on A, B and C.

## Decision

No production policy is selected from this Phase 3 fixture alone. Phase 2 tested ten
geometry and label configurations: observation A/C lowered signal error in
nine, but was unavailable for the overlapping-Halo case (design rank 2 of 4).
The one noisy end-to-end replay supports its promise, but does not establish
behavior across noisy overlap, boundary, and sparse-donor cases or production
scale. The tested observation method also uses a dense matrix, which cannot be
carried into production under the memory constraints. Geometry kNN has bounded
lookup cost and no label reads, but the ten-case sweep did not consistently
lower worst-case error.

The project gives no absolute acceptance threshold for Halo signal, Peeling,
or downstream bias. Phase 4 selected geometry kNN-4 using the broader paired
sweep, where mean, median, and worst signal RMSE were slightly below the
GroupKey baseline and all donors were available. In this three-Target case,
kNN-4 matches the baseline because both use the same three donors. This is a
baseline-level operational choice, not a certification of general scientific
accuracy. Observation-driven reconstruction still needs a bounded-memory
design and an agreed science threshold before production use.

## Limitations

- The end-to-end Local and Group comparison covers one three-Target, two-Halo
  geometry. The broad ten-case sweep in Phase 2 scores Halo prediction but not
  the full Local/Group pipeline.
- Noise is deterministic and small; it does not model a map-specific noise
  distribution.
- Group alpha training and UMAP inclusion/exclusion are not exercised in this
  three-member group case.
- The density and rank behavior of the observation fit at production scale is
  unknown.
