# Fixed-state factor / normal-action acceptance

Production still defaults to `SearchMethod::LegacyCompact`. Implementation and
bounded numerical validation are complete. These measurements establish observed
fixed-state/fixed-step improvements, not full-search promotion or scalability.

## Validation and provenance

The baseline is `ec24f3056427408cb1978cff004ed477908fb0d0` with the archived
instrumentation-only patch. A uses its factor and composed action; B uses the new
factor with composed action; C uses the new factor and ApplyNormal. Both backends
load the same baseline-Eigen frozen beta/eta, free face and parent scale/rank
context. All production-source, wrapper, input, state and binary fingerprints
are retained. The measured production source matches the final production source.

Both backend CTest runs passed all nine selected groups, including joint core,
endpoint/reference, observable/partial selection, uncertainty, persistence,
physical runtime and CLI coverage. Joint core: SPQR 113 passed; Eigen 110 passed
with three pre-existing SPQR-only skips. The dedicated fixed comparator suite
passes five tests; the strengthened search comparator suite passes nine tests.
The latter includes post-measurement tightening of missing-evidence rejection;
its shared finite-value helper is unchanged. Post-measurement timing bounds do not alter numerical comparisons.

The one campaign completed 72/72 processes and 24/24 three-repetition A/B, B/C,
A/C comparisons. Every action, spectrum/rank, same-state control and all three
fixed-preconditioner step comparisons passed their declared tolerances. Normal
actions used exactly two Q/Q' calls versus six for the composed control. There
were zero reference solves or derivative preparations in the measured paths.

Elapsed campaign time: 1190.516 seconds (19.84 minutes), below the 30-minute limit. No process exceeded 300 seconds or 4 GiB; no restart or budget extension occurred.

## Measured medians

Seconds; three independent processes per cell. Recorded step subtotals include
preparation, metric, applicable partition/preconditioner construction, solve and
prediction. The original diagnostic omitted gradient setup from these subtotals.
They also exclude input reconstruction and audit serialization. Exact complete
step timing is therefore unverified in this campaign. The final driver separately
times gradient setup and includes it; no campaign was rerun. Timers are nested;
phase medians must not be summed to derive totals.

| Case | Backend | Prepare A → C | Normal A → C | Schwarz subtotal A → C | Subtotal reduction |
|---|---|---:|---:|---:|---:|
| single-128 | EIGEN | 0.6908 → 0.6872 | 0.0508 → 0.0168 | 0.9970 → 0.8103 | 18.73% |
| single-128 | SPQR | 0.1942 → 0.0475 | 0.0358 → 0.0145 | 0.4055 → 0.1567 | 61.35% |
| heterogeneous-168 | EIGEN | 2.7679 → 2.6632 | 0.1938 → 0.0644 | 4.5748 → 3.3660 | 26.42% |
| heterogeneous-168 | SPQR | 0.8913 → 0.1761 | 0.1385 → 0.0568 | 2.1377 → 0.8038 | 62.40% |
| single-512 | EIGEN | 40.2868 → 39.8260 | 0.7630 → 0.2489 | 45.3320 → 41.7601 | 7.88% |
| single-512 | SPQR | 18.8407 → 1.7208 | 1.1113 → 0.4517 | 26.0317 → 5.3369 | 79.50% |

For every displayed case/backend, normal actions and all identity/diagonal/
Schwarz subtotal comparisons have a lower median and all three paired times are
faster. Preparation alone does not consistently improve for Eigen single-128
or single-512. Very small chain-8 timing is retained as a diagnostic control,
not used to claim representative workload speed.

## Post-measurement timing correction

The raw receipts are unchanged. For the sum of tracked numerical phases plus
the omitted gradient **adjoint action**, the recorded subtotal is a lower bound;
adding all recorded adjoint action time is an upper bound. In C that aggregate
contains the audit adjoint and gradient adjoint. This deliberately overcounts
candidate work. It does not bound vector-expression or other uninstrumented
setup overhead and is not complete wall-clock step timing.

For Schwarz, all three candidate upper bounds are below their paired baseline
lower bounds; the same median direction holds:

| Case | Backend | Median A lower bound (s) | Median C upper bound (s) | Consistent improvement of tracked phases |
|---|---|---:|---:|---|
| single-128 | EIGEN | 0.9970 | 0.8607 | yes |
| single-128 | SPQR | 0.4055 | 0.1998 | yes |
| heterogeneous-168 | EIGEN | 4.5748 | 3.5587 | yes |
| heterogeneous-168 | SPQR | 2.1377 | 0.9667 | yes |
| single-512 | EIGEN | 45.3320 | 42.4880 | yes |
| single-512 | SPQR | 26.0317 | 6.6913 | yes |

These are conservative bounds, not corrected point estimates. Per-preconditioner
bounds and inconclusive factor-only comparisons are retained in comparison.json.
Exact gradient-inclusive complete step totals require a future campaign; this
one was closed within its original budget. The corrected driver was rebuilt
and functionally checked separately, without adding performance repetitions.

## Attribution and resource tradeoff

- SPQR single-128: factor-only preparation median 0.1942 → 0.0426 s; normal-action median on the new factor 0.0427 → 0.0145 s.
- SPQR heterogeneous-168: factor-only preparation median 0.8913 → 0.1721 s; normal-action median on the new factor 0.1650 → 0.0568 s.
- SPQR single-512: factor-only preparation median 18.8407 → 1.7793 s; normal-action median on the new factor 1.3720 → 0.4517 s.

SPQR single-512 compact extraction alone fell from 17.313357 to 0.000560 seconds
(A/C medians); that former cost was primarily extraction, not SVD. The public
Householder representation makes composed actions slower in some B controls,
so factor acquisition and normal-action changes must be assessed separately.
Eigen single-512 preparation remains about 40 seconds and dominates its total.

The exported Householder representation has a memory tradeoff. Sampled
single-512 process-tree RSS ranges (GiB) are:

| Backend | A | B | C |
|---|---:|---:|---:|
| EIGEN | 1.233–1.596 | 1.168–1.378 | 1.246–1.796 |
| SPQR | 0.787–0.957 | 1.469–1.621 | 1.274–1.641 |

The maximum sampled RSS across all measured runs was 1.796 GiB. OS peaks and
sampling gaps are retained separately. Exported-array byte counts exclude the
retained design and construction scratch; R fill is not total factor storage.
The p-by-p compact, endpoint dense reductions and full uncertainty remain.
The evidence supports proceeding with bounded factor/rank work; it does not
justify a larger automatic operator cutoff or 2k/5k/10k workflow claims.

## Reproduction and retained evidence

See the [implementation/reproduction guide](joint-fixed-actions.md),
[compact summary](figures/joint-fixed-actions/summary.json), and
[full comparison decisions](figures/joint-fixed-actions/comparison.json).
The [archive manifest](figures/joint-fixed-actions/manifest.json) identifies the
[raw receipts and logs](figures/joint-fixed-actions/receipts.tar.xz).
The archive contains campaign receipts, frozen states, per-process output,
resource records, measured wrappers/source patch and regression logs. Original
PR0-PR3 receipts, cancellation records and failed promotion gates are unchanged.
