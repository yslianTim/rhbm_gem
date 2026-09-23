# Observable-halo experiment conclusions

The 2026-09-23–24 experiment started from `d634fa0d`. Singleton-halo profiling
removed the diagnosed initial A/C obstruction, but **6Z6U did not pass complete
correctness or performance acceptance**. Testing was stopped at the user's
request. This page preserves conclusions only; experiment-specific runners,
diagnostic extensions, additional timers, raw measurements, source archives and
temporary builds were removed afterward. The cleanup was not rebuilt or retested.
See the [production contract](joint-observable-halo.md) for retained behavior.

## 6Z6U findings

The original model and 0.50 A simulated map matched the previously recorded
hashes. The immutable problem contained 1,559 targets, 2,192 contributors,
262,801 rows and 866,216 memberships, in components of 2,167 and 25 atoms.

The large component had the previously diagnosed 20 singleton halos. Across
both components there were **27 singleton halos on 26 distinct rows**, leaving
2,165 FullABC atoms and 262,775 informative rows. Halos sharing a row were grouped
without assigning individual A/B/C estimates.

Both initial constrained A/C solves had valid KKT certificates: free-face ranks
were **4,194 and 33**, with 100 and 3 active A boundaries respectively. All 2,165
FullABC seeds were valid first-stage fits; median fallback was not needed for this
dataset. Reduced halo seeds were recorded as not required.

There were also **17 remaining FullABC halos with exactly two informative rows**;
seven had positive A at the initial solution. When their A/C columns span those
two rows, their width derivative lies in the same span and cannot be identified
separately. A valid initial A/C solve therefore does not establish ABC
identifiability. Active A boundaries impose an additional uncertainty limitation.
This does not prove that every selected target is individually unidentifiable,
or establish the rank of an unobserved final endpoint. These atoms were not
recursively downgraded and no rank threshold was relaxed.

## Resource outcomes

SPQR Release used one numerical worker with the existing selection,
normalization and solver budgets. The correctness budget covered the complete
analysis-to-export flow: **3,600 seconds / 8 GiB**. The run stopped at 3,600.02
seconds, with sampled peak RSS **4.29 GiB**, inside outer search. Diagnostic stacks
included derivative preparation, SPQR Q operations, cancellation fallback tiled QR
and compact SVD. No final endpoint, runtime-convergence acceptance, uncertainty,
peeling, SQLite save or export was completed. The later reviewed/optimized source
was not subjected to another complete 60-minute run.

The reviewed performance baseline hit the **600-second / 4-GiB** gate at 600.03
seconds, with sampled peak RSS 3.70 GiB. Its subsequent repetitions were cancelled.
The optimized 6Z6U attempt was interrupted by the user before completion. Its last
printed progress was 180 seconds / 3,198 MiB; these are incomplete, rounded
observations, not final resource measurements. No 6Z6U median or speedup is claimed.

## Completed performance comparisons

Each version completed three serial fresh-process analysis-to-export runs for
Single 128 and Single 512, under 600 seconds / 4 GiB per run.

| Case | Reviewed stage-3 median | Optimized median | Speedup | Baseline / optimized maximum RSS |
| --- | ---: | ---: | ---: | ---: |
| Single 128 | 4.68 s | 3.34 s | 1.40x | 0.608 / 0.551 GiB |
| Single 512 | 289.79 s | 122.93 s | 2.36x | 2.964 / 3.020 GiB |

All completed runs passed runtime convergence, scalar replay, persistence and
export checks. Paired comparisons preserved ranks, active faces, availability and
convergence states. A/C/B, prediction, peeling and normalized objective were
identical in these cases. The largest covariance relative difference was
`7.98e-13`, below `1e-10`. These are historical measurements made before removal
of experiment instrumentation, not new measurements of the cleaned source.

Single 512 endpoint assessment fell from about **67.01 to 17.54 seconds**, and
uncertainty from **142.39 to 25.01 seconds**. Search remained about 78.5 seconds;
derivative preparation remained about 61.6 seconds and Jacobian reduction about
29 seconds. These numerical sub-timings overlap workflow phases and cannot be
added to them. Persistence took about 0.16 seconds and export 0.21 seconds.
The main measured gain was endpoint/uncertainty SVD work; complete-flow peak
memory did not decrease for Single 512. Planned isolated SVD/QR measurements
were not run, so individual optimization qualification remained incomplete.

## Regression scope and unresolved checks

Before cleanup, SPQR passed 25/27 CTest entries and EIGEN 26/27. The new 11-case
observable-profile suite and two numerical equivalence cases passed. The failures
were independently confirmed on the original baseline: missing historical
sparse-acceptance documentation attachments and, for SPQR, the existing
`active-a / first-stage-double` normalized-spectrum discrepancy of
`1.2990450737934428e-09`. Rank agreed; no tolerance was relaxed.

The Python validation/watchdog suite passed 19 tests before its experiment-only
comparison test was removed. Python bindings were not built. The final installed
consumer smoke expectation was updated to JSON v4 but not executed. Permanent
small regression tests remain, including the two-row identifiability
counterexample. Full 6Z6U correctness, performance and the unfinished checks must
not be reported as complete.
