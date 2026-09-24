# PR2 / PR3 validation receipt

The production default remains `SearchMethod::LegacyCompact`. The baseline
campaign reached its original 80-minute deadline before completing every
required comparison. This alone prevents promotion, independently of later
candidate timing or numerical results. No budget was increased or campaign
restarted.

## Frozen implementations and tests

The baseline production source is pristine
`c6c869cd4f4939104dfde96984ae17f736a8ce4a`, built separately for Eigen and SPQR.
The recorded environment is Darwin arm64, Eigen 5.0.1, SuiteSparse 7.14.1,
SPQR 4.3.6 and CHOLMOD 5.3.5, with one numerical thread. Complete compiler,
linker, cache and binary fingerprints accompany the receipts.
The candidate adds operator Guarded LM, identity/diagonal controls, deterministic
128-core one-ring partitions, local Schwarz and resource diagnostics. Public
API, A/C, reference, endpoint evidence and uncertainty remain unchanged.

Before candidate measurement, review tightened the radius guard to require a
step inside the radius, capped damping overrides at 20, and added per-attempt
regularization records and mapping checks. The candidate was rebuilt; eight
new operator-search tests and all eight selected CTest regression groups passed
on both backends. Six Python comparison/preflight tests passed. Candidate
measurement did not start with stale binaries. The earlier baseline measurement
wrapper differs only in code excluded by `PR23_BASELINE_DRIVER`.

The joint core group ran 111 tests: SPQR passed all 111; Eigen passed 108 with
the three existing SPQR-only workspace/derivative tests skipped. All eight new
operator-search tests passed on both backends.

Regression coverage includes joint numerics, PR1 operator parity/cancellation,
observable profile, partial selection, runtime fixtures, endpoint certification,
uncertainty, persistence and CLI round-trip. Kernel tests establish correctness
on small systems; they do not establish full-search endpoint parity at every
benchmark size.

## Baseline campaign

Each process/pipeline used one numerical thread, at most 600 seconds and 4 GiB
sampled process-tree RSS. The campaign deadline was 80 minutes. There were 41
completed numerical runs, one campaign-deadline termination, and 16 fixture
entries not run after the deadline. Both 6Z6U controls were not run for the same
reason. Every completed baseline endpoint passed its existing convergence
checks.

| Case | Eigen legacy search median (s) | SPQR legacy search median (s) |
|---|---:|---:|
| single-128 | 5.826 | 2.968 |
| heterogeneous-168 | 20.757 | 12.296 |
| single-512 | 229.232 | 217.933 |
| chain-128 | 8.331 | 1.944 |
| cube-128 | 4.140 | 3.489 |
| chain-512 | 371.783 | 49.669 |
| cube-512 | 213.037 | Incomplete: two of three |

Only three completed independent processes produce a reported median. Search
time excludes endpoint assessment; full process time and RSS include it. The
SPQR cube-512 third process was terminated by the remaining campaign time,
not by a fresh 600-second allowance.

## Candidate and large-local campaigns

The user stopped all remaining measurements on 2026-09-24. The SPQR campaign
was interrupted during its first single-512 diagonal process; its three Schwarz
repetitions are complete. Large-local never started. These are user-cancelled
or not-run results, not numerical passes or resource failures. The Eigen campaign ended at
its original deadline with 11 completed processes, one deadline termination
(the third heterogeneous-168 Schwarz repetition), and 41 not-run-budget case
entries. Its 6Z6U control was also not run. No continuation campaign was opened.

Eigen single-512 Schwarz completed three independent
processes and passed endpoint parity. Its search median is 419.562 seconds,
1.830 times the legacy median, so the 1.10 performance gate fails.
Eigen diagonal also completed three processes with passed endpoint convergence;
its median is 460.460 seconds. Schwarz is 8.88% faster on this case, which has
four blocks and a maximum of 216 atoms per block (128-core partitions therefore
have actual structural overlap). This does not offset the failed legacy gate.
Identity completed three processes with passed endpoint parity and a median of
456.175 seconds. All nine single-512 candidate processes passed scientific
parity against the corresponding legacy repetitions. The two completed
heterogeneous-168 Schwarz endpoints passed convergence, but two repetitions do
not establish the required median/comparison.

SPQR single-512 Schwarz also completed all three repetitions with passed
scientific parity. Its median is 350.104 seconds, 1.606 times legacy: the 1.10
performance gate therefore fails on both backends. The first SPQR repetition
spent 215.365 seconds in PCG and 105.821 in operator preparation, including
98.667 in the retained compact rank check. Local matrices and shifted factors
together took 0.186 seconds. These inclusive phase measurements must not be
summed with their child J/J-adjoint or rank timers.

No further campaign was launched. SPQR diagonal/identity comparisons and the
remaining search cases are incomplete. A Schwarz advantage over diagonal on
the same overlapping case across both backends has therefore not been established.
The 2k/5k/10k supplied-state local-build campaign was not run; no validation of
those sizes or global search capability is claimed.

## Final acceptance status

| Deliverable | Status |
|---|---|
| Operator search | Small-system kernels and regressions pass on both backends; single-512 Schwarz endpoint parity passes on both. Full search acceptance matrix is incomplete. |
| Schwarz | Small-system topology, mapping, regularization, SPD and lifetime tests pass. Large-local validation and the cross-backend diagonal performance comparison are incomplete. |
| Production default | **Not switched.** Both single-512 legacy cost gates fail, and required comparisons remain incomplete. |
| 128 / 168 search | Baseline measurements exist; candidate coverage is incomplete. Eigen heterogeneous-168 has two completed Schwarz repetitions. |
| 512 search | Single-512 Schwarz passes endpoint parity on both backends. Eigen identity/diagonal also pass; SPQR controls are incomplete. Frozen chain/cube candidate searches are not validated. |
| 2k / 5k / 10k | Large-local campaign not run at user request; global search remains outside this batch's scope. |
| 6Z6U | Not run in this batch; earlier resource termination is not a numerical pass. |

The [audit archive](figures/joint-operator-search/receipts.tar.gz) preserves raw
campaign receipts, process logs, frozen source patches, environment fingerprints,
test logs and a separate cancellation receipt. Original partial campaign files
are preserved without marking them complete. The adjacent
[archive manifest](figures/joint-operator-search/manifest.json) records its SHA-256.

The inherited `factor_nonzeros_upper_bound` probe describes R fill, not complete
Q/Householder storage. Exact global factor bytes are unavailable from that
probe. Process RSS includes those factors; local Schwarz storage bounds are
reported separately. Neither the probe nor sparse storage alone establishes
linear global factor memory growth.

See [implementation and reproduction](joint-operator-search.md) for the exact
contracts, commands and promotion gates.
