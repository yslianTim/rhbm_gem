# Bounded search and rank validation results

The implementation builds on `5f61bb6e`. Production remains `LegacyCompact`;
all measured searches retain the dense rank backend. The SPQR rank prototype is
an explicitly selected offline path. Only the textual results of this experiment
are retained; reproduction instructions and experiment artifacts were removed.

## Regression evidence

Both backends were built in Release with Eigen 5.0.1, SuiteSparse 7.14.1,
SPQR 4.3.6 and CHOLMOD 5.3.5 on Darwin arm64. The numerical thread environment
was fixed to one.

The joint core group contains 119 tests: SPQR passes all 119; Eigen passes 115
with four SPQR-only tests skipped. New coverage includes rank bounds, equality
and overrides, structural deficiency, resource budgets, factor identity, and
prototype operator parity with an active A face. New Python tests cover real
comparison CLI exit codes, cancellation, scheduling, fingerprints, wall-timer
completeness, rank evidence and tiny actual diagnostic-driver executions.

The selected regression groups cover joint numerics/search, runtime fixtures,
physical smoke, persistence/schema and CLI workflow. Eigen passes these groups.
SPQR has one retained runtime-fixture failure, described below; the remaining
selected groups pass. The resource watchdog test initially could not invoke
`ps` inside the sandbox and passed when run with process-monitoring permission.
The permission failure is not a numerical result.

### Pre-existing SPQR failure: active-a / first-stage-double

The same-state normalized spectrum comparison at the actual endpoint is
`1.2990450737934428e-9`, exceeding its existing `1e-10` tolerance. Rank agreement,
gradient, derivative action and local-correction checks agree. A separate pristine
`5f61bb6e` SPQR build, with the same configuration and one numerical
thread, reproduces the exact same failure and numerical values. No tolerance
was changed and no production endpoint code was modified.

This remains a failed regression gate, not a newly introduced regression and
not a waived pass. The independent bounded campaign may provide diagnostic
search/rank evidence, but this increment does not claim all acceptance gates
pass or qualify a production promotion.

## Campaign

The declared limits are five minutes fixed-step, forty-five minutes search and
ten minutes rank, with no budget transfer, restart or extension. Inputs for
single-128/512 match the prior fixed-actions campaign hashes exactly.

All 128 and heterogeneous-168 method/backend numerical comparisons pass, each
with one diagnostic repetition. They are not three-repetition performance claims.
The first 512 Schwarz attempt exceeds the 4 GiB process-tree RSS limit on both
backends, during search: Eigen after 206.709 seconds at a sampled 4,299,882,496
bytes; SPQR after 13.111 seconds at 4,418,551,808 bytes. The watchdog terminates
both, and their remaining repetitions are not run. These peaks include sampling
overshoot and are not completed-process memory totals; OS peaks are unavailable
for the terminated process groups.

There is no completed 512 Schwarz endpoint, search total, or final search-work
counter snapshot. Therefore this campaign cannot establish 512 numerical parity,
a search median ratio, or a performance-gate pass. The independent legacy controls
continue according to the frozen schedule. Resource termination is not relabeled
as numerical deficiency or convergence failure.

The single campaign finished in **1,358.614 seconds (22.64 minutes)** with exit
code **3 (incomplete)**. Fixed-step consumed 3.142 seconds, search 1,350.145 seconds,
and rank 4.809 seconds, including their preparation. There were 38 completed
processes, two RSS terminations, four skipped Schwarz repetitions and eight
skipped larger-rank/oracle processes. No budget was increased or reused.

| Backend, single-512 | Qualified legacy repetitions | Legacy search median (s) | Legacy assessment median (s) | Schwarz repetitions completed | Schwarz / legacy ratio |
|---|---:|---:|---:|---:|---|
| Eigen | 3 | 165.799679 | 51.955452 | 0 | unavailable |
| SPQR | 3 | 74.616436 | 16.812436 | 0 | unavailable |

The 1.10 performance gate remains unassessed. No cross-backend Schwarz-versus-
diagonal superiority or production promotion is claimed.

## Rank prototype evidence

Small deterministic tests obtain both certified full-rank and deficient decisions
and compare them with the dense oracle. The real tiny SPQR diagnostic driver also
obtains a full-rank certificate. However, **neither representative 128-atom
campaign prototype obtains a decision**: both exhaust the configured 100 million
charged sparse-entry work units before finishing the streamed reconstruction
certificate. This is an incomplete scalability acceptance, not a pass based on
all-unavailable results.

| SPQR prototype | Rows / free columns | Decision | Rank wall, including factor (s) | Charged work units | Sampled process-tree peak RSS (bytes) | Exported factor bytes |
|---|---|---|---:|---:|---:|---:|
| chain-128 | 58,300 / 256 | unavailable: rank-work-budget | 0.161702 | 99,995,408 | 64,143,360 | 10,610,096 |
| cube-128 | 47,792 / 256 | unavailable: rank-work-budget | 0.121975 | 99,987,038 | 133,464,064 | 55,111,568 |

Separate Eigen and SPQR dense-oracle processes report rank 256 for both inputs.
Prototype instrumentation records zero compact extractions and zero free-design
SVDs. The estimated additional rank workspaces are 2,847,552 and 2,343,168 bytes;
these exclude the design and sparse factor. Instrumentation is not an allocation
trace. The incomplete reconstruction error was unavailable, not zero.

Per the stopping rule, neither topology advances to 512 or 2,000 atoms. The next
rank work package must address certificate cost and calibrate a new, explicitly
budgeted experiment; these results do not establish linear scaling or 2k
capability. The next search work package should first locate the 512 Schwarz
memory growth before proposing another performance campaign.

## Post-measurement guard correction

A final boundary regression reproduced a structural-column scan returning a
decision after its rank time budget. The implementation now checks the deadline
per input column and before returning any decision. The new regression fails
before the fix and passes afterward. Both backends were rebuilt, and the joint
core, fixed-runner and bounded-runner groups passed again. The 119-test counts
above describe this final version.

The campaign was **not repeated** after this correction. Its measurements
therefore describe the version before the final deadline guard fix. The
correction does not change dense search, rank formulas or resource limits.
