# Complete-command resource envelope

## Measured environment and inputs

The measured host is an Apple M1 iMac, 8 cores, 16 GiB RAM, macOS 26.6.2 arm64.
Release builds use one numerical worker. Every analysis has
`--exclude-hydrogen=true --asymmetry=false --only-backbone=false
--map-normalization=false`. Map normalization divisor is verified as 1 in exported
outcomes. Analysis and export run serially in fresh processes and databases.

The user model is `/Users/yslian/data/6z6u.cif`, SHA-256
`349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa`.
The map is `sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map`, SHA-256
`6d7771a9fd584e562550da473d624307b22d20396f11a2a79a064bff2551682a`.
Both match its simulation sidecar. It is float32 storage, 258³ voxels at 0.5 Å,
origin (53.5, 53.5, 53.5) Å. The model contains 70,190 loaded atoms, of which
37,406 are non-hydrogen catalogue atoms. `--asymmetry=false` retains the first
chain of each entity (A, Y, BA, YA), giving 1,559 targets.

The production builder confirms 633 halo atoms, 262,801 rows, 866,216 memberships,
and two components of 2,167 and 25 contributors. The model copy includes hydrogen
atoms even though they are excluded from contributors, so its full-model cost
must not be estimated from the non-hydrogen count alone.

## Controlled scale cases and measurement

Catalogue cases keep two target atoms and 681 rows / 1,025 memberships fixed,
while remote non-hydrogen catalogue atoms grow to 1,000 / 10,000 / 37,406. All
maps have the same 258³ grid. The connected-size series has 128 / 512 / 2,167
atoms on a 3.5 Å lattice. The multiple-component control has 2,167 total atoms,
at most 32 per block, with 20 Å block translations. The actual row/membership
counts and partition, not a center-distance approximation, are verified by the
production builder before each case.

Wall time includes the real `potential_analysis` command and subsequent
`result_dump`: import, input hashing, Q-score/reference Gaussian, construction,
initialization, fitting, assessment, assembly, SQLite storage and JSON/CSV export.
Input generation and census are separate from command timings but charged to the
resource-stage budget. Existing phase costs are retained; search-reference time
is already a subset of search time and is not added twice. The remainder outside
construction/initialization/search/assessment/assembly is reported as a combined
command overhead, not falsely assigned to a single preprocessing step.

The 4 GiB watchdog samples the process tree at a nominal 100 ms interval and
terminates it when the sampled threshold is crossed. Actual sampling gaps and
OS-reported per-process peak RSS are retained. Sampling is not a hard OS memory
limit and can miss short peaks; detected post-exit peak violations also disqualify
a run. The analysis/export pair shares one 600-second deadline. A process limit,
unavailable numerical state and a completed-but-nonconverged result are different
outcomes. File-cache state is uncontrolled; these are local sequential runs, not
cold-cache or cross-platform benchmarks.

## Preflight findings retained

The [preflight records](figures/joint-validation/resource-preflight.json) retain
three catalogue-case failures caused by the old SQLite identity check, which
required results to contain every non-hydrogen atom rather than just actual
contributors. Unit and CLI regressions reproduced the failure before the minimal
subset-membership fix and pass afterwards. Unknown/hydrogen/duplicate IDs still
fail, saved results round-trip, and invalid replacement rolls back. Neither the
solver nor the JSON/SQLite schema changed.

A preflight 1.5 Å lattice with 512 atoms failed initialization and did no search;
its short runtime is not solver-scale evidence. That outcome remains recorded.
The fixed 3.5 Å lattice series is the final factorization workload. The preflight
also exposed a census-label bug (counting all loaded atoms as eligible); the
corrected census reports both counts. Its interrupted first attempt consumed
28.51 seconds, retained in the same 80-minute C allocation. Final records are
kept separately, not substituted for the preflight evidence.

## Results

See the [machine-readable measurements](figures/joint-validation/resources.json)
for every repetition, input/build fingerprints, actual census, phase costs,
resource stops and numerical statuses. Only complete analysis **and export**
constitute an end-to-end completed run. Such completion does not imply a usable
state or passed runtime convergence.

| Case | Rows / memberships | Component sizes | Completed / attempted | Wall seconds | Maximum sampled tree RSS, MiB | Numerical result |
| --- | --- | --- | --- | --- | --- | --- |
| Catalogue 1,000; two contributors | 681 / 1,025 | 2 | 3 / 3 | 0.872–0.947 | 223 | passed, all three |
| Catalogue 10,000; two contributors | 681 / 1,025 | 2 | 3 / 3 | 1.201–1.262 | 232 | passed, all three |
| Catalogue 37,406; two contributors | 681 / 1,025 | 2 | 3 / 3 | 2.181–3.140 | 282 | passed, all three |
| Single 128 | 51,939 / 65,920 | 128 | 3 / 3 | 17.693–19.740 | 368 | passed, all three |
| Single 512 | 186,694 / 263,680 | 512 | 0 / 1 | 600.030, time limit | 759 | no exported result |
| Single 2,167 | unmeasured | unmeasured | 0 / 0 | skipped after 512 limit | unmeasured | unmeasured |
| Multiple, total 2,167 | 856,379 / 1,116,005 | 67 × 32 + 23 | 0 / 1 | 600.032, time limit | 1,570 | no exported result |
| User 6Z6U, 2,192 contributors | 262,801 / 866,216 | 2,167 + 25 | 0 / 1 | 600.037, time limit | 1,956 | no exported result |

The three time-limited cases produced no database or exported numerical outcome;
their numerical convergence is unknown, not a numerical failure or success
inferred from termination. None crossed the sampled memory limit. Killing the
process group prevented `/usr/bin/time` from returning its final OS peak for
these cases; only the observed tree peak is available. Sampling gaps, rather
than an assumed exact 100 ms cadence, are recorded per process.

The largest **completed single component** here is 128, with 51,939 rows. The
37,406-atom catalogue result applies to a two-contributor domain; it does not
establish completion for the same catalogue's full selected domain. The supplied
6Z6U workload is outside the demonstrated 10-minute envelope. Multiple small
components also do not establish full-command scalability at 2,167 total atoms.
Those cases have different row/membership counts and are not a controlled
estimate of component-count effects alone.

For the completed 128 case, search took 13.45–15.16 seconds, assessment
3.01–3.15 seconds and initialization 0.18–0.19 seconds. For the largest catalogue
control, construction took 0.216–0.227 seconds; fitting remained about 0.002
seconds, while 1.94–2.89 seconds lay outside estimator phase costs. Output sizes
were approximately 0.23 / 0.99 / 3.35 MB for the catalogue controls and 3.10 MB
for the 128 case. Three local repetitions show the observed range only, not
stable general performance or a maximum supported problem size.

## Bounded bottleneck sampling

Three separate diagnostic executions used the same input hashes and CLI binary,
with one-second macOS stack samples at 5 / 15 / 25 seconds and an intentional
stop at 32 seconds. Their [receipts and selected frames](figures/joint-validation/diagnostic-profiles.json)
retain hashes of the complete local stack files. They are excluded from the
uninstrumented timing samples and charged to C's allocation.

- Single 512: samples at 5 and 15 seconds are in reference QR during trial trust
  checks; the 25-second sample is in the associated Jacobi SVD. Independent
  reference factorization is an observed cost, not just sparse primary fitting.
- User 6Z6U: the five-second window includes production alpha/MDPDE initialization;
  15 and 25 seconds are in sparse QR within the initial profiled linear solve of
  component search. The large connected solve is an observed bottleneck.
- Multiple 2,167: all three early windows show component search and assessment,
  including QR, derivative and spectrum work. They do **not** establish where
  the rest of the 600-second run was spent. The code also performs global
  assembly/assessment, but assigning the timeout primarily to that stage would
  require later-phase evidence absent from these bounded samples.

Thus neither a universal complexity law nor a full phase breakdown is claimed
for a terminated command. Larger single-component cases were skipped as planned;
none of the capped cases was retried to completion, given extra solver budget,
split into artificial production components, or run with relaxed checks.

C consumed 2,002.37 seconds (33.37 minutes), including preflight and 96.26 seconds
of separate profiling, within its 80-minute allocation. Unspent time did not
expand the measurement matrix. Time-capped cases establish a limit of the tested
workflow on this host, not that a longer or differently engineered run can never
finish.
