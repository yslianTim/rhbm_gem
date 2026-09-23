# Endpoint-reference incremental acceptance

This experiment compares `e32919f3` (compact SVD with reference checks during
search) against `c9e0f8c9` (search replay with endpoint reference certification).
It does not change production algorithms, rank policies or convergence gates.

## Reproduction

Use a separate SPQR Release build of `e32919f3`, and SPQR/EIGEN Release builds
of `c9e0f8c9`. Match SYSTEM dependencies, OpenMP ON, ROOT/UMAP/Python OFF,
and enable joint extended/offline tests. Build `tests_all`, `rhbm_gem_cli` and
`joint_sparse_benchmark` for candidates; the baseline needs the latter two.
Complete regression tests and a testing-disabled build before measurement.

The runner reuses the original immutable compact campaign. It verifies its
source fingerprints against Git archives, its harness against `e32919f3`, input
hashes, raw fixed samples, audits, matrix replay and completed SQLite/JSON/CSV
exports. Historical executable hashes describe the binaries used then; current
build files are not substituted for them. The retained historical source is
`f1ac45c36945e8e89e588af5c305158a7723b268a21c5c7ab0b22a3ae676a755`;
the new candidate is
`62d178d17643a4a8195aceac17a5e4ced5af59cfb79a7e288f87bdd5d82c2de1`.

```sh
python3 tests/integration/joint_reference_validation.py \
  --prior-work-dir build/joint-compact-measurements-final \
  --baseline-spqr build/joint-reference-baseline-spqr \
  --spqr build/joint-compact-spqr --eigen build/joint-compact-eigen \
  --work-dir build/joint-reference-measurements
```

A fresh output directory is required. All workers are single-threaded and no
build or other benchmark may overlap measurement. Initialization is checked
against frozen widths before command timing. Each completed 128/512 case has
three fresh processes per version, alternating baseline and candidate. Each
analysis-to-export pipeline shares 600 seconds and a sampled 4 GiB process-tree
limit. A resource stop ends repetitions for that version/case. The entire new
experiment, including initialization checks and audits, is bounded to 80 minutes.
The monitor is sampled, not an OS hard memory limit; OS peak RSS is retained too.

After command timing, both backends audit all three compact modes on 128,
heterogeneous 168 and 512. Full spectra and derivative arrays are compared
against latest legacy mode and the corresponding historical audit. No new
fixed-state performance sample is inferred from these audit runs.

The report retains each command's phase costs, search counts, termination,
initialization, parameter/rank/active-face/evidence differences, exports and
fingerprints. Search-reference time is nested inside search time and is not
added to the phase totals. An active A is strictly positive, matching the
canonical free face. Exported `JointState.objective` is already normalized:
comparisons use its value directly, without another division by observation
scale. Historical 128 exports are rechecked with this corrected comparison.

Only complete three-sample groups produce medians and speedup ratios. Latest
128 requires coefficient/width agreement at 1e-10 and normalized objective
agreement at 1e-12. Latest 512 requires runtime convergence, persistence/export,
and resource compliance; trajectory differences are reported separately.
Candidate search-reference counts and time must be zero. Existing compact
30%/10% thresholds remain unchanged; no new 30% endpoint-speedup gate is imposed.
Missing evidence, source mismatches, failed checks and incomplete samples never
produce a passing latest gate. Preserve a failed attempt before repairing or
repeating only the affected case.

The receipt contains a hash manifest of the historical and new evidence.
Reaggregation works after relocation and needs neither the original binaries
nor Git:

```sh
python3 tests/integration/joint_reference_validation.py \
  --work-dir /path/to/extracted/evidence --report-only
```

## Results

All three gates passed: verified historical compact acceptance, latest numerical
controls, and latest Single 512 end-to-end completion. The new experiment took
2628.80 seconds (43.81 minutes), within the 80-minute budget. No new case hit a
time or memory stop. All 18 fixed-state audits agreed with latest legacy and
historical results, including spectra, ranks, active faces, coefficients,
correction, complete derivatives, gradients, KKT and trust.

Each row below uses three fresh processes per version. Time is the median of
analysis plus export; RSS is the maximum of all samples and of the OS/process-tree
measurements. RSS describes observed resident memory, not allocated bytes.

| Case | Baseline seconds | Candidate seconds | Time reduction | Baseline peak GiB | Candidate peak GiB |
|---|---:|---:|---:|---:|---:|
| Single 128 | 5.229 | 4.912 | 6.07% | 0.451 | 0.449 |
| Single 512 | 307.805 | 300.564 | 2.35% | 2.500 | 1.872 |

Every command completed SQLite persistence and JSON/CSV export and passed
runtime convergence. Across all nine baseline/candidate pairs for each size,
A/C, B and normalized objective were identical; identities, ranks, active face
and certification statuses also agreed. Each search used seven profile
evaluations and five accepted updates. Search reference evaluations decreased
from six to zero, and candidate search-reference time was zero in every sample.

| Single 512 phase (median seconds) | Baseline | Candidate |
|---|---:|---:|
| Search, inclusive | 201.865 | 195.431 |
| Search reference, nested within search | 10.114 | 0.000 |
| Endpoint assessment | 102.330 | 103.211 |
| Assembly | 1.419 | 1.391 |

The observed total reduction is 7.241 seconds; peak RSS decreased by 25.11%.
The measured improvement is incremental. Search and endpoint assessment remain
the dominant costs, at about 195 and 103 seconds respectively. Phase medians
are reported independently and must not be summed into an exact total; the
nested reference row must not be added to search. No solver or assessment
algorithm was changed by this acceptance work.

Both EIGEN and SPQR passed all 27 CTests, plus 18 runner tests and repository
guards. The testing-disabled SPQR build completed the 128 workflow with runtime
convergence and zero search-reference time. Runtime controls cover replay-only
search, endpoint rejection, saved-coefficient fallback, exhausted fallback and
preservation of existing stop reasons; frozen failure expectations were unchanged.

The [evidence index](figures/joint-reference-acceptance/README.md) links the full
historical and incremental records, exports, source proofs, test logs and six
archives. The [machine-readable summary](figures/joint-reference-acceptance/summary.json)
retains every numerical comparison. See the
[compact acceptance record](joint-component-compact-svd-acceptance.md) for the
separate original compact thresholds and historical command results.
