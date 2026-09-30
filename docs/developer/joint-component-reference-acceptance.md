# Endpoint-reference incremental acceptance

Status: Historical. `joint_reference_validation.py` was retired. Current search
measurements are documented in [`joint-benchmark.md`](joint-benchmark.md), while
independent certification remains in the offline reference tests.

This experiment compares `e32919f3` (compact SVD with reference checks during
search) against `c9e0f8c9` (search replay with endpoint reference certification).
It does not change production algorithms, rank policies or convergence gates.

## Historical provenance

The original campaign compared `e32919f3` with `c9e0f8c9`; its retained source
and harness fingerprints are recorded in the historical report. The archive
member layout, former paths, source and retrieval commits, compressed sizes, and
SHA-256 values are in the [Joint evidence archive inventory](figures/joint-reference-acceptance/archives.json).
Retrieve a blob with `git show <last_present_commit>:<former_path>`. Historical
report reaggregation requires the runner from Git history; current measurements
use the benchmark profile documented in `joint-benchmark.md`.

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

The [Joint evidence guide](joint-component-evidence.md) summarizes retained
comparisons and retrieval instructions. The [machine-readable summary](figures/joint-reference-acceptance/summary.json)
retains every numerical comparison. See the
[compact acceptance record](joint-component-compact-svd-acceptance.md) for the
separate original compact thresholds and historical command results.
