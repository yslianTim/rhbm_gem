# Measured acceptance — PR0 / PR1 / PR3 contracts

Measured on macOS arm64, Release, one numerical thread, against pristine `15e71e3fc230205e6efd7f60c8a915eb1989733b`. Full receipts (including source, input and binary hashes, build configuration, matrix probes, factor fill, process limits and statuses) are preserved beside this report as compressed JSON.

## Verdict

- All 18 fixed-state audits passed; primary/reference coefficients and objectives were unchanged. Both physical runtime comparisons passed.
- Joint tests: Eigen 100 passed / 3 existing backend-specific skips; SPQR 103 passed. All 11 new operator/contract tests passed on both backends. Existing evidence, runtime runner, runtime fixtures, physical and CLI round-trip CTest groups passed; runner comparison tests passed (4 cases). Logs are retained beside this report.
- All 20 preparation-only workloads completed, with identical input hashes across backends, one connected component, complete 515-membership supports per atom, and zero factorization, reference or SVD work counters.
- 6Z6U hit the sampled RSS watchdog in both pristine and candidate runs on both backends. It has **no numerical verdict**. The overall comparison correctly exits 2 (incomplete real-data control), not success.
- The original fixed-state driver cannot decode null halo widths in 6Z6U; those cases are explicitly not run and are covered by the public runtime control. No replacement widths were invented.

## Operator medians

Three independent processes per row. Times include dedicated factor preparation; rank time is a subset of prepare. RSS covers the combined legacy oracle/operator audit, not isolated operator storage.

| Backend / case | Prepare s | Rank s | Apply s | Adjoint s | Max relative action error | Audit peak GiB |
|---|---:|---:|---:|---:|---:|---:|
| eigen/single-128 | 0.693 | 0.012 | 0.0257 | 0.0256 | 1.32e-13 | 0.419 |
| eigen/heterogeneous-168 | 2.597 | 0.019 | 0.0964 | 0.0974 | 2.52e-13 | 0.641 |
| eigen/single-512 | 37.784 | 0.280 | 0.3656 | 0.3616 | 5.31e-13 | 2.085 |
| spqr/single-128 | 0.174 | 0.154 | 0.0171 | 0.0170 | 4.42e-15 | 0.313 |
| spqr/heterogeneous-168 | 0.829 | 0.720 | 0.0677 | 0.0669 | 6.44e-15 | 0.473 |
| spqr/single-512 | 17.270 | 16.170 | 0.5384 | 0.5332 | 5.54e-15 | 1.772 |

## 10k preparation only

| Backend / layout | Rows | Memberships | Process wall s | OS peak GiB |
|---|---:|---:|---:|---:|
| eigen/chain-10000 | 4,550,060 | 5,150,000 | 3.869 | 1.664 |
| eigen/cube-10000 | 3,434,036 | 5,150,000 | 3.341 | 1.396 |
| spqr/chain-10000 | 4,550,060 | 5,150,000 | 3.882 | 1.669 |
| spqr/cube-10000 | 3,434,036 | 5,150,000 | 3.345 | 1.400 |

No large A/C solve, compact rank, profile operator, search, reference, assessment or uncertainty was executed. Sparse basis derivative entries are preparation data, not the full profile Jacobian.

## Real-data resource stops

| Version / backend | Status | Wall s | Sampled tree peak GiB |
|---|---|---:|---:|
| baseline/eigen | rss-limit | 569.35 | 4.608 |
| baseline/spqr | rss-limit | 121.52 | 4.139 |
| candidate/eigen | rss-limit | 122.07 | 4.069 |
| candidate/spqr | rss-limit | 126.85 | 4.037 |

The 4 GiB watchdog samples process-tree RSS and can overshoot between samples; it is not an OS hard cap. Missing output and missing post-kill OS peaks remain missing. Budgets were not increased.

## Scope and provenance

The candidate cube-128 workflow completed for each backend and records basis, A/C, derivative, search, reference, assessment, assembly and uncertainty. Inclusive times overlap and must not be summed. Shape probes are known allocations, not a complete allocator trace.

The successful candidate campaign followed an interrupted preflight with redundant Eigen Q-transpose compact construction. Final measurements use extraction of the existing permuted R. The original baseline and interrupted preflight remain under `/private/tmp/joint-pr01-final`; final raw files are under `/private/tmp/joint-pr01-final-r2`. `provenance.json.gz` identifies the reused pristine baseline. Temporary build/input trees are not committed; recipes and fingerprints are retained.

This validates the full fixed-state operator and preconditioner interface foundation. The transient p-by-p rank check remains. It does not establish 10k search or full-workflow scalability; the 6Z6U numerical gate remains unverified.
