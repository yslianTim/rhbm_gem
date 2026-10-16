# 6Z6U Local A/C KKT Cost Profile

Date: 2026-10-16

Repository revision: `e1b8cb7c25c2043f29e2246c914c0959444e2cc3`

Branch: `develop`

## Question and method

This study checks whether the production Local A/C KKT limit of `1e-10`
triggers extra work during the Global Practical 6Z6U chain-A run. It uses the
same 12-atom forward Gauss-Seidel blocks, one trusted update per visit,
FirstStage A/B fixed halos with C=0, 200-sweep limit, normalized map, and
Global Practical criteria. The diagnostic adds aggregate counters only; it
does not retain per-trial matrices or coefficients.

Input SHA-256 values:

| Input | SHA-256 |
| --- | --- |
| `/Users/yslian/data/6z6u.cif` | `349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa` |
| `sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map` | `6d7771a9fd584e562550da473d624307b22d20396f11a2a79a064bff2551682a` |
| Matching `.simulation.json` | `33a3c221ec3f2fff18dcd2a0ea6484f68334971fc688f7cd100e220ad5848931` |

The map normalization divisor was `0.44031545732594374`. The run used the
repository's `joint_fixed_halo_validation` research executable with
`production practical --cli-normalization`.

Reproduction command (write the output JSON to a local scratch path):

```sh
build/bin/joint_fixed_halo_validation \
  /Users/yslian/data/6z6u.cif \
  /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map \
  /private/tmp/joint-local-strict-p0.json production practical --cli-normalization
```

## P0 results

| Measure | Observed |
| --- | ---: |
| Observations / targets / fixed halos | 229,783 / 1,418 / 659 |
| Sweeps / blocks per sweep / total block visits | 79 / 142 / 11,218 |
| Profile evaluations / accepted updates / failed blocks | 22,439 / 10,785 / 0 |
| Local linear solve calls / valid / invalid | 22,439 / 22,439 / 0 |
| Rank-deficient solves / active-set releases | 0 / 0 |
| Initial KKT mean / maximum | `2.05148e-17` / `9.02056e-17` |
| Initial KKT at or below `1e-10` / `1e-8` / `1e-7` | 22,439 / 22,439 / 22,439 |
| Final KKT at or below `1e-10` / `1e-8` / `1e-7` | 22,439 / 22,439 / 22,439 |
| Bounded refinement eligible / triggered | 22,439 / 0 |
| Refinement iterations / successes / failures | 0 / 0 / 0 |
| Refinement QR factorizations / wall time | 0 / 0 s |
| Local profile KKT rejections / replay KKT rejections | 0 / 0 |
| Local replay disagreements / invalid gradients / width safeguard rejections | 0 / 0 / 0 |
| Local solver time / total search time | 11.31 s / 339.23 s |
| Total run wall time / peak RSS | 340.48 s / 633,339,904 bytes (604 MiB) |

The strict local threshold did not reject any profile or cause a refinement.
Relaxing it to `1e-8` or `1e-7` would therefore save **0 seconds of refinement
time on this observed trajectory**. Local linear solve time was about 3.3% of
search time; the measurable cost attributable to the strict KKT threshold was
zero.

The final global result remained qualified: 1,418/1,418 target estimates,
Global A/C KKT `9.197902e-8`, width-gradient infinity norm `4.231091e-10`,
complete-sweep log-width change `5.666148e-7`, independent replay passed, and
full-observation replay passed. Candidate objective was `0.00410208212`.

All 142 reported local free A/C block designs were full rank (24/24) and
classified `WellConditioned`; the largest reported free A/C condition
indicator was `7.01366`.

## Local LM stopping

Of 11,218 block visits, 10,040 stopped at the one-accepted-update budget and
1,178 stopped in Eigen LM. The native LM status counts were:

| Status | Count | Interpretation |
| --- | ---: | --- |
| `RelativeReductionTooSmall` (code 1) | 1,171 | `ftol` reduction test |
| `RelativeErrorAndReductionTooSmall` (code 3) | 7 | Combined `ftol` and `xtol` test |
| `CosinusTooSmall` (code 4) | 0 | `gtol` test |

No block stopped on `gtol`; most visits stopped at the existing accepted-update
budget. This run gives no evidence to loosen LM tolerances. `ftol` ended some
searches, but this task did not test a change because local KKT qualification
was not limiting the run.

## Promotion decision

Keep production Local A/C KKT at **`1e-10`**. Do not promote `1e-8` or `1e-7`,
and do not change LM stopping tolerances. P1/P2 threshold plumbing and six
additional production/Oracle runs are not justified by this baseline: every
observed local solve and replay already met the strict threshold, so a looser
gate cannot improve the measured path.

There is also a relevant stress case outside this final trajectory:
`tests/fixtures/joint_local_ac_kkt/block-21-sweep-4.json` records a full-rank,
ill-conditioned local design (condition indicator `6.6235e10`) with projected
KKT `4.6566e-10`. The existing strict regression keeps that candidate rejected.
It was not encountered by this 79-sweep production run, and this P0 study does
not establish that accepting such a candidate under `1e-8` or `1e-7` preserves
parameter stability. There is no measured speed benefit to justify taking that
additional risk.

The Global Practical criteria, replay consistency limits, and production
estimator remain unchanged. The report makes no claim about Fixed-Halo C=0
model error or recovery accuracy beyond the recorded run metrics.
