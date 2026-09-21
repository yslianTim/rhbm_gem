# Joint estimator capabilities and limitations

The supported estimator is fixed-position, equal-weight joint LS with A ≥ 0,
B > 0, signed C and fixed 2.5 Å support. Selected non-hydrogen atoms fix the
observation domain; every eligible catalogue contributor intersecting that domain
is estimated once. Closure is complete for the supplied catalogue and this support
model, not for unknown atoms or physical contributions outside the truncation.

## What an outcome establishes

The production initializer operates on a full model copy and can read beyond the
fitting domain. Global command preprocessing also reads the full map/model.
Search completion, a usable state, runtime convergence and an offline certificate
are separate facts. A failed check stays failed. Missing components prevent a
complete prediction/objective; other component states can still be retained.
An available estimate or a small RSS does not establish parameter accuracy,
statistical uncertainty, a global optimum or real-data scientific validity.

JSON schema 3 / SQLite v17 persist all contributors and their target/halo roles.
Contributor IDs must belong to non-hydrogen atoms in the saved model; the result
may be a strict subset of the full catalogue. Foreign, hydrogen and duplicate IDs
remain invalid. Export reads the saved outcome without rerunning the estimator.
A command exit code of zero means saving succeeded, not that numerical convergence
passed. A/C are input-map scale only after applying the recorded normalization
divisor; B is in Å. No absolute charge calibration is implied.

## Evidence and applicable limits

- [Weak-halo attribution](joint-weak-halo-attribution.md): one frozen problem,
  three starts, actual-state comparisons, bounded restarts and independent
  precision controls. Weak sensitivity and nonconvergence remain explicit.
- [Noise and position mismatch](joint-noise-mismatch-validation.md): two geometries,
  20 fixed noise seeds per condition, separate roles and initialization modes.
  Available estimates frequently fail current convergence checks; failed results
  are included in the error analysis. There is no general noise-tolerance claim.
- [Complete-command resource envelope](joint-command-resource-envelope.md): local
  Apple M1 / 16 GiB, one worker, 10-minute / 4-GiB per-pipeline limits, and an
  80-minute resource-stage allocation within a two-hour experimental budget.
  An observed limit is not a declared maximum supported size on other machines.

The input 6Z6U map is a hash-verified simulation. It is not experimental-map pilot
evidence. A future real-data claim requires a named dataset, fixed selection and
geometry, initialization/normalization provenance, and independent-map or half-map
comparisons where available. Agreement alone is not atom-wise truth, and shared
initialization can make apparently separate fits dependent.

## Reproduce the bounded campaign

The default installed library/API does not link multiprecision audit code.
Configure a Release test build with `RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON`, then
build `joint_validation`. Python requires NumPy; resource monitoring uses `ps`
and `/usr/bin/time`. On macOS, sandboxed invocations need permission to inspect
child-process RSS. No external inputs are downloaded.

```sh
cmake --build build/joint-v1-offline --target joint_validation -j 2
python3 tests/integration/joint_validation.py --stage a \
  --work-dir build/joint-validation-new --executable build/joint-v1-offline/bin/joint_validation
python3 tests/integration/joint_validation.py --stage b \
  --work-dir build/joint-validation-new --executable build/joint-v1-offline/bin/joint_validation
python3 tests/integration/joint_validation.py --stage c \
  --work-dir build/joint-validation-new --executable build/joint-v1-offline/bin/joint_validation \
  --cli build/joint-v1-default/bin/RHBM-GEM \
  --model /Users/yslian/data/6z6u.cif \
  --map /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map
python3 tests/integration/joint_validation_report.py \
  --work-dir build/joint-validation-new --output build/joint-validation-new/report
```

Use a fresh work directory. The ledger refuses silent stage reruns; the explicit
`--resume-stage --stage-dir NEW_DIRECTORY` recovery option requires a fresh directory, preserves the prior
attempt and subtracts its consumed time. It is intended for diagnosed harness or
implementation defects, not retries of cases exceeding their limits. Large inputs,
raw outputs, databases and full traces remain under the work directory; compact
scientific records and their [hash manifest](figures/joint-validation/manifest.json)
are checked in here. The retained campaign includes an explicitly recorded
preflight before the persistence and census fixes; a fresh campaign uses the
corrected tools directly.

## Validation and provenance

The baseline is commit `6b5b020721e64ac3d99333635dace3d193de0ea6`.
The original A/B attempt ran before the persistence fix (production source
SHA-256 `e346827a…`); final A/B and command measurements use the fixed build
(`89f96f22…`). The
[campaign ledger](figures/joint-validation/campaign.json) and individual outcomes
retain full source, configuration, executable and build hashes rather than
pretending that every stage used one binary. The only production-source change
is contributor-subset validation in SQLite storage. Public APIs, schema 3,
SQLite v17, objective, kernels, optimizer and convergence policies are unchanged.

The [10 related runtime/data/CLI checks](figures/joint-validation/related-tests.txt)
and [four offline/extended checks](figures/joint-validation/offline-tests.txt)
passed, including frozen regressions and high-precision controls. The added
persistence unit and expanded CLI smoke both failed before the fix and passed
afterwards; their receipts and logs are in the evidence manifest. These test
passes verify contracts and regression expectations, including expected numerical
failures; they do not change the failed outcomes in the experimental reports.

[Eight experiment-runner tests](figures/joint-validation/runner-tests.txt) also passed: fixed seeds and matrix, variance
matching, fixed selection under mismatch, missing values and statistical
denominators, failure intervals, expired budgets, and watchdog termination.
The report verifies frozen weak-input hashes and exact equality between each
audited A/C/log-B state and its formal outcome. Report generation does not modify
any runtime outcome. Archived executed tool sources remain in the raw output
directory. An exact-state assertion found a one-ULP coefficient change in the
experiment wrapper's default JSON parser. Production decoding already used
precise parsing. The wrapper now also uses precise parsing and verifies A/C and
log-B against the actual runtime state before saving. Original A/B outputs remain
preserved; A/B were repeated in separate directories after this tool fix, charged
to their original allocations. The runner's explicit recovery option was extended
to these stages; no experimental condition, solver policy or tolerance changed.
The [recovery comparison](figures/joint-validation/serialization-recovery.json)
confirms identical input hashes and all 450 numerical statuses, with at most
one ULP (`4.44e-16`) difference in parameter representation. The final
[joint/CLI smoke rerun](figures/joint-validation/final-joint-cli-tests.txt) also passed.

## Delivery decision

| Work item | Status | Meaning |
| --- | --- | --- |
| P1-A | Complete with limitations | Two available weak endpoints support coexisting search, precision and sensitivity explanations. The third start has no trusted state; its finer cause remains undistinguished. No solver defect or converged restart was established. |
| P1-B | Complete with limitations | All 450 fits recorded; 65 passed convergence, 385 failed with explicit evidence. Noise/mismatch behavior is interpretable for this matrix, not a general accuracy guarantee. |
| P1-C | Complete with limitations | Three successful repetitions for the catalogue controls and single 128. Single 512, multiple 2,167 and the user case hit the time limit; single 2,167 was deliberately unmeasured. Late-phase attribution of the multiple-component timeout remains unavailable. |
| P2-D | Complete | Reproduction, evidence, regression validation, applicable limits and this delivery decision are recorded. No real cryo-EM accuracy claim is made. |

**Deliver the bounded validation and persistence fix; do not declare a generally
complete or noise-qualified estimator, or support for this full 6Z6U workload
within the requested resource cap.** Explicit numerical/resource failures remain
limitations, not new success conditions. Any future general-noise, global-optimum,
large-domain performance or real-data claim is blocked until it has its own
evidence. Including recovery runs, the campaign used 76.37 / 101.57 / 2,002.37
seconds for A / B / C, respectively: 36.34 minutes in total. The saved regression
checks add approximately 5.6 minutes; even counting them, the total is below
43 minutes and the two-hour budget. Compilation and documentation are excluded.
