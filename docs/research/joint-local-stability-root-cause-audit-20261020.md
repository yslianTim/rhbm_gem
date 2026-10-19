# Joint local stability root-cause audit

## Scope and source state

This audit starts from `develop` at `a8f5eac81142231ddd36408427b933bc9ba18e4d`.
The checkout and its cached `origin/develop` ref were equal (`0 0`). A live
GitHub HEAD check could not resolve `github.com` in this environment, so the
remote's current HEAD is unverified.

The audit adds measurements only. It does not change any solver, global
qualification, bound, replay, or convergence threshold. The saved bounded
runner exercises the same FirstStage calibration and bounded FixedNeighbor
entry used by the production candidate. Conditional A/C reprofiling fixes the
accepted global log-widths and calls the existing constrained coefficient
solver; its state is diagnostic and is never copied into the accepted Joint
state.

## Criteria found in the current source

`CertifyLinear` marks a local A/C state strictly qualified at projected KKT
`<= 1e-10`. Local Parameter Stability separately requires its Width gradient
`<= 1e-12`. The accepted bounded Joint endpoint uses practical A/C KKT
`<= 1e-7` and projected Width KKT `<= 1e-9`, with its existing parameter
change checks. The local Width value is computed by the same bound-aware
`ProfileGradientInfinityNorm(evaluation, context)` function, so at a bound it
is the projected KKT value rather than the unconstrained gradient.

The local diagnostics are not runtime solver gates. `JointFitResult::RuntimeConvergence`
and `TargetRuntimeConvergence` consume numerical evidence and target parameter
validity, not `JointBlockParameterStability`. An unavailable stability
diagnostic marks affected `point_quality` as `Warning`. The downstream
`QualifiedTarget` check requires passed runtime statuses, valid parameters and
adequate local conditioning; it does not check parameter stability. Therefore
this diagnostic status does not disqualify a target from the current
downstream path.

## 6Z6U simulation control

The Release bounded rerun selected block `1/block-117` from 138 blocks. It
reported:

| Measurement | Original global slice | Fixed-width local A/C reprofile |
| --- | ---: | ---: |
| Local A/C projected KKT | 2.3282681338443467e-11 | 2.42861286636753e-17 |
| Projected Width KKT | 1.432187701766452e-13 | 3.808064974464287e-14 |
| Local objective | 4.389101300502815e-18 | difference 3.037728304414033e-22 |
| Scaled A/C coefficient difference | — | 3.778658854419188e-11 |
| Relative prediction difference | — | 2.834674209554505e-10 |
| Width upper-bound slack in log-width | 1.5155294861319157 | same widths |
| Active A constraints | 0 | — |
| Independent replay | Passed | Passed |

The selected block was strictly qualified before reprofiling. The local
stability assessment remains `Stable` with six perturbation solves and 24
profile evaluations (0.04931 process CPU seconds for that block). The run
retained the original candidate values: objective `1.2330854370811315e-13`,
global A/C KKT `9.434574186999167e-8`, projected Width KKT
`8.776510069988752e-10`, 83 sweeps, 22,908 profile evaluations, and
`RuntimeConvergence=Passed`. Independent replay passed. Process Peak RSS was
`637,943,808 B`.

## 6Z6U unsharpened failure block

The Release bounded rerun selected the previously reported failure, block
`1/block-5` of 140. Its first failing check was `baseline.certificate.kkt_passed`,
which tests local A/C KKT at `<= 1e-10`:

| Measurement | Original global slice | Fixed-width local A/C reprofile |
| --- | ---: | ---: |
| Local A/C projected KKT | 3.5208046721135755e-10 | 1.3877787807814457e-17 |
| Projected Width KKT | 5.4610760358286825e-12 | 7.245815059064853e-12 |
| Local objective | 2.863055030028584e-4 | difference 4.1257784491995526e-19 |
| Scaled A/C coefficient difference | — | 1.2363583329392307e-9 |
| Relative prediction difference | — | 1.1157094118308093e-8 |
| Minimum Width bound slack in the block | 0 | 0 |
| Active A constraints | 0 | — |
| Independent replay under practical criteria | Passed | Passed |

The baseline A/C KKT is 3.52 times the local strict threshold, while remaining
284 times below the practical `1e-7` threshold. The baseline projected Width
KKT is 5.46 times the local `1e-12` threshold, while remaining 183 times below
the practical `1e-9` threshold. The block contains at least one Width-bound
active target (zero minimum slack); the complete endpoint has two active-bound
targets. Exact local active-bound multiplicity was not emitted by this audit
record.

The reported failure is therefore not a constrained-solver failure: it is an
accepted global slice failing a stricter post-fit local A/C certificate, and
that check runs before the similarly over-strict local Width test. The
conditional reprofile improves local A/C stationarity. It changes normalized
objective by `4.13e-19`; its relative prediction difference `1.12e-8` is close
to the existing local warning level. The formal state remains the original
global slice, so the conditional optimum is retained only as diagnostic
evidence. The test report below confirms the global endpoint remains
qualified: objective `0.018417539368903946`, global A/C KKT
`8.174157334426418e-8`, projected Width KKT `9.402726930396454e-10`, two
active-bound targets, independent replay passed, and
`RuntimeConvergence=Passed`. The component completed 124 sweeps and 34,848
profile evaluations; process Peak RSS was `1,426,980,864 B` and wall time was
`435.395722208 s`.

Both 6Z6U maps support using the original global slice as the diagnostic
baseline and the same practical local criteria as global qualification.
Conditional reprofile remains useful as an audit comparison, but using it as
the stability baseline would measure a nearby conditional optimum rather than
the accepted Joint estimate.

## Reproduction and test evidence

The bounded no-instrumentation runner was built with
`cmake --build build --target rhbm_tests joint_width_repro_no_instrumentation --parallel 4`.
The focused CTest command used one internal numerical thread and
`--parallel 4`; all 8 selected tests passed, covering local numerical tests,
bounded synthetic, persistence, workflow, cleanup, and regression contracts.
The fitted state remains unchanged by the added local audit calculations; the
aggregate candidate values match the retained B4 records to reported precision.
