# 6Z6U FullABC width trajectory baseline

## Reproduction

The checked-in `joint_width_trajectory` diagnostic runs the production
FixedNeighbor search on the selected Chain A observations and records the
requested atoms' profile trials. It is test instrumentation: trajectory-only
rank, projected-width spectrum, and A/C cancellation calculations do not run
in ordinary production fits.

With the local inputs used for this capture, run:

```bash
build/bin/joint_width_trajectory \
  /Users/yslian/data/6z6u.cif \
  /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map \
  /private/tmp/6z6u-width-trajectory \
  11721,1,11650 1024
```

The output directory must not already exist. The adjacent map metadata records
an element-scaled truth width: 0.50 Å for carbon and other elements, 0.45 Å for
nitrogen, and 0.40 Å for oxygen. Atom 11721 is a Chain D / HIS 174 / C halo with
77 observation-support rows. Atom 1 is the first selected Chain A atom and is
a normal target control; atom 11650 is a second weak-response halo control.
The tool writes `width-trajectory.jsonl` for local profile trials and
`sweeps.jsonl` for full-component objective, A/C KKT, width-gradient, eta
change, and replay checks after each complete sweep.

## Baseline at source `7ec349289e3524978d36e43d3914bd96b7ca01b3`

The selected problem has 1,418 targets, 626 eligible halos, 2,044 contributors,
229,783 observation voxels, and one structural component. The trajectory for
halo 11721 was:

| Event | Width (Å) | Evidence |
| --- | ---: | --- |
| First-stage seed | 0.5107510641 | Initial estimate |
| Sweep 1, Block 21 | 0.4942385633 | Accepted; local design condition 4.75 |
| Sweep 2, Block 21 | 0.4942385633 | No width change |
| Sweep 3, Block 21 | 85,164.2929583 | Accepted `Δlog B = 12.0570745` |
| Sweep 4, Block 21 | 85,164.2929583 | Current-width A/C candidate rejected |

The Sweep 3 profile loss fell from `3.4421267642e-5` to `2.5189071623e-5`.
The LM actual and predicted reductions were `0.2682119704` and
`0.2554559517` (ratio `1.0499343181`). The full objective fell from
`0.00476824477` to `0.00475234217`, a decrease of `1.59025998e-5` (0.334% of
the pre-block objective). The accepted endpoint had local projected-width
gradient `2.1100e-17`, full A/C design rank 24/24, and condition indicator
`6.62347e10`. Its Gaussian/charge basis correlation rounded to 1.0. This is a
good LM reduction but weak width information: the ratio alone did not expose
the direction's numerical sensitivity.

At Sweep 4 / Block 21 the 930 × 24 A/C design remained rank 24/24, with free
face rank 23 and condition about `6.62347e10`. The feasible double candidate's
projected KKT was `4.6566128731e-10`, above the unchanged `1e-10` limit. The
50-digit fixed-design reference and the difference between rounded-double and
functional-kernel certification are documented in
[`joint-fullabc-width-degeneracy.md`](joint-fullabc-width-degeneracy.md).

The same conditional profile study supplies a well-conditioned moderate-width
control (`B=0.4942385633`, A/C pair condition about 1.53) and a weak small-width
halo control (atom 11650). It is not a substitute for a fitted normal-target
comparison; atom 1 is included in the reproduction command for that purpose.
The existing `block-82.json` local fixture contains selected Chain A carbon
atom 1762. Its conditional target profile can be scanned with the same bounded
high-precision reference path:

```bash
build/bin/joint_fullabc_width_study \
  tests/fixtures/joint_local_ac_kkt/block-82.json \
  /private/tmp/6z6u-target-width-profile --target 1762
```

The fixture's true carbon width is 0.50 Å in the source simulation metadata.
This scan is local to the captured Block 82 response and holds the other 11
local widths fixed; it does not qualify global target accuracy.
The bounded `joint_fixed_width_oracle` tool provides the fixed-width controls:

```bash
cmake --build build --target joint_fixed_width_oracle --parallel 4
build/bin/joint_fixed_width_oracle \
  /Users/yslian/data/6z6u.cif \
  /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map \
  /private/tmp/6z6u-fixed-width-oracle
```

It runs (1) all-contributor oracle widths with A/C block-coordinate fitting,
(2) oracle halo widths with target A/C/B fitting, and (3) first-stage halo
widths with target A/C/B fitting. It uses the sidecar's element-scaled width,
charge, and atomic-number truth. These restricted-model endpoints report their
own free-parameter stationarity and independent replay; they are never marked
as formal unconstrained Joint convergence. Per-atom estimates and the summary
are written under the fresh output directory.

## Capture fields

When enabled, each recorded profile trial includes its proposed and accepted
eta changes, local objective reduction, LM actual/predicted reductions and
ratio, local A/C rank and condition, projected-width singular values, active
nonnegative-A state, and A/C contribution cancellation ratio. Capture remains
bounded by the explicit event, atom, and row limits.

## Weak-width trial safeguard

The production profile search now performs a detailed width-Jacobian check
only when an LM proposal changes a log-width by more than `ln(10)`. It projects
the proposed step onto the right-singular subspace below the existing local
rank threshold. The step is rejected only when at least 95% lies in that weak
subspace and its local objective decrease is no larger than the existing
objective-replay enclosure (`1e-12 + 2e-12 |F_before|`). LM then shrinks its
trust-region radius and retries; the candidate does not replace the accepted
state. There is no width bound or objective change.

This rule deliberately does not veto the captured Sweep 3 update: its local
objective decrease was about `9.23e-6`, far above the replay enclosure, even
though the endpoint's projected-width information was weak. The trajectory now
records the proposed step's weak-subspace fraction and whether this narrow
replay-scale rejection fired. This preserves material descent while stopping
only large moves whose apparent improvement is within replay uncertainty.
