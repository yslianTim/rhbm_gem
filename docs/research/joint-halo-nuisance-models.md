# Joint halo nuisance parameterization study

## Controlled models

`tests/experiments/joint_fixed_width_oracle.cpp` compares four fits on the
same selected Chain A observation rows and the same 2,044-contributor problem:

| Run | Target widths | Halo widths | Role |
| --- | --- | --- | --- |
| `full-abc` | Estimated | Estimated | Model A baseline, initialized from first-stage widths |
| `oracle-halo-widths` | Estimated | Fixed to simulation truth | Model B oracle control |
| `seed-halo-widths` | Estimated | Fixed to first-stage estimates | Model B deployable-seed control |
| `all-widths-oracle` | Fixed to simulation truth | Fixed to simulation truth | Upper-bound A/C fitting control |

The Model B runs start target widths from the same first-stage initialization
as the FullABC baseline. The oracle control takes each contributor's effective
width, charge, and atomic number from the map's `.simulation.json` sidecar;
it does not assume one width for all elements. All runs keep the selected
observation domain, objective, nonnegative A constraint, signed C parameter,
core size 12, forward serial block order, and 200-sweep cap unchanged.

The diagnostic disables endpoint certification only to bound the experiment
to the fit comparison. It separately computes the full-model A/C KKT, the
width gradient over free width parameters, and independent replay. Its
`search_stationarity_qualified` field describes that restricted experimental
parameter space; `formal_joint_converged` is always false. Per-atom CSVs and
`summary.json` contain target and halo A/C/B errors, target and halo prediction
contribution norms, objective, replay, block activity, recovery counts,
estimated peak memory, and wall time. Process peak RSS is a shared high-water
mark across sequential cases, so it is not an isolated per-model measurement.
QR factorization counts are unavailable in this diagnostic.

## Adaptive Model C

No adaptive halo parameterization is promoted or fitted in this change. The
current local conditioning and projected-width spectra can identify weak
directions for a particular fitted block, but the repository has no validated
threshold that classifies halo widths across components before fitting. An
arbitrary threshold would make a data-dependent model choice and could shift
target A/C/B estimates. A future Model C study needs a predeclared rule,
independent simulation cases for calibration, explicit free-parameter
stationarity, and target-bias measurements against the same fixed observation
domain.

## Omitted halo signal

The controlled pruning study in
[`joint-halo-pruning-validation.md`](joint-halo-pruning-validation.md)
already varies omitted weak, strong, and overlapping low-support halo signals
without changing observation rows. It shows that omitted signal can be
absorbed as target error or leave the fit unqualified, depending on its size
and overlap. These examples do not establish omitted halo signal as the cause
of the 6Z6U width runaway. The archived threshold sweep predates the bounded
block-recovery change; final validation reruns the executable on the updated
source.

## Promotion criteria

Fixed-width halo models remain experimental. A production proposal would need
to improve target width and A/C accuracy as well as target prediction error,
while reporting the full observation objective, free-space stationarity,
independent replay, failure rate, runtime, and memory. Lower halo parameter
error alone is not sufficient evidence to change the default model.
