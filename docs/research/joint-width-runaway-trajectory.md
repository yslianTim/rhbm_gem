# Joint width runaway trajectory

## Scope and method

This report reconstructs the supplied 6Z6U run through the first FixedNeighbor
failure. It uses the original CIF and 0.50 Å map, chain A, FullABC
initialization, core size 12, forward block order, and the production 200-sweep
limit. The run stops at the existing component 1 / sweep 4 / block 21 KKT
failure. No width-search policy or accepted state was changed for this study.

The opt-in `joint_width_trajectory` audit executable records selected atom IDs
or all atoms (`all`) into a new output directory. Trial history is limited to
the selected atoms, a configurable event cap, 2,048 atoms, and 4,096 local rows.
The diagnostic is compiled only with joint test instrumentation. It evaluates
the bounded local state for recorded trials and does not alter the solver state.
For example:

```sh
build/joint-validation/tests/joint_width_trajectory \
  /Users/yslian/data/6z6u.cif \
  /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map \
  /private/tmp/joint-width-trajectory-run 11721,11650 1024
```

`profile_objective_reduction` is the change between successive profile
evaluations. The first profile evaluation includes the block-entry A/C
re-profile. `accepted_local_objective_improvement` measures the accepted
profile candidate against the local block-entry response, while
`accepted_global_objective_improvement` is the original global objective
change across the whole accepted block update. The latter can include both
the entry A/C re-profile and the LM width update. Null fields mean that no
finite qualified endpoint was available.

## Run result

The all-atom accepted-update capture contained all 620 accepted updates across
the run; it reached neither the 2,048 event cap nor the 2,048 atom cap. The
component had 2,044 input contributors. Search stopped after three completed
sweeps and 5,056 profile evaluations, with 620 accepted updates, at the same
sweep 4 / block 21 failure. The bounded trace does not reconstruct rejected
trials for all 2,044 atoms; the focused trace records those trials for atoms
11721 and 11650.

### Atom 11721

| Event | Width before → after (Å) | Log-width step | Evidence |
|---|---:|---:|---|
| First-stage seed | 0.5107510641 | — | Normal initialization; no extreme input width |
| Sweep 1, block 21, accepted update | 0.5107510641 → 0.4942385633 | −0.0329 | Local profile reduction 9.25765e−8; global block reduction 8.07065e−4 (mostly the entry A/C re-profile) |
| Sweep 2, block 21 | unchanged | — | No width update accepted; trial search explored invalid or rank-deficient proposals and ended with an LM rejection |
| Sweep 3, block 21, accepted update | 0.4942385633 → 85,164.2929583 | +12.0571 | One large accepted LM step; profile objective reduction 9.23220e−6; global block objective reduction 1.59026e−5 |
| Sweep 4, block 21 | remains 85,164.2929583 | — | Local A/C endpoint failed KKT at 4.65661e−10; no new width update was accepted |

The runaway therefore did not come from the first-stage seed or from many
small updates. A single sweep 3 / block 21 LM proposal caused the first
extreme width. Its LM ratio was 1.0499343 (actual decrease 0.268212 versus
predicted decrease 0.255456 in the LM residual convention). At the accepted
endpoint, local A/C KKT was 1.38778e−17 and the width-gradient magnitude was
2.11003e−17. The constrained profile objective decreased by 9.23220e−6; this
was a finite, material reduction in the original objective, not a change
indistinguishable from the recorded binary64 objective. The separate global
block objective decreased by 1.59026e−5, which also includes entry A/C
re-profiling.

At that accepted endpoint, A was zero and C was 4,672.2875; the normalized
charge coefficient was 0.384112, so the accepted step itself did not have
large A/C cancellation. By sweep 4, with the same width, the block's local
conditional A/C solution changed to A=1.0233649e24 and C=−1.1228081e13. Their
column-norm-scaled contributions were +9.23067447e8 and −9.23067446e8. The
large cancellation thus emerged after later neighboring block updates changed
the conditional problem; the trace does not support attributing it to the
initial width proposal alone.

### Other atoms and neighborhood effects

Accepted large-width updates occurred at multiple locations and sweeps, so
11721 is not unique. The first accepted jump for atom 16409 was 0.457613 Å to
251,475.687 Å at sweep 1 / block 58 (Δlog B=13.2168). Atom 17247 moved from
0.085638 Å to 239,833.713 Å at sweep 2 / block 78 (Δlog B=14.8453). Other
examples include atom 11483 (0.287373 Å to 44,291.579 Å, sweep 2 / block 19),
atom 55975 (0.501566 Å to 23,357.994 Å, sweep 1 / block 82), and atom 16610
(0.442268 Å to 18,768.106 Å, sweep 1 / block 63).

Some accepted wide states already had severe normalized A/C cancellation.
At the sweep 1 / block 63 endpoint for atom 16610, the normalized A/C terms
were +4.04948e7 and −4.04948e7; at sweep 2 / block 78 for atom 17247, they
were +1.15177e10 and −1.15177e10. Across the component, accepted large-width
states and coefficient cancellation appeared at different blocks and
interleaved across sweeps. This is consistent with the conditional
landscape changing as neighboring A/C/B values are updated, but a full
counterfactual replay with frozen neighbors was not run here.

At large B, the local Gaussian/charge correlation rounded to 1.0 for these
examples. The local full-design condition indicators ranged from about
5.3e9 to 2.5e12 for the examples above (11721: 6.62347e10 at the accepted
endpoint). These are evidence of weak practical separation, not proof of
mathematical rank deficiency. Detailed fixed-design, functional-kernel, and
profile-slice qualification is deferred to
[`joint-fullabc-width-degeneracy.md`](joint-fullabc-width-degeneracy.md).

### Atom 11650

Atom 11650 also started from an ordinary first-stage width, 0.410172309 Å.
The first accepted block 21 update reduced it to 0.123236457 Å. At that
state, the Gaussian column norm was 2.10279e−38 while the charge column norm
was 2.44415; A was on its nonnegative active boundary at zero. The width
gradient magnitude was about 1.92e−45. The Gaussian response is nonzero in
the floating-point evaluation, but its amplitude and width information are
negligible at this support. This is a small-width practical-identifiability
issue distinct from the large-width A/C cancellation.

## Conclusions and limits

1. Atom 11721's extreme width was produced by one accepted update at sweep 3,
   block 21, not by input data or first-stage initialization.
2. The LM accepted it because the local constrained profile objective
   decreased materially and the local endpoint passed its local trust,
   A/C KKT, and width-gradient checks. The reduction was reliable at the
   recorded binary64 objective scale.
3. A/C cancellation was not present at that accepted endpoint. It became
   severe by the next sweep's block 21 endpoint after intervening block
   updates, so neighborhood state changes are implicated.
4. Similar large-width states occur in many atoms; some already exhibit
   enormous opposing normalized A/C contributions. This supports a broad
   practical-identifiability concern rather than an isolated corrupted atom.
5. This trace by itself does not justify a hard width cap or a changed
   objective. It also does not prove that an alternative local step would
   produce a scientifically better estimate. The next study must quantify
   profile flatness, normalized spectra, and constrained precision references
   before changing production search policy.

The run had the supplied CIF and map available, so the trajectory is a
reproducible full-input diagnostic through the first numerical qualification
failure. It is not a successful full 6Z6U qualification.
