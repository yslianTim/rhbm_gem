# Width recovery and identifiability closeout

## Scope and evidence

Audit date: 2026-10-30. Initial `develop` HEAD:
`f3f3da85711152aeac85fdf924a0a3c6795e4d2f` (clean tree).
No branch/worktree creation, push, new dependency, or Notion update.
Optimization qualification, Width identifiability, and eligibility for a specific
downstream use are separate questions. Stationarity does not imply a narrow
confidence interval, and an active amplitude constraint does not imply a zero
Width derivative.

R1 regression entry points are
`WidthFailureIsScopedToOneBlockVisit`,
`TransientWidthFailureRecoversAndAssemblesQualifiedComponent`, and
`WidthRecoveryExhaustionAndUnqualifiedTerminalStatesRejectAssembly`.
The scoped callback reaches the existing LM failure hook **after** initial A/C
profiling. Search records `width_failed`; the block rejects the entire visit,
leaves its accepted state intact, and retries on the next sweep. The integration
test calls `FitFixedNeighborComponents`, including production component evidence
construction and assembly, rather than the diagnostic facade. Actual endpoint
verification and independent replay pass. The retained historical LM check is
Failed (`failed=1`), while component/global RuntimeConvergence pass. Negative
tests exhaust real recovery, stop after recovery with nonstationary KKT, corrupt
a candidate, and invalidate/remove required terminal evidence. No test fabricates
Passed evidence for the successful recovery chain. R1-A: 51 focused tests passed;
R1-B: 62 focused tests passed.

## Code audit (R2-A)

`TiledDerivative.cpp` forms each raw Width column from **both** derivative
columns, `A*dG/deta + C*dQ/deta`, then projects out the free A/C design. The
complete profile derivative also includes its residual-dependent correction.
The amplitude column is free only for A>0; the charge column remains free.

`LocalConditioning.cpp` assesses the free A/C, projected Width, and complete
profile spectra with column normalization. RankPolicy uses machine epsilon times
the larger of row/column counts, with explicit rank-boundary sensitivity checks.
Zero columns/rank loss yield Singular; roundoff amplification or uncertain rank
boundaries yield IllConditioned; invalid states, resource limits, and failed
derivatives can yield Unavailable/partial. A=0 is explicitly reported as
`active_amplitude_count`, `active-face-*`, and `active-set-boundary`: even an
informative charge-only direction is conservatively ActiveFace. Classification
is per block, so otherwise informative neighbors can share its warning.

`TargetEvidence.cpp` maps only assessed WellConditioned blocks to Adequate;
ActiveFace, Singular and IllConditioned become Warning. Target runtime requires
convergence and valid parameters, independently of this warning.
**Adequate is a numerical local geometry check, not an absolute sensitivity,
noise calibration, global uniqueness, or statistical precision certificate.**
The regression `NormalizedConditioningDoesNotMeasureAbsoluteWidthInformation`
scales A/C and observations by 1e-12: normalized conditioning stays
WellConditioned while absolute projected sensitivity shrinks by over 1e10.
`ZeroAmplitudeWidthInformationDependsOnCharge` checks zero and informative
charge, actual KKT, full-model derivatives and constrained A/C re-profiling.
All 9 LocalConditioning tests passed.

ParameterStability was deliberately removed before this task (schema 22 / point
validation v2). It is **not assessed**, and was not restored. No perturbation
audit below is a persisted ParameterStability or confidence-interval claim.

`QualifiedJointTargetDonors` requires runtime/target convergence, valid A/B/C,
and Adequate conditioning. It does not examine A alone. The current conditioning
policy nevertheless excludes all A=0 active-face blocks, including informative
ones. This is an existing conservative policy, not proof that their Widths lack
information. There is no separate absolute Width-information gate.

`PostFitPeeling.cpp` uses each qualified Joint target's entire A*G+C*Q model.
Its Local response subtracts the tricubically interpolated C*Q at the **Joint**
Width. Local does not require Adequate target conditioning; it fits a new Gaussian
A/B and checks positive parameters, MDPDE qualification, working covariance and
positive-definite information. `PostJointLocalValidation.cpp` enforces current
source/run/sample/fixed-charge pairing. `JointGroupFitting.cpp` accepts only this
qualified Local information and Local Gaussian point. It never takes Joint B
directly as a Gaussian measurement, but charge correction can affect it indirectly.

## 7A6A evidence (R2-A)

Original CIF/map were available. The current Release CLI ran with chain A,
normalization, and four-thread limits, writing a **new**
`/tmp/rhbmgem-jce-closeout/final.sqlite`, key `7a6a-closeout`.
The older successful `/tmp/rhbmgem-jce-fix/final.sqlite` was read-only.
Input hashes and the fitted parameters match that baseline. Current run ID:
`3edd06cc-6ed1-48c7-a138-58dfd71565fb`.

All 1418 targets are optimization-qualified. A/C KKT = 6.86763424e-8;
projected Width KKT = 9.57221247e-10; original B upper bound =
3.0660309206224183 angstrom. The 111-sweep component, endpoint, replay and global
assembly pass, with zero unresolved or failed LM subproblems.

The following norms use the original component observation scale
S=2098.24637605. Raw is `||A*dG+C*dQ||/S`; projected is the per-column norm after
removing the free local A/C span, also divided by S. All 136 production blocks
and 1418 columns were audited using their actual effective responses and fixed
outside neighbors, rather than a Gaussian-only derivative.

| Atom ID / chain A | A | B | C | Raw norm | Projected norm | Width KKT | Conditioning / Halo pool | Local / Group |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 393 LYS50 NZ (N) | 0 | .01396310 | -.14040946 | 1.20e-71 | 1.15e-71 | 1.11e-75 | ActiveFace / excluded | invalid-local-parameters / excluded |
| 1421 GLU163 OE2 (O) | 0 | .03382466 | -.17016283 | 3.95e-19 | 3.85e-19 | 4.68e-24 | ActiveFace / excluded | Available / included |
| 1415 GLU163 CG (C) | 0 | 3.06603092 | 1.31542754 | 2.87e-3 | 1.74e-4 | 0 (upper-bound projected) | ActiveFace / excluded | invalid-local-parameters / excluded |
| 544 GLU65 CD (C) | .10713355 | .17397750 | -.23234023 | 5.31e-4 | 2.80e-4 | 2.30e-16 | Adequate / eligible | Available / included |
| 435 LYS54 NZ (N) | 3.61328020 | .74776986 | .16548032 | 2.53e-3 | 8.56e-4 | 6.68e-15 | Adequate / eligible | invalid-local-parameters / excluded |
| 590 LYS69 NZ (N) | 1.22395538 | .57385054 | .07863686 | 1.25e-3 | 5.70e-4 | 8.64e-16 | Adequate / eligible | invalid-local-parameters / excluded |
| 787 ASP90 O (O) | 29.03248427 | .55836896 | -1.15064642 | 2.62e-2 | 1.40e-2 | 5.08e-13 | Adequate / eligible | local-mdpde-unqualified / excluded |
| 1 MET1 N (N) | 5.16328999 | .73858733 | -.16032554 | 2.99e-3 | 1.29e-3 | 7.50e-16 | Adequate / eligible | Available / included |

These are **all three exact A=0 atoms**. ID 544 is the smallest positive A and
the smallest raw Width norm among A>0 targets; it remains measurably informative,
so neither small A nor small B alone justifies exclusion. The next-smallest raw
norm is ID 590 at 1.25e-3. NZ393 and OE2(1421) are many orders below that gap.
Their complete profile norms are 1.16e-71 and 3.85e-19 respectively; CG1415's
is 1.74e-4. Its raw gradient is -4.55916209e-8, correctly projected to zero at
the original upper bound. It has informative charge and is not a flat direction.

The 18 Warning atoms are block 51 IDs
1391,1392,1393,1394,1395,1396,1397,1398,1399,1405,1415,1423;
block 55 ID 1421; block 96 IDs 385,387,389,391,393.
All three blocks are ActiveFace, not an ill-conditioned normalized spectrum.
There are 1400 **eligible** Halo donors, of which 713 distinct IDs actually appear
in the four-neighbor lists for 734 halos. Eligible and selected donor membership
must not be confused. Neither nearly flat target is in either donor set.

The four invalid-local-parameters atoms are IDs 393,435,590,1415; the sole
local-mdpde-unqualified atom is 787. Thus two of four invalid parameter results
have A=0, but the other two have informative Width and Adequate conditioning.
OE2's very weak Width still produces a qualified Local fit. **These facts do not
establish causality for the Local failures.**

The five singular-group-covariance groups are the CYS N, CA, C, CB and O groups,
each with two eligible members (CYS103/CYS131): IDs (895,1131), (896,1132),
(897,1133), (899,1135), (898,1134). None overlaps A=0 or conditioning warnings.
The remaining 21 unavailable groups are single-member. NZ393 and CG1415 are
excluded from their otherwise available groups; OE2 is included in an available
group. Counts alone do not identify a statistical cause.

### Controlled profile sensitivity

For IDs 1,393,435,544,590,787,1415,1421, multiply B by .5,.8,1.25,2 within the
original upper bound; keep other Widths and outside-block contributions fixed,
and re-solve constrained A/C for the block. Compare with the same block's
re-profiled baseline. No priors or lower bounds were added. Objective differences
are normalized by S^2; prediction norms by S. Charge differences below use the
actual 100 stored Local sample geometries and production tricubic interpolation.

| Case | Profile objective difference | Prediction difference | Maximum fixed-C Local correction difference |
| --- | --- | --- | --- |
| NZ393, B*.8, *1.25, *2 | 0 at double resolution | 0 | 0 |
| NZ393, B*.5 | full-design profile rejects rank loss; A=0 face only gives 0 | 6.76e-19 on fixed-neighbor face | 0 |
| OE2(1421), B*.5 or *.8 | 0 | 0 | 0 |
| OE2(1421), B*1.25 | 3.87e-20 | 2.90e-15 | 5.59e-12 |
| OE2(1421), B*2 | 9.06e-14 | 5.89e-9 | 1.14e-5 |
| CG1415, B*.8 / *.5 | 1.49e-8 / 1.18e-7 | 4.78e-5 / 2.30e-4 | .0856 / .342 |
| low-A CD544, B*.8 / *1.25 | 1.39e-9 / 2.88e-9 | 5.20e-5 / 6.86e-5 | .133 / .133 |
| informative N1, B*.8 / *1.25 | 5.29e-8 / 4.04e-8 | 2.96e-4 / 2.62e-4 | .0423 / .0343 |

NZ393's half-width Gaussian column underflows; that failed full A/C profile is
reported rather than counted as a successful solve. The additional A=0-face
experiment fixes its neighbors and solves only its free charge coefficient under
the same constraints. Smaller negative perturbations and both positive
directions have valid full profiles. Re-profiled C changes do not materially
alter the flat-target conclusion (NZ correction <=2.22e-16; OE2*2 <=1.08e-5).

This establishes local, sample-dependent flatness for NZ393/OE2, not global
nonidentifiability over the entire Width domain. Larger Width changes can become
observable. A normalized objective difference alone is not a likelihood ratio
or confidence interval; no calibrated map noise model is available here.
Informative CG1415 and CD544 can change Local responses if their Widths change,
but this is not evidence that their fitted Widths are unreliable.

The audit executable/source, full per-atom TSV, logs, SQLite and timing data are
temporary under `/tmp/rhbmgem-jce-closeout/`; no binaries or large data are tracked.
R2-A does not change production numerical or downstream eligibility policies.
