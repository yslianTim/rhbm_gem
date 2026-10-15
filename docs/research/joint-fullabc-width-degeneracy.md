# FullABC width degeneracy and Block 21 qualification

## Scope

This is a conditional local study of the failure-only fixture
tests/fixtures/joint_local_ac_kkt/block-21-sweep-4.json, captured from
component 1, sweep 4, block 21 of the 6Z6U FixedNeighbor run. The fixture
contains 930 local observation rows and 12 local contributors (24 A/C columns),
the exact binary64 support distances and effective response, current log
widths, numerical scale, and production A/C candidate. It does not contain
the full component state.

The study leaves the objective

\[
L(A,C,B)=\frac{1}{2s^2}\|X(B)(A,C)-y\|_2^2,\quad A\ge 0,\ C\in\mathbb R
\]

unchanged. The constrained profile re-solves A/C whenever a width changes.
All other local widths stay at the captured Block 21 values. The width sweep
is therefore a conditional profile slice, not a global uniqueness claim.

The bounded offline program is joint_fullabc_width_study. It writes
summary.json and profile-slices.csv to a fresh output directory:

    cmake --build build/joint-validation --target joint_fullabc_width_study --parallel 4
    build/joint-validation/bin/joint_fullabc_width_study \
      tests/fixtures/joint_local_ac_kkt/block-21-sweep-4.json \
      /private/tmp/joint-fullabc-width-study

The two slices contain 39 widths for atom 11721 (0.1 Å to approximately
1e6 Å, log step 0.5, plus the captured, pre-run and centered-difference widths)
and 24 widths for atom 11650 (0.01 Å to 0.8 Å, log step 0.25, plus its
first-stage and centered-difference widths).
Each point records the production binary64 constrained profile and a bounded
50-digit, column-normalized Householder-QR constrained reference for both
the rounded double design and the functional kernels. The reference uses an
active-set loop and verifies free-gradient and active-boundary KKT conditions;
it does not use normal equations. The maximum local reference shape is
930 × 24.

## Kernel dependence and rank evidence

For \(B\gg r\), the kernels have the expansions

\[
G(r;B)=\frac{1}{(2\pi)^{3/2}B^3}
  \left(1-\frac{r^2}{2B^2}+O(B^{-4})\right),
\qquad
K(r;B)=\sqrt{\frac{2}{\pi}}\frac{1}{B}
  \left(1-\frac{r^2}{6B^2}+O(B^{-4})\right).
\]

On a fixed finite support, each normalized column approaches the same
constant direction. Their angle shrinks as width increases. An A/C pair can
then retain a modest prediction through very large opposing coefficients.
This is a weakly identified parameter direction even while its finite-width
design remains mathematically full rank.

At the captured width of atom 11721, \(B=85,164.2929583\) Å:

| Measure | Value |
|---|---:|
| Gaussian column norm | \(9.01992495\times10^{-16}\) |
| Charge column norm | \(8.22106115\times10^{-5}\) |
| Normalized G/K correlation | 1.0 in binary64 |
| Smallest singular value of the normalized two-column pair | \(3.52256\times10^{-11}\) |
| Normalized pair condition indicator | \(4.01473\times10^{10}\) |
| Local full normalized-design rank | 24 / 24 |
| Local full normalized-design singular-value range | \(2.10476345\) to \(3.17773444\times10^{-11}\) |
| Local full normalized-design condition indicator | \(6.62347182\times10^{10}\) |
| Free-face rank and condition | 23 / 23; \(6.61051736\times10^{10}\) |

The normalized pair singular value is computed from the difference of the two
unit columns, rather than \(\sqrt{1-\rho}\), which loses all digits when the
dot product rounds to 1. The normalized correlation is consequently a
rounded diagnostic; the direct pair singular value retains useful evidence.

At this width the production coefficients are
\(A=1.02336489\times10^{24}\) and
\(C=-1.12280815\times10^{13}\). Their column-norm-scaled terms are
\(+9.23067447\times10^8\) and \(-9.23067446\times10^8\). The predicted
contribution is small relative to either term. The local design is full rank,
but the large coefficients and small pair singular value make the estimate
highly sensitive to perturbations in basis evaluation and residual arithmetic.

For comparison, the conditional slice at \(B=0.4942385633\) Å has G/K
correlation 0.401525, pair \(\sigma_{\min}=0.773612\), pair condition 1.5303,
and local normalized-design condition 11.8254. Its normalized A/C terms are
0.0266853 and 1.0576861. This is a well-conditioned local point on the same
fixed-neighbor slice.

## Small-width response

For \(B\ll r_{\min}\),

\[
G(r;B)\propto B^{-3}\exp\left(-\frac{r^2}{2B^2}\right).
\]

The exponential dominates \(B^{-3}\), so the Gaussian response can be
mathematically nonzero yet carry no usable amplitude or width information.
For atom 11650 in this fixture, \(r_{\min}=1.65595229\) Å and
\(B=0.1232364573\) Å:

| Measure | Value |
|---|---:|
| Gaussian column norm | \(2.10279168\times10^{-38}\) |
| Charge column norm | 2.44414642 |
| Gaussian norm divided by observation scale \(s\) | \(6.06369\times10^{-41}\) |
| Pair correlation / pair \(\sigma_{\min}\) | 0.247073 / 0.867714 |
| A/C profile state | \(A=0\) on its active nonnegative boundary; \(C=-0.938758\) |
| 50-digit constrained profile width gradient | \(1.15465\times10^{-45}\) |

Here the issue is not G/K collinearity: the pair directions are well separated.
The Gaussian column is simply too small to measure, and the constrained
profile places A at zero. Changing B then has essentially no effect on the
Gaussian prediction.

The sweep separates response underflow from norm underflow. At \(B=0.057546\)
Å, 14 double Gaussian entries are nonzero and the stable column norm is
\(5.13376\times10^{-178}\), but the production-style sum-of-squares norm
rounds to zero. At \(B=0.044817\) Å, one entry remains nonzero
(\(2.44414\times10^{-294}\)) while that norm still rounds to zero. At and below
about 0.034903 Å, all Gaussian entries underflow to zero in binary64 on this
support. The 50-digit functional-kernel profile still evaluates nonzero
Gaussian responses; the corresponding profile gradients fall below about
\(10^{-494}\) and then \(10^{-5960}\). These are numerical loss-of-information
states, not evidence that the mathematical kernel has a zero response.

## Local A/C solver accuracy and independent references

The current fixture is feasible and the 24-column double design has rank 24.
The active A face fixes local atom index 5 (atom 11650); the remaining 23
columns are free. Multiple independent binary64 methods all miss the existing
\(10^{-10}\) KKT threshold:

| Method | Scaled KKT | Outcome |
|---|---:|---|
| Production SparseQR | \(4.65661\times10^{-10}\) | Reject |
| Dense rank-revealing QR | \(3.02325\times10^{-10}\) | Reject |
| Dense SVD | \(2.32831\times10^{-9}\) | Reject |
| Independent active-face enumeration | \(6.25375\times10^{-10}\) | Reject |
| Dense free-face QR | \(1.39698\times10^{-9}\) | Reject |
| Host long-double free-face QR | \(1.86265\times10^{-9}\) | Reject |

On this host long double has 53 binary digits, the same as double; it is
not a higher-precision oracle.

The separate 50-digit constrained Householder QR reference solves the rounded
double design with the same active A constraint. Its KKT is
\(2.99887\times10^{-66}\). The production coefficients, evaluated in 50-digit
arithmetic on those exact rounded double columns, have KKT
\(4.84606\times10^{-10}\). Their raw half-RSS objective is
2.01592388428234613 versus 2.01592388428231799 at the reference optimum.
After division by \(s^2\), that objective difference is only about
\(2.34\times10^{-19}\), despite a maximum column-norm-scaled coefficient
difference of about \(3.63\times10^3\). The objective is stable while the
coefficient vector is not.

This is a full-rank but ill-conditioned solve with a real binary64
stationarity error. It is not mathematical rank deficiency. It is also not
fixed by agreement among binary64 solvers: QR, SVD, enumeration, and the
production solve all miss the threshold.

### Rounded double design versus functional kernels

The second reference recomputes the Gaussian and charge functions in 50-digit
arithmetic from the fixture's rounded squared distances and widths. It retains
the \(A\ge0\) constraints and reaches active atom 5 with KKT
\(2.67478\times10^{-66}\). For this functional design:

| Measure | Value |
|---|---:|
| Production beta high-precision KKT | \(3.89990\times10^{-10}\) |
| Rounded-design HP solution checked on functional kernels | \(8.74593\times10^{-10}\) |
| Maximum kernel absolute difference from rounded double design | \(3.48913\times10^{-16}\) |
| Maximum relative kernel difference | \(1.74027\times10^{-14}\) |
| Maximum production prediction difference | \(6.65085\times10^{-8}\) |
| Functional versus rounded-design HP objective difference, divided by \(s^2\) | \(7.26\times10^{-14}\) |

The kernel entries differ by only a few ulps, but large opposing A/C terms
amplify those differences into observable prediction and KKT changes. The
fixed-double-design reference answers the linear-algebra question for exactly
the columns production formed. The functional-kernel reference answers the
separate question about rounding in kernel evaluation.

There is also a near-threshold certification discrepancy. Rounding the
50-digit fixed-design optimum to binary64 gives a double certificate KKT of
\(5.47387\times10^{-11}\), which passes the current gate, while evaluating
those exact rounded coefficients and columns in 50-digit arithmetic gives
\(1.11200\times10^{-10}\), just above it. Binary64 KKT alone can therefore
understate the true residual stationarity at this scale. This is why a
binary64/replay agreement must not be treated as independent accuracy proof.

## Width Jacobian and practical information

For the captured production coefficients, the local free-face design \(X_F\)
has rank 23. For each local atom, the width column is

\[
D_i=A_i\,\partial_{\log B_i}G_i+C_i\,\partial_{\log B_i}K_i,
\]

and the study projects \(D\) through \(I-P_{X_F}\). It reports the raw column
scales as well as column-normalized singular values:

| Atom | \(\|D_i\|_2\) | \(\|(I-P_{X_F})D_i\|_2\) | Projected fraction | State width-gradient magnitude |
|---|---:|---:|---:|---:|
| 11721 | \(1.84613\times10^9\) | \(6.74154\times10^{-6}\) | \(3.65170\times10^{-15}\) | \(2.57985\times10^{-3}\) |
| 11650 | \(3.76737\times10^{-39}\) | \(3.61111\times10^{-39}\) | 0.958524 | \(1.15465\times10^{-45}\) |

Atom 11721's raw width derivative is almost entirely explainable by changes in
the free A/C columns. The projected derivative is only \(1.94\times10^{-8}\)
after dividing by \(s\). Its profile-width direction is practically flat
after A/C are re-estimated. Atom 11650's projected fraction is large, but both
the derivative and its projection are vanishingly small in the observation
scale; a good angle does not mean useful information.

Under the project's relative rank scale \(\epsilon\max(m,n)\), the raw
projected-width matrix has numerical rank 11 of 12 and condition indicator
\(4.38\times10^{39}\). Normalizing each projected column raises that rank to
12 of 12 and condition 2.80, which hides the atom 11650 scale loss. The
column-normalized combined Jacobian \([X_F,D]\) has rank 34 of 35 and
condition \(2.60\times10^{16}\); including the active A direction gives rank
34 of 36 and condition \(3.19\times10^{16}\). These are local Jacobian
diagnostics, not a replacement for full-component structural-rank,
projected-width-rank, or target-estimability qualification.

## Conditional profile scans

The width scans keep the other 11 local widths fixed and re-solve all local
A/C coefficients at every grid point. The 50-digit references satisfy
constrained KKT at every recorded point. The binary64 result remains separately
marked: for atom 11721, 27 of 39 points pass the current inner KKT check; for
atom 11650, only 2 of 24 pass because the same local A/C problem still contains
the ill-conditioned atom 11721 direction.

Selected atom 11721 points:

The table reports the binary64 profile candidate objective. Rows that fail
KKT are retained as diagnostic candidates and are not treated as qualified
profile minima; the 50-digit constrained losses are stored separately in the
CSV and summary output.

| B (Å) | Binary64 candidate \(L\) | Pair \(\sigma_{\min}\) | Local condition | Normalized A/C terms | Binary64 KKT |
|---:|---:|---:|---:|---:|---:|
| 0.49424 | \(1.67618674922\times10^{-5}\) | 0.773612 | 11.83 | \(+0.02669,\ +1.05769\) | \(6.94\times10^{-18}\) |
| 1.21825 | \(1.67604174724\times10^{-5}\) | 0.208724 | 11.84 | \(0,\ +1.05419\) | \(3.47\times10^{-18}\) |
| 16,275.48 | \(1.67631701549\times10^{-5}\) | \(9.65\times10^{-10}\) | \(2.42\times10^9\) | \(+3.37118\times10^7,\ -3.37118\times10^7\) | \(2.91\times10^{-11}\) |
| 85,164.29 | \(1.67631700845\times10^{-5}\) | \(3.52\times10^{-11}\) | \(6.62\times10^{10}\) | \(+9.23067\times10^8,\ -9.23067\times10^8\) | \(4.66\times10^{-10}\) |
| 120,260.43 | \(1.67631701314\times10^{-5}\) | \(1.77\times10^{-11}\) | \(1.32\times10^{11}\) | \(+1.84061\times10^9,\ -1.84061\times10^9\) | \(9.31\times10^{-10}\) |
| 888,611.05 | \(1.67631769517\times10^{-5}\) | \(3.24\times10^{-13}\) | \(7.21\times10^{12}\) | \(+1.00417\times10^{11},\ -1.00417\times10^{11}\) | \(4.34\times10^{-8}\) |

The 50-digit functional profile gradient at 85,164 Å is
\(4.39\times10^{-18}\), while the fixed-double-design derivative is
\(-3.19\times10^{-13}\). Both are below the existing \(10^{-12}\) width
stationarity requirement. The production candidate's raw width-gradient
magnitude is \(2.58\times10^{-3}\), but that candidate fails A/C KKT; the
constrained high-precision profile places the same width direction on a very
flat plateau. Thus the production gradient is not reliable until its A/C
subproblem is solved accurately.

The centered finite-difference slice at that width, with a 0.25 log-width
half-step, gives functional-kernel slope (4.57\times10^{-18}) and curvature
(-8.96\times10^{-18}). On the independently re-evaluated rounded-double
designs, the same finite difference gives slope (-1.46\times10^{-13}) and
curvature (1.56\times10^{-12}). The larger fixed-design finite difference
reflects kernel-evaluation rounding under cancellation, not a different
functional-kernel profile. At atom 11650, the functional analytic gradient
at the captured width is (1.15\times10^{-45}); over a 0.125 log-width
half-step, the functional centered slope and curvature are only
\(1.35\times10^{-38}\) and \(2.17\times10^{-37}\). These are finite-step
flatness indicators, not exact Hessian eigenvalues.

The atom 11721 profile loss spans about \(4.16\times10^{-9}\) across the
bounded grid. In the large-width tail, the objective changes by only a few
\(10^{-12}\) while normalized A/C terms grow from \(10^7\) through \(10^{11}\).
At the far endpoint, the rounded double design drifts from the high-precision
functional profile as the cancellation grows. This is both practical
non-identifiability and precision sensitivity, not a unique mathematical
minimum lost solely to the nonlinear optimizer.

For atom 11650, the high-precision conditional profile holds A at zero across
the scan. Between 0.01 and 0.2 Å the loss is essentially flat at the printed
precision, and the width gradient at the captured 0.123 Å is
\(1.15\times10^{-45}\). At the smallest widths, the double design can have an
exactly zero Gaussian column from exponential underflow, or nonzero entries
whose squared norm underflows. Those are distinct from the captured-width
state, where the 28 Gaussian entries remain nonzero but their amplitude
information is negligible.

## Classification and answers

| Observation | Classification | Evidence |
|---|---|---|
| Atom 11721 near 0.494 Å | Well-conditioned identifiable locally | Pair condition 1.53; local condition 11.83; KKT passes |
| Atom 11721 at tens of thousands of Å | Weakly identifiable | G/K directions approach one; profile objective plateaus; normalized A/C terms grow together |
| Block 21 fixed double design | Numerically degenerate, but full rank | Local rank 24; condition \(6.62\times10^{10}\); all tested double solvers miss KKT |
| Block 21 production coefficients | Inaccurate solver result | 50-digit fixed-design KKT \(4.85\times10^{-10}\); constrained QR reference KKT \(3.0\times10^{-66}\) |
| Atom 11650 below about 0.058 Å | Binary64 information loss | Stable norm is nonzero while production sum-of-squares norm rounds to zero; below about 0.035 Å all sampled basis values underflow |
| Full-component qualification | Unavailable | The full-data workflow stops before endpoint, rank, and target-estimability qualification |

Block 21's local A/C KKT failure is chiefly a numerical solve/certificate
problem in an ill-conditioned, near-cancelling G/K direction. The same
direction also makes the parameters practically weakly identifiable: very
different coefficients yield an objective change near \(10^{-19}\) in the
fixed-design constrained solution. The local design is not rank deficient at
the captured state. This separates solver accuracy, numerical rank, and
parameter accuracy.

No arbitrary width cap, fixed-width substitution, prior, penalty, or parameter
constraint change is supported. The sweep 3 accepted step in the separate
trajectory report materially decreased the original objective and passed its
local KKT and width-gradient checks at that time. The later Block 21 state
retains a flat, constrained high-precision profile, while the production A/C
candidate is not accurate enough to qualify. A conditioning-only width veto
would conflate weak parameter uniqueness with an unreliable accepted
objective decrease. The appropriate current behavior is to preserve the
original objective and constraints, keep the failed A/C state rejected, and
retain the width and numerical provenance diagnostics.

Local evidence is sufficient to explain why Block 21 can fail without
asserting that every extreme width is invalid. It does not establish scientific
accuracy of any component A/C/B estimate or provide the full-component rank,
projected-width-rank, target-estimability, endpoint, or global-assembly
evidence required for downstream inference.

## Permanent regression coverage

The Release CTest suite includes the bounded offline
`joint_fullabc_width_precision_regression`. It checks a well-conditioned
conditional profile against constrained 50-digit QR, the captured large-width
and small-width controls, full-rank Block 21 KKT rejection, severe A/C
cancellation, and the near-threshold binary64/50-digit discrepancy. The
`joint_local_ac_kkt_sweep4_block21_rejection` fixture keeps the production
block rejected. The `joint_local_replay_fixture` regression requires the
captured Block 82 trial's `CheckReplay()` to pass while independently reporting
the scalar-kernel A/C KKT failure, and retains wrong-gradient, wrong-KKT and
wrong-coefficient rejection checks. The `AccurateProfileTrialCanRemainWidthNonstationary`
test uses a separate well-conditioned synthetic problem where both binary64
replay paths pass with a nonstationary width derivative. Resource-limited and
threshold-sensitive rank tests continue to require `Unavailable` rather than
converting missing evidence to rank deficiency or success.
