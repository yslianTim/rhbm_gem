# Operator search and Schwarz preconditioners

Current solve measurements use the
[Joint benchmark guide](joint-benchmark.md). Historical promotion results and
incomplete campaigns are indexed in the
[canonical historical evidence](joint-component-evidence.md).

Production still defaults to SearchMethod::LegacyCompact. The operator path is
available through the internal EvaluationContext search policy and
FitWithSearchPolicy; it does not extend the installed API or serialized schema.
Search completion is independent of the returned state's convergence evidence.
The current solver has no two-level Schwarz or coarse correction. The
`joint_scaling_analysis.py` report only classifies whether the measured
one-level PCG iteration trend warrants further coarse-space investigation; it
does not alter search policy or launch a solver experiment.
Returned-state assessment uses guarded compact stacked QR when its rank
decision is outside the one-ULP boundary enclosure and falls back to the exact
observation-scale Jacobian TSQR when the decision is ambiguous.

## Search and ownership

The operator path solves

    (J'J + mu Gamma'Gamma) step = -J'(r/s)

using the full profile operator, including residual correction. It constructs no
global normal matrix, reduced Jacobian, or dense LM QR. Operator rank selection
is explicit in `SearchPolicy`: `auto` resolves to bounded SPQR rank when the
SPQR backend is active and to dense rank on Eigen; `dense` selects the dense
control on either backend. `spqr-bounds` is unavailable on Eigen. The SPQR
bounded route uses the free-design factor's rank view and does not extract a
global p-by-p compact rank matrix or run a dense rank SVD. An unavailable or
unresolved bounded result stops the operator search with its rank reason and
never switches to dense. Rank checks use the policy's finite time, entry, and
workspace budgets. A/C, reference, assessment, assembly and uncertainty retain
their current paths.

The production OperatorPcg path reuses the accepted profile factor when its
normalized free design and column map match the operator design. For bounded
SPQR rank, reuse is accepted only when the existing no-factor rank check
certifies FullRank; otherwise the dedicated fixed-factor path supplies the
same bounded rank view. `Unavailable` remains unavailable. The operator holds
the accepted factor read-only, and the trial workspace copy-on-writes before
its next factorization, so a later evaluation cannot invalidate the operator.
Rejected trials retain the accepted linearization. Accepting a trial rebuilds
state-dependent numerical data. Only structural block topology survives
accepted updates. There is no state/face cache.

The stored profile factor is for normalized free columns `Z = X_f N^-1`, where
`N` is the diagonal of raw free-column norms. Operator action tests compare
reuse against a dedicated normalized factor and against the raw-design scaling
identities for Apply, ApplyAdjoint, and ApplyNormal. They check the adjoint and
normal identities, bounded rank status/certificate/reason, and workspace
mutation lifetime.

Raw-width column norms d_j = norm(D_j)/s provide an approximate search metric.
Positive floors are sqrt(epsilon)*max(d) (all-zero norms use ones), and
accepted updates take coordinatewise maxima. All preconditioner controls use
the same metric. This approximation changes trajectory, not objective, KKT,
rank, or endpoint stationarity thresholds.

PCG starts from zero and checks the true metric-scaled residual at relative
1e-10, every 32 iterations and before success. Residual replacement restarts
the recurrence. The iteration limit is min(4m,1000); internal test overrides
may only reduce it. Damping, metric, operator and inverse action stay fixed
throughout each solve. Curvature, nonfinite action, stale context and exhaustion
have explicit failure reasons, with no dense or diagonal fallback.

The initial radius is 0.1*norm(Gamma eta), or 0.1 if zero. Initial damping is
1e-3, floored at 1e-12. Bracket search returns a step inside the radius, with at
most 20 damping solves per accepted linearization; internal overrides may only
reduce this limit. Invalid or untrusted trials quarter the radius and quadruple
damping. Acceptance uses the full trial objective and predicted reduction
    -g.dot(step) - 0.5*norm(J step)^2
requires ratio at least 1e-4, and passes scalar replay. Radius adaptation uses
0.25/0.75 ratio boundaries. Reduction and step tolerances remain 1e-14 and
1e-12; gradient stopping uses the full gradient infinity norm 1e-12. These are
search stops, not endpoint evidence. Profile/update budgets remain 200/100.
Search performs zero reference solves. Endpoint certification and newest-first
fallback independently check the saved coefficients.

## Schwarz preconditioner

Structural partitions and coordinate identity are owned by the
[profile-operator guide](joint-profile-operator.md). The search uses informative
frozen-support incidence, parent atom/row indices, and separate mappings for
Width and canonical free-A/C coordinates. Preconditioner factors are not
reference factors or covariance factors.

Schwarz approximates local Gauss-Newton curvature after local A/C elimination.
It uses the global state's beta, active face and parent scale without local
nonlinear fits. For normalized local free design Z and
E = D*Gamma_local^-1/s, it forms sparse products A=Z'Z, B=Z'E, C=E'E, then

    S = C - B'(A + lambda I)^-1 B

where lambda = sqrt(epsilon)*max(1,norm(A)_infinity) affects only this
preconditioner approximation. It does not replace the global
residual-corrected Jacobian.

The symmetrized S is factored with damping and a local shift:
S + (mu + tau) I. Tau starts at
64*epsilon*max(1,norm(S)_infinity,mu) and may grow tenfold for at most six LLT
attempts. Local rank deficiency is not a global rank verdict. Both restriction
and scatter use weights 1/sqrt(coordinate membership), and local solves scale
on both sides by Gamma_local^-1. The resulting additive inverse action is fixed
and SPD during PCG. Identity and raw-width diagonal inverses are controls.

`SchwarzPolicy` records `core_atoms`, `overlap_hops`, `max_block_atoms`,
`storage_bytes`, and `scratch_bytes`. Cores are deterministic disjoint groups
up to the requested core size; each block expands its core by the requested
number of hops on the structural shared-row graph. Zero hops is core-only,
one hop is the complete one-ring overlap, and two hops includes the second
structural ring. The default is 128 core atoms and one overlap hop, with a 512
atom maximum block, 512 MiB storage limit, and 256 MiB construction scratch
limit. Resource-limit failures return unavailable without truncating overlap
or changing the solver. Search retains its last accepted state for endpoint
certification. These settings affect only the preconditioner; the full global
operator and objective remain unchanged.

## Validation and promotion prerequisites

Permanent tests cover the damped dense oracle, inverse symmetry/SPD,
permutation, coverage, parent/nuisance mappings, local rank regularization,
resource failure, stale contexts, budgets, absence of global derivative
reduction, replay-only search, and runtime uncertainty/persistence. Operator
tests also cover cancellation, active-face changes, factor lifetime and tiled
parity.

Current bounded solve measurements use an explicit preconditioner profile:

    python3 tests/integration/joint_benchmark.py \
      --profile solve --case chain-8 --preconditioner schwarz \
      --operator-rank auto --schwarz-core-atoms 32 --schwarz-overlap-hops 1 \
      --build-dir build/joint-eigen --output build/search-schwarz.json

A future production-default change requires both backend regressions, agreement
on returned-state availability and required convergence evidence, complete
resource-bounded comparisons, and measured benefit that includes setup and
rebuild cost. Incomplete comparisons cannot justify promotion. The production
default remains LegacyCompact; the historical evidence index records why
previous operator/Schwarz campaigns did not qualify.

The rigorous bounded `LocalSupport` rank certificate has removed the rank-work
blocker under the existing 120-second, 100,000,000-entry, 256-MiB budget.
Current evidence is reported on six distinct scalability axes: rank
certification, Krylov iteration scaling, search operator throughput, search
factor memory, assessment derivative reduction, and end-to-end returned-state
assessment. The `search` benchmark profile runs `SearchProfile(...)` and exits before
`AssessComponentSearch(...)`; its search-only scope does not establish runtime
convergence or endpoint qualification.

The [cube multi-block gate](figures/joint-search-scaling-r2/search-scaling-analysis.json)
uses actual block counts 2/4/7/9 at sizes 256/512/768/1024. It reports four
valid points, slope 0.0841, endpoint growth 1.167, and `stable` under the
unchanged 0.25 slope / 1.5 growth thresholds. The chain multi-block gate is
also stable. Current synthetic evidence does not warrant a two-level Schwarz
investigation; the cube-2048 RSS failures concern search memory and do not
reopen the Krylov decision.

### Cube search memory and SPQR ordering

The [cube memory campaign](figures/cube-search-memory-r1/cube-memory-analysis.json)
adds process-tree RSS samples and active/last-completed search-stage snapshots.
With production COLAMD, cube-2048 diagonal and Schwarz runs both stopped in
`spqr-fixed-factor` at about 4.02 GiB. Diagonal had completed
`profile-evaluation`; Schwarz had completed `schwarz-partition`, but neither
Schwarz local-build nor PCG had begun. This locates the shared baseline blocker
at the global fixed factor, not in Schwarz local blocks or Krylov convergence.

The benchmark-only COLAMD/DEFAULT/BEST/METIS comparison completed cube-512 and
cube-1024 with diagonal preconditioning. DEFAULT matched COLAMD fill. At
cube-1024, BEST reduced fixed-factor nnz from 481,584 to 255,930 and exported
factor storage from 1,472 MB to 744 MB; sampled RSS fell from 3.860 GiB to
3.377 GiB. METIS produced 257,141 fixed-factor nnz, 750 MB storage, and
3.203 GiB RSS. The candidates agreed on rank certificate, accepted updates,
profile evaluations, search stop reason, PCG iterations, objective, gradient,
residual, and returned state within floating-point tolerance. Despite the
material fill reduction, none completed cube-2048 within the original 4 GiB
envelope. BEST/diagonal under 4 GiB reached PCG/trial evaluation but later hit
the limit in `rank-certificate`; other recorded candidates hit
`spqr-fixed-factor`. Production ordering therefore remains COLAMD.

A single diagnostic-only BEST/diagonal run under 6 GiB completed at 5.672 GiB
in 395.7 seconds with mean 7.8 PCG iterations per solve. It is excluded from
the formal 4 GiB scaling gate and provides no Schwarz 2048 point. The
preselected cube-768 point completed the four-point multi-block gate. No coarse
correction was implemented, and the RSS failures are not reported as a Schwarz
solver failure.

### Factor residency and representation

The [factor residency analysis](figures/joint-factor-residency-r2/residency-analysis.json)
records the dedicated-fixed control with two overlapping global factors: an
older profile/trial workspace and the next operator factor under construction.
Production reuses the accepted profile factor. During a trial factorization,
the accepted profile and new trial factor can still overlap, so the current
route can also reach two factors. The owned-factor byte estimate is distinct
from sampled process-tree RSS.

The [operator representation frontier](figures/joint-operator-factor-frontier-r1/factor-frontier-analysis.json)
compared the benchmark-only native SuiteSparseQR state with the dedicated
exported `Fixed` control. On matched completed runs, it preserved the
stop reason, accepted updates, rank certificate, all per-solve PCG counts, and
returned search state within floating-point differences. At cube-512/COLAMD,
native QR lowered search time from 75.58 s to 64.55 s and sampled RSS from
2.04 GiB to 1.89 GiB. At cube-1024/BEST, search time changed from 153.53 s to
132.78 s and RSS from 2.94 GiB to 2.62 GiB. BEST is benchmark-only. For
cube-1024/COLAMD, the native run finished only about 29 MiB below the 4 GiB cap;
the current run varied between a small RSS-limit failure and an earlier
completion near 4.12 GB. All six cube-2048 retries stopped at the 4 GiB limit
before PCG. Production ordering remains COLAMD; native QR remains
benchmark-only. The dedicated fallback retains exported `Fixed`.

The [accepted-factor ownership campaign](figures/joint-factor-ownership-r1/factor-ownership-analysis.json)
compares dedicated factors, copy-on-write reuse, and a handoff prototype. At
cube-512/COLAMD, COW preserved the full search trace while search time changed
from 75.72 s to 62.93 s and peak RSS from 2.21 GiB to 1.06 GiB. At
cube-1024/COLAMD, reuse matched the completed dedicated-fixed reference trace;
search time changed from 273.61 s to 220.79 s and peak RSS from 3.83 GiB to
2.32 GiB. A fresh dedicated-fixed retry stopped at `spqr-fixed-factor` under
the same 4-GiB limit, while COW completed with about 1.68 GiB of RSS headroom.
The raw references, exact policy, and analyzer output are listed in the
[campaign manifest](figures/joint-factor-ownership-r1/campaign-manifest.json).

At cube-2048, both COW search-only profiles remained `rss_limit` under the
original 4-GiB / 600-s envelope: Diagonal stopped at 592.33 s and Schwarz at
570.46 s, with sampled peaks of 4.03 / 4.01 GiB, each during trial
`spqr-numeric` after 28 completed PCG solves.
Both had last completed stage `linear-symbolic`, peak factor count two, and an
owned-byte estimate of 4.60 GiB. The factor-reuse route clears the former
operator-factor construction stage but does not complete a 2048 search; no
2048 search trajectory parity is claimed.

### Assessment derivative reduction

Under the current compact-Jacobian assessment route, 512-atom production
assessments complete in 454.88 s for chain and 410.41 s for cube. Projected QR
uses 304.15 s and 235.94 s respectively; compact Jacobian QR uses about
1.71 s in each case. The [preselected 640 frontier](figures/joint-post-compact-frontier-r1/frontier-analysis.json)
times out in `derivative-projected-qr` for both topologies under the same
600-second / 4-GiB envelope. This makes observation-scale projected-width QR
the first end-to-end assessment blocker after compact-Jacobian promotion.

The [structured projected-width candidate](figures/joint-projected-width-r1/projected-width-analysis.json)
matches full Assessment and returned state at 128/256, with 1.61x–2.23x total
assessment speedups and much lower projected-reduction time. It exceeds the
4-GiB limit at 512 on both topologies, so it remains benchmark-only. Production
continues to use observation-tiled projected QR.

The [compact identity diagnostics](figures/joint-compact-jacobian-r1/compact-jacobian-analysis.json)
measured a maximum `Z^T P` relative residual of 1.40e-15 and full Gram
reconstruction relative Frobenius error of 1.51e-15. The stacked QR
`[R_P; -R_Z C]` matched rank across chain/cube 8/32/128/256; maximum
singular-value relative difference was 3.05e-15 and correction difference
4.20e-15. Boundary tests cover clearly-below, one-ULP-below, exact-boundary,
one-ULP-above, and clearly-above cases plus cancellation, active-face,
row-permutation, scale, and applicable Eigen/SPQR cases. Rank-ambiguous cases
fall back to exact observation-scale TSQR.

The [compact acceptance campaign](figures/joint-compact-acceptance-r2/acceptance-analysis.json)
passed full Assessment and returned-state parity for the required endpoints
and 128/256 acceptance cases. Median assessment speedups were 1.41x–1.57x,
with peak RSS changes no larger than 3.4%. After production promotion,
chain-512 completed assessment in 454.88 s at 627 MB peak RSS; cube-512 took
410.41 s at 2.56 GB, within the original 600-second / 4-GiB envelope. Production
uses guarded compact Jacobian reduction and retains current observation-scale
TSQR as its exact rank-boundary fallback. Rank, objective, and convergence
thresholds are unchanged. The [benchmark guide](joint-benchmark.md) links the
projected-width method and end-to-end frontier evidence.
