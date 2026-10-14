# Joint FullABC halo eligibility

## Production contributor model

Selected non-hydrogen atoms define the fixed observation domain
`V = union(support(target))`. The builder reads each unselected non-hydrogen
atom's support on those same rows and treats it as a halo candidate. Selected
targets are always retained, including targets with zero or insufficient
support. A halo is included in `JointProblemInput` only when it passes the
geometry-based FullABC eligibility rule:

```text
m_h = number of halo support rows on V
q_h = number of exactly distinct squared distances on those rows
eligible(h) = (m_h >= 3) and (q_h >= 3)
```

Squared-distance equality is exact and deterministic. No tolerance or
numerical-distinguishability claim is used; adjacent representable distances
count as distinct. The production minimum is three observations, the smallest
possible for the three FullABC parameters `(A, C, log(B))`. The threshold is not
a measured optimum and does not assert identifiability.

All retained contributors use

```text
f_i(r) = A_i G(r; B_i) + C_i K(r; B_i)
A_i >= 0, B_i > 0, C_i unrestricted
```

`C` multiplies the charge kernel and is not a constant background offset.
Production has no singleton-halo nuisance amplitude, nuisance-only component,
or alternate profiling route.

## Fixed rows and changed model

Halo eligibility only selects contributors. It does not rebuild the selected
row set, remove rows, or change observations, row IDs, radii, or observation
scale. A row with no retained contributor remains in the objective as a
constant residual. Structural components and bridge relationships are built
after pruning from the retained contributors.

Pruning changes the fitted additive model. The observations still contain any
excluded halo signal, while the fitted prediction contains only selected
targets and eligible retained halos. An excluded signal can therefore increase
residuals or shift selected-target estimates. Production adds no background
term, regularizer, pseudo-observation, or redistribution of that signal.

## Numerical qualification

The gate uses geometry before first-stage width seeds exist. `m_h >= 3` and
`q_h >= 3` are necessary local conditions only. Repeated radii can limit the
local Jacobian rank to `q_h`; passing the gate does not prove that the halo's
fixed-width A/C columns are independent or that its full `(A, C, B)` parameters
are identifiable in the complete component.

`DiagnoseHaloSupport()` remains a diagnostic. Its fixed-width A/C numerical
rank and condition indicator are distinct from the full FullABC Jacobian and
from the component rank checks. Existing endpoint replay, constrained A/C KKT,
width stationarity, full design rank, projected-width rank, and target
estimability requirements remain unchanged. A retained but rank-deficient
component can still fail numerical qualification.

## Exclusion provenance

New Map/Model problems record the contributor policy
`fullabc-eligible-halo-pruning-v1` in `JointSelectionDomain`. The provenance
stores the minimum observation and distinct-distance counts, candidate and
retained halo counts, plus each excluded atom's identity, support-row count,
distinct squared-distance count, and reason (`insufficient-observations` or
`insufficient-radial-information`). Excluded atoms are not fitted failures and
receive no A/C/B estimate or covariance.

Joint JSON schema 11 stores this provenance. Schemas 3 through 10 remain
readable, including saved singleton-halo layouts and nuisance amplitudes; old
results retain their original parameterization when re-encoded. Current CSV
atom rows contain retained contributors. The nuisance contribution CSV is
written only when exporting a historical result that has nuisance groups.

## Regression and scientific sensitivity

Builder regressions cover fixed rows, selected-target retention, exclusion,
first-stage contributor mapping, and component rebuilding after a bridge halo
is removed. SparseQR rank regression includes a design with more rows than
columns that remains rank deficient. The deterministic weak, strong low-support,
overlapping-halo cases and the 3/4/6/8/12 observation sweep are reported in the
[halo pruning validation study](../research/joint-halo-pruning-validation.md).
Those measurements separate numerical convergence from selected-target
parameter error; they do not establish a universally best threshold.
