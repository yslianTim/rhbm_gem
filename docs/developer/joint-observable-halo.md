# Joint fixed-halo contributor model

## Production contract

Selected non-hydrogen atoms define the unchanged observation domain
`V = union(support(target))`. Every unselected non-hydrogen atom with at least
one support row on `V` is retained as a fixed halo. Atoms with no support on
`V` are omitted because they cannot change the target objective. Observation
support count and radial diversity do not prune a halo.

Each retained halo uses a validated FirstStage local MDPDE snapshot:

```text
(A_h, B_h, C_h) = (A_h_FirstStage, B_h_FirstStage, 0)
A_h is finite and >= 0; B_h is finite and > 0
```

The snapshot records contributor identity, FirstStage sample count and status,
fit qualification, and the source contract. Joint-domain support count remains
available separately in the contributor geometry. A failed or invalid
FirstStage fit does not silently become a zero halo or median-width fallback;
the Joint fit reports `fixed-halo-initialization-unavailable`.

Only selected targets have free `(A, B, C)` parameters. Structural components,
FixedNeighbor cores, local conditioning, parameter stability, and numerical
convergence use targets only. A halo that overlaps observations from several
disconnected target components remains one problem-level fixed snapshot and
does not connect those components.

## Objective and complete prediction

For a halo support row `v`, the estimator evaluates the existing Gaussian
kernel only on stored support entries:

```text
H_v = sum_h A_h_FirstStage * G(r_vh; B_h_FirstStage)
y*_v = y_v - H_v
```

Targets are fit to `y*`, while the normalization scale stays
`max(1, ||y_V||_2)` using the original observations. The complete prediction is
`H_v + prediction_targets`; its residual and objective must replay the
target-only residual and objective. Fixed halo gradients do not enter target
KKT or width stationarity checks.

The selected observation rows, row IDs, and scale do not change with the number
of halos. Rows without a free target contribution remain in the objective as
constant rows. The full contributor representation is retained for prediction
replay and PostFitPeeling.

## Saved results and downstream use

New results use contributor policy `first-stage-fixed-halo-v1`, parameterization
contract `target-only-fixed-halo-v1`, and Joint JSON schema 14. The snapshot
stores fixed halo parameters independently from target component state. Target
stage estimates are marked as Joint targets; fixed halo estimates preserve the
FirstStage source and are not labeled as Joint-optimized convergence. Schemas
3 through 13 remain readable with their original free-halo or pruned-halo
meaning; decoding an older result does not reinterpret it as fixed-halo state.

PostFitPeeling reconstructs each target component plus the problem-level fixed
halos, subtracts other target and halo contributions once, and preserves the
target's existing self-Gaussian/self-charge behavior. Post-Joint local MDPDE and
group fitting remain target-only. Fixed halos are contributors to the Joint
prediction, not members of downstream target inference.

## Historical pruning research

`fullabc-eligible-halo-pruning-v1` and its support/radial thresholds describe
older FullABC results only. The former policy used minimum observation and
distinct-distance counts to remove halos before fitting. It is retained for
historical JSON decoding and old research records, but is not a production
selection route. The [halo pruning validation study](../research/joint-halo-pruning-validation.md)
records the measurements and assumptions of that retired policy.
