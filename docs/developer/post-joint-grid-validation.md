# Post-Joint grid-domain production decision

## Adopted method

Post-Joint processing uses Original-grid Peeling and direct own-charge
correction. Candidate collection is inclusive at **2.0 Å**; qualified Local
fitting is fixed at **1.0 Å**. The completed 1.0/1.5/2.0 Å study did not justify
promoting either wider fitting range. These are separate collection and
estimation responsibilities.

For Target i, collect original map voxels V_i with ||x_v-r_i|| <= 2 Å from
its Original support, in deterministic Original-row order. Each observation
retains the Original row, global grid index, physical position, squared
distance and raw voxel response. The Fitting-domain mask does not reduce this
candidate set.

```text
y_peel(i,v) = y_v - sum[j != i, r_jv <= 2.5 Å]
                       (A_j G(r_jv;tau_j) + C_j K(r_jv;tau_j))
z(i,v)      = y_peel(i,v) - C_hat_i K(r_iv;tau_hat_i)
```

Selected Targets use qualified Joint endpoints; Halos use qualified
geometry-local kNN-4 Target median reconstructions. Sparse Original support
contributions are evaluated once and non-self contributions are summed
directly. Missing contributor state makes affected responses unavailable.
Own charge uses the existing stable kernel, including its analytic origin
limit. No Post-Joint map or prediction interpolation is used.

The finite positive z observations within 1 Å form one selected dataset for
Local alpha training, OLS/MDPDE, QR rank qualification, working covariance and
measurement information. Finite nonpositive corrections are valid observations
but are excluded from Gaussian solver rows. Empty datasets retain zero actual
rows. FirstStage and Two-Stage keep their existing sampling and fitting rules.

Group fitting retains its current alpha training, MDPDE and WEB kernels and
accepts only qualified Local evidence with matching source/run identity.
Gaussian posterior and fixed Joint charge remain separate. The
[stage integration guide](joint-stage-integration.md) owns provenance,
qualification, unavailable semantics and invalidation. Persistence remains
Stage JSON 5 / SQLite 21 / Joint JSON 24, with strict rejection of unsupported
formats; the Local fitting-radius field records the estimation contract.

## Why fitting remains 1 Å

The completed paired synthetic study used 48 configurations with 16 trials
each: spacing 0.35/0.7 Å, aligned/off-grid centers, base tau 0.4/0.7/0.9 Å,
noise SD 0.5%/3% of Gaussian peak, independent or spatially correlated noise,
and signed known charge. The same innovations were paired across fitting
radii. The study conditioned on true endpoints; it did not propagate Joint or
Halo estimation uncertainty. Representative aggregate results were:

| Metric | Grid 1 Å | Grid 1.5 Å | Grid 2 Å |
| --- | ---: | ---: | ---: |
| Qualified Local / 768 | 767 | 766 | 768 |
| A bias | +0.012441 | -0.030434 | -0.063734 |
| A RMSE | 0.098797 | 0.152310 | 0.501155 |
| tau bias (Å) | +0.000476 | +0.023456 | +0.115998 |
| tau RMSE (Å) | 0.012598 | 0.057009 | 0.214714 |
| Positive retained / radius candidates | 99.77% | 93.42% | 86.49% |
| Mean information condition | 68.35 | 24.60 | 19.59 |
| Nominal 95% posterior ellipse coverage | 44.59% | 27.81% | 13.67% |

At base tau=0.4 Å, the 2 Å mean width bias was +0.3416 Å with independent
noise and +0.3343 Å with correlated noise, versus +0.00083/+0.00407 Å at 1 Å.
Broader peaks sometimes benefited from more radial leverage, but improved
conditioning and numerical qualification did not compensate for narrow-peak
tail selection bias and deteriorating coverage. No qualification threshold
was relaxed. At 1 Å, posterior coverage was 62.24% for independent noise and
26.89% for correlated noise, still far below nominal 95%.

An independent 16-trial overlapping-Target study with signed charges confirmed
consistent subtraction with exact Halo medians. The narrow Target width RMSE
was 0.00143/0.02388/0.23306 Å at 1/1.5/2 Å. Deliberately misspecified Halo
profiles also biased 1 Å estimates, demonstrating a remaining conditional-model
limitation rather than a reason to widen the production range.

The statistical evidence is historical, not an active runtime fitting-policy
selector. Raw per-Target/per-Group research CSV/JSON outputs are retired; their
history remains in Git at `98b101cd6670451e1f4f7f37abbd51ef4eb96d56`.

## Representative adoption performance

The prior paired benchmark used 6Z6U / EMD-11103, 70,190 model atoms,
181 selected non-H Targets (chain A, sequence IDs 30–50), a 480³ map with
0.492 Å spacing, Release Apple Clang/Eigen/ROOT/OpenMP, and one numerical
thread. Three independent processes per implementation compared the former
interpolated path with Original-grid processing, both fitting within 1 Å.
These are adoption measurements, not before/after measurements of subsequent
production cleanup.

| Median phase | Former path (s) | Original grid (s) |
| --- | ---: | ---: |
| Peeling | 0.605913083 | 0.030749875 |
| Own charge / preparation | 0.030867829 | 0.007354877 |
| Local alpha | 0.571445703 | 0.346706747 |
| Local fit | 0.005666874 | 0.002600571 |
| Group WEB | 0.000129289 | 0.000132541 |
| Total postprocessing | 1.222804708 | 0.410749625 |
| Process peak RSS (B) | 1,498,513,408 | 1,498,923,008 |
| Estimated postprocessing payload (B) | 19,451,926 | 8,658,344 |

Joint A/tau/C bit patterns and objective 0.01566223279243328 were unchanged.
All 50,936 grid Peeling observations were available. Both paths had 178/181
qualified Locals and 25/85 available Groups, although member identities and
estimates changed with the observation domain. The 66.4% postprocessing time
reduction and +0.03% process RSS difference apply to this selection. Real-map
truth was unknown; these observations establish neither accuracy nor universal
scaling. The sparse implementation does not allocate a dense voxel-by-atom
matrix.

## Regression authority and remaining limits

Normal CTest retains independent Gaussian/erf references, inclusive 2 Å grid
selection and 2.5 Å support, overlap/Halo subtraction, own charge, missing and
nonfinite contributors, nonpositive selection, coarse-grid rank failure,
endpoint qualification, source/run/grid identity, map/geometry/snapshot
invalidation, Group eligibility and JSON/SQLite round-trips. FirstStage
tricubic sampler coverage and malformed/legacy-format rejection remain active.
See the [test guide](../../tests/README.md) and
[consolidation audit](post-joint-production-audit.md).

Low-SNR positive selection bias, correlated voxel noise, Halo reconstruction
uncertainty, Local width identifiability and Group posterior covariance
calibration remain separate statistical research questions. Rank/SPD
qualification certifies legal numerical evidence, not unbiased truth recovery
or nominal coverage. Joint charge and Halo profiles remain conditional inputs;
this production decision does not redesign those statistical models.
