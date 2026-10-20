# Observation-driven Halo FullABC — P0 geometry census

**Date:** 2026-10-20  
**Branch:** `develop`  
**Runner:** `tests/experiments/joint_halo_geometry_census.cpp`  
**Raw result:** [joint-halo-geometry-census-6z6u.json](joint-halo-geometry-census-6z6u.json)

## Inputs and domain contract

The census used the local 6Z6U model at `/Users/yslian/data/6z6u.cif`, the
map used by the prior 6Z6U production run at
`/Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map`,
and selected non-hydrogen targets in chain A. It constructs the domain through
`BuildPureHaloSupportExclusionProblem` and reads the saved
`JointProblem::OriginalInput()`; it does not substitute the Halo-excluded fit
rows. The original domain contains 1,418 Targets and 229,783 observation rows.
Its row IDs are unique global map voxel IDs. Existing Pure Exclusion tests
cover the fitting-row to original-row mapping and preservation of the original
domain.

The Halo candidate set for observation counts is every non-hydrogen model atom
not selected as a Target, including atoms with no support in the original
domain. That yields 35,988 candidates. The runner deduplicates support row IDs
before counting observations. Of these candidates, 659 intersect the original
domain and affect 49,305 distinct observation rows; the remaining 35,329 have
`n_h=0` in this domain and remain potential reconstruction-only contributors.

## Observation counts

| Distinct observations `n_h` | Halo count |
| ---: | ---: |
| 0 | 35,329 |
| 1 | 20 |
| 2 | 13 |
| 3 | 9 |
| 4 | 12 |
| 5 or more | 605 |

Under the requested preclassification, 35,362/35,988 (98.2605%) are fixed by
`n_h < 3`, and 626/35,988 (1.7395%) enter the free set. Restricting the
denominator to the 659 supported Halos, 33 are fixed and 626 are free. The JSON
artifact contains the complete count-frequency distribution.

## Atomic centers and support overlap

Distances below are Euclidean center distances. “All Halo candidates” means
the 35,988 non-target heavy atoms; “supported Halos” means the 659 atoms with
non-empty intersection against `V_original`.

| Center-distance summary | All candidates | Supported Halos |
| --- | ---: | ---: |
| Nearest-neighbor minimum | 1.1782 Å | 1.1790 Å |
| 5th percentile | 1.2114 Å | 1.2166 Å |
| Median | 1.3640 Å | 1.4923 Å |
| 95th percentile | 2.7044 Å | 3.1388 Å |
| Maximum nearest-neighbor distance | 6.3383 Å | 6.3383 Å |
| Halo–Halo pairs below 1.5 Å | 21,740 | 250 |

No candidate pair is closer than 1.0 Å, and no pair has coincident centers.
The closest pairs in the source coordinates are same-chain peptide carbonyl
C–O pairs (for example, atom IDs 18297/18298 in chain G at 1.1790 Å). Thus the
real model does contain legitimate covalently bonded Halo pairs below 1.5 Å;
this is not evidence of coincident or duplicated atom centers. The raw JSON
keeps the first ten close-pair identities for inspection.

Kernel-support overlap is measured independently by intersecting the
distinct `OriginalInput` row sets of supported Halos. Among 216,811 supported
Halo pairs, 3,354 have center distance below 5 Å and 1,783 share at least one
voxel. Shared-voxel counts per overlapping pair are: minimum 1, median 8,
P95 49.9, maximum 239. There are no shared rows for center distances at or
above 5 Å in this discrete grid. Some pairs below 5 Å have no shared voxel,
which reflects the actual lattice sampling rather than a sphere-volume
approximation.

| Halo center distance | Candidate pairs | Pairs sharing voxels | Sum of pairwise shared voxels |
| --- | ---: | ---: | ---: |
| 0–1.5 Å | 250 | 238 | 8,367 |
| 1.5–2 Å | 123 | 119 | 2,111 |
| 2–3 Å | 702 | 604 | 11,094 |
| 3–4 Å | 868 | 531 | 4,571 |
| 4–5 Å | 1,411 | 291 | 1,020 |
| 5 Å or more | 213,457 | 0 | 0 |

The sums count a voxel once per overlapping pair, so a voxel shared by three
Halos contributes to three pairs. Atomic-center proximity and discrete kernel
support overlap are therefore separate measurements.

For context, the median nearest Halo-to-Target distance across all candidate
Halos is 54.108 Å (P05 9.435 Å, P95 87.782 Å); this set includes atoms far
outside the selected chain's observation domain.

## Existing synthetic geometry audit

The historical `joint_label_free_halo_policy` experiment is not a physical
geometry benchmark. Its `overlapping-halos` case places both Halo centers at
exactly `(2, 0, 0)`. Its dense donor case uses coordinate increments of 0.10 Å,
and other donor clusters also contain 0.1–0.3 Å separations. These cases remain
useful as adversarial numerical diagnostics, but cannot establish normal
Halo-Halo identifiability or support a Production promotion decision. P3 adds
separately identified physical scenarios.

## Reproduction and validation

```sh
cmake --build build/debug --target joint_halo_geometry_census --parallel 4
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  ./build/debug/bin/joint_halo_geometry_census \
  /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map \
  /Users/yslian/data/6z6u.cif A
```

The 6Z6U runner completed in 30.5 seconds on this host. `jq` checks verified
that the observation-count histogram sums to all Halo candidates, the
supported count equals nonzero observations, and fixed plus free counts equal
the candidate count. Peak RSS was not captured for this census. The runner
suppresses file-parser chatter so its stdout is valid JSON. Production defaults
and Target Joint fitting are unchanged.
