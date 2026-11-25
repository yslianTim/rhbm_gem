Buffer Expansion Validation research tools

The joint_buffer_expansion target uses the existing Joint solver, Halo
reconstruction, local fitting, grouping and UMAP functions. Configure a
Release build with BUILD_TESTING=ON and build the target. An OpenMP-off
build makes serial solver and initialization execution explicit.

Data layout
  A research directory contains snapshot/SHA256.json and the immutable input
  copies named in that manifest. Use the dataset keys in the manifest.
  Run joint_buffer_expansion geometry MODEL MAP CHAINS OUTPUT_DIRECTORY
  for each dataset. Place exports under geometry/<dataset>. The binary
  incidence records contain a native-endian uint32 row and float64 squared
  radius without padding. rows.tsv and geometry.json supply record order,
  identity and grid metadata.

Typical sequence
  export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=4
  python3 tests/research/buffer/geometry.py STUDY
  python3 tests/research/buffer/geometry.py STUDY --scales
  python3 tests/research/buffer/bundle_selection.py STUDY
  python3 tests/research/buffer/run.py STUDY RUNNER pilot
  python3 tests/research/buffer/check_pilot.py STUDY
  python3 tests/research/buffer/run.py STUDY RUNNER matrix
  python3 tests/research/buffer/compare.py STUDY TAG archive g2
  python3 tests/research/buffer/audit_umap.py STUDY RUNNER
  RUNNER truth STUDY/truth-fixtures.json STUDY/truth-results.jsonl
  python3 tests/research/buffer/truth.py STUDY summarize
  python3 tests/research/buffer/truth.py STUDY workflow
  python3 -m unittest discover -s tests/research/buffer -v

run.py serializes full fits. For a fixed-domain run the output calibration
argument derives both width calibration and the reference norm from that
run's own input and selected outputs. compare.py scores supported MRC layouts
in raw units through the reported divisor. Archived baseline projections are
recomputed on recorded geometry from certified state and Halo medians.

The C++ truth mode applies the original controlled perturbations to a
four-contributor fixture and uses the existing native constrained solver
with the documented width interval. Its oracle prediction metric supplies
known external Halo truth after fitting. The 6Z6U simulation experiment
checks the complete workflow under the standard donor qualification rule.
Use the frozen truth-fixtures.json from the E3–E4 research archive with seed
2907 and SHA-256 6b0d8a589e09bf10446ce285a6747badcc020b73c8275feb1862e76b37a9d9c5.
The manifest and historical reports record the source of the numerical draws.
To inspect uncertified endpoints, copy the desired fixture conditions into
a new file and add "diagnose": true to those conditions. The truth runner
then saves the terminal candidate, its component and block identity, the
native replay checks, and independent 50/100-digit reference calculations.
These diagnostics do not publish trusted parameters.

confirm_truth.py fixtures FROZEN_JSON NEW_JSON prepares independent train
and scoring perturbations for the narrow, coarse and wide cases.
confirm_truth.py workflow TRUE_MAP ORIGINAL_ROWS_TSV SEED NEW_MAP prepares
a noisy complete volume and records its seed, scale and source digests.
The full-donor-audit research stage compares fitted-target donors with
output-only donors on the same certified fitted state. It saves paired
Peeling observations and postprocessing in a separate output directory.

For E5-W, run confirm_truth.py widths E5_SYNTHETIC_FIXTURES STUDY
with --geometry-only. Run the native truth mode on geometry-fixtures.json to
export the incidence, then repeat the widths command without the flag.
Run truth on W2–W5 and confirm fixture files in the generated manifest.
The reports use true Gaussian width bins through truth.stratify_matrix.
An empty selected bundle reuses the certified or failed G0 result.
Forced coarse-grid fixtures are controls of the added contributor, separate
from the selector. Each arm uses the same pinned B1 interval and noise draw.

Fitted state and donor records are research artifacts. The general
production selection path, CLI behavior and SQLite schema need no setup for
these tests. Use independent result directories and the committed study
scripts to reproduce the report.
