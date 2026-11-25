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
  python3 tests/research/buffer/half_validate.py STUDY RUNNER TAG
  python3 tests/research/buffer/truth.py STUDY prepare
  RUNNER truth STUDY/truth-fixtures.json STUDY/truth-results.jsonl
  python3 tests/research/buffer/truth.py STUDY summarize
  python3 tests/research/buffer/truth.py STUDY workflow
  python3 -m unittest discover -s tests/research/buffer -v

run.py serializes full fits. For a fixed-domain run the output calibration
argument derives both width calibration and the reference norm from that
run's own input and selected outputs. The training runner does not open
the validation map. compare.py scores supported MRC layouts using raw units
through the reported training divisor. Archived baseline projections are
recomputed on recorded geometry from certified state and Halo medians.

The C++ truth mode applies the original controlled perturbations to a
four-contributor fixture and uses the existing native constrained solver
with the documented width interval. Its oracle prediction metric supplies
known external Halo truth after fitting. The 6Z6U simulation experiment
checks the complete workflow under the standard donor qualification rule.

Fitted state and donor records are research artifacts. The general
production selection path, CLI behavior and SQLite schema need no setup for
these tests. Use independent result directories and the committed study
scripts to reproduce the report.
