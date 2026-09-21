# Joint component estimator v1 acceptance and freeze

## Frozen baseline

**Accepted code commit: `7c0b7a4543fdc43390ca59b3aa3b5d9267fbaac7`.** All four lanes below ran from this same
clean tracked tree. The later documentation-only commit records acceptance; it
is not represented as another fully tested code revision. The runtime source
fingerprint is `e601fd9b727d59ec9668e08a1a3949ddef2365eca859c688ae4737767435c450` in all four builds.

The joint estimator remains opt-in. This freezes the fixed-position,
all-non-hydrogen, structural 2.5 Angstrom support, equal-weight Guarded joint-LS
contract, not the project package version (which remains 2.0.0). No tag or package
release was created. There were no solver, objective, convergence-threshold or
frozen-fixture changes relative to `64127cca`.

## Validation

| Lane | Configuration and outcome |
| --- | --- |
| Default Release | All 25 CTests passed, including C++ groups, CLI, Python bindings and eight frozen runtime cases; repository guards passed. |
| Extended/offline Release | All 4 CTests passed: three offline tests and one extended regression covering three additional cases, including heterogeneous-168. |
| Testing-disabled install | BUILD_TESTING=OFF, offline/extended OFF, Python bindings OFF; installed library consumer compiled, executed and exported schema-2 JSON/CSV successfully. |
| Source package without Git | Archive of the accepted commit extracted outside the checkout, with no .git and no discoverable parent repository; testing-disabled build/install and consumer execution passed with the same source fingerprint. |

The host was macOS 26.6.2 arm64, AppleClang 21,
CMake 4.4.3 and Python 3.14.7. System dependencies included
Eigen 5.0.1, Boost 1.92.0, SQLite 3.51.0, CLI11 2.7.2, pybind11 3.1.0 and GTest
1.18.0; OpenMP 5.1 and ROOT 6.40.04 were enabled. UMAP and experimental commands
were enabled in all lanes. The same installed dependencies and cached pinned UMAP
sources were used for the source-package check; this is not a claim of a vendored,
dependency-free distribution. Some build/test phases overlapped; elapsed test
times are diagnostics, not performance measurements.

The [machine-readable receipt](joint-component-v1-acceptance.json) records exact
commands, resolved options, per-lane source/configuration/build fingerprints,
test names, log hashes, source archive hash and the immutable fixture hashes.
Compact logs: [default CTests](figures/joint-component-v1-acceptance/default-tests.txt),
[extended/offline CTests](figures/joint-component-v1-acceptance/offline-tests.txt),
[repository guards](figures/joint-component-v1-acceptance/default-guards.txt).
Full local build and installation logs remain under `build/joint-v1-acceptance/`.

To reproduce, use the recorded commands with paths adjusted for the new checkout
and dependency cache. Configure separate Release directories for default tests,
extended/offline tests, and testing-disabled installation. Build `tests_all`, run
all default CTests, then run `ctest -L 'joint:extended|joint:offline'` in the offline
build. Run `lint_repo` in the checkout and `lint_install_smoke` in both installation
builds. The latter now executes the consumer and fails on a nonzero exit code.
Create the no-Git source copy using `git archive` of the accepted commit, not a
copy of a potentially modified working tree.

## Metadata and unit acceptance

Production joint JSON is schema 2 inside unchanged SQLite v17. Older joint JSON
is rejected without migration; non-joint records are unaffected. C++ consumers
must rebuild. The [command contract](commands/potential_analysis.md#provenance-and-map-units-joint-json-schema-2)
defines the saved fields and conversion rules.

The new regression coverage verifies actual normalization divisors, including
sample SD=1, nonzero means, constant maps, disabled normalization and simulation
bypass. The operation remains division by sample standard deviation with no mean
subtraction. Unknown in-memory provenance stays null. The objective's parent
observation scale is separate.

Raw file SHA-256 values identify inputs, and capture stamps the existing library
version/source/configuration/build identity. Tests independently compare file
hashes, change bytes at the same input path, and preserve stored identity during
round trips. CLI, C++ request and Python request checks delete their temporary
model/map files before export. Invalid metadata rolls back storage; old document
reads fail without changing database bytes. Missing/unconverged states and offline
`NotRun` remain unchanged.

For fitted map unit U, A is U*Angstrom^3, C is U*Angstrom and B is Angstrom. A
runtime-converged synthetic scaling control verifies that multiplying fitted A/C by the
normalization divisor restores input-map scale while preserving B and prediction.
CSV retains fitted values and its existing columns; its companion JSON provides
units and provenance. No physical map calibration is inferred.

## Numerical limits retained

Existing frozen-case acceptance remains scaled 1e-10 for converged identifiable
parameters and 1e-12 for normalized objectives, with same-state backend checks
where applicable. No historical all-start matrix or oracle regeneration was run.
Weak-1e-4 retains its offline derivative resolution limit; active-a retains
local-correction failure; zero-signal retains width-identification failure;
duplicate retains no usable state. A changed active face for a nonconverged
endpoint is not promoted to a new certificate. Runtime convergence is distinct
from search completion, usable state and offline regular certification.

This evidence supports engineering acceptance on the recorded macOS Release
configuration only. It does not establish cross-platform behavior, real-data
accuracy, global optimality or end-to-end resource limits. Further runtime,
test or build changes require a new acceptance baseline.
