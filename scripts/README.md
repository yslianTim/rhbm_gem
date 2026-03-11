# Scripts Directory Policy

- Keep reusable project scripts under `scripts/`.
- Put machine-specific one-off scripts under `scripts/local/` (ignored by git).
- Local scripts must use environment variables for machine-specific paths.
- Do not hardcode absolute user paths in scripts.
- Use `scripts/run_local_template.sh` as the baseline template.
- Do not place scripts at repository root.

## Built-in Command Automation

- `scripts/command_scaffold.py`: generate a new command/binding/test/doc scaffold.
  - add `--wire` to also update built-in registration + manifest files.
  - add `--wire --strict` to fail-fast on anchor drift with repair hints.
- `scripts/check_builtin_command_sync.py`: validate `BuiltInCommandList.def` sync with
  headers/sources/bindings/tests/CMake manifests, docs placeholders, and order consistency.
- `scripts/run_clang_format_check.sh`: check formatting for painter/parser directories (check-only).
- `scripts/run_clang_tidy_check.sh`: run clang-tidy for painter/parser source files (check-only),
  and optionally compare warning counts with `scripts/clang_tidy_baseline.json`.
- `scripts/check_clang_tidy_baseline.py`: compare/update clang-tidy warning baseline.
- `scripts/run_ctest_with_classification.sh`: run ctest and print contract/runtime failure split on failure.
- `scripts/classify_ctest_failures.py`: classify failed CTest entries by `intent:*` labels.
