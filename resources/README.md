# Resources Directory

This directory groups non-core project assets under a single top-level path.

- User-facing examples live under `resources/examples/`.
- User execution templates live under `resources/tools/user/`.
- Developer workflow automation lives under `resources/tools/developer/`.

## User Examples

- `resources/examples/python/00_quickstart.py`: smoke-test the Python bindings.
- `resources/examples/python/01_end_to_end_from_test_data.py`: run the bundled end-to-end Python demo.

## User Scripts

- `resources/tools/user/run_local_template.sh`: baseline local execution template for `potential_analysis`.

## Developer Scripts

- `resources/tools/developer/command_scaffold.py`: generate a new command/binding/test/doc scaffold.
  - add `--wire` to also update `CommandList.def`.
  - add `--wire --strict` to fail-fast if manifest update fails.
- `resources/tools/developer/run_clang_format_check.sh`: check formatting for painter/parser directories (check-only).
- `resources/tools/developer/run_clang_tidy_check.sh`: run clang-tidy for painter/parser source files (check-only),
  and optionally compare warning counts with `resources/tools/developer/clang_tidy_baseline.json`.
- `resources/tools/developer/check_clang_tidy_baseline.py`: compare/update clang-tidy warning baseline.
- `resources/tools/developer/run_ctest_with_classification.sh`: run ctest and print contract/runtime failure split on failure.
- `resources/tools/developer/classify_ctest_failures.py`: classify failed CTest entries by `intent:*` labels.
