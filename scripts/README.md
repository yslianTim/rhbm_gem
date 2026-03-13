# Scripts Directory Policy

- User-facing execution templates under `scripts/user/`.
- Developer workflow automation under `scripts/developer/`.

## User Scripts

- `scripts/user/run_local_template.sh`: baseline local execution template for `potential_analysis`.


## Developer Scripts

- `scripts/developer/generate_command_artifacts.py`: generate all command derived artifacts from `CommandList.def`
  (generated `CommandId` entries, catalog include/entry fragments, binding-name traits,
  command CMake source lists, docs generated blocks).
- `scripts/developer/command_scaffold.py`: generate a new command/binding/test/doc scaffold.
  - add `--wire` to also update `CommandList.def` and regenerate derived artifacts.
  - add `--wire --strict` to fail-fast if manifest generation fails.
- `scripts/developer/check_command_sync.py`: validate `CommandList.def` sync with
  generated artifacts plus required headers/sources/bindings/tests surface files.
- `scripts/developer/run_clang_format_check.sh`: check formatting for painter/parser directories (check-only).
- `scripts/developer/run_clang_tidy_check.sh`: run clang-tidy for painter/parser source files (check-only),
  and optionally compare warning counts with `scripts/developer/clang_tidy_baseline.json`.
- `scripts/developer/check_clang_tidy_baseline.py`: compare/update clang-tidy warning baseline.
- `scripts/developer/run_ctest_with_classification.sh`: run ctest and print contract/runtime failure split on failure.
- `scripts/developer/classify_ctest_failures.py`: classify failed CTest entries by `intent:*` labels.
