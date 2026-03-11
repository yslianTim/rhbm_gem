# PotentialAnalysisCommand

Scaffold generated for CLI command `potential_analysis`.

## Registration Checklist

1. Add `PotentialAnalysisCommand` into `src/core/internal/BuiltInCommandList.def`.
2. Include `PotentialAnalysisCommand.hpp` in `src/core/command/BuiltInCommandCatalog.cpp`.
3. Add source files to `src/rhbm_gem_sources.cmake` and `bindings/CMakeLists.txt`.
4. Wire binding function in `bindings/BindingHelpers.hpp` and `bindings/CoreBindings.cpp`.
5. Add tests to `tests/cmake/core_tests.cmake`.
