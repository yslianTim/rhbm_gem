#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="${1:-build}"

if ! command -v clang-tidy >/dev/null 2>&1; then
    echo "clang-tidy is not available on PATH" >&2
    exit 1
fi

if [[ ! -f "${BUILD_DIR}/compile_commands.json" ]]; then
    echo "compile_commands.json was not found under ${BUILD_DIR}" >&2
    echo "Run CMake with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON first." >&2
    exit 1
fi

mapfile -t FILES < <(
    find src/core/painter src/data/io/file -type f -name '*.cpp' | sort
)

if [[ "${#FILES[@]}" -eq 0 ]]; then
    echo "No source files found for clang-tidy check." >&2
    exit 1
fi

clang-tidy -p "${BUILD_DIR}" "${FILES[@]}"
echo "clang-tidy check passed for ${#FILES[@]} files."
