#!/usr/bin/env bash
set -euo pipefail

if ! command -v clang-format >/dev/null 2>&1; then
    echo "clang-format is not available on PATH" >&2
    exit 1
fi

TARGET_DIRS=(
    "src/core/painter"
    "src/core/internal/plot"
    "src/data/io/file"
    "src/data/internal/io/file"
)

FILES=()
while IFS= read -r file; do
    FILES+=("${file}")
done < <(
    find "${TARGET_DIRS[@]}" -type f \( -name '*.cpp' -o -name '*.hpp' \) | sort
)

for file in "${FILES[@]}"; do
    clang-format --dry-run --Werror "${file}"
done

echo "clang-format check passed for ${#FILES[@]} files."
