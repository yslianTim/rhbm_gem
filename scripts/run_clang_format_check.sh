#!/usr/bin/env bash
set -euo pipefail

if ! command -v clang-format >/dev/null 2>&1; then
    echo "clang-format is not available on PATH" >&2
    exit 1
fi

FILES=(
    "src/core/painter/PotentialPlotBuilder.cpp"
    "src/core/painter/PotentialPlotBuilderSequenceGraphs.cpp"
    "src/core/painter/PotentialPlotBuilderPlotFunctions.cpp"
    "src/core/painter/RootPlotRenderBackend.cpp"
    "src/core/internal/plot/IPlotDataBuilder.hpp"
    "src/core/internal/plot/QueryPlotDataBuilder.hpp"
    "src/core/internal/plot/IPlotRenderBackend.hpp"
    "src/core/internal/plot/RootPlotRenderBackend.hpp"
    "src/data/io/file/CifFormat.cpp"
    "src/data/io/file/ICifCategoryParser.cpp"
    "src/data/internal/io/file/CifFormat.hpp"
    "src/data/internal/io/file/ICifCategoryParser.hpp"
)

for file in "${FILES[@]}"; do
    clang-format --dry-run --Werror "${file}"
done

echo "clang-format check passed for ${#FILES[@]} targeted files."
