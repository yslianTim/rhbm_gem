#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="build"
if [[ $# -gt 0 && "${1}" != --* ]]; then
    BUILD_DIR="${1}"
    shift
fi

BASELINE_FILE=""
UPDATE_BASELINE=0
while [[ $# -gt 0 ]]; do
    case "${1}" in
        --baseline)
            if [[ $# -lt 2 ]]; then
                echo "--baseline requires a file path" >&2
                exit 1
            fi
            BASELINE_FILE="${2}"
            shift 2
            ;;
        --update-baseline)
            UPDATE_BASELINE=1
            shift
            ;;
        *)
            echo "Unknown argument: ${1}" >&2
            exit 1
            ;;
    esac
done

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

REPORT_FILE="$(mktemp -t rhbm_clang_tidy.XXXXXX)"
cleanup() {
    rm -f "${REPORT_FILE}"
}
trap cleanup EXIT

set +e
clang-tidy -p "${BUILD_DIR}" "${FILES[@]}" 2>&1 | tee "${REPORT_FILE}"
CLANG_TIDY_STATUS="${PIPESTATUS[0]}"
set -e

if [[ "${CLANG_TIDY_STATUS}" -ne 0 ]]; then
    exit "${CLANG_TIDY_STATUS}"
fi

if [[ -n "${BASELINE_FILE}" ]]; then
    BASELINE_ARGS=(--report "${REPORT_FILE}" --baseline "${BASELINE_FILE}")
    if [[ "${UPDATE_BASELINE}" -eq 1 ]]; then
        BASELINE_ARGS+=(--update-baseline)
    fi
    python3 scripts/check_clang_tidy_baseline.py "${BASELINE_ARGS[@]}"
fi

echo "clang-tidy check passed for ${#FILES[@]} files."
