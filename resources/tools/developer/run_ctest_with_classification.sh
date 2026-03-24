#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="${1:-build}"
if [[ $# -gt 0 ]]; then
    shift
fi

CTEST_ARGS=("$@")
if [[ "${#CTEST_ARGS[@]}" -eq 0 ]]; then
    CTEST_ARGS=(-j2)
fi

LOG_FILE="$(mktemp -t rhbm_ctest_log.XXXXXX)"
cleanup() {
    rm -f "${LOG_FILE}"
}
trap cleanup EXIT

set +e
ctest --test-dir "${BUILD_DIR}" --output-on-failure "${CTEST_ARGS[@]}" 2>&1 | tee "${LOG_FILE}"
CTEST_STATUS="${PIPESTATUS[0]}"
set -e

if [[ "${CTEST_STATUS}" -eq 0 ]]; then
    echo "ctest passed."
    exit 0
fi

python3 resources/tools/developer/classify_ctest_failures.py \
    --build-dir "${BUILD_DIR}" \
    --ctest-log "${LOG_FILE}"
exit "${CTEST_STATUS}"
