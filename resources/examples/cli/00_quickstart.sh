#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
#
# Script:
#   00_quickstart.sh
#
# Purpose:
#   Download one example model/map pair (PDB-6Z6U/EMD-11103),
#   run `potential_analysis`, and then export the estimation results with `result_dump`.
#
# Audience / Platform:
#   New RHBM-GEM CLI users on macOS or Linux.
#   This script is not intended for Windows environments.
#
# Prerequisites:
#   Resolve `RHBM-GEM` via `--executable`, `RHBM_GEM_EXECUTABLE`, or `PATH`.
#   Require `curl` or `wget` for downloads.
#   Use `--skip-download` for offline reruns after staging files under `<workdir>/data`.
#
# Inputs / Outputs:
#   Download and read inputs from `<workdir>/data`.
#   Write generated outputs under `<workdir>`.
#
# References:
#   docs/user/getting-started.md (CLI Quickstart)
#   resources/README.md
#   resources/examples/cli/common.sh
#   https://www.rcsb.org/structure/6Z6U
set -euo pipefail

script_dir="$(
    cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 &&
        pwd -P
)"
# shellcheck source=resources/examples/cli/common.sh
source "${script_dir}/common.sh"

usage() {
    cat <<'EOF'
Usage:
  bash resources/examples/cli/00_quickstart.sh <work_dir> [--executable /path/to/RHBM-GEM] [--skip-download]
  bash resources/examples/cli/00_quickstart.sh --workdir <work_dir> [--executable /path/to/RHBM-GEM] [--skip-download]

Run the RHBM-GEM CLI quickstart on one example model/map pair.
EOF
    rhbm_gem_usage_example
}

require_file() {
    local path="$1"

    if [ ! -f "$path" ]; then
        echo "Error: Required input file not found: $path" >&2
        echo "Prepare the file manually or rerun without --skip-download." >&2
        exit 1
    fi
}

work_dir=""
explicit_executable=""
skip_download=false

while [ $# -gt 0 ]; do
    case "$1" in
        --workdir)
            if [ $# -lt 2 ]; then
                echo "Error: --workdir requires a path." >&2
                usage >&2
                exit 1
            fi
            work_dir="$2"
            shift 2
            ;;
        --executable)
            if [ $# -lt 2 ]; then
                echo "Error: --executable requires a path." >&2
                usage >&2
                exit 1
            fi
            explicit_executable="$2"
            shift 2
            ;;
        --skip-download)
            skip_download=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "Error: Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
        *)
            if [ -n "$work_dir" ]; then
                echo "Error: work_dir was already provided as '$work_dir'." >&2
                usage >&2
                exit 1
            fi
            work_dir="$1"
            shift
            ;;
    esac
done

if [ -z "$work_dir" ]; then
    usage >&2
    exit 1
fi

executable="$(rhbm_gem_resolve_executable "$explicit_executable" "${BASH_SOURCE[0]}")"

mkdir -p "$work_dir"
work_dir="$(
    cd "$work_dir" >/dev/null 2>&1 &&
        pwd -P
)"

echo "Working directory: $work_dir"
echo "Using executable: $executable"

data_dir="$work_dir/data"
mkdir -p "$data_dir"

model_name="$data_dir/6Z6U.cif"
map_name="$data_dir/emd_11103_additional.map"
compressed_map_name="${map_name}.gz"

if [ "$skip_download" = false ]; then
    echo "Downloading data files under $data_dir"
    (
        cd "$data_dir"
        rhbm_gem_download_file "https://files.rcsb.org/download/6Z6U.cif"
        curl -L "https://drive.usercontent.google.com/download?id=1MA7jQRZAKc0hGMLom_I2RAPcEqDc6ZSh&confirm=xxx" -o $map_name
        #rhbm_gem_download_file "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-11103/other/emd_11103_additional.map.gz"
        #gunzip -f "$compressed_map_name"
    )
else
    echo "Skipping downloads; expecting existing files under $data_dir"
fi

require_file "$model_name"
require_file "$map_name"

database="$data_dir/database.sqlite"
alpha_train=false
jobs=4
key_tag="6z6u"

"$executable" potential_analysis \
    -d "$database" \
    -m "$map_name" \
    -a "$model_name" \
    -k "$key_tag" \
    -j "$jobs" \
    --training-alpha "$alpha_train"

echo "Output gaussian estimation results of following model: $key_tag"
"$executable" result_dump \
    -o "$work_dir" \
    -d "$database" \
    -k "$key_tag" \
    -p 2
