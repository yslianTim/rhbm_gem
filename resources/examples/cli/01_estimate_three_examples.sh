#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
#
# Script:
#   01_estimate_three_examples.sh
#
# Purpose:
#   Download three example model/map pairs (PDB-6Z6U/EMD-11103, PDB-6DRV/EMD-8908, PDB-9GXM/EMD-51667),
#   run `potential_analysis`, and then export three estimation results with `result_dump`.
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
#   resources/README.md
#   resources/examples/cli/common.sh
#   https://www.rcsb.org/structure/6Z6U
#   https://www.rcsb.org/structure/6DRV
#   https://www.rcsb.org/structure/9GXM
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
  bash resources/examples/cli/01_estimate_three_examples.sh <work_dir> [--executable /path/to/RHBM-GEM]
  bash resources/examples/cli/01_estimate_three_examples.sh --workdir <work_dir> [--executable /path/to/RHBM-GEM]
EOF
}

work_dir=""
explicit_executable=""

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
echo "Downloading data files under $data_dir"
(
    cd "$data_dir"
    rhbm_gem_download_file "https://files.rcsb.org/download/6Z6U.cif"
    rhbm_gem_download_file "https://files.rcsb.org/download/6DRV.cif"
    rhbm_gem_download_file "https://files.rcsb.org/download/9GXM.cif"
    curl -L "https://drive.usercontent.google.com/download?id=1MA7jQRZAKc0hGMLom_I2RAPcEqDc6ZSh&confirm=xxx" -o $map_name
    #rhbm_gem_download_file "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-11103/other/emd_11103_additional.map.gz"
    rhbm_gem_download_file "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-8908/other/emd_8908_additional.map.gz"
    rhbm_gem_download_file "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-51667/other/emd_51667_additional_1.map.gz"
    #gunzip -f emd_11103_additional.map.gz
    gunzip -f emd_8908_additional.map.gz
    gunzip -f emd_51667_additional_1.map.gz
)

database="$data_dir/database.sqlite"
alpha_train=false
jobs=4

key_tag_params=("6z6u" "6drv" "9gxm")
model_params=("6Z6U.cif" "6DRV.cif" "9GXM.cif")
map_params=("emd_11103_additional.map" "emd_8908_additional.map" "emd_51667_additional_1.map")

for i in "${!model_params[@]}"; do
    map_name="$data_dir/${map_params[i]}"
    model_name="$data_dir/${model_params[i]}"
    "$executable" potential_analysis \
        -d "$database" \
        -m "$map_name" \
        -a "$model_name" \
        -k "${key_tag_params[i]}" \
        -j "$jobs" \
        --training-alpha "$alpha_train"
done

key_tag_list="$(IFS=,; echo "${key_tag_params[*]}")"
echo "Output gaussian estimation results of following models: $key_tag_list"
"$executable" result_dump \
    -o "$work_dir" \
    -d "$database" \
    -k "$key_tag_list" \
    -p 2
