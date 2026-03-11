#!/usr/bin/env bash
set -euo pipefail

: "${RHBM_GEM_EXECUTABLE:=./build/bin/RHBM-GEM}"
: "${RHBM_GEM_DATA_DIR:=${HOME}/.rhbmgem/data}"
: "${RHBM_GEM_DATABASE_PATH:=${RHBM_GEM_DATA_DIR}/database.sqlite}"
: "${RHBM_GEM_MODEL_PATH:?set RHBM_GEM_MODEL_PATH to a model file (.cif/.pdb)}"
: "${RHBM_GEM_MAP_PATH:?set RHBM_GEM_MAP_PATH to a map file (.map/.mrc/.ccp4)}"
: "${RHBM_GEM_KEY_TAG:=model}"
: "${RHBM_GEM_JOBS:=4}"
: "${RHBM_GEM_TRAINING_ALPHA:=false}"

"${RHBM_GEM_EXECUTABLE}" potential_analysis \
  -d "${RHBM_GEM_DATABASE_PATH}" \
  -m "${RHBM_GEM_MAP_PATH}" \
  -a "${RHBM_GEM_MODEL_PATH}" \
  -k "${RHBM_GEM_KEY_TAG}" \
  -j "${RHBM_GEM_JOBS}" \
  --training-alpha "${RHBM_GEM_TRAINING_ALPHA}"
