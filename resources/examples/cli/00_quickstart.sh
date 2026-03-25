#!/bin/bash

# Find executable
executable="./build/bin/RHBM-GEM"

# Define working directory (Require user input)
work_dir=$1 # Check if work_dir is provided
if [ -z "$work_dir" ]; then
    echo "Usage: $0 <work_dir>"
    exit 1
fi

mkdir -p "$work_dir"
echo "Working directory: $work_dir"
ln -s $executable "$work_dir/RHBM-GEM"
cd "$work_dir" || { echo "Error: Failed to enter directory $work_dir"; exit 1; }

# Create data/ folder under work_dir and download required model/map files
data_dir="$work_dir/data"
mkdir -p "$data_dir"
echo "Downloading data files under $data_dir"
cd "$data_dir" || { echo "Error: Failed to enter directory $data_dir"; exit 1; }

wget https://files.rcsb.org/download/6Z6U.cif
wget https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-11103/other/emd_11103_additional.map.gz
gunzip emd_11103_additional.map.gz

# Run RHBM-GEM
cd "$work_dir" || { echo "Error: Failed to enter directory $work_dir"; exit 1; }

database="$data_dir/database.sqlite"
alpha_train=false
jobs=4
key_tag="6z6u"
map_name="$data_dir/emd_11103_additional.map"
model_name="$data_dir/6Z6U.cif"
$work_dir/RHBM-GEM potential_analysis -d $database -m $map_name -a $model_name -k $key_tag -j $jobs --training-alpha $alpha_train

# Run Gaus Estimation Result Dumping
echo "Output gaussian estimation results of following model: $key_tag"
$work_dir/RHBM-GEM result_dump -o $work_dir -d $database -k $key_tag -p 2
