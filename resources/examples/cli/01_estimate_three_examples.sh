#!/bin/bash

# Define working directory (Require user input)
work_dir=$1 # Check if work_dir is provided
if [ -z "$work_dir" ]; then
    echo "Usage: $0 <work_dir>"
    exit 1
fi

mkdir -p "$work_dir"
echo "Working directory: $work_dir"
cd "$work_dir" || { echo "Error: Failed to enter directory $work_dir"; exit 1; }

# Create data/ folder under work_dir and download required model/map files
data_dir="$work_dir/data"
mkdir -p "$data_dir"
echo "Downloading data files under $data_dir"
cd "$data_dir" || { echo "Error: Failed to enter directory $data_dir"; exit 1; }

wget https://files.rcsb.org/download/6Z6U.cif
wget https://files.rcsb.org/download/6DRV.cif
wget https://files.rcsb.org/download/9GXM.cif
wget https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-11103/other/emd_11103_additional.map.gz
wget https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-8908/other/emd_8908_additional.map.gz
wget https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-51667/other/emd_51667_additional_1.map.gz

gunzip emd_11103_additional.map.gz
gunzip emd_8908_additional.map.gz
gunzip emd_51667_additional_1.map.gz

# Run RHBM-GEM
cd "$work_dir" || { echo "Error: Failed to enter directory $work_dir"; exit 1; }
executable="./build/bin/RHBM-GEM"
database="$data_dir/database.sqlite"
alpha_train=false
jobs=4

key_tag_params=("6z6u" "6drv" "9gxm")
model_params=("6Z6U.cif" "6DRV.cif" "9GXM.cif")
map_params=("emd_11103_additional.map" "emd_8908_additional.map" "emd_51667_additional_1.map")

for i in "${!model_params[@]}"; do
    map_name="$data_dir/${map_params[i]}"
    model_name="$data_dir/${model_params[i]}"
    $executable potential_analysis -d $database -m $map_name -a $model_name -k ${key_tag_params[i]} -j $jobs --training-alpha $alpha_train
done


# Run Gaus Estimation Result Dumping
key_tag_list=$(IFS=,; echo "${key_tag_params[*]}")
echo "Output gaussian estimation results of following models: $key_tag_list"
$executable result_dump -o "$work_dir" -d "$database" -k "$key_tag_list" -p 2
