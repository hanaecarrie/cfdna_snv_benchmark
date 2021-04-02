#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default

export patient=$1
export germline_vcf_name=$2
export snp_database=$3
export chunk_start=$4
export chunk_end=$5
export path_data=$6

# STEP 1: preprare reads to remove

python3 ~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.py $patient $germline_vcf_name $snp_database $chunk_start $chunk_end $path_data


