#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default

export patientid_date=$1
export snp_database=$2
export chunk_start=$3
export chunk_end=$4
export path_data=$5

# STEP 1: preprare reads to remove

python3 ~/cfdna_snv_benchmark/analysis/preprocess_preparereadstoremove.py $patientid_date $snp_database $chunk_start $chunk_end $path_data


