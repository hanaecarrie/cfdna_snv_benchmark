#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default

export patientid_date=$1
export path_data=$2

# STEP 2: preprare bamsurgeon input

python3 ~/cfdna_snv_benchmark/analysis/preprocess_preparebamsurgeon.py $patientid_date $path_data