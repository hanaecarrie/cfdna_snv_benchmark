#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default

export patient=$1
export chunk=$2
export ncpus=$3
export path_data=$4
echo $patient
echo $chunk
echo $ncpus

python3 /mnt/projects/carriehc/cfDNA/utils/bamsurgeon/bin/addsnv.py \
  -v $path_data/data/prepare_pooled_healthy/varfile_snv_${patient}_${chunk}.bed \
  -f $path_data/data/pooled_healthy/healthy_chr22_merged-ready_filter_${patient}.bam \
  -r /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa \
  -o $path_data/data/pooled_healthy/healthy_chr22_merged-ready_${patient}_filter_snv.bam \
  --mindepth 0 --maxdepth 200000 --ignoresnp --ignoreref --force --tagreads  \
  --picardjar /mnt/projects/carriehc/cfDNA/utils/picard.jar --aligner mem \
  --tmpdir $path_data/data/pooled_healthy/addsnv_${patient}.tmp --skipmerge -p $ncpus --seed 1


