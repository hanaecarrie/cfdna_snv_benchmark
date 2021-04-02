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

cd ~/ || exit

python3 ~/bamsurgeon/bin/addsnv.py \
  -v $path_data/varfile_snv_${patient}_total.bed \
  -f $path_data/healthy_chr22_merged-ready_${patient}_filter.bam \
  -r $path_data/GRCh37/GRCh37.fa \
  -o $path_data/healthy_chr22_merged-ready_${patient}_filter_snv.bam \
  --mindepth 0 --maxdepth 100000 --ignoresnp --ignoreref --force --insane --covdiff 1 --tagreads  \
  --picardjar $path_data/picard.jar --aligner mem \
  --tmpdir addsnv_${patient}.tmp -p $ncpus --seed 1

samtools index $path_data/healthy_chr22_merged-ready_${patient}_filter_snv.bam
