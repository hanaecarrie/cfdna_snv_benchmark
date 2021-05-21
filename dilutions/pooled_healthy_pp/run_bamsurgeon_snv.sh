#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default

export patient=$1
export chunk=$2
export ncpus=$3
echo $patient
echo $chunk
echo $ncpus

cd /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-${patient}-chunks

python3 /mnt/projects/carriehc/cfDNA/utils/bamsurgeon/bin/addsnv.py -v /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-${patient}-chunks/varfile_snv_${patient}_dbsnp_${chunk}.bed -f /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-${patient}/healthy_chr22_merged-ready_filter.bam -r /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -o /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-${patient}-chunks/healthy_chr22_merged-ready_${patient}_filter_snv.bam --mindepth 0 --maxdepth 20000 --ignoresnp --ignoreref --force --tagreads  --picardjar /mnt/projects/carriehc/cfDNA/utils/picard.jar --aligner mem --tmpdir /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-${patient}-chunks/addsnv.tmp --skipmerge -p $ncpus --seed 1


