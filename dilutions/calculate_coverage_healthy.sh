#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default

export cov=$(samtools depth -r 22 -a  /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/healthy/healthy_chr22_merged-ready.bam | awk '{sum+=$3} END {print sum/NR}')
echo $cov
echo $cov >> /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/healthy/coverage_healthy_chr22.txt
