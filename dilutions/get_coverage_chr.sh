#!/bin/bash

# def path and environment
source /mnt/projects/carriehc/anaconda3/etc/profile.d/conda.sh
conda activate default

export bamfile=$1
export chr=$2
export outputfile=$3

export cov=$(samtools depth -r $chr -a  $bamfile | awk '{sum+=$3} END {print sum/NR}')
echo $cov
echo $cov >> $outputfile
