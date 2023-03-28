#!/bin/bash

# def path and environment
source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
conda activate default

export bamfile=$1
export outputfile=$2
echo $bamfile
echo $outputfile

export cov=$(samtools depth -a $bamfile | awk '{sum+=$3} END {print sum/NR}')
echo $cov
if ! -f $outputfile ; then echo $cov >> $outputfile ; fi  

