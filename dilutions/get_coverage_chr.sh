#!/bin/bash

# def path and environment
source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
conda activate default

export bamfile=$1
export chr=$2
export outputfile=$3

echo $bamfile
echo $chr
echo $outputfile

if [ ! -f ${bamfile}.bai ] ; then samtools index $bamfile ; fi

export cov=$(samtools depth -r $chr -a $bamfile | awk '{sum+=$3} END {print sum/NR}')
echo $cov
echo $cov >> $outputfile
