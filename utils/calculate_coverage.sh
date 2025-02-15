#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input a bamfile, a path to ichorCNA extdata folder, a condapath, and a chromosome.
# step 1: create read count WIG file using readCounter
# step 2: run ichorCNA while restricting training and analysis to the given chromosome.

if [ $# == 0 ]; then
    echo "Usage: $0 bamfile condapath"
    echo "* bamfile: string. full path to the cfDNA bam file with bam.bai indexfile located in the same directory."
    echo "* condapath: string. full path to the conda.sh file to source."
    echo "Example:"
    echo "$ bash $0 /home/users/astar/gis/carriehc/scratch-LOCCG/carriehc/data/mixtures/mixtures_chr1/mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T/mixture_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_80x/mixture_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_80x.bam \
     /apps/anaconda3-individual-edition/2020.11/etc/profile.d/conda.sh"
    exit 1
fi

export bamfile=$1
export condapath=$2
echo "Bamfile: "$bamfile
echo "Path to conda file: " $condapath

# activate conda environment
source $condapath
conda activate default

export cov=$(samtools depth -a $bamfile | awk '{sum+=$3} END {print sum/NR}')
echo $cov
if [ ! -f $(dirname $bamfile)/coverage_$(basename $bamfile).txt ] ; then
  echo $cov >> $(dirname $bamfile)/coverage_$(basename $bamfile).txt ; 
fi

