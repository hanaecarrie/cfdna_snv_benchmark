#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input a hg38 bedile, a liftover chain from hg38 to hg19,
# a path to output hg19 bedfile and to output unlifted bedfile, and a condapath.
# lifts over the coordinates of the hg38 bedfile to hg19 reference human genome.

if [ $# == 0 ]; then
    echo "Usage: $0 bedfilehg38 bedfilehg19 hg18tohg19overchain unlifted condapath"
    echo "* bedfilehg38: string. full path to the input bedfile in hg38."
    echo "* bedfilehg19: string. full path to the output bedfile in hg19."
    echo "* hg18tohg19overchain: string. full path to the liftover chain."
    echo "* unlifted: string. full path to the output unlifted bed."
    echo "* condapath: string. full path to the conda.sh file to source."
    echo "Example:"
    echo "$ bash $0 Cancer226-targets_hg38.bed /mnt/projects/zhug/cfDNA/vaf-calling-crc/hg38ToHg19.over.chain.gz \
    Cancer226-targets_hg19.bed unlifted226gene.bed /apps/anaconda3-individual-edition/2020.11/etc/profile.d/conda.sh"
    exit 1
fi

export bedfilehg38=$1
export bedfilehg19=$2
export hg38tohg19overchain=$3
export unlifted=$4
export condapath=$5
echo "Input bedfile in hg38: " $bedfilehg38
echo "Output bedfile in hg19: " $bedfilehg19
echo "Liftover chain hg38tohg19: " $hg38tohg19overchain
echo "Output unlifted bedfile: " $unlifted
echo "Path to conda file: " $condapath

# activate conda environment
source $condapath
conda activate default

liftOver $bedfilehg38 $hg38tohg19overchain $bedfilehg19 $unlifted

