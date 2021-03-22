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

if [ ! -d $path_data/CRC-${patient} ] ; then
  mkdir $path_data/CRC-${patient}
fi

cd  $path_data/CRC-${patient} || exit

python3 /mnt/projects/carriehc/cfDNA/utils/bamsurgeon/bin/addsnv.py \
  -v $path_data/CRC-${patient}/varfile_snv_${patient}_dbsnp_${chunk}.bed \
  -f $path_data/CRC-${patient}/healthy_chr22_merged-ready_filter.bam \
  -r /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa \
  -o $path_data/CRC-${patient}-chunks/healthy_chr22_merged-ready_${patient}_filter_snv.bam \
  --mindepth 0 --maxdepth 200000 --ignoresnp --ignoreref --force --tagreads  \
  --picardjar /mnt/projects/carriehc/cfDNA/utils/picard.jar --aligner mem \
  --tmpdir $path_data/CRC-${patient}-chunks/addsnv.tmp --skipmerge -p $ncpus --seed 1


