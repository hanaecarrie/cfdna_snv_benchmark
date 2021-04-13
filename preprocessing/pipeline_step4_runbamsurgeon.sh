#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default

export type=$1
export patient=$2
export chunk=$3
export ncpus=$4
export path_data=$5
echo $type
echo $patient
echo $chunk
echo $ncpus

cd ~/ || exit

if [ $type == 'snv' ]

then python3 ~/bamsurgeon/bin/addsnv.py \
  -v $path_data/varfile_snv_${patient}_total.bed \
  -f $path_data/healthy_chr22_merged-ready_${patient}_filter.bam \
  -r $path_data/GRCh37/GRCh37.fa \
  -o $path_data/healthy_chr22_merged-ready_${patient}_filter_snv.bam \
  --mindepth 0 --maxdepth 100000 --ignoresnp --ignoreref --force --insane --coverdiff 1 --tagreads  \
  --picardjar $path_data/picard.jar --aligner mem \
  --tmpdir addsnv_${patient}.tmp -p $ncpus --seed 1

samtools index $path_data/healthy_chr22_merged-ready_${patient}_filter_snv.bam

fi


if [ $type == 'indel' ]
then python3 ~/bamsurgeon/bin/addindel.py \
  -v $path_data/varfile_indel_${patient}_total.bed \
  -f $path_data/healthy_chr22_merged-ready_${patient}_filter_snv.bam \
  -r $path_data/GRCh37/GRCh37.fa \
  -o $path_data/healthy_chr22_merged-ready_${patient}_filter_snv_indel.bam \
  --mindepth 0 --maxdepth 100000 --ignoresnp --ignoreref --force --insane --coverdiff 1 --tagreads  \
  --picardjar $path_data/picard.jar --aligner mem \
  --tmpdir addsnv_${patient}.tmp -p $ncpus --seed 1

samtools index $path_data/healthy_chr22_merged-ready_${patient}_filter_snv_indel.bam

fi
