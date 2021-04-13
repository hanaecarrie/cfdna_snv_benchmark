#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default

export bamfile=$1
export outbam=$2
export mutsfolder=$3
export tagreads=$4
export seed=$5
echo $bamfile
echo $outbam
echo $mutsfolder
echo $tagreads
echo $seed

# STEP 5: merge bamsurgeon after addsnv

cd /mnt/projects/carriehc/cfDNA/utils/bamsurgeon

if [ $tagreads == 'tagreads' ] ;
then python3 ~/cfdna_snv_benchmark/analysis/preprocess_mergebamsurgeon.py -f $bamfile -o $outbam -m $mutsfolder --tagreads --seed $seed
else python3 ~/cfdna_snv_benchmark/analysis/preprocess_mergebamsurgeon.py -f $bamfile -o $outbam -m $mutsfolder --seed $seed

