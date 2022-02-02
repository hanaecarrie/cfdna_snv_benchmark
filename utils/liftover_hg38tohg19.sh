#!/bin/bash

# def path and environment
source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
conda activate default

cd /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/utils

/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/liftOver Cancer226-targets_hg38.bed /mnt/projects/zhug/cfDNA/vaf-calling-crc/hg38ToHg19.over.chain.gz Cancer226-targets_hg19.bed unlifted226gene.bed

