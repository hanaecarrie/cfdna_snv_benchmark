#!/bin/bash

cd  /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/

export patient=$1

if [ ! -f ]; then
cat  /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/data/prepare_pooled_healthy/readfile_${$patient}*_genomad_*.txt  /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/data/prepare_pooled_healthy/readfile_${$patient}*_dbsnp.txt >>  /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/readfile_${$patient}_total.txt
fi

/mnt/projects/simngl/wgs/tools/anaconda3/bin/java -jar /mnt/projects/carriehc/cfDNA/utils/picard.jar FilterSamReads -I /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/healthy_chr22_merged-ready.bam -O /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-${$patient}/healthy_chr22_merged-ready_filter_${$patient}.bam --READ_LIST_FILE /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-${$patient}/readfile_${$patient}_total.txt --FILTER excludeReadList

/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-${$patient}/healthy_chr22_merged-ready_filter_${$patient}.bam
