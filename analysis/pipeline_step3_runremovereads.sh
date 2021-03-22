#!/bin/bash

export patient=$1
export path_data=$2

cd  $path_data || exit

if [ ! -f ${path_data}/data/prepare_pooled_healthy/readfile_${patient}_total.txt ]; then
cat  $path_data/data/prepare_pooled_healthy/readfile_${patient}*_genomad_*.txt  $path_data/data/prepare_pooled_healthy/readfile_${patient}*_dbsnp.txt >>  $path_data/data/prepare_pooled_healthy/readfile_${patient}_total.txt
fi

/mnt/projects/simngl/wgs/tools/anaconda3/bin/java -jar /mnt/projects/carriehc/cfDNA/utils/picard.jar FilterSamReads -I $path_data/data/healthy_chr22_merged-ready.bam -O $path_data/data/pooled_healthy/healthy_chr22_merged-ready_filter_${patient}.bam --READ_LIST_FILE $path_data/data/prepare_pooled_healthy/readfile_${patient}_total.txt --FILTER excludeReadList

/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $path_data/data/pooled_healthy/healthy_chr22_merged-ready_filter_${patient}.bam
