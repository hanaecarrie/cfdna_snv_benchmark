#!/bin/bash

export patient=$1
export path_data=$2

cd  $path_data || exit

if [ ! -f ${path_data}/data/prepare_pooled_healthy/readfile_${patient}_total.txt ]; then
cat  $path_data/data/prepare_pooled_healthy/readfile_${patient}*_genomad_*.txt  $path_data/data/prepare_pooled_healthy/readfile_${patient}*_dbsnp.txt >>  $path_data/data/prepare_pooled_healthy/readfile_${patient}_total.txt
fi

/mnt/projects/simngl/wgs/tools/anaconda3/bin/java -jar /mnt/projects/carriehc/cfDNA/utils/picard.jar FilterSamReads -I $path_data/data/healthy_chr22_merged-ready.bam -O $path_data/data/pooled_healthy/healthy_chr22_merged-ready_${patient}_filter.bam --READ_LIST_FILE $path_data/data/prepare_pooled_healthy/readfile_${patient}_total.txt --FILTER excludeReadList

/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $path_data/data/pooled_healthy/healthy_chr22_merged-ready_${patient}_filter.bam

# test the program worked
echo "start testing..."
read_id="$(shuf -n 1 $path_data/data/prepare_pooled_healthy/readfile_${patient}_total.txt)"
export read_id
echo "pick a random read to filter out in the readfile"
echo $read_id
echo "test read ID is present in original bam"
samtools view $path_data/data/healthy_chr22_merged-ready.bam | grep $read_id
echo "test read ID is not present in processed bam"
samtools view $path_data/data/pooled_healthy/healthy_chr22_merged-ready_${patient}_filter.bam | grep $read_id
echo "...done testing"

