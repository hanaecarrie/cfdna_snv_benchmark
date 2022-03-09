#!/bin/bash

# def path and environment
source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
conda activate default

export currentdate=$(date +"%y-%m-%d_%H-%M-%S")
echo $currentdate
if [ -f /mnt/projects/skanderupamj/wgs/data/cfdna.crc/benchmark/patients_wgs.tsv ] ; then mv /mnt/projects/skanderupamj/wgs/data/cfdna.crc/benchmark/patients_wgs.tsv /mnt/projects/skanderupamj/wgs/data/cfdna.crc/benchmark/samplestsv/patients_wgs_${currentdate}.tsv ; fi

cp /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/bcbio/bcbio_config.tsv /mnt/projects/skanderupamj/wgs/data/cfdna.crc/benchmark/patients_wgs.tsv
chmod 715  /mnt/projects/skanderupamj/wgs/data/cfdna.crc/benchmark/patients_wgs.tsv

cd /mnt/projects/skanderupamj/wgs/data/cfdna.crc/benchmark

python /home/simngl/tools/skandlab_pipeline/wgs_pipeline/scripts/start_skandlab_bcbio.py -t patients_wgs.tsv -d /mnt/projects/skanderupamj/wgs/data/cfdna.crc/benchmark

