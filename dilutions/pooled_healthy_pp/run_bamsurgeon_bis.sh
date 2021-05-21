#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default

# patient 986

python3 /mnt/projects/carriehc/cfDNA/utils/bamsurgeon/bin/addsnv.py -v /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-986/varfile_snv_986_dbsnp.bed -f /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-986/healthy_chr22_merged-ready_filter.bam -r /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -o /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-986/healthy_chr22_merged-ready_986_filter_snv_bis.bam -p 16 --mindepth 0 --maxdepth 200000 --ignoresnp --ignoreref --force --seed 1 --tagreads --tmpdir /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-986/addsnv_bis.tmp

#/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-986/healthy_chr22_merged-ready_986_filter_snv.bam

#python3 /mnt/projects/carriehc/cfDNA/utils/bamsurgeon/bin/addindel.py -v /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-986/varfile_indel_986_dbsnp.bed -f /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-986/healthy_chr22_merged-ready_986_filter_snv.bam -r /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -o /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-986/healthy_chr22_merged-ready_986_filter_snvindel.bam -p 16 --mindepth 0 --maxdepth 200000 --ignorepileup --force --seed 1 --tagreads --tmpdir /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-986/addindel.tmp

#/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/CRC-986/healthy_chr22_merged-ready_986_filter_snvindel.bam
