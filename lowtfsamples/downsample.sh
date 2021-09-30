#!/bin/bash

# dilution factor 45
#/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b -s 0.0222222 -@ 8 -o /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/raw_data/NCC_CRC-986-300316-CW-T_0.0222222.bam /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-986_300316-CW/bcbio_work/align/NCC_CRC-986_300316-CW-T/NCC_CRC-986_300316-CW-T-sort.bam

# dilution factor 30
/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b -s 0.0333333 -@ 8 -o /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/raw_data/NCC_CRC-1014-090516-CW-T_0.0333333.bam /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-1014_090516-CW/bcbio_work/align/NCC_CRC-1014_090516-CW-T/NCC_CRC-1014_090516-CW-T-sort.bam
/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index -@ 8 /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/raw_data/NCC_CRC-1014-090516-CW-T_0.0333333.bam
