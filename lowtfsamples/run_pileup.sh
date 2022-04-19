#!/bin/bash

source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
# def path and environment
conda activate default

cd /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/lowtfsamples/


echo "######## PATIENT 986 #########"

echo "## low tb sample 300316 ##"

#bash check_pileup.sh /mnt/projects/skanderupamj/wgs/data/clonal.hema/bams/WHC432/986_300316_P-T-sort.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg19.bed /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-986_300316-CW-T_deepWGS_pileup.txt /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa 

#bash check_pileup.sh /mnt/projects/huangwt/wgs/hanae/targeted/CCG_226_986.300316.reordered.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg38.bed /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-986_300316-CW-T_targeted_pileup.txt  /mnt/projects/ngbhs/cirqseq/references/hg38.fa 

echo "## high tb sample 100215 ##"

bash check_pileup.sh  /mnt/projects/skanderupamj/wgs/data/training/ready.bams/cfdna/NCC_CRC-986_100215-CW-T-ready.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg19.bed  /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-986_100215-CW-T_deepWGS_pileup.txt /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa  &

bash check_pileup.sh /mnt/projects/huangwt/wgs/hanae/targeted/CCG_226_986.110215.P.reordered.bam  /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg38.bed /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-986_100215-CW-T_targeted_pileup.txt /mnt/projects/ngbhs/cirqseq/references/hg38.fa &

echo "## high tb sample 261016 ##"

bash check_pileup.sh /mnt/projects/skanderupamj/wgs/data/training/ready.bams/cfdna/NCC_CRC-986_261016-CW-T-ready.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg19.bed  /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-986_261016-CW-T_deepWGS_pileup.txt /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa &

bash check_pileup.sh /mnt/projects/huangwt/wgs/hanae/targeted/CCG_226_986.261016.P.reordered.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg38.bed  /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-986_261016-CW-T_targeted_pileup.txt /mnt/projects/ngbhs/cirqseq/references/hg38.fa &



echo "######## PATIENT 1014 #########"

echo "## low tb sample 090516 ##"

bash check_pileup.sh /mnt/projects/skanderupamj/wgs/data/clonal.hema/bams/WHC433/1014_090516_P-T-sort.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg19.bed /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-1014_090516-CW-T_deepWGS_pileup.txt  /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa &

bash check_pileup.sh /mnt/projects/huangwt/wgs/hanae/targeted/CCG_226_1014.090516.reordered.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg38.bed /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-1014_090516-CW-T_targeted_pileup.txt /mnt/projects/ngbhs/cirqseq/references/hg38.fa &

echo "## high tb sample 180816 ##"

bash check_pileup.sh /mnt/projects/skanderupamj/wgs/data/training/ready.bams/cfdna/NCC_CRC-1014_180816-CW-T-ready.bam  /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg19.bed  /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-1014_180816-CW-T_deepWGS_pileup.txt /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa &

bash check_pileup.sh /mnt/projects/zhug/cfDNA/vaf-calling-crc/sortpanel-1014_180816.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg38.bed  /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-1014_180816-CW-T_targeted_pileup.txt  /mnt/projects/ngbhs/cirqseq/references/hg38.fa &

echo "## high tb sample 110116 ##"

bash check_pileup.sh /mnt/projects/skanderupamj/wgs/data/training/ready.bams/cfdna/NCC_CRC-1014_110116-CW-T-ready.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg19.bed /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-1014_110116-CW-T_deepWGS_pileup.txt /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa &

bash check_pileup.sh /mnt/projects/huangwt/wgs/hanae/targeted/CCG_226_1014.110116.reordered.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg38.bed  /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-1014_110116-CW-T_targeted_pileup.txt /mnt/projects/ngbhs/cirqseq/references/hg38.fa &



echo "######## PATIENT 123 #########"

echo "## low tb sample 121115 ##"

bash check_pileup.sh  /mnt/projects/skanderupamj/wgs/data/cfdna.crc/cfDNA.20211019/NCC_CRC-123_121115_P/bcbio_work/align/NCC_CRC-123_121115_P-T/NCC_CRC-123_121115_P-T-sort.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg19.bed NCC_CRC-123_121115-CW-T_deepWGS_pileup.txt  /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa &

bash check_pileup.sh /mnt/projects/huangwt/wgs/hanae/targeted/CCG_226_123.121115.reordered.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg38.bed  /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-123_121115-CW-T_targeted_pileup.txt /mnt/projects/ngbhs/cirqseq/references/hg38.fa &

echo "## high tb sample 310715 ##"

bash check_pileup.sh /mnt/projects/skanderupamj/wgs/data/cfdna.crc/cfDNA.20211019/NCC_CRC-123_310715_P/bcbio_work/align/NCC_CRC-123_310715_P-T/NCC_CRC-123_310715_P-T-sort.bam  /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg19.bed  /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-123_310715-CW-T_deepWGS_pileup.txt  /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa &

bash check_pileup.sh /mnt/projects/huangwt/wgs/hanae/targeted/CCG_226_123.310715.reordered.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg38.bed  /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-123_310715-CW-T_targeted_pileup.txt /mnt/projects/ngbhs/cirqseq/references/hg38.fa &

echo "## high tb sample 070116 ##"

bash check_pileup.sh /mnt/projects/skanderupamj/wgs/data/cfdna.crc/cfDNA.20211019/NCC_CRC-123_070116_P/bcbio_work/align/NCC_CRC-123_070116_P-T/NCC_CRC-123_070116_P-T-sort.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg19.bed /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-123_070116-CW-T_deepWGS_pileup.txt /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa &

bash check_pileup.sh /mnt/projects/huangwt/wgs/hanae/targeted/CCG_226_123.0701.reordered.bam /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/extdata/Cancer226-targets_hg38.bed /mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-123_070116-CW-T_targeted_pileup.txt /mnt/projects/ngbhs/cirqseq/references/hg38.fa 

