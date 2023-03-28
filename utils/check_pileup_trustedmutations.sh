#!/bin/bash

source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
# def path and environment
conda activate default


echo "#########################"
echo "##### PATIENT 986 #######"
echo "#########################"

echo "######## EPHB2 ########"
echo " T -> G "
echo "chr1 pos 22895544 (hg38)"
samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r 1:23222037-23222037 /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-986_300316-CW/bcbio_work/align/NCC_CRC-986_300316-CW-T/NCC_CRC-986_300316-CW-T-sort.bam

echo "######## TP53 ########"
echo " T -> C "
echo "chr17 pos 7675076 (hg38)"
samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r 17:7578394-7578394 /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-986_300316-CW/bcbio_work/align/NCC_CRC-986_300316-CW-T/NCC_CRC-986_300316-CW-T-sort.bam

echo "######## SOX9 ########"
echo " C -> CGA "
echo "chr17 pos 72123617 to 72123619 (hg38)"
samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r 17:70119758-70119760 /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-986_300316-CW/bcbio_work/align/NCC_CRC-986_300316-CW-T/NCC_CRC-986_300316-CW-T-sort.bam

echo "######## PIK3CA ########"
echo " A -> G "
echo "chr3 pos 179218304 (hg38)"
samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r 3:178936092-178936092 /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-986_300316-CW/bcbio_work/align/NCC_CRC-986_300316-CW-T/NCC_CRC-986_300316-CW-T-sort.bam

echo "######## APC ########"
echo " C -> T "
echo "chr5 pos 112792446 (hg38)"
samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r 5:112128143-112128143 /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-986_300316-CW/bcbio_work/align/NCC_CRC-986_300316-CW-T/NCC_CRC-986_300316-CW-T-sort.bam

echo "#########################"
echo "##### PATIENT 1014 ######"
echo "#########################"

echo "######## AKT1 ########"
echo " C -> A "
echo "chr14 pos 104773077 (hg38)"
samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r 14:105239414-105239414 /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-1014_090516-CW/bcbio_work/align/NCC_CRC-1014_090516-CW-T/NCC_CRC-1014_090516-CW-T-sort.bam

echo "######## TP53 ########"
echo " CTG -> C "
echo "chr17 pos 7675111 (hg38)"
samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r 17:7578429-7578431 /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-1014_090516-CW/bcbio_work/align/NCC_CRC-1014_090516-CW-T/NCC_CRC-1014_090516-CW-T-sort.bam

echo "######## SOX9 ########"
echo " A -> T "
echo "chr17 pos 72122793 (hg38)"
samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r 17:70118934-70118934 /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-1014_090516-CW/bcbio_work/align/NCC_CRC-1014_090516-CW-T/NCC_CRC-1014_090516-CW-T-sort.bam

echo "######## PIK3CA ########"
echo " A -> G "
echo "chr3 pos 179234297 (hg38)"
samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r 3:178952085-178952085 /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-1014_090516-CW/bcbio_work/align/NCC_CRC-1014_090516-CW-T/NCC_CRC-1014_090516-CW-T-sort.bam

echo "######## APC ########"
echo " G -> T "
echo "chr5 pos 112839783 (hg38)"
samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r 5:112175480-112175480 /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-1014_090516-CW/bcbio_work/align/NCC_CRC-1014_090516-CW-T/NCC_CRC-1014_090516-CW-T-sort.bam
