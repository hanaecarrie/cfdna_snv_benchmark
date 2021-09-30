#!/bin/bash

source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
# def path and environment
conda activate default

cd /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/utils

echo "######## PATIENT 986 #########"
touch NCC_CRC-986_300316-CW-T_pileup.txt
while IFS="" read -r p || [ -n "$p" ]
do
  printf '%s\n' "$p"
  printf '%s\n' "$p" >> NCC_CRC-986_300316-CW-T_pileup.txt
  chr="$( cut -f1 <<< "$p" )"; chr=${chr#"chr"} ; echo "$chr"
  pos="$( cut -f2 <<< "$p" )"; echo "$pos"
  name="$( cut -f3 <<< "$p" )"; echo "$name"
  samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r $chr:$pos-$pos /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-986_300316-CW/bcbio_work/align/NCC_CRC-986_300316-CW-T/NCC_CRC-986_300316-CW-T-sort.bam >> NCC_CRC-986_300316-CW-T_pileup.txt
done < /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/utils/Cancer226-targets_hg19.bed

#echo "######## PATIENT 1014 #########"
#touch NCC_CRC-1014_090516-CW-T_pileup.txt
#while IFS="" read -r p || [ -n "$p" ]
#do
#  printf '%s\n' "$p"
#  printf '%s\n' "$p" >> NCC_CRC-1014_090516-CW-T_pileup.txt
#  chr="$( cut -f1 <<< "$p" )"; chr=${chr#"chr"} ; echo "$chr"
#  pos="$( cut -f2 <<< "$p" )"; echo "$pos"
#  name="$( cut -f3 <<< "$p" )"; echo "$name"
#  samtools mpileup -f /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa -r $chr:$pos-$pos /mnt/projects/skanderupamj/wgs/data/cfdna.crc/NCC_CRC-1014_090516-CW/bcbio_work/align/NCC_CRC-1014_090516-CW-T/NCC_CRC-1014_090516-CW-T-sort.bam >> NCC_CRC-1014_090516-CW-T_pileup.txt
#done < /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/utils/Cancer226-targets_hg19.bed


