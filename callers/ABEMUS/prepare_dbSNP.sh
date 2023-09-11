#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input XXX

if [ $# == 0 ]; then
    echo "Usage: $0 extdata repopath"
    echo "* extdata: XXX"
    echo "* repopath: XXX"
    echo "Example:"
    echo "$ bash $0 XXX XXX"
    exit 1
fi

# Parse input to get initial configuration file
export extdata=$1
export repopath=$2
echo "Extdata location: " $1
echo "Repositiory location: " $2

# 1. Download it in /data/extdata folder
cd ${extdata}/dbsnp_vcf
if [ ! -f GRCh37_latest_dbSNP_all.vcf.gz ] ; then
  wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_dbSNP_all.vcf.gz
fi
tabix GRCh37_latest_dbSNP_all.vcf.gz
# to get correspondance refs and chromosome: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/
#TODO recover convert_chrom_names.txt
bcftools annotate --rename-chrs ${extdata}/GRCh37/convert_chrom_names.txt -o dbSNP.vcf GRCh37_latest_dbSNP_all.vcf.gz
bcftools view dbSNP.vcf -Oz -o dbSNP.vcf.gz
bcftools index dbSNP.vcf.gz

for chr in {1..22} ; do

  echo "Process chromosome ${chr}"

  # 2. Split per chromosome
  if [ ! -f ${extdata}/dbsnp_vcf/dbSNP_hg19_chr${chr}.vcf ] ; then
    bcftools filter -r ${chr} ${extdata}/dbsnp_vcf/dbSNP.vcf.gz -o ${extdata}/dbsnp_vcf/dbSNP_hg19_chr${chr}.vcf ;
  fi

  # 3. Edit vcf to keep only SNPs with one Alt base
  echo ${extdata}/dbsnp_vcf/dbSNP_hg19_chr${chr}_edited.vcf
  if [ ! -f ${extdata}/dbsnp_vcf/dbSNP_hg19_chr${chr}_edited.vcf ] ; then
    python ${repopath}/callers/ABEMUS/edit_vcf.py --vcfpath $extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}.vcf ;
  fi

done

# 4. concat all vcf edited for WES

# need to add header before concatenate
for chr in {1..22} ; do echo $chr ; cp dbSNP_hg19_chr${chr}_edited.vcf dbSNP_hg19_chr${chr}_edited_header.vcf ; done
# TODO automatize display and copy manually header chrom1
cat dbSNP_hg19_chr1.vcf | head -350
# TODO manually add to _edited_header.vcf files
# then concat
bcftools concat -Ov -o dbSNP_hg19_chrall_edited.vcf dbSNP_hg19_chr{1..22}_edited_header.vcf