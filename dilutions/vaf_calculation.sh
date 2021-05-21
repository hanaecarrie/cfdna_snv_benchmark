#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default

export outputdir=/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/results_vcf

for input_vcf in $outputdir/*vcf;
do vcftools --vcf $input_vcf --freq --out $outputdir/$(basename $input_vcf .vcf)_outputvaf.vcf ;
done

