#!/bin/bash

# def path and environment
source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
conda activate default

export outputfolder="/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/dilutions_series_chr22/"
export samplename_tumor="CRC-809_110914"
export samplename_healthy="pooledhealthy"
export chr=22
export dilutionfactor_tumor=0.75
export dilutionfactor_healthy=0.765

export outputdir=$outputfolder/dilutions_${samplename_tumor}/dilution_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}

export outputfile=dilution_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.sorted.bam

export cov=$(samtools depth -a $outputdir/$outputfile | awk '{sum+=$3} END {print sum/NR}')
echo $cov
#if [ ! -f $outputdir/coverage_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt ] ; then echo $cov >> $outputdir/coverage_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt ; fi
echo $cov > $outputdir/coverage_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt 

