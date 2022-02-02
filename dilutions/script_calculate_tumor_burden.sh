#!/bin/bash

# def path and environment
source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
conda activate default

export outputfolder="/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/dilutions_series_chr22/"
export samplename_tumor="CRC-809_030915"
export samplename_healthy="pooledhealthy"
export chr=22
export dilutionfactor_tumor=1
export dilutionfactor_healthy=0
export tffile='/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/dilutions/ctDNA_burden_estimation_deepwgs_plasma_CRC_samples.csv'
export covhealthy=400.814


export outputdir=$outputfolder/dilutions_${samplename_tumor}/dilution_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}
echo $tffile

while read line ; do export A="$(cut -d',' -f1 <<<"$line")" ;
if [ $A == $samplename_tumor ] ; then echo $A ; export median_tumor_burden="$(cut -d',' -f3 <<<"$line")" ; export cov_tumor="$(cut -d',' -f2 <<<"$line")" ; fi ; done < $tffile
echo $median_tumor_burden
echo $cov_tumor
export cov_healthy=$covhealthy
echo $cov_healthy

cov_t=$(echo "$median_tumor_burden * $cov_tumor * $dilutionfactor_tumor" | bc)
echo $cov_t
cov_tot=$(echo "$cov_healthy * $dilutionfactor_healthy + $cov_tumor * $dilutionfactor_tumor" | bc)
echo $cov_tot
mixed_sample_tf=$(echo "$cov_t / $cov_tot" | bc -l)
echo $mixed_sample_tf
#if [ ! -f $outputdir/estimated_tf_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt ] ; then echo $mixed_sample_tf >> $outputdir/estimated_tf_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt ; fi  
echo $mixed_sample_tf > $outputdir/estimated_tf_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt
