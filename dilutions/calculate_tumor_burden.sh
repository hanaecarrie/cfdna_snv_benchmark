#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default


# function to parse config file
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}


while getopts ":c:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

#echo $sample_tumor
#echo $sample_healthy
#echo $sample_tumor
#echo $samplename_tumor
#echo $samplename_healthy
#echo $dilutionfactor_tumor
#echo $dilutionfactor_healthy
#echo $outputdir
export outputdir=$outputdir/dilutions_${samplename_tumor}/dilution_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}
#echo $outputdir


while read line ; do export A="$(cut -d',' -f1 <<<"$line")" ;
if [ $A == $samplename_tumor ] ; then echo $A ; export median_tumor_burden="$(cut -d',' -f3 <<<"$line")" ; export cov_tumor="$(cut -d',' -f2 <<<"$line")" ; fi ; done < /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/ctDNA_burden_estimation_deepwgs_plasma_CRC_samples.csv ;
echo $median_tumor_burden
echo $cov_tumor
export cov_healthy=187.6145
echo $cov_healthy

cov_t=$(echo "$median_tumor_burden * $cov_tumor * $dilutionfactor_tumor" | bc)
echo $cov_t
cov_tot=$(echo "$cov_healthy * $dilutionfactor_healthy + $cov_tumor * $dilutionfactor_tumor" | bc)
echo $cov_tot
mixed_sample_tf=$(echo "$cov_t / $cov_tot" | bc -l)
echo $mixed_sample_tf
if [ ! -f $outputdir/estimated_tf.txt ] ; then echo $mixed_sample_tf >> $outputdir/estimated_tf.txt ; fi  

