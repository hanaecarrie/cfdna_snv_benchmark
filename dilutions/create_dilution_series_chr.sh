#!/bin/bash

# def path and environment
source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
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

echo $config_file
echo $sample_healthy
echo $sample_tumor
echo $sample_buffycoat 

echo $samplename_healthy
echo $samplename_tumor
echo $samplename_buffycoat

echo $dilutionfactors
echo $chr
echo $outputfolder
echo $tffile

if [ ! -d $outputfolder ] ; then mkdir $outputfolder ; fi

for dilutionfactor in $dilutionfactors ;

do qsub -pe OpenMP 4 -l mem_free=48G,h_rt=24:00:00 -o /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/logs/ -e /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/logs/ /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/dilutions/create_dilution_chr.sh -c $config_file -d $dilutionfactor

done

