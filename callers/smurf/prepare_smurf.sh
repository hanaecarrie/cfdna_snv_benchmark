#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input a configuration file indicating a mixture series to prepare the tree and files
# in order to run SMuRF.

if [ $# == 0 ]; then
    echo "Usage: $0 -c [config_file]"
    echo "* config_file: string. full path to the configuration .yaml file."
    echo "Example:"
    echo "$ cd cfdna_snv_benchmark/callers/smurf"
    echo "$ bash $0 -c config/config_smurf_template.yml"
    exit 1
fi

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


echo $datadir
echo $bcbiodir

cd $bcbiodir

for a in ${datadir}/mixtures/${regexpattern} ; do
        export b=$(basename $a) ;
	if [ ! -f $b/bcbio_final/smurf/time.txt ] ; then
                mkdir $b/bcbio_work ;
		mkdir $b/bcbio_work/smurf_workdir ;
                cp run_smurf_sampleid.sh $b/bcbio_work/smurf_workdir/run_smurf_${b}.sh ;
                sed -i "s/XXX/${b}/g" $b/bcbio_work/smurf_workdir/run_smurf_${b}.sh ;
                sed -i "s/bcbiodir/${bcbiodir}/g" $b/bcbio_work/smurf_workdir/run_smurf_${b}.sh ;
                sed -i "s/smurfdir/${smurfdir}/g" $b/bcbio_work/smurf_workdir/run_smurf_${b}.sh ;
                sed -i "s/aparam/${aparam}/g" $b/bcbio_work/smurf_workdir/run_smurf_${b}.sh ;
                sed -i "s/sparam/${sparam}/g" $b/bcbio_work/smurf_workdir/run_smurf_${b}.sh ;
        else
                echo $b
  fi
done

