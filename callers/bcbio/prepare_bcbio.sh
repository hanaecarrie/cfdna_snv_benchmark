#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input a configuration file indicating a mixture series to prepare the tree and files
# in order to run 5 callers included in the bcbio-nextgen pipeline in tumor-normal mode without realignment.

if [ $# == 0 ]; then
    echo "Usage: $0 -c [config_file]"
    echo "* config_file: string. full path to the configuration .yaml file."
    echo "Example:"
    echo "$ cd cfdna_snv_benchmark/callers/bcbio"
    echo "$ bash $0 -c config/config_bcbio_template.yml"
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
echo $extdatadir
echo $bcbiodir
echo $mode

cd $bcbiodir

for a in ${datadir}/mixtures/${regexpattern} ; do
  export b=$(basename $a) ;
  if [ ! -d $b ] ; then
    # create bcbio_work directory
    mkdir $b ;
    mkdir $b/bcbio_work ;
    # copy configuration file and changes names
    cp config.yaml $b/config.yaml ;
    export chr=$(echo $b | cut -d '_' -f 2) ; echo $chr ;
    export patient=$(echo $b | cut -d '-' -f 2 | cut -d '_' -f 1) ; echo $patient ;
    if [ mode == 'exome' ] ; then sed -i "s!exomebed!${bcbiodir}/$b/exome_${chr}.bed!g" $b/config.yaml ; fi
    sed -i "s/sampleid/${b}/g" $b/config.yaml ;
    export tumorbamfile=$a/$b.bam
    export normalbamfile=${datadir}/initialsamples/initialsamples_${chr}/initialsamples_${chr}_CRC-${patient}-CW-N/CRC-${patient}-CW-N_${chr}.bam
    echo $tumorbamfile ; echo $normalbamfile ;
    sed -i "s!tumorbamfile!${tumorbamfile}!g" $b/config.yaml ;
    sed  -i "s!normalbamfile!${normalbamfile}!g" $b/config.yaml ;
    # copy run script and change names
    cp run_sampleid.sh $b/bcbio_work/run_${b}.sh ;
    sed -i "s/sampleid/${b}/g" $b/bcbio_work/run_${b}.sh ;
    export t=$(echo $b | cut -d '_' -f 5) ; echo $t ;
    export n=$(echo $b | cut -d '_' -f 8) ; echo $n ;
    export jobid=${t}${n}_${chr}_${patient} ; echo $jobid ;
    sed -i "s/jobid/${jobid}/g" $b/bcbio_work/run_${b}.sh ;
    sed -i "s/bcbionextgenpath/${bcbionextgenpath}/g" $b/bcbio_work/run_${b}.sh ;
    sed -i "s/bcbiodir/${bcbiodir}/g" $b/bcbio_work/run_${b}.sh ;
    # if WGS change job type to long
    if [ mode == 'wholegenome' ] ; then
      sed -i "s/-t 3-00:00:00/-t 14-00:00:00/g" $b/bcbio_work/run_${b}.sh ;
      sed -i "s/-p normal/-p long/g" $b/bcbio_work/run_${b}.sh ;
    fi
    # copy exome bed
    if [ mode == 'exome' ] ; then cp ${extdatadir}/exome_bed/exome_${chr}.bed $b/exome_${chr}.bed ; fi
  else
    echo $b
fi
done
