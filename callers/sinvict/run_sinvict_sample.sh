#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input a configuration file indicating a mixture series and a specific plasma bam file of this series
# to apply SiNVICT caller on this sample.

if [ $# == 0 ]; then
    echo "Usage: $0 -c [config_file] -p [plasma]"
    echo "* config_file: string. full path to the configuration .yaml file."
    echo "* plasma: string. full path to the plasma bam file of the mixture series."
    echo "Example:"
    echo "$ cd cfdna_snv_benchmark/callers/sinvict"
    echo "$ bash $0 -c config_sinvict_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml -p XXX"
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

while getopts ":c:p:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
    p) plasma=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

source $condapath
conda activate default

cd ${repopath}/sinvict

echo $config_file
echo $dilutionseriesfolder
echo $buffycoatbam
echo $chr
echo $extdata
echo $outdir

echo $plasma

export outdirplasma=$outdir/$(basename $plasma .bam)
echo $outdirplasma
if [ ! -d $outdirplasma ] ; then mkdir $outdirplasma ; fi

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>$outdirplasma/log.out 2>&1
# Everything below will go to the file 'log.out'

echo 'Start sinvict per sample'
echo $plasma

export nchunk=$(ls $extdata/${mode}_bed/${mode}_hg19_chr${chr}_*.bed | wc -l)
echo $nchunk

export outdirplasma=$outdir/$(basename $plasma .bam)
echo $outdirplasma
if [ ! -d $outdirplasma ] ; then mkdir $outdirplasma ; fi

for n in $(seq -f "%02g" 0 $(($nchunk - 1))) ; do
	echo "run sinvict on chunk ${n}"
        export npid=$((${n#0} + 1))
	bash ${repopath}/sinvict/run_sinvict_i.sh -c $config_file -i $n -p $plasma &  pids[${npid}]=$!
	sleep 2
done
# wait for all pids
for pid in ${pids[*]}; do
	wait $pid
done

#### SINVICT ###
if [ ! -d $outdirplasma/results ] ; then mkdir $outdirplasma/results ; fi
${repopath}/sinvict/sinvict --min-depth $mindepth -t ${outdirplasma}/bam-readcount -o ${outdirplasma}/results
