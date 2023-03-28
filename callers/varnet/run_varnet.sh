#!/bin/bash

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


source ${condapath}
echo $PATH
conda activate varnet

echo $normalbam
echo $dilutionseriesfolder
echo $includeallelefrequency
echo $outdir
echo $py
echo $chr

cd $repopath

####### generate info file #######
if [ ! -d $outdir ] ; then mkdir $outdir ; fi
cp $config_file $outdir/

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>$outdir/log.out 2>&1
# Everything below will go to the file 'log.out'

# Start logging the RAM usage and CPU usage
bash ${repopath}/log_mem_cpu.sh $outdir/log_mem_cpu.out  & export logpid=$!

# mark and remove duplicates in normal bam
if [ ! -f $(dirname $normalbam)/$(basename $normalbam .bam)_rmd.bam ] ; then  
	$java -jar $picard MarkDuplicates \
		-I $normalbam \
		-O $(dirname $normalbam)/$(basename $normalbam .bam)_rmd.bam \
		-M $(dirname $normalbam)/$(basename $normalbam .bam)_rmd_metrics.txt \
		--REMOVE_DUPLICATES true
fi
# index rmd normal bam
if [ ! -f $(dirname $normalbam)/$(basename $normalbam .bam)_rmd.bam.bai ] ; then
	$samtools index $(dirname $normalbam)/$(basename $normalbam .bam)_rmd.bam
fi

for tumorbam in ${dilutionseriesfolder}/*/*[Tx].bam ; do 
	${repopath}/run_varnet_sample.sh -c $config_file -t $tumorbam 
done



 
