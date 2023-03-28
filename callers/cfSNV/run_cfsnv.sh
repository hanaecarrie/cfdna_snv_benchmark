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

#parse config file
eval $(parse_yaml $config_file)

source $condapath
conda activate cfSNV

cd ${repopath}/cfSNV


echo $config_file
echo $outdir
echo $tmpdir
echo $dilutionseriesfolder

if [ ! -d $outdir ]; then mkdir $outdir ; fi

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>$outdir/log.out 2>&1
# Everything below will go to the file 'log.out':

# Start logging the RAM usage and CPU usage
bash ${repopath}/ABEMUS/log_mem_cpu.sh $outdir/log_mem_cpu.out  & export logpid=$!

#echo 'Clean tmp folder'
#echo $tmpdir
#rm -r $tmpdir

touch $outdir/logtime.out


###### prepare bedfile ######
# Download reference genome (here all with hg19)
# Index it
# Create dictionary
# Create bed file chr
# Split bed file chr by chunks of 5,000 lines
if [ ! -f $extdata/exome_bed/exome_hg19_chr${chr}_00.bed ] ; then split -l 5000 --numeric-suffixes --additional-suffix='.bed' $extdata/exome_bed/exome_hg19_chr${chr}.bed $extdata/exome_bed/exome_hg19_chr${chr}_ ; fi


echo "dilfolder ${dilutionseriesfolder}"

export normal=$buffycoatbam
if [ ! -d $outdir/$(basename $normal .bam) ] ; then mkdir $outdir/$(basename $normal .bam) ; fi

export normalfastqoutdir=$(dirname $normal)
echo $normalfastqoutdir
if [ ! -f $(dirname $normal)/$(basename $normal .bam)_R1.fastq.gz ] ; then python ${repopath}/cfSNV/generate_fastqs_yaml.py -i $normal -t $(dirname $normal)/tmp -o $normalfastqoutdir -c $config_file ; fi

# get bam for plasma and normal as well as notcombined and extendedfrags bams for plasma
# run parameter recommend function
export normalid=$(basename $normal .bam)
echo $normalid
Rscript run_getbamalign_buffycoat.R --config_file $config_file


export npid=0
for plasma in ${dilutionseriesfolder}/*/*[Tx].bam ; do
        echo "plasma ${plasma}" ;
	export npid=$((npid+1))
	echo "n PID ${npid}" 
        bash ${repopath}/cfSNV/run_cfsnv_sample.sh -c $config_file -p $plasma  &  pids[${npid}]=$!
done

# wait for all pids
for pid in ${pids[*]}; do
        wait $pid
done


kill -9 $logpid

