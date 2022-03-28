#!/bin/bash

source /home/ubuntu/anaconda3/etc/profile.d/conda.sh
conda activate cfSNV

cd /home/ubuntu/cfSNV

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
bash /home/ubuntu/ABEMUS/log_mem_cpu.sh $outdir/log_mem_cpu.out  & export logpid=$!

echo 'Clean tmp folder'
echo $tmpdir
rm -r $tmpdir

touch $outdir/logtime.out

echo "dilfolder ${dilutionseriesfolder}"

export normal=$buffycoatbam
if [ ! -d $outdir/$(basename $normal .bam) ] ; then mkdir $outdir/$(basename $normal .bam) ; fi

export normalfastqoutdir=$(dirname $normal)
echo $normalfastqoutdir
if [ ! -f $(dirname $normal)/$(basename $normal .bam)_R1.fastq.gz ] ; then python /home/ubuntu/cfSNV/generate_fastqs_yaml.py -i $normal -t $(dirname $normal)/tmp -o $normalfastqoutdir -c $config_file ; fi

export npid=0
for plasma in ${dilutionseriesfolder}/*/*[Tx].bam ; do
        echo "plasma ${plasma}" ;
	export npid=$((npid+1))
	echo "n PID ${npid}" 
        bash /home/ubuntu/cfSNV/run_cfsnv_sample.sh -c $config_file -p $plasma  &  pids[${npid}]=$!
done

# wait for all pids
for pid in ${pids[*]}; do
        wait $pid
done


kill -9 $logpid

