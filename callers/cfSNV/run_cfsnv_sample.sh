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

while getopts ":c:p:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
    p) plasma=${OPTARG} ;;
  esac
done

#parse config file
eval $(parse_yaml $config_file)

echo $config_file
echo $outdir
echo $tmpdir
echo $dilutionseriesfolder

if [ ! -d $outdir ]; then mkdir $outdir ; fi

export outdirplasma=$outdir/$(basename $plasma .bam)
echo $outdirplasma
if [ ! -d $outdirplasma ] ; then mkdir $outdirplasma ; fi

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>$outdirplasma/log.out 2>&1
# Everything below will go to the file 'log.out':

if [ ! -d $outdir ]; then mkdir $outdir ; fi

echo $tmpdir

echo "dilfolder ${dilutionseriesfolder}"

export normal=$buffycoatbam
if [ ! -d $outdir/$(basename $normal .bam) ] ; then mkdir $outdir/$(basename $normal .bam) ; fi

echo "plasma ${plasma}"

### Convert BAM to FASTQ ###
export plasmafastqoutdir=$(dirname $plasma)
echo $plasmafastqoutdir
if [ ! -f $(dirname $plasma)/$(basename $plasma .bam)_R1.fastq.gz ] ; then python /home/ubuntu/cfSNV/generate_fastqs_yaml.py  -i $plasma -t $(dirname $plasma)/tmp -o $plasmafastqoutdir -c $config_file ; fi
export normalfastqoutdir=$(dirname $normal)
echo $normalfastqoutdir
if [ ! -f $(dirname $normal)/$(basename $normal .bam)_R1.fastq.gz ] ; then python /home/ubuntu/cfSNV/generate_fastqs_yaml.py -i $normal -t $(dirname $normal)/tmp -o $normalfastqoutdir -c $config_file ; fi

### Run cfSNV pipeline ###

startcfsnv=$(date +%s)

# get bam for plasma and normal as well as notcombined and extendedfrags bams for plasma
# run parameter recommend function
export plasmaid=$(basename $plasma .sorted.bam)
echo $plasmaid
Rscript run_getbamalign.R --config_file $config_file --plasmaid $plasmaid

if [ ! -f $outdirplasma/parameter.txt ] ; then Rscript run_parameter.R --config_file $config_file --plasmaid $plasmaid ;
tail -n 13 $outdirplasma/log.out > $outdirplasma/parameter.txt ;
fi

# apply cfSNV mutation calling function per batch
if [ ! -d $outdirplasma/results ] ; then mkdir $outdirplasma/results ; fi
for targetbed in $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}_*.bed ; do 
	echo $targetbed 
	Rscript run_variantcalling.R --config_file $config_file --plasmaid $plasmaid --targetbed $targetbed
done

####### copy results to common folder ########
echo $finaloutdir
if [ ! -d $finaloutdir ] ; then mkdir $finaloutdir ; fi
echo "plasma ${plasma}" ;
if [ ! -d $finaloutdir/$(basename $plasma .bam) ] ; then mkdir $finaloutdir/$(basename $plasma .bam) ; fi
if [ ! -d $finaloutdir/$(basename $plasma .bam)/cfsnv ] ; then mkdir $finaloutdir/$(basename $plasma .bam)/cfsnv ; fi
scp $outdirplasma/results/* $finaloutdir/$(basename $plasma .bam)/cfsnv/

endcfsnv=$(date +%s)
timecfsnv=$(($endcfsnv-$startcfsnv))
hmscfsnv=$(printf '%02dh:%02dm:%02ds\n' $((timecfsnv/3600)) $((timecfsnv%3600/60)) $((timecfsnv%60)))
        echo "Elapsed Time cfSNV on ${sampleid} for chr${chr}: ${hmscfsnv}"
        echo "Elapsed Time cfSNV on ${sampleid} for chr${chr}: ${hmscfsnv}" >> $outdir/logtime.out

