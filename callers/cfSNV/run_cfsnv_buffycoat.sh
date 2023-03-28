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

while getopts ":c:p:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
    p) plasma=${OPTARG} ;;
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

echo $tmpdir

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


### Convert BAM to FASTQ ###
export normalfastqoutdir=$(dirname $normal)
echo $normalfastqoutdir
if [ ! -f $(dirname $normal)/$(basename $normal .bam)_R1.fastq.gz ] ; then python ${repopath}/cfSNV/generate_fastqs_yaml.py -i $normal -t $(dirname $normal)/tmp -o $normalfastqoutdir -c $config_file ; fi

### Run cfSNV pipeline ###

# get bam for normal
Rscript run_getbamalign_buffycoat.R --config_file $config_file 


