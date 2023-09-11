#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input a configuration file of a mixture series to:
# step 1: create the a mixture bam file for a given chromosome sorted and indexed,
# step 2: submit jobs for individual mixture samples creation,
# step 3: create also the corresponding buffycoat bam file for the given chromosome.

if [ $# == 0 ]; then
    echo "Usage: $0 -c [config_file]"
    echo "* config_file: string. full path to the configuration .yaml file"
    echo "Example:"
    echo "$ cd cfdna_snv_benchmark/mixtures/"
    echo "$ bash $0 -c config/config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml"
    echo "Remark 1: In practice, for step 1 and 3, 24Gb RAM and > 8 CPUs recommended. Takes several hours to complete."
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

# activate conda environment
source $condapath
conda activate default

# print information from config file
echo $config_file
echo "path to samtools: "$samtools
echo "computer machine type: Aquila (qsub), Ronin (nohup) or Pluto (slurm): " $machine
echo "initial high tumor fraction (TF) cfDNA sample: " $sample_tumor
echo "initial low TF cfDNA sample: " $sample_healthy
echo "buffycoat sample as matched normal for germline subtraction: " $sample_buffycoat
echo "ID of high TF sample: " $samplename_tumor
echo "ID of low TF sample: " $samplename_healthy
echo "ID of buffycoat sample: " $samplename_buffycoat
echo "chromosome: " $chr
echo "main output folder: " $outputfolder
echo "file with estimated intital TF: " $tffile
echo "dilution factors list (x reads of highTFcfDNA _ x reads of lowTFcfDNA): " $dilutionfactors
echo "folder of high TF cfDNA sample: " $tumordir
echo "folder of low TF cfDNA sample: "$healthydir
echo "folder of buffycoat sample: "$buffycoatdir
# define derived useful variables
export sample_tumor_chr=$tumordir/${samplename_tumor}_chr${chr}.bam
export sample_healthy_chr=$healthydir/${samplename_healthy}_chr${chr}.bam
export sample_buffycoat_chr=$buffycoatdir/${samplename_buffycoat}_chr${chr}.bam
export tumor_chr_coverage=$tumordir/coverage_${samplename_tumor}_chr${chr}.txt
export healthy_chr_coverage=$healthydir/coverage_${samplename_healthy}_chr${chr}.txt
export buffycoat_chr_coverage=$buffycoatdir/coverage_${samplename_buffycoat}_chr${chr}.txt
# print created variables
echo "path to bam of the high TF cfDNA sample restricted to the studied chromosome: " $sample_tumor_chr
echo "path to bam of the low TF cfDNA sample restricted to the studied chromosome: " $sample_healthy_chr
echo "path to coverage file of the high TF cfDNA sample restricted to the studied chromosome: " $tumor_chr_coverage
echo "path to coverage file of the low TF cfDNA sample restricted to the studied chromosome: " $healthy_chr_coverage
echo "path to coverage file of the buffycoat sample restricted to the studied chromosome: " $buffycoat_chr_coverage
# create necessary folders if needed
if [ ! -d $outputfolder/mixtures_chr${chr}_${samplename_tumor}_${samplename_healthy} ] ; then mkdir -p $outputfolder/mixtures_chr${chr}_${samplename_tumor}_${samplename_healthy} ; fi
if [ ! -d $tumordir ] ; then mkdir -p $tumordir ; fi
if [ ! -d $healthydir ] ; then mkdir -p $healthydir ; fi
if [ ! -d $buffycoatdir ] ; then mkdir -p $buffycoatdir ; fi
# copy configuration file to output directory for reference
cp $config_file $outputfolder/mixtures_chr${chr}_${samplename_tumor}_${samplename_healthy}/ 

#######################################################################################################################
echo "Step 1: Select chr ${chr} only for the high TF and the low TF cfDNA samples"
#######################################################################################################################

if [ ! -d $outputdir ] ; then mkdir -p $outputdir ; fi
# for high TF cfDNA sample
if [ ! -d $tumordir ] ; then mkdir $tumordir ; fi
  if [ "$chr" == 'all' ] ; then
    # if no chr selection is needed, simply copy the whole bam file
    if [ ! -f $sample_tumor_chr ] ; then cp $sample_tumor $sample_tumor_chr ; fi
  else
    # if chromosome selection is needed, use samtools view
	  if [ ! -f $sample_tumor_chr ] ; then $samtools view -b $sample_tumor $chr > $sample_tumor_chr ; fi
fi
# index corresponding bam file
if [ ! -f ${sample_tumor_chr}.bai ] ; then $samtools index $sample_tumor_chr ; fi
# calculate coverage
if  [ ! -f $tumor_chr_coverage ] ; then export tumor_cov=$($samtools depth -a $sample_tumor_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $tumor_cov >  $tumor_chr_coverage ; else export tumor_cov=$(cat $tumor_chr_coverage) ; fi
# for low TF cfDNA sample
if [ ! -d $healthydir ] ; then mkdir $healthydir ; fi
if [ "$chr" == 'all' ] ; then
  # if no chromosome selection is needed, simply copy the whole bam file
	if [ ! -f $sample_healthy_chr ] ; then cp $sample_healthy $sample_healthy_chr ; fi
else
  # if chromosome selection is needed, use samtools view
	if [ ! -f $sample_healthy_chr ] ; then $samtools view -b $sample_healthy $chr > $sample_healthy_chr ; fi
fi
# index corresponding bam file
if [ ! -f ${sample_healthy_chr}.bai ] ; then $samtools index $sample_healthy_chr ; fi
# calculate coverage
if  [ ! -f $healthy_chr_coverage ] ; then export healthy_cov=$($samtools depth -a $sample_healthy_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $healthy_cov >  $healthy_chr_coverage ; else export healthy_cov=$(cat $healthy_chr_coverage) ; fi
# print obtained effective coverage
echo "effective coverage of high TF cfDNA sample: " $tumor_cov
echo "effective coverage of low TF cfDNA sample: " $healthy_cov

#######################################################################################################################
echo "Step 2: submit jobs for mixture samples creation"
#######################################################################################################################

for dilutionfactor in $dilutionfactors ; do
  echo "Mixture with dilution factor: " $dilutionfactor
  export dilutionfactor_tumor=$(echo $dilutionfactor | cut -f1 -d-)
  export dilutionfactor_healthy=$(echo $dilutionfactor | cut -f2 -d-)
  export outputdir=$outputfolder/mixtures_chr${chr}_${samplename_tumor}_${samplename_healthy}/mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x
  if [ ! -d $outputdir/logs ] ; then mkdir -p  $outputdir/logs ; fi
  if [ $machine == 'Aquila' ] ; then # on HPC with qsub job submission
    /opt/uge-8.1.7p3/bin/lx-amd64/qsub -pe OpenMP 4 -l mem_free=24G,h_rt=24:00:00 -o $outputdir/logs/ -e $outputdir/logs/ $repopath/cfdna_snv_benchmark/mixtures/create_mixtures_chr.sh -c $config_file -d $dilutionfactor ;
  elif [ $machine == 'Ronin' ] ; then # on AWS virtual machine, nohup job submission
    nohup $repopath/cfdna_snv_benchmark/mixtures/create_mixtures_chr.sh -c $config_file -d $dilutionfactor > $outputdir/logs/log_mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.out &
  elif [ $machine == 'Pluto' ] ; then # on HPC with slurm job submission
    sbatch -p normal -J chr${chr}_mix -t 24:00:00 -N 1 --mem 24G --output $outputdir/logs/log_mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.o --error $outputdir/logs/log_mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.e $repopath/cfdna_snv_benchmark/mixtures/create_mixtures_chr.sh -c $config_file -d $dilutionfactor
  else
    echo "Computer machine type handled are 'Aquila' (qsub), 'Ronin' (nohup) and 'Pluto' (slurm)."
    echo "For other job submission format, a new category should be added."
  fi
done


#######################################################################################################################
echo "Step 3: Select chr ${chr} only for the buffycoat sample if it does not exists yet"
#######################################################################################################################

if [ ! -d $buffycoatdir ] ; then mkdir $buffycoatdir ; fi
if [[ $chr == 'all' ]] ; then
	if [ ! -f $sample_buffycoat_chr ] ; then cp $sample_buffycoat $sample_buffycoat_chr ; fi
else
	if [ ! -f $sample_buffycoat_chr ] ; then $samtools view -b $sample_buffycoat $chr > $sample_buffycoat_chr ; fi
fi
if [ ! -f ${sample_buffycoat_chr}.bai ] ; then  $samtools index $sample_buffycoat_chr ; fi
if [[ $bedfile == 'wgs' ]] ; then #TODO
	if  [ ! -f $buffycoat_chr_coverage ] ; then export buffycoat_cov=$($samtools depth -a $sample_buffycoat_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $buffycoat_cov >  $buffycoat_chr_coverage ; else export buffycoat_cov=$(cat $buffycoat_chr_coverage) ; fi
else
	if  [ ! -f $buffycoat_chr_coverage ] ; then export buffycoat_cov=$($samtools depth -a -b $bedfile $sample_buffycoat_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $buffycoat_cov >  $buffycoat_chr_coverage ; else export buffycoat_cov=$(cat $buffycoat_chr_coverage) ; fi
fi

