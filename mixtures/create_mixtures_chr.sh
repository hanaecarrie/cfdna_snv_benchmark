#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input a configuration file of a mixture series and a dilutionfactor to:
# step 1 a-e: create the a mixture bam file for a given chromosome sorted and indexed,
# step 2: create also the corresponding buffycoat bam file for the given chromosome
# step 3 a-c: create companion information i.e. coverage, tumor fraction (TF) rough estimate and ichorCNA TF estimate.

if [ $# == 0 ]; then
    echo "Usage: $0 -c [config_file] -d [dilutionfactor]"
    echo "* config_file: string. full path to the configuration .yaml file"
    echo "* dilutionfactor: string 'D1-D2'. with D1x depth of coverage to extract from high TF cfDNA sample"
    echo "                                   and D2x depth of coverage to extract from low TF cfDNA sample."
    echo "Example:"
    echo "$ cd cfdna_snv_benchmark/mixtures/"
    echo "$ bash $0 -c config/config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml -d '70-80'"
    echo "Remark 1: In practice, on 100x WGS, 24Gb RAM and > 8 CPUs recommended. Takes several hours to complete."
    echo "Remark 2: The obtained mixture would be about 70+80=150x cov."
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

# get inputs
while getopts ":c:d:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
    d) dilutionfactor=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

# activate conda environment
source $condapath
conda activate default

# print information from config file
echo "initial high tumor fraction (TF) cfDNA sample: " $sample_tumor
echo "initial low TF cfDNA sample: " $sample_healthy
echo "buffycoat sample as matched normal for germline subtraction: " $sample_buffycoat
echo "ID of high TF sample: " $samplename_tumor
echo "ID of low TF sample: " $samplename_healthy
echo "ID of buffycoat sample: " $samplename_buffycoat
echo "chromosome: " $chr
echo "main output folder: " $outputfolder
echo "file with estimated intital TF: " $tffile
echo "dilution factor (x reads of highTFcfDNA _ x reads of lowTFcfDNA): " $dilutionfactor
echo "folder of high TF cfDNA sample: " $tumordir
echo "folder of low TF cfDNA sample: "$healthydir
echo "folder of buffycoat sample: "$buffycoatdir
# define derived useful variables
export dilutionfactor_tumor=$(echo $dilutionfactor | cut -f1 -d-) # get nb reads to extract from high TF cfDNA sample
export dilutionfactor_healthy=$(echo $dilutionfactor | cut -f2 -d-) # get nb reads to extract from low TF cfDNA sample
export outputdir=$outputfolder/mixtures_chr${chr}_${samplename_tumor}_${samplename_healthy}/mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x
export sample_tumor_chr=$tumordir/${samplename_tumor}_chr${chr}.bam
export sample_healthy_chr=$healthydir/${samplename_healthy}_chr${chr}.bam
export sample_buffycoat_chr=$buffycoatdir/${samplename_buffycoat}_chr${chr}.bam
export tumor_chr_coverage=$tumordir/coverage_${samplename_tumor}_chr${chr}.txt
export healthy_chr_coverage=$healthydir/coverage_${samplename_healthy}_chr${chr}.txt
export buffycoat_chr_coverage=$buffycoatdir/coverage_${samplename_buffycoat}_chr${chr}.txt
# print created variables
echo $dilutionfactor_tumor "x depth of coverage to extract from high TF cfDNA sample"
echo $dilutionfactor_healthy "x depth of coverage to extract from low TF cfDNA sample"
echo "specific output folder: " $outputdir
echo "path to bam of the high TF cfDNA sample restricted to the studied chromosome: " $sample_tumor_chr
echo "path to bam of the low TF cfDNA sample restricted to the studied chromosome: " $sample_healthy_chr
echo "path to coverage file of the high TF cfDNA sample restricted to the studied chromosome: " $tumor_chr_coverage
echo "path to coverage file of the low TF cfDNA sample restricted to the studied chromosome: " $healthy_chr_coverage
echo "path to coverage file of the buffycoat sample restricted to the studied chromosome: " $buffycoat_chr_coverage

echo "START MIXTURE CREATION"

#######################################################################################################################
echo "Step 1a: Select chr ${chr} only for the high TF and the low TF cfDNA samples"
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
echo "Step 1b: Subsample high and low TF cfDNA samples to obtain the desired coverage"
#######################################################################################################################

# calculate the fraction of coverage used for each sample
export dilutionfraction_tumor=$(echo "scale=4; $dilutionfactor_tumor / $tumor_cov" | bc -l)
# if fraction >= 1 then set it to 1
if [[ $(echo "$dilutionfraction_tumor>1" | bc) == 1 ]]; then dilutionfraction_tumor=1 ; fi
export dilutionfraction_healthy=$(echo "scale=4; $dilutionfactor_healthy / $healthy_cov" | bc -l)
# if fraction >= 1 then set it to 1
if [[ $(echo "$dilutionfraction_healthy>1" | bc) == 1 ]]; then dilutionfraction_healthy=1 ; fi
export sample_tumor_chr_downsample=$tumordir/${samplename_tumor}_chr${chr}_${dilutionfactor_tumor}x.bam
export sample_healthy_chr_downsample=$healthydir/${samplename_healthy}_chr${chr}_${dilutionfactor_healthy}x.bam
export dilutionname=mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x
export rgname=rg_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.txt
if [[ $dilutionfraction_tumor != 1 ]] ; then
  # fraction < 1 then use sametools view -s seed.fraction to subsample. Here seed is always set to 0.
  if [ ! -f $sample_tumor_chr_downsample ] ; then $samtools view -b -s $dilutionfraction_tumor -o $sample_tumor_chr_downsample $sample_tumor_chr ; fi
else
  # fraction = 1 then just copy whole file without subsampling.
  if [ ! -f $sample_tumor_chr_downsample ] ; then cp $sample_tumor_chr $sample_tumor_chr_downsample ; fi
fi
if [[ $dilutionfraction_healthy == 0 ]] ; then # if case of D2 = 0 i.e. no mixing with low TF cfDNA sample, simply copy
  if [ ! -f $outputdir/${dilutionname}.bam ] ; then  cp  $sample_tumor_chr_downsample $outputdir/${dilutionname}.bam ; fi
else # D2 != 0 mixing of high and low TF cfDNA samples
  if [[ $dilutionfraction_healthy != 1 ]] ; then
    # fraction < 1 then use sametools view -s seed.fraction to subsample. Here seed is always set to 0.
    if [ ! -f $sample_healthy_chr_downsample ] ; then $samtools view -b -s $dilutionfraction_healthy -o $sample_healthy_chr_downsample $sample_healthy_chr ; fi
  else
    # fraction = 1 then just copy whole file without subsampling.
    if [ ! -f $sample_healthy_chr_downsample ] ; then cp $sample_healthy_chr $sample_healthy_chr_downsample ; fi
  fi
  # label reads coming from the high and low TF cfDNA samples
  if [ ! -f $outputdir/$rgname ] ; then perl -e 'print "@RG\\tID:${samplename_tumor}\\tSM:hs\\tLB:${samplename_tumor}\\tPL:Illumina\\n@RG\\tID:${samplename_healthy}\\tSM:hs\\tLB:${samplename_healthy}\\tPL:Illumina\\n"' > $outputdir/$rgname ; fi
fi
# print obtained information and paths
echo 'coverage fraction used from high TF cfDNA sample: ' $dilutionfraction_tumor
echo 'coverage fraction low TF cfDNA sample:' $dilutionfraction_healthy
echo "path to downsampled chrom high TF cfDNA sample: " $sample_tumor_chr_downsample
echo "path to downsampled chrom low TF cfDNA sample: " $sample_healthy_chr_downsample
echo "bamfile name without extension: " $dilutionname
echo "path to read group name file: "$rgname

#######################################################################################################################
echo "Step 1c: Merge both subsampled cfDNA samples into a single mixture sample"
#######################################################################################################################

if [[ $dilutionfraction_healthy != 0 ]] ; then # D2 != 0 mixing of both samples, need to merge
if [ ! -f $outputdir/${dilutionname}.bam ] ; then $samtools merge -rh $outputdir/$rgname $outputdir/${dilutionname}.unsorted.bam $sample_tumor_chr_downsample $sample_healthy_chr_downsample ; fi
fi

#######################################################################################################################
echo "Step 1d: Sort mixture"
#######################################################################################################################

# using 4 processes and 2Gb of memory
if [ ! -f $outputdir/${dilutionname}.bam ] ; then $samtools sort -o $outputdir/${dilutionname}.bam -@ 4 -m 2G $outputdir/${dilutionname}.unsorted.bam ; fi

#######################################################################################################################
echo "Step 1e: Index mixture"
#######################################################################################################################

# using 4 processes
if [ ! -f $outputdir/${dilutionname}.bam.bai ] ; then $samtools index -@ 4 $outputdir/${dilutionname}.bam ; fi
# clean unsorted bam for space
if [  -f $outputdir/${dilutionname}.unsorted.bam ] ; then rm $outputdir/${dilutionname}.unsorted.bam ; fi

#######################################################################################################################
echo "Step 2: Select chr ${chr} only for the buffycoat sample if it does not exists yet"
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

#######################################################################################################################
echo "Step 3a: Estimate tumor burden by rough calculation"
#######################################################################################################################

# the tffile should contain tissue purity estimates of the initial high TF cfDNA sample
# here the low TF cfDNA sample is supposed to be cancer free
if [ ! -f $outputdir/estimated_tf_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.txt ] ; then
while read line ; do export A="$(cut -d',' -f1 <<<"$line")" ;
if [[ "$A" == *${samplename_tumor}* ]] ; then echo $A ; export median_tumor_burden="$(cut -d ',' -f3 <<<"$line")" ;  fi ; done < $tffile
# rough estimate of nb tumor reads per position
cov_t=$(echo "$median_tumor_burden * $tumor_cov * $dilutionfraction_tumor" | bc)
# rough estimate of overall coverage per position
cov_tot=$(echo "$healthy_cov * $dilutionfraction_healthy + $tumor_cov * $dilutionfraction_tumor" | bc)
mixed_sample_tf=$(echo "$cov_t / $cov_tot" | bc -l) # rough TF estimate
echo $mixed_sample_tf > $outputdir/estimated_tf_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.txt
echo "file paths with tissue purity estimates of the initial high TF cfDNA sample: " $tffile
echo "median of purity estimate methods on the initial high TF cfDNA sample: "$median_tumor_burden
echo "about cov_t=" $cov_t "x tumor reads"
echo "about cov_tot=" $cov_tot "x total reads"
echo "rough estimation of mixture TF: cov_t/cov_tot=" $mixed_sample_tf
fi

#######################################################################################################################
echo "Step 3b: Calculate mixture coverage"
#######################################################################################################################

if [ ! -f $outputdir/coverage_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.txt ] ; then
if [[ $bedfile == 'wgs' ]] ; then #TODO
	export cov=$($samtools depth -a $outputdir/${dilutionname}.bam | awk '{sum+=$3} END {print sum/NR}')
echo $cov > $outputdir/coverage_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.txt 
else
	export cov=$($samtools depth -b $bedfile -a $outputdir/${dilutionname}.bam | awk '{sum+=$3} END {print sum/NR}')
echo $cov > $outputdir/coverage_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.txt
fi
fi

#######################################################################################################################
echo "Step 3c: Estimate mixture tumor fraction (TF) by ichorCNA"
#######################################################################################################################

bash ${repopath}/cfdna_snv_benchmark/utils/run_ichorcna_chr.sh $outputdir/${dilutionname}.bam $ichorcnaextdata $condapath $chr

echo "DONE MIXTURE CREATION"
