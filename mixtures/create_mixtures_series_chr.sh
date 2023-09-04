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

echo $config_file
echo $sample_healthy
echo $sample_tumor
echo $sample_buffycoat 

echo $samplename_healthy
echo $samplename_tumor
echo $samplename_buffycoat

echo $dilutionfactors
echo $chr
echo $outputfolder
echo $tffile
echo $samtools

# def path and environment
source $condapath
conda activate default

if [ ! -d $outputfolder ] ; then mkdir $outputfolder ; fi
if [ ! -d $outputfolder/mixtures_chr${chr}_${samplename_tumor}_${samplename_healthy} ] ; then mkdir $outputfolder/mixtures_chr${chr}_${samplename_tumor}_${samplename_healthy} ; fi
echo $tumordir
echo $healthydir
echo $buffycoatdir

if [ ! -d $tumordir ] ; then mkdir $tumordir ; fi
if [ ! -d $healthydir ] ; then mkdir $healthydir ; fi
if [ ! -d $buffycoatdir ] ; then mkdir $buffycoatdir ; fi
export sample_tumor_chr=$tumordir/${samplename_tumor}_chr${chr}.bam
export sample_healthy_chr=$healthydir/${samplename_healthy}_chr${chr}.bam
export sample_buffycoat_chr=$buffycoatdir/${samplename_buffycoat}_chr${chr}.bam
echo $sample_tumor_chr
echo $sample_healthy_chr

export tumor_chr_coverage=$tumordir/coverage_${samplename_tumor}_chr${chr}.txt
export healthy_chr_coverage=$healthydir/coverage_${samplename_healthy}_chr${chr}.txt
export buffycoat_chr_coverage=$buffycoatdir/coverage_${samplename_buffycoat}_chr${chr}.txt
echo $tumor_chr_coverage
echo $healthy_chr_coverage
echo $buffycoat_chr_coverage
echo $machine

cp $config_file $outputfolder/mixtures_chr${chr}_${samplename_tumor}_${samplename_healthy}/ 

# Select chr only
echo "Select chr ${chr} only for the tumor and the healthy sample..."
if [ ! -d $healthydir ] ; then mkdir $healthydir ; fi
if [[ $chr == 'all' ]] ; then
	if [ ! -f $sample_healthy_chr ] ; then cp $sample_healthy $sample_healthy_chr ; fi
else
	if [ ! -f $sample_healthy_chr ] ; then $samtools view -b $sample_healthy $chr > $sample_healthy_chr ; fi
fi
if [ ! -f ${sample_healthy_chr}.bai ] ; then $samtools index $sample_healthy_chr ; fi
if [[ $bedfile == 'wgs' ]] ; then
if  [ ! -f $healthy_chr_coverage ] ; then export healthy_cov=$($samtools depth -a $sample_healthy_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $healthy_cov >  $healthy_chr_coverage ; else export healthy_cov=$(cat $healthy_chr_coverage) ; fi
else 
if  [ ! -f $healthy_chr_coverage ] ; then export healthy_cov=$($samtools depth -a -b $bedfile $sample_healthy_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $healthy_cov >  $healthy_chr_coverage ; else export healthy_cov=$(cat $healthy_chr_coverage) ; fi
fi


if [ ! -d $tumordir ] ; then mkdir $tumordir ; fi
if [[ $chr == 'all' ]] ; then
	if [ ! -f $sample_tumor_chr ] ; then cp $sample_tumor $sample_tumor_chr ; fi
else
	if [ ! -f $sample_tumor_chr ] ; then $samtools view -b $sample_tumor $chr > $sample_tumor_chr ; fi
fi
if [ ! -f ${sample_tumor_chr}.bai ] ; then $samtools index $sample_tumor_chr ; fi
if [[ $bedfile == 'wgs' ]] ; then
if  [ ! -f $tumor_chr_coverage ] ; then export tumor_cov=$($samtools depth -a $sample_tumor_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $tumor_cov >  $tumor_chr_coverage ; else export tumor_cov=$(cat $tumor_chr_coverage) ; fi
else
if  [ ! -f $tumor_chr_coverage ] ; then export tumor_cov=$($samtools depth -a -b $bedfile $sample_tumor_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $tumor_cov >  $tumor_chr_coverage ; else export tumor_cov=$(cat $tumor_chr_coverage) ; fi
fi


for dilutionfactor in $dilutionfactors ;

do echo $dilutionfactor
export dilutionfactor_tumor=$(echo $dilutionfactor | cut -f1 -d-)
export dilutionfactor_healthy=$(echo $dilutionfactor | cut -f2 -d-)
export outputdir=$outputfolder/mixtures_chr${chr}_${samplename_tumor}_${samplename_healthy}/mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x
if [ ! -d $outputdir ] ; then mkdir $outputdir ; fi
if [ ! -d $outputdir/logs ] ; then mkdir $outputdir/logs ; fi

if [ $machine == 'Aquila' ] ; then /opt/uge-8.1.7p3/bin/lx-amd64/qsub -pe OpenMP 4 -l mem_free=24G,h_rt=24:00:00 -o $outputdir/logs/ -e $outputdir/logs/ $repopath/cfdna_snv_benchmark/mixtures/create_mixtures_chr.sh -c $config_file -d $dilutionfactor ; fi
if [ $machine == 'Ronin' ] ; then
	nohup $repopath/cfdna_snv_benchmark/mixtures/create_mixtures_chr.sh -c $config_file -d $dilutionfactor > $outputdir/logs/log_mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.out &
fi
if [ $machine == 'Pluto' ] ; then
	sbatch -p normal -J chr${chr}_mix -t 24:00:00 -N 1 --mem 24G --output $outputdir/logs/log_mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.o --error $outputdir/logs/log_mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.e $repopath/cfdna_snv_benchmark/mixtures/create_mixtures_chr.sh -c $config_file -d $dilutionfactor
fi

done

# check buffy coat select chr exists
echo "buffy coat..."
if [ ! -d $buffycoatdir ] ; then mkdir $buffycoatdir ; fi
if [[ $chr == 'all' ]] ; then
	if [ ! -f $sample_buffycoat_chr ] ; then cp $sample_buffycoat $sample_buffycoat_chr ; fi
else
	if [ ! -f $sample_buffycoat_chr ] ; then $samtools view -b $sample_buffycoat $chr > $sample_buffycoat_chr ; fi
fi
if [ ! -f ${sample_buffycoat_chr}.bai ] ; then  $samtools index $sample_buffycoat_chr ; fi
if [[ $bedfile == 'wgs' ]] ; then
	if  [ ! -f $buffycoat_chr_coverage ] ; then export buffycoat_cov=$($samtools depth -a $sample_buffycoat_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $buffycoat_cov >  $buffycoat_chr_coverage ; else export buffycoat_cov=$(cat $buffycoat_chr_coverage) ; fi
else
	if  [ ! -f $buffycoat_chr_coverage ] ; then export buffycoat_cov=$($samtools depth -a -b $bedfile $sample_buffycoat_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $buffycoat_cov >  $buffycoat_chr_coverage ; else export buffycoat_cov=$(cat $buffycoat_chr_coverage) ; fi
fi

