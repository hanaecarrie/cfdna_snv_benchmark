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

if [ ! -d $outputfolder ] ; then mkdir $outputfolder ; fi

export outputdir=$outputfolder/dilutions_${samplename_tumor}/dilution_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}
if [ ! -d $outputfolder/dilutions_${samplename_tumor} ] ; then mkdir $outputfolder/dilutions_${samplename_tumor} ; fi
echo $outputdir
echo $tumordir
echo $healthydir
echo $buffycoatdir

export sample_tumor_chr=$tumordir/${samplename_tumor}_chr${chr}.bam
export sample_healthy_chr=$healthydir/${samplename_healthy}_chr${chr}.bam
export sample_buffycoat_chr=$buffycoatdir/${samplename_buffycoat}_chr${chr}.bam
echo $sample_tumor_chr
echo $sample_healthy_chr

export tumor_chr_coverage=$tumordir/coverage_${samplename_tumor}_chr${chr}.txt
export healthy_chr_coverage=$healthydir/coverage_${samplename_healthy}_chr${chr}.txt
echo $tumor_chr_coverage
echo $healthy_chr_coverage

cp $config_file $outputfolder/

# Select chr only
echo "Select chr ${chr} only for the tumor and the healthy sample..."
if [ ! -d $tumordir ] ; then mkdir $tumordir ; fi

if [ ! -f $sample_tumor_chr ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $sample_tumor $chr > $sample_tumor_chr ; fi
if [ ! -f ${sample_tumor_chr}.bai ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $sample_tumor_chr ; fi
if  [ ! -f $tumor_chr_coverage ] ; then export tumor_cov=$(/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools depth -a $sample_tumor_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $tumor_cov >  $tumor_chr_coverage ; else export tumor_cov=$(cat $tumor_chr_coverage) ; fi
if [ ! -d $healthydir ] ; then mkdir $healthydir ; fi
if [ ! -f $sample_healthy_chr ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $sample_healthy $chr > $sample_healthy_chr ; fi
if [ ! -f ${sample_healthy_chr}.bai ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $sample_healthy_chr ; fi
if  [ ! -f $healthy_chr_coverage ] ; then export healthy_cov=$(/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools depth -a $sample_healthy_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $healthy_cov >  $healthy_chr_coverage ; else export healthy_cov=$(cat $healthy_chr_coverage) ; fi

# check buffy coat select chr exists
echo "buffy coat..."
if [ ! -d $buffycoatdir ] ; then mkdir $buffycoatdir ; fi
if [ ! -f $sample_buffycoat_chr ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $sample_buffycoat $chr > $sample_buffycoat_chr ; fi
if [ ! -f ${sample_buffycoat_chr}.bai ] ; then  /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $sample_buffycoat_chr ; fi

for dilutionfactor in $dilutionfactors ;

do qsub -pe OpenMP 4 -l mem_free=48G,h_rt=24:00:00 -o /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/logs/ -e /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/logs/ /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/dilutions/create_dilution_chr.sh -c $config_file -d $dilutionfactor

done

