#!/bin/bash

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


while getopts ":c:d:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
    d) dilutionfactor=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

echo $sample_tumor
echo $sample_healthy
echo $sample_buffycoat
echo $samplename_tumor
echo $samplename_healthy
echo $samplename_buffycoat
echo $chr
echo $outputdir
echo $tffile

echo $dilutionfactor
export dilutionfactor_tumor=$(echo $dilutionfactor | cut -f1 -d-)
export dilutionfactor_healthy=$(echo $dilutionfactor | cut -f2 -d-)
export outputdir=$outputfolder/mixtures_chr${chr}_${samplename_tumor}/mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x
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
export buffycoat_chr_coverage=$buffycoatdir/coverage_${samplename_buffycoat}_chr${chr}.txt
echo $tumor_chr_coverage
echo $healthy_chr_coverage
echo $buffycoat_chr_coverage


# Select chr only
echo "Select chr ${chr} only for the tumor and the healthy sample..."
if  [ ! -f $tumor_chr_coverage ] ; then export tumor_cov=$($samtools depth -a $sample_tumor_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $tumor_cov >  $tumor_chr_coverage ; else export tumor_cov=$(cat $tumor_chr_coverage) ; fi
if  [ ! -f $healthy_chr_coverage ] ; then export healthy_cov=$($samtools depth -a $sample_healthy_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $healthy_cov >  $healthy_chr_coverage ; else export healthy_cov=$(cat $healthy_chr_coverage) ; fi
echo $tumor_cov
echo $healthy_cov
echo $dilutionfactor_tumor
echo $dilutionfactor_healthy
                                   
export dilutionfraction_tumor=$(echo "scale=4; $dilutionfactor_tumor / $tumor_cov" | bc -l)
if [ $(echo "$dilutionfraction_tumor>1" | bc) == 1 ]; then dilutionfraction_tumor=1 ; fi
export dilutionfraction_healthy=$(echo "scale=4; $dilutionfactor_healthy / $healthy_cov" | bc -l)
if [ $(echo "$dilutionfraction_healthy>1" | bc) == 1 ]; then dilutionfraction_healthy=1 ; fi
echo 'Dilution fraction tumor'
echo $dilutionfraction_tumor
echo 'Dilution fractio healthy'
echo $dilutionfraction_healthy

export sample_tumor_chr_downsample=$tumordir/${samplename_tumor}_chr${chr}_${dilutionfactor_tumor}x.bam
export sample_healthy_chr_downsample=$healthydir/${samplename_healthy}_chr${chr}_${dilutionfactor_healthy}x.bam
echo $sample_tumor_chr_downsample
echo $sample_healthy_chr_downsample
export dilutionname=dilution_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x
export rgname=rg_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.txt
echo $dilutionname
echo $rgname

echo "estimate tumor burden by calculation..."
#if [ ! -f $outputdir/estimated_tf_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.txt ] ; then
echo $tffile
while read line ; do export A="$(cut -d',' -f1 <<<"$line")" ;
if [[ "$A" == *${samplename_tumor}* ]] ; then echo $A ; export median_tumor_burden="$(cut -d ',' -f3 <<<"$line")" ;  fi ; done < $tffile
echo $median_tumor_burden
cov_t=$(echo "$median_tumor_burden * $tumor_cov * $dilutionfraction_tumor" | bc)
echo $cov_t
cov_tot=$(echo "$healthy_cov * $dilutionfraction_healthy + $tumor_cov * $dilutionfraction_tumor" | bc)
echo $cov_tot
mixed_sample_tf=$(echo "$cov_t / $cov_tot" | bc -l)
echo $mixed_sample_tf
echo $mixed_sample_tf > $outputdir/estimated_tf_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}x_${samplename_healthy}_${dilutionfactor_healthy}x.txt
#fi
