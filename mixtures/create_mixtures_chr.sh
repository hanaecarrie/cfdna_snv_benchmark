#!/bin/bash

# def path and environment
source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
conda activate default


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
export outputdir=$outputfolder/mixtures_${samplename_tumor}/mixture_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}
if [ ! -d $outputfolder/mixtures_${samplename_tumor} ] ; then mkdir $outputfolder/mixtures_${samplename_tumor} ; fi
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

if [ ! -d $outputdir ] ; then mkdir $outputdir ; fi

# Select chr only
echo "Select chr ${chr} only for the tumor and the healthy sample..."
if [ ! -d $tumordir ] ; then mkdir $tumordir ; fi

if [ ! -f $sample_tumor_chr ] ; then $samtools view -b $sample_tumor $chr > $sample_tumor_chr ; fi
if [ ! -f ${sample_tumor_chr}.bai ] ; then $samtools index $sample_tumor_chr ; fi
if  [ ! -f $tumor_chr_coverage ] ; then export tumor_cov=$(/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools depth -a $sample_tumor_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $tumor_cov >  $tumor_chr_coverage ; else export tumor_cov=$(cat $tumor_chr_coverage) ; fi
if [ ! -d $healthydir ] ; then mkdir $healthydir ; fi
if [ ! -f $sample_healthy_chr ] ; then $samtools view -b $sample_healthy $chr > $sample_healthy_chr ; fi
if [ ! -f ${sample_healthy_chr}.bai ] ; then $samtools index $sample_healthy_chr ; fi
if  [ ! -f $healthy_chr_coverage ] ; then export healthy_cov=$(/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools depth -a $sample_healthy_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $healthy_cov >  $healthy_chr_coverage ; else export healthy_cov=$(cat $healthy_chr_coverage) ; fi
echo $tumor_cov
echo $healthy_cov
echo $dilutionfactor_tumor
echo $dilutionfactor_healthy
                                   
export dilutionfraction_tumor=$(echo "scale=4; $dilutionfactor_tumor / $tumor_cov" | bc -l)
if [ $(echo "$dilutionfraction_tumor>1" | bc) == 1 ]; then dilutionfraction_tumor=1 ; fi
export dilutionfraction_healthy=$(echo "scale=4; $dilutionfactor_healthy / $healthy_cov" | bc -l)
if [ $(echo "$dilutionfraction_healthy>1" | bc) == 1 ]; then dilutionfraction_healthy=1 ; fi
echo $dilutionfraction_tumor
echo $dilutionfraction_healthy

export sample_tumor_chr_downsample=$tumordir/${samplename_tumor}_chr${chr}_downsampled_${dilutionfactor_tumor}.bam
export sample_healthy_chr_downsample=$healthydir/${samplename_healthy}_chr${chr}_downsampled_${dilutionfactor_healthy}.bam
echo $sample_tumor_chr_downsample
echo $sample_healthy_chr_downsample
export dilutionname=dilution_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}
export rgname=rg_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt
echo $dilutionname
echo $rgname

echo "subsample tumor and healthy samples as desired..."

if [ $dilutionfraction_tumor != 1 ] ;
then if [ ! -f $sample_tumor_chr_downsample ] ; then $samtools view -b -s $dilutionfraction_tumor -o $sample_tumor_chr_downsample $sample_tumor_chr ; fi
else if [ ! -f $sample_tumor_chr_downsample ] ; then cp $sample_tumor_chr $sample_tumor_chr_downsample ; fi
fi


if [ $dilutionfraction_healthy != 0 ] ;

then if [ $dilutionfraction_healthy != 1 ] ;
then if [ ! -f $sample_healthy_chr_downsample ] ; then $samtools view -b -s $dilutionfraction_healthy -o $sample_healthy_chr_downsample $sample_healthy_chr ; fi

else if [ ! -f $sample_healthy_chr_downsample ] ; then cp $sample_healthy_chr $sample_healthy_chr_downsample ; fi
fi

echo "label reads coming from the tumor and ones from the merged healthy..."
if [ ! -f $outputdir/$rgname ] ; then perl -e 'print "@RG\\tID:tumor\\tSM:hs\\tLB:tumor\\tPL:Illumina\\n@RG\\tID:healthy\\tSM:hs\\tLB:healthy\\tPL:Illumina\\n"' > $outputdir/$rgname ; fi

echo "merge helthy and tumor reads into a silico mixture sample..."
if [ ! -f $outputdir/${dilutionname}.sorted.bam ] ; then $samtools merge -rh $outputdir/$rgname $outputdir/${dilutionname}.bam $sample_tumor_chr_downsample $sample_healthy_chr_downsample ; fi

# dilutionfactor healthy == 0
else  cp  $sample_tumor_chr_downsample $outputdir/${dilutionname}.bam ;
fi

# sort
echo "sort mixture..."
if [ ! -f $outputdir/${dilutionname}.sorted.bam ] ; then $samtools sort -o $outputdir/${dilutionname}.sorted.bam -@ 4 $outputdir/${dilutionname}.bam ; fi
# index
echo "index mixture..."
if [ ! -f $outputdir/${dilutionname}.sorted.bam.bai ] ; then $samtools index -@ 4 $outputdir/${dilutionname}.sorted.bam ; fi

if [  -f $outputdir/${dilutionname}.bam ] ; then rm $outputdir/${dilutionname}.bam ; fi

# check buffy coat select chr exists
echo "buffy coat..."
if [ ! -d $buffycoatdir ] ; then mkdir $buffycoatdir ; fi
if [ ! -f $sample_buffycoat_chr ] ; then $samtools view -b $sample_buffycoat $chr > $sample_buffycoat_chr ; fi
if [ ! -f ${sample_buffycoat_chr}.bai ] ; then  $samtools index $sample_buffycoat_chr ; fi

echo "estimate tumor burden by calculation..."
if [ ! -f $outputdir/estimated_tf_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt ] ; then
echo $tffile
while read line ; do export A="$(cut -d',' -f1 <<<"$line")" ;
echo $A
if [ $A == $samplename_tumor ] ; then echo $A ; export median_tumor_burden="$(cut -d ',' -f3 <<<"$line")" ; export cov_tumor="$(cut -d ',' -f2 <<<"$line")" ; fi ; done < $tffile
echo $median_tumor_burden
cov_t=$(echo "$median_tumor_burden * $tumor_cov * $dilutionfraction_tumor" | bc)
echo $cov_t
cov_tot=$(echo "$healthy_cov * $dilutionfraction_healthy + $tumor_cov * $dilutionfraction_tumor" | bc)
echo $cov_tot
mixed_sample_tf=$(echo "$cov_t / $cov_tot" | bc -l)
echo $mixed_sample_tf
echo $mixed_sample_tf > $outputdir/estimated_tf_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt
fi

echo "calculate coverage..."
if [ ! -f $outputdir/coverage_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt ] ; then
export cov=$(/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools depth -a $outputdir/${dilutionname}.sorted.bam | awk '{sum+=$3} END {print sum/NR}')
echo $cov > $outputdir/coverage_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt  ; fi

echo "estimate tumor burden by ichorCNA..."
bash /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/mixtures/run_ichorcna_chr.sh $outputdir/${dilutionname}.sorted.bam $ichorcnaextdata $chr

