#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
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


while getopts ":c:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

export sample_healthy=$samples_healthy
echo $sample_healthy
echo $samples_tumor
echo  $samples_buffycoat 

export samplename_healthy=$samplenames_healthy
echo $samplename_healthy
echo $samplenames_tumor
echo $samplenames_buffycoat

echo $dilutionfactors
echo $chr
echo $outputfolder
if [ ! -d $outputfolder ] ; then mkdir $outputfolder ; fi

# counter
export c=1
#echo "${samplenames_tumor[0]}"
#echo "${samplenames_tumor[$c]}"

for sample_tumor in $samples_tumor ; do

for dilutionfactor in $dilutionfactors ; do

export samplename_tumor=$(echo $samplename_tumor | cut -f $c -d ' ')
export tumordir=$(echo $tumordir | cut -f $c -d ' ')
export buffycoatdir=$(echo $buffycaotdirs | cut -f $c -d ' ')
export dilutionfactor_tumor=$(echo $dilutionfactor | cut -f1 -d-)
export dilutionfactor_healthy=$(echo $dilutionfactor | cut -f2 -d-)
export outputdir=$outputfolder/dilutions_${samplename_tumor}/dilution_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}
echo $outputdir
echo $tumordir
echo $healthydir
echo $buffycoatdir

export sample_tumor_chr=$tumordir/${samplename_tumor}_chr${chr}.bam
export sample_healthy_chr=$healthydir/${samplename_healthy}_chr${chr}.bam
export sample_buffycoat_chr=$buffycoatdir/${samplename_buffycoat}_chr${chr}.bam
export sample_tumor_chr_downsample=$tumordir/${samplename_tumor}_chr${chr}_downsampled_${dilutionfactor_tumor}.bam
export sample_healthy_chr_downsample=$healthydir/${samplename_healthy}_chr${chr}_downsampled_${dilutionfactor_healthy}.bam
echo $sample_tumor_chr
echo $sample_healthy_chr
echo $sample_tumor_chr_downsample
echo $sample_healthy_chr_downsample
export dilutionname=dilution_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}
export rgname=rg_chr${chr}_${samplename_tumor}_${dilutionfactor_tumor}_${samplename_healthy}_${dilutionfactor_healthy}.txt
echo $dilutionname
echo $rgname


'''

if [ ! -d $outputdir ] ; then mkdir $outputdir ; fi

# Select chr only
echo "Select chr ${chr} only for the tumor and the healthy sample..."
if [ ! -f $sample_tumor_chr ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $sample_tumor $chr > $sample_tumor_chr ; fi
if [ $dilutionfactor_healthy != 0 ] ; 
then if [ ! -f $sample_healthy_chr ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $sample_healthy $chr > $sample_healthy_chr ; fi
fi 

# subsample tumor and healthy samples are desired
echo "subsample tumor and healthy samples as desired..."

if [ $dilutionfactor_tumor != 1 ] ; 
then if [ ! -f $sample_tumor_chr_downsample ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b -s $dilutionfactor_tumor -o $sample_tumor_chr_downsample $sample_tumor_chr ; fi
else if [ ! -f $sample_tumor_chr_downsample ] ; then cp $sample_tumor_chr $sample_tumor_chr_downsample ; fi
fi


if [ $dilutionfactor_healthy != 0 ] ;

then if [ $dilutionfactor_healthy != 1 ] ; 
then if [ ! -f $sample_healthy_chr_downsample ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b -s $dilutionfactor_healthy -o $sample_healthy_chr_downsample $sample_healthy_chr ; fi

else if [ ! -f $sample_healthy_chr_downsample ] ; then cp $sample_healthy_chr $sample_healthy_chr_downsample ; fi
fi

# label reads coming from the tumor and ones from the merged healthy
echo "label reads coming from the tumor and ones from the merged healthy..."
if [ ! -f $outputdir/$rgname ] ; then perl -e 'print "@RG\\tID:tumor\\tSM:hs\\tLB:tumor\\tPL:Illumina\\n@RG\\tID:healthy\\tSM:hs\\tLB:healthy\\tPL:Illumina\\n"' > $outputdir/$rgname ; fi
# merge helthy and tumor reads into a silico mixture sample
echo "merge helthy and tumor reads into a silico mixture sample..."
if [ ! -f $outputdir/${dilutionname}.sorted.bam ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools merge -rh $outputdir/$rgname $outputdir/${dilutionname}.bam $sample_tumor_chr_downsample $sample_healthy_chr_downsample ; fi

# dilutionfactor healthy == 0
else  cp  $sample_tumor_chr_downsample $outputdir/${dilutionname}.bam ;
fi

# sort
echo "sort mixture..."
if [ ! -f $outputdir/${dilutionname}.sorted.bam ] ; then  /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools sort -o $outputdir/${dilutionname}.sorted.bam $outputdir/${dilutionname}.bam ; fi
if [  -f $outputdir/${dilutionname}.bam ] ; then rm $outputdir/${dilutionname}.bam ; fi
# index
echo "index mixture..."
if [ ! -f $outputdir/${dilutionname}.sorted.bam.bai ] ; then  /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $outputdir/${dilutionname}.sorted.bam ; fi

# check buffy coat select chr exists
if [ ! -f $sample_buffycoat_chr ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $sample_buffycoat $chr > $sample_buffycoat_chr ; fi
if [ ! -f ${sample_buffycoat_chr}.bai ] ; then  /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $sample_buffycoat_chr ; fi

# filter for mutations of the healthy samples
if [ ! -f $outputdir/${dilutionname}.filtered.sorted.bam ] ; then python3 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/filter_vcf_positions.py $vcffile $outputdir/${dilutionname}.sorted.bam $outputdir/${dilutionname}.filtered.sorted.bam ; fi
if [ ! -f $outputdir/${dilutionname}.filtered.sorted.bam.bai ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $outputdir/${dilutionname}.filtered.sorted.bam ; fi
if [ ! -f $buffycoatdir/${samplename_buffycoat}_chr${chr}.filtered.bam ] ; then python3 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/filter_vcf_positions.py $vcffile $buffycoatdir/${samplename_buffycoat}_chr${chr}.bam $buffycoatdir/${samplename_buffycoat}_chr${chr}.filtered.bam ; fi
if [ ! -f $buffycoatdir/${samplename_buffycoat}_chr${chr}.filtered.bam.bai ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $buffycoatdir/${samplename_buffycoat}_chr${chr}.filtered.bam ; fi
#rm  $outputdir/${dilutionname}.sorted.bam*

# calculate tf and coverage of obtained diluted file
if [ ! -f $outputdir/estimated_tf.txt ] ; then bash /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/calculate_tumor_burden.sh -c $config_file ; fi
if [ ! -f $outputdir/coverage.txt ] ; then bash /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/calculate_coverage.sh -c $config_file ; fi

'''
done
c=$((c+1))
done


