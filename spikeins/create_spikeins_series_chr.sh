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

while getopts ":c:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

echo $sample_pseudohealthy
echo $sample_buffycoat
echo $samplename_pseudohealthy
echo $samplename_buffycoat
echo $chr
echo $outputfolder
echo $mutsbed
echo $type

export outputdir=$outputfolder/spikeins_${samplename_pseudohealthy}
if [ ! -d $outputfolder ] ; then mkdir $outputfolder ; fi
echo $outputdir
echo $pseudohealthydir
echo $buffycoatdir

export sample_pseudohealthy_chr=$pseudohealthydir/${samplename_pseudohealthy}_chr${chr}.bam
export sample_buffycoat_chr=$buffycoatdir/${samplename_buffycoat}_chr${chr}.bam
echo $sample_pseudohealthy_chr

if [ ! -d $outputdir ] ; then mkdir $outputdir ; fi

cp $config_file $outputdir/

# Select chr only
echo "Select chr ${chr} only for the pseudohealthy sample..."
if [ ! -d $pseudohealthydir ] ; then mkdir $pseudohealthydir ; fi

if [ ! -f $sample_pseudohealthy_chr ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $sample_pseudohealthy $chr > $sample_pseudohealthy_chr ; fi
if [ ! -f ${sample_pseudohealthy_chr}.bai ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $sample_pseudohealthy_chr ; fi

# check buffy coat select chr exists
echo "buffy coat..."
if [ ! -d $buffycoatdir ] ; then mkdir $buffycoatdir ; fi
if [ ! -f $sample_buffycoat_chr ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $sample_buffycoat $chr > $sample_buffycoat_chr ; fi
if [ ! -f ${sample_buffycoat_chr}.bai ] ; then  /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $sample_buffycoat_chr ; fi

for vaf in $vafs ;

do echo $vaf
export outputdirnew=$outputfolder/spikeins_${samplename_pseudohealthy}/spikeins_chr${chr}_${samplename_pseudohealthy}_${vaf}
if [ ! -d $outputdirnew/logs ] ; then mkdir $outputdirnew/logs ; fi
qsub -pe OpenMP 4 -l mem_free=48G,h_rt=24:00:00 -o $outputdirnew/logs/ -e $outputdirnew/logs/ /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/spikeins/create_spikeins_chr.sh -c $config_file -v $vaf

done

