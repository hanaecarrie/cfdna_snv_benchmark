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

while getopts ":c:t:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
    t) tumorbam=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

source ${condapath}
echo $PATH
conda activate varnet

cd $repopath

echo $normalbam
echo $tumorbam

export samplename=$(basename $tumorbam .bam)
export outputdir=$outdir
export ref=${extdata}/GRCh37/GRCh37.fa
export bed=${extdata}/${mode}_bed/${mode}_hg19_chr${chr}.bed

# mark and remove duplicates in tumor bam
if [ ! -f $(dirname $tumorbam)/$(basename $tumorbam .bam)_rmd.bam ] ; then  
	$java -jar $picard MarkDuplicates \
		-I $tumorbam \
		-O $(dirname $tumorbam)/$(basename $tumorbam .bam)_rmd.bam \
		-M $(dirname $tumorbam)/$(basename $tumorbam .bam)_rmd_metrics.txt \
		--REMOVE_DUPLICATES true \
		--TMP_DIR $(dirname $tumorbam)/tmp
fi
# index rmd tumor bam
if [ ! -f $(dirname $tumorbam)/$(basename $tumorbam .bam)_rmd.bam.bai ] ; then
        $samtools index $(dirname $tumorbam)/$(basename $tumorbam .bam)_rmd.bam
fi

# run Varnet filter.py
$py ${tooldir}/filter.py \
	--sample_name $samplename \
	--normal_bam $normalbam \
	--tumor_bam $tumorbam \
	--processes $processesfilter \
	--output_dir $outputdir \
	--reference $ref \
	--region_bed $bed 


# run Varnet predict.py
$py ${tooldir}/predict.py \
        --include_allele_frequency $includeallelefrequency \
	--sample_name $samplename \
        --normal_bam $normalbam \
        --tumor_bam $tumorbam \
        --processes $processespredict \
        --output_dir $outdir \
        --reference $ref
 
