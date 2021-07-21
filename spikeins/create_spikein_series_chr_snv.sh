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

echo $inputsample_name
echo $inputsample_bamfile
echo $commonmutationdir
echo $tfs
echo $chr
echo $outputdir

if [ ! -d $outputdir ] ; then mkdir $outputdir ; fi

# prepare input healthy sample
echo 'prepare input healthy sample chr...'
if [ ! -f $inputsample_dsfile ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $inputsample_bamfile $chr > $inputsample_dsfile ; fi
if [ ! -f ${inputsample_dsfile}.bai ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $inputsample_dsfile ; fi

# for each subtype
for st in $commonmutationdir/*/ ;

do echo "$st"
export subtype=$(basename "$st")
echo $subtype
if [ ! -d $outputdir/$subtype ] ; then mkdir $outputdir/$subtype ; fi

# for each tumor fraction
for tf in $tfs ; 

do if [ ! -d $outputdir/$subtype/$tf ] ; then mkdir $outputdir/$subtype/$tf ; fi

# prepare bedfile SNV
export bedfile_snv=$outputdir/$subtype/$tf/${subtype}_${tf}_snv.bed
cp $commonmutationdir/$subtype/chr${chr}_${subtype}_SNV.bed $bedfile_snv
awk -v tfv=$tf '$4=tfv*$4' $bedfile_snv > $outputdir/$subtype/$tf/tmp_snv.bed 
sed -e 's/    /\t/g'  $outputdir/$subtype/$tf/tmp_snv.bed  > $bedfile_snv
rm $outputdir/$subtype/$tf/tmp_snv.bed

# run bamsurgeon with snv
echo 'run bamsurgeon SNV...'
if [ ! -f $outputdir/$subtype/$tf/${inputsample_name}_chr${chr}_${subtype}_${tf}_snv.bam ] && [ ! -f $outputdir/$subtype/$tf/${inputsample_name}_chr${chr}_${subtype}_${tf}_snv.sorted.bam ];
then python3 /mnt/projects/carriehc/cfDNA/utils/bamsurgeon/bin/addsnv.py \
    -v $bedfile_snv \
    -f $inputsample_dsfile \
    -r $refgenome \
    -o $outputdir/$subtype/$tf/${inputsample_name}_chr${chr}_${subtype}_${tf}_snv.bam \
    --picardjar $picardjarfile --aligner $aligner --mindepth $mindepth --maxdepth $maxdepth --tagreads \
    --tmpdir $outputdir/$subtype/$tf/addsnv.tmp -p 4 --seed 1 ;
fi
echo 'index output file...'
if [ ! -f $outputdir/$subtype/$tf/${inputsample_name}_chr${chr}_${subtype}_${tf}_snv.sorted.bam ]; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools sort -o $outputdir/$subtype/$tf/${inputsample_name}_chr${chr}_${subtype}_${tf}_snv.sorted.bam -@ 4 $outputdir/$subtype/$tf/${inputsample_name}_chr${chr}_${subtype}_${tf}_snv.bam ; fi
if [ ! -f $outputdir/$subtype/$tf/${inputsample_name}_chr${chr}_${subtype}_${tf}_snv.sorted.bam.bai ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index -@ 4 $outputdir/$subtype/$tf/${inputsample_name}_chr${chr}_${subtype}_${tf}_snv.sorted.bam ; fi
rm -r $outputdir/$subtype/$tf/addsnv.tmp
rm $outputdir/$subtype/$tf/${inputsample_name}_chr${chr}_${subtype}_${tf}_snv.bam
done
done
