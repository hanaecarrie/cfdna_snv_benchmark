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


while getopts ":c:v:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
    v) vaf=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

source $condapath
conda activate default

echo $sample_pseudohealthy
echo $sample_buffycoat
echo $samplename_pseudohealthy
echo $samplename_buffycoat
echo $chr
echo $vaf
echo $mutsbed_snv
echo $mutsbed_indel

export outputdir=$outputfolder/spikeins_chr${chr}_${samplename_knowntumormuts}_${samplename_pseudohealthy}/spikeins_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}
if [ ! -d $outputfolder/spikeins_chr${chr}_${samplename_knowntumormuts}_${samplename_pseudohealthy} ] ; then mkdir $outputfolder/spikeins_chr${chr}_${samplename_knowntumormuts}_${samplename_pseudohealthy} ; fi
echo $outputdir
echo $pseudohealthydir
echo $buffycoatdir

export sample_pseudohealthy_chr=$pseudohealthydir/${samplename_pseudohealthy}_chr${chr}.bam
export sample_buffycoat_chr=$buffycoatdir/${samplename_buffycoat}_chr${chr}.bam
echo $sample_pseudohealthy_chr

if [ ! -d $outputdir ] ; then mkdir $outputdir ; fi

# Select chr only
echo "Select chr ${chr} only for the pseudohealthy sample..."
if [ ! -d $pseudohealthydir ] ; then mkdir $pseudohealthydir ; fi

if [ ! -f $sample_pseudohealthy_chr ] ; then $samtools view -b $sample_pseudohealthy $chr > $sample_pseudohealthy_chr ; fi
if [ ! -f ${sample_pseudohealthy_chr}.bai ] ; then $samtools index $sample_pseudohealthy_chr ; fi

####################

export pseudohealthy_chr_coverage=$pseudohealthydir/coverage_${samplename_pseudohealthy}_chr${chr}.txt
if  [ ! -f $pseudohealthy_chr_coverage ] ; then export pseudohealthy_cov=$($samtools depth -a $sample_pseudohealthy_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $pseudohealthy_cov >  $pseudohealthy_chr_coverage ; else export pseudohealthy_cov=$(cat $pseudohealthy_chr_coverage) ; fi
echo $pseudohealthy_cov
echo $dilutionfactor_pseudohealthy # 150x

export dilutionfraction_pseudohealthy=$(echo "scale=4; $dilutionfactor_pseudohealthy / $pseudohealthy_cov" | bc -l)
if [ $(echo "$dilutionfraction_pseudohealthy>1" | bc) == 1 ]; then dilutionfraction_pseudohealthy=1 ; fi
echo 'Dilution fraction pseudohealthy'
echo $dilutionfraction_pseudohealthy

export sample_pseudohealthy_chr_downsample=$pseudohealthydir/${samplename_pseudohealthy}_chr${chr}_${dilutionfactor_pseudohealthy}x.bam
echo $sample_pseudohealthy_chr_downsample

echo "subsample pseudohealthy sample as desired..."

if [ $dilutionfraction_pseudohealthy != 1 ] ;
then if [ ! -f $sample_pseudohealthy_chr_downsample ] ; then $samtools view -b -s $dilutionfraction_pseudohealthy -o $sample_pseudohealthy_chr_downsample $sample_pseudohealthy_chr ; fi
else if [ ! -f $sample_pseudohealthy_chr_downsample ] ; then cp $sample_pseudohealthy_chr $sample_pseudohealthy_chr_downsample ; fi
fi
if [ ! -f ${sample_pseudohealthy_chr_downsample}.bai ] ; then $samtools index $sample_pseudohealthy_chr_downsample ; fi

####################

export sample_pseudohealthy_chr=$sample_pseudohealthy_chr_downsample

if [ $vaf == 0 ] ; then cp $sample_pseudohealthy_chr $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.sorted.bam 
cp ${sample_pseudohealthy_chr}.bai  $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.sorted.bam.bai ; fi

if [ $vaf != 0 ] ; then echo $vaf ; 
 
# prepare bedfile SNV
cp $mutsbed_snv $outputdir/$(basename $mutsbed_snv "1.bed")${vaf}.bed
cp $mutsbed_indel $outputdir/$(basename $mutsbed_indel "1.bed")${vaf}.bed
awk -v tfv=$vaf '$4=tfv' $mutsbed_snv > $outputdir/tmp.bed
sed -e 's/ /\t/g'  $outputdir/tmp.bed  > $outputdir/$(basename $mutsbed_snv "1.bed")${vaf}.bed
rm $outputdir/tmp.bed
export mutsbed_snv=$outputdir/$(basename $mutsbed_snv "1.bed")${vaf}.bed
awk -v tfv=$vaf '$4=tfv' $mutsbed_indel > $outputdir/tmp.bed
sed -e 's/ /\t/g'  $outputdir/tmp.bed  > $outputdir/$(basename $mutsbed_indel "1.bed")${vaf}.bed
rm $outputdir/tmp.bed
export mutsbed_indel=$outputdir/$(basename $mutsbed_indel "1.bed")${vaf}.bed

# change working directory
cd $outputdir

# run bamsurgeon with snv
echo 'run bamsurgeon SNV...'
if [ ! -f $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.bam ] && [ ! -f $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.sorted.bam ];
then python3 $bamsurgeonpath/bamsurgeon/bin/addsnv.py \
    -v $mutsbed_snv \
    -f $sample_pseudohealthy_chr \
    -r $refgenome \
    -o $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.bam \
    --picardjar $picardjarfile --aligner $aligner --mindepth $mindepth --maxdepth $maxdepth --tagreads \
    --tmpdir $outputdir/addsnv.tmp -p 4 --seed 1 ;
fi
echo 'index output file...'
if [ ! -f $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.sorted.bam ]; then $samtools sort -o $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.sorted.bam -@ 4 $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.bam ; fi
if [ ! -f $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.sorted.bam.bai ] ; then $samtools index -@ 4 $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.sorted.bam ; fi
rm -r $outputdir/addsnv.tmp
rm $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.bam

# run bamsurgeon with indel
echo 'run bamsurgeon INDEL...'
if [ ! -f $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.bam ] && [ ! -f $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.sorted.bam ];
then python3 $bamsurgeonpath/bamsurgeon/bin/addindel.py \
    -v $mutsbed_indel \
    -f $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.sorted.bam \
    -r $refgenome \
    -o $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.bam \
    --picardjar $picardjarfile --aligner $aligner --mindepth $mindepth --maxdepth $maxdepth --tagreads \
    --tmpdir $outputdir/addindel.tmp -p 4 --seed 1 ;
fi
echo 'sort and index output file...'
if  [ ! -f $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.sorted.bam ] ; then $samtools sort -o $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.sorted.bam -@ 4 $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.bam ; fi
if [ ! -f $outputdir/spikein_chr${chr}_${samplename_pseudohealthy}_vaf${vaf}_snv_indel.sorted.bam.bai ] ; then $samtools index -@ 4 $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.sorted.bam ; fi
rm -r $outputdir/addindel.tmp
rm $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.bam

rm $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.sorted.bam
rm $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv.sorted.bam.bai

fi

if [ -f  $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.sorted.bam ] ; then mv  $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.sorted.bam  $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}.bam ; fi
if [ -f  $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.sorted.bam.bai ] ; then mv $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}_snv_indel.sorted.bam.bai  $outputdir/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}.bam.bai ; fi

# check buffy coat select chr exists
echo "buffy coat..."
if [ ! -d $buffycoatdir ] ; then mkdir $buffycoatdir ; fi
if [ ! -f $sample_buffycoat_chr ] ; then $samtools view -b $sample_buffycoat $chr > $sample_buffycoat_chr ; fi
if [ ! -f ${sample_buffycoat_chr}.bai ] ; then $samtools index $sample_buffycoat_chr ; fi

