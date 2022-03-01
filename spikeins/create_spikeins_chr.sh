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


while getopts ":c:v:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
    v) vaf=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

echo $sample_pseudohealthy
echo $sample_buffycoat
echo $samplename_pseudohealthy
echo $samplename_buffycoat
echo $chr
echo $vaf
echo $mutsbed
echo $type

export outputdir=$outputfolder/spikeins_${samplename_pseudohealthy}/spikeins_chr${chr}_${samplename_pseudohealthy}_${vaf}
if [ ! -d $outputfolder/spikeins_${samplename_pseudohealthy} ] ; then mkdir $outputfolder/spikeins_${samplename_pseudohealthy} ; fi
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

if [ ! -f $sample_pseudohealthy_chr ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $sample_pseudohealthy $chr > $sample_pseudohealthy_chr ; fi
if [ ! -f ${sample_pseudohealthy_chr}.bai ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $sample_pseudohealthy_chr ; fi

                         
# prepare bedfile SNV
cp $mutsbed $outputdir/$(basename $mutsbed "1.bed")${vaf}.bed
awk -v tfv=$vaf '$4=tfv' $mutsbed > $outputdir/tmp.bed
sed -e 's/ /\t/g'  $outputdir/tmp.bed  > $outputdir/$(basename $mutsbed "1.bed")${vaf}.bed
rm $outputdir/tmp.bed
export mutsbed=$outputdir/$(basename $mutsbed "1.bed")${vaf}.bed

# change working directory
cd $outputdir

if [ $type = 'SNV' ]; then 

# run bamsurgeon with snv
echo 'run bamsurgeon SNV...'
if [ ! -f $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_snv.bam ] && [ ! -f $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_snv.sorted.bam ];
then python3 /mnt/projects/carriehc/cfDNA/utils/bamsurgeon/bin/addsnv.py \
    -v $mutsbed \
    -f $sample_pseudohealthy_chr \
    -r $refgenome \
    -o $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_snv.bam \
    --picardjar $picardjarfile --aligner $aligner --mindepth $mindepth --maxdepth $maxdepth --tagreads \
    --tmpdir $outputdir/addsnv.tmp -p 4 --seed 1 ;
fi
echo 'index output file...'
if [ ! -f $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_snv.sorted.bam ]; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools sort -o $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_snv.sorted.bam -@ 4 $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_snv.bam ; fi
if [ ! -f $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_snv.sorted.bam.bai ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index -@ 4 $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_snv.sorted.bam ; fi
rm -r $outputdir/addsnv.tmp
rm $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_snv.bam

elif [ $type = 'INDEL' ] ; then

# run bamsurgeon with indel
echo 'run bamsurgeon INDEL...'
if [ ! -f $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_indel.bam ] && [ ! -f $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_indel.sorted.bam ];
then python3 /mnt/projects/carriehc/cfDNA/utils/bamsurgeon/bin/addindel.py \
    -v $mutsbed \
    -f $sample_pseudohealthy_chr \
    -r $refgenome \
    -o $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_indel.bam \
    --picardjar $picardjarfile --aligner $aligner --mindepth $mindepth --maxdepth $maxdepth --tagreads \
    --tmpdir $outputdir/addindel.tmp -p 4 --seed 1 ;
fi
echo 'sort and index output file...'
if  [ ! -f $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_indel.sorted.bam ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools sort -o $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_indel.sorted.bam -@ 4 $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_indel.bam ; fi
if [ ! -f $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_indel.sorted.bam.bai ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index -@ 4 $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_indel.sorted.bam ; fi
rm -r $outputdir/addindel.tmp
rm $outputdir/${samplename_pseudohealthy}_chr${chr}_vaf${vaf}_indel.bam
fi

# check buffy coat select chr exists
echo "buffy coat..."
if [ ! -d $buffycoatdir ] ; then mkdir $buffycoatdir ; fi
if [ ! -f $sample_buffycoat_chr ] ; then /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $sample_buffycoat $chr > $sample_buffycoat_chr ; fi
if [ ! -f ${sample_buffycoat_chr}.bai ] ; then  /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index $sample_buffycoat_chr ; fi

