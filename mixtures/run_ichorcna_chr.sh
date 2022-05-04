#!/bin/bash

#source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
#conda activate ichorcna

export bamfile=$1
export outputdir=$(dirname $bamfile)
export saveid=$(basename $bamfile .bam)
export ichorcna_extdata=$2
export chr=$3
echo $bamfile
echo $saveid
echo $outputdir
echo $ichorcna_extdata  # /mnt/projects/carriehc/cfDNA/anaconda3/envs/ichorcna/lib/R/library/ichorCNA/extdata/

# def path and environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate default

export outputdir=$outputdir/ichorcna/

if [ ! -f $outputdir ] ; then mkdir $outputdir ; fi

if [ ! -f "$outputdir/${saveid}.wig" ] ; then
echo "Create WIG file $outputdir/${saveid}.wig"
# Create WIG Files
#/mnt/projects/carriehc/cfDNA/utils/hmmcopy_utils/bin/readCounter --window 50000 --quality 20 --chromosome "$chr" $bamfile > $outputdir/${saveid}.wig
/home/ubuntu/anaconda3/envs/default/bin/readCounter --window 50000 --quality 20 --chromosome "$chr" $bamfile > $outputdir/${saveid}.wig
echo "Done WIG file"
fi


# ichorCNA run
echo "ichorCNA processing $outputdir/${saveid}.wig"
if [ ! -f "$outputdir/${saveid}.cna.seg" ] ; then
#/mnt/projects/carriehc/cfDNA/anaconda3/envs/ichorcna/bin/Rscript /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/utils/runIchorCNA.R \
/home/ubuntu/anaconda3/envs/default/bin/Rscript /home/ubuntu/cfdna_snv_benchmark/utils/runIchorCNA.R \
--id $saveid \
--WIG $outputdir/${saveid}.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
--gcWig $ichorcna_extdata/gc_hg19_50kb.wig \
--mapWig $ichorcna_extdata/map_hg19_50kb.wig \
--centromere $ichorcna_extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
--includeHOMD False --chrs "$chr" --chrTrain "$chr" --chrNormalize "$chr" \
--estimateNormal True --estimatePloidy True --estimateScPrevalence True \
--scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 \
--outDir $outputdir
echo "Done ichorCNA"
fi

