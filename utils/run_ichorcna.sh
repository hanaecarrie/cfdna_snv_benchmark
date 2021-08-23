#!/bin/bash

source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
# def path and environment
conda activate ichorcna

export bamfile=$1
export outputdir=$(dirname $bamfile)
export saveid=$(basename $bamfile .bam)
export ichorcna_extdata=$3
echo $bamfile
echo $saveid
echo $outputdir
echo $ichorcna_extdata  # /mnt/projects/carriehc/cfDNA/anaconda3/envs/ichorcna/lib/R/library/ichorCNA/extdata/

if [ ! -f "$outputdir/${saveid}.wig" ] ; then
echo "Create WIG file $outputdir/${saveid}.wig"
# Create WIG Files
/mnt/projects/carriehc/cfDNA/utils/hmmcopy_utils/bin/readCounter --window 1000000 --quality 20 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" $bamfile > $outputdir/${saveid}.wig
echo "Done WIG file"
fi

# ichorCNA run
echo "ichorCNA processing $outputdir/${saveid}.wig"
if [ ! -f "$outputdir/${saveid}.cna.seg" ] ; then
/mnt/projects/carriehc/cfDNA/anaconda3/envs/ichorcna/bin/Rscript /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/utils/runIchorCNA.R \
--id $saveid \
--WIG $outputdir/${saveid}.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
--gcWig $ichorcna_extdata/gc_hg19_1000kb.wig \
--mapWig $ichorcna_extdata/map_hg19_1000kb.wig \
--centromere $ichorcna_extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
--normalPanel $ichorcna_extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
--includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
--estimateNormal True --estimatePloidy True --estimateScPrevalence True \
--scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 \
--outDir $outputdir
echo "Done ichorCNA"
fi

