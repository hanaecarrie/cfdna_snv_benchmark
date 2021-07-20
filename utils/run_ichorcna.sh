#!/bin/bash

# def path and environment
conda activate default

export bamfile=$1
export outputdir=$(dirname $bamfile)
export saveid=$(basename $bamfile .bam)
export ichorcna_extdata=$2
echo $bamfile
echo $saveid
echo $outputdir
echo $ichorcna_extdata  # /mnt/projects/carriehc/cfDNA/anaconda3/envs/default/lib/R/library/ichorCNA/extdata

if [ ! -f "$outputdir/${saveid}.wig" ] ; then
echo "Create WIG file $output_dir/${saveid}.wig"
# Create WIG Files
/home/carriehc/programs/hmmcopy_utils/bin/readCounter \
--window 1000000 --quality 20 \
--chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
$bamfile > $outputdir/${saveid}.wig
echo "Done WIG file"
fi

# ichorCNA run
echo "ichorCNA processing $outputdir/${saveid}.wig"
if [ ! -f "$outputdir/${saveid}.cna.seg" ] ; then
/mnt/projects/carriehc/cfDNA/anaconda3/envs/default/bin/Rscript /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/utils/runIchorCNA.R \
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

