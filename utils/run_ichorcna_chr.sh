#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input a bamfile, a path to ichorCNA extdata folder, a condapath, and a chromosome.
# step 1: create read count WIG file using readCounter
# step 2: run ichorCNA while restricting training and analysis to the given chromosome.

if [ $# == 0 ]; then
    echo "Usage: $0 bamfile ichorcnaextdata condapath chrom [downsampleto5x]"
    echo "* bamfile: string. full path to the cfDNA bam file with bam.bai indexfile located in the same directory."
    echo "* ichorcnaextdata: string. full path to the folder containing ichorCNA extdata."
    echo "* condapath: string. full path to the conda.sh file to source."
    echo "* chrom: string. chrom number from '1' to '22', 'X' or 'Y' or 'all' to say all autosomes from chrom 1 to 22."
    echo "* downsampleto5x: bool. optional. default false. if true, downsample input bam to 5x before running ichorCNA."
    echo "Example:"
    echo "$ cd cfdna_snv_benchmark/utils/"
    echo "$ bash $0 /home/users/astar/gis/carriehc/scratch-LOCCG/carriehc/data/mixtures/mixtures_chr1/mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T/mixture_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_80x/mixture_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_80x.bam \
    /home/users/astar/gis/carriehc/.conda/envs/ichorCNA/lib/R/library/ichorCNA/extdata/ /apps/anaconda3-individual-edition/2020.11/etc/profile.d/conda.sh '1'"
    exit 1
fi

export bamfile=$1
export outputdir=$(dirname $bamfile)
export saveid=$(basename $bamfile .bam)
export ichorcna_extdata=$2
export condapath=$3
export chr=$4
export downsampleto5x=$5
if [ "$#" -lt 5 ]; then
  export downsampleto5x=false
else
  export downsampleto5x=$5
fi
echo "Bamfile: "$bamfile
echo "File ID: " $saveid
echo "Output directory: " $outputdir
echo "Path to ichorCNA extdata folder: " $ichorcna_extdata
echo "Path to conda file: " $condapath
echo "Downsample to 5x: " $downsampleto5x

# activate conda environment
source $condapath
conda activate ichorCNA

export outputdir=$outputdir/ichorcna/
if [ ! -f $outputdir ] ; then mkdir $outputdir ; fi

if [ $downsampleto5x ] ; then
  #####################################################################################################################
  echo "Downsample input bam to 5x cov"
  #####################################################################################################################
  # calculate cov
  export cov=$(samtools depth -a $bamfile | awk '{sum+=$3} END {print sum/NR}')
  echo $cov > $outputdir/coverage_${saveid}.txt
  echo "Coverage initial bam: ${cov}x"
  # derive downsample factor
  export factor=$(echo "scale=4; 5 / $cov" | bc -l)
  echo "Dilution factor to apply: " $factor
  # downsample using samtools view
  samtools view -b -s $factor -@ 4 -o $outputdir/${saveid}_5x.bam $bamfile
  # use new bam to run ichorCNA
  export bamfile=$outputdir/${saveid}_5x.bam
  export saveid=$(basename $bamfile .bam)
  echo "New bamfile: "$bamfile
  echo "New file ID: " $saveid
  export newcov=$(samtools depth -a $bamfile | awk '{sum+=$3} END {print sum/NR}')
  echo $newcov > $outputdir/coverage_${saveid}.txt
  echo "Coverage initial bam: ${newcov}x"
fi

#######################################################################################################################
echo "Start create read count WIG file $outputdir/${saveid}.wig"
#######################################################################################################################

if [ ! -f "$outputdir/${saveid}.wig" ] ; then
  if [[ $chr == 'all' ]] ; then
    readCounter --window 1000000 --quality 20 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22" $bamfile > $outputdir/${saveid}.wig
  else
    readCounter --window 50000 --quality 20 --chromosome "$chr" $bamfile > $outputdir/${saveid}.wig
  fi
fi
echo "Done WIG file"

#######################################################################################################################
echo "ichorCNA processing $outputdir/${saveid}.wig"
#######################################################################################################################

if [[ $chr == 'all' ]] ; then
  export chr="c(1:22)"
  export binsize='1000'
else
  export binsize='50'
fi
echo "binsize: ${binsize}kb"

if [ ! -f "$outputdir/${saveid}.params.txt" ] ; then
  Rscript runIchorCNA.R \
  --id $saveid \
  --WIG $outputdir/${saveid}.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --gcWig $ichorcna_extdata/gc_hg19_${binsize}kb.wig \
  --mapWig $ichorcna_extdata/map_hg19_${binsize}kb.wig \
  --centromere $ichorcna_extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
  --includeHOMD False --chrs "$chr" --chrTrain "$chr" --chrNormalize "$chr" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 \
  --outDir $outputdir
fi
echo "Done ichorCNA"

