#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input a bamfile, a bedfile, a path to the output file to be created,
# a path to the reference human genome, and a path to the conda.sh file to source
# to extract the pileup output of each locus indicated in the bedfile.

if [ $# == 0 ]; then
    echo "Usage: $0 bamfile bedfile outputfile refgenome"
    echo "* bamfile: string. path to bamfile to analyse. should be sorted and indexed."
    echo "* bedfile: string. path to bedfile with genome coordinates to look for. The coordinate system should match the ref human genome used."
    echo "* outputfile: string. path to the output text file to be created."
    echo "* refgenome: string. path to the reference human genome used i.e. hg19 or hg38."
    echo "* condapath: string. path to the conda.sh file to source."
    echo "Example:"
    echo "$ bash $0 XXX XXX XXX XXX XXX"
    exit 1
fi

export bamfile=$1
export bedfile=$2
export outputfile=$3
export refgenome=$4
export condapath=$5
echo "Bam file: " $bamfile
echo "Bed file: " $bedfile
echo "Output file: " $outputfile
echo "Reference human genome: " $refgenome
echo "Conda path: " $condapath

# activate conda environment
source $condapath
conda activate default

touch $outputfile
while IFS="" read -r p || [ -n "$p" ]
do
  printf '%s\n' "$p"
  printf '%s\n' "$p" >> $outputfile
  chr="$( cut -f1 <<< "$p" )"; if [ $(basename $refgenome) ==  'GRCh37.fa' ] ; then chr=${chr#"chr"} ; fi ; echo "$chr"
  startpos="$( cut -f2 <<< "$p" )"; echo "$startpos"
  endpos="$( cut -f3 <<< "$p" )"; echo "$endpos"
  name="$( cut -f4 <<< "$p" )"; echo "$name"
  samtools mpileup -f $refgenome -r $chr:$startpos-$endpos $bamfile >> $outputfile
done < $bedfile

