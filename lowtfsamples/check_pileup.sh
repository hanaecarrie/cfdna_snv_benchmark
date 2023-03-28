#!/bin/bash

source /mnt/projects/carriehc/cfDNA/anaconda3/etc/profile.d/conda.sh
# def path and environment
conda activate default

cd /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/lowtfsamples/

export bamfile=$1
export bedfile=$2
export outputfile=$3
export refgenome=$4



echo "######## ${outputfile} #########"
touch $outputfile
while IFS="" read -r p || [ -n "$p" ]
do
  printf '%s\n' "$p"
  printf '%s\n' "$p" >> $outputfile
  chr="$( cut -f1 <<< "$p" )"; if [ $refgenome ==  /mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa ] ; then chr=${chr#"chr"} ; fi ; echo "$chr"
  startpos="$( cut -f2 <<< "$p" )"; echo "$startpos"
  endpos="$( cut -f3 <<< "$p" )"; echo "$endpos"
  name="$( cut -f4 <<< "$p" )"; echo "$name"
  samtools mpileup -f $refgenome -r $chr:$startpos-$endpos $bamfile >> $outputfile
done < $bedfile

