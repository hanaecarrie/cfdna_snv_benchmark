#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input XXX

if [ $# == 0 ]; then
    echo "Usage: $0 buffycoatbam patientid datadir extdatadir condapath mode"
    echo "* buffycoatbam: XXX"
    echo "* patientid: cancer type accronym - patient id - N ex. BRA-693-N "
    echo "* datadir: XXX"
    echo "* extdatadir: XXX"
    echo "* condapath: XXX"
    echo "* mode: 'WGS' or 'WES'"
    echo "Example:"
    echo "$ bash $0 XXX XXX XXX XXX"
    exit 1
fi

# Parse input to get initial configuration file
export buffycoatbam=$1
export patientid=$2
export datadir=$3
export extdatadir=$4
export condapath=$5
echo

# activate conda environment
source $condapath
conda activate default

cd ${datadir}/PoNbuffycoat

if [[ $mode == 'WGS' ]] ; then

  for chr in {1..22} ; do
    echo "Processing chromosome " $chr
    samtools view -b ${buffycoatbam} ${chr} > ${datadir}/PoNbuffycoat/PoNbuffycoat_chr${chr}/${patientid}_chr${chr}.bam
    samtools index  ${datadir}/PoNbuffycoat/PoNbuffycoat_chr${chr}/${patientid}_chr${chr}.bam
  done

elif [[ $mode == 'WES' ]] ;

  export datawesdir=$(dirname $datadir)/data_wes
  echo $datawesdir
  for chr in {1..22} ; do
    echo "Processing chromosome " $chr
    samtools view -b -L ${extdatadir}/exome_bed/exome_hg19_chr${chr}.bed \
    ${datadir}/PoNbuffycoat/PoNbuffycoat_chr${chr}/${patientid}_chr${chr}.bam \
     > ${datawesdir}/PoNbuffycoat/PoNbuffycoat_chrall/${patientid}_chr${chr}.bam
    samtools index ${datawesdir}/PoNbuffycoat/PoNbuffycoat_chrall/${patientid}_chr${chr}.bam
  done

  samtools merge ${datawesdir}/PoNbuffycoat/PoNbuffycoat_chrall/${patientid}_chrall.bam \
  ${datawesdir}/PoNbuffycoat/PoNbuffycoat_chrall/${patientid}_chr*.bam
  samtools index  ${datawesdir}/PoNbuffycoat/PoNbuffycoat_chrall/${patientid}_chrall.bam

  if [ -f ${datawesdir}/PoNbuffycoat/PoNbuffycoat_chrall/${patientid}_chrall.bam.bai ]; then
    rm ${datawesdir}/PoNbuffycoat/PoNbuffycoat_chrall/${patientid}_chr[0-9]*.bam
  fi

else
  print "Supported modes should be 'WES' or 'WGS' but here is ${mode}"
fi


