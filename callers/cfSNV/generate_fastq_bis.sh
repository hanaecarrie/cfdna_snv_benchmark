#!/bin/bash

source /apps/anaconda3-individual-edition/2020.11/etc/profile.d/conda.sh
conda activate default

export dilutionseriesfolder=$1
export normal=$2
export config_file=/home/users/astar/gis/carriehc/cfdna_snv_benchmark/callers/cfSNV/config_generate_fastq.yaml

for plasma in ${dilutionseriesfolder}/*/*[Tx].bam ; do

echo $plasma

export plasmafastqoutdir=$(dirname $plasma)
export plasmafastqoutdir=$(echo "${plasmafastqoutdir/\/rfs-storageservice\/GIS\/Projects\/LOCCG\/carriehc/\/home\/users\/astar\/gis\/carriehc\/scratch\/fastq}")  
echo $plasmafastqoutdir
if [ ! -d $plasmafastqoutdir ] ; then mkdir -p $plasmafastqoutdir ; fi
if [ ! -f $plasmafastqoutdir/$(basename $plasma .bam)_R1.fastq.gz ] ; then
	python /home/users/astar/gis/carriehc/cfdna_snv_benchmark/callers/cfSNV/generate_fastqs_yaml.py  -i $plasma -t $plasmafastqoutdir/tmp -o $plasmafastqoutdir -c $config_file  
fi

done

export normalfastqoutdir=$(dirname $normal)
export normalfastqoutdir=$(echo "${normalfastqoutdir/\/rfs-storageservice\/GIS\/Projects\/LOCCG\/carriehc/\/home\/users\/astar\/gis\/carriehc\/scratch\/fastq}") 
echo $normalfastqoutdir
if [ ! -d $normalfastqoutdir ] ; then mkdir -p $normalfastqoutdir ; fi
if [ ! -f $normalfastqoutdir/$(basename $normal .bam)_R1.fastq.gz ] ; then
	python /home/users/astar/gis/carriehc/cfdna_snv_benchmark/callers/cfSNV/generate_fastqs_yaml.py -i $normal -t $normalfastqoutdir/tmp -o $normalfastqoutdir -c $config_file 
fi
