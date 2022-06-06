#!/bin/bash

source /home/ubuntu/anaconda3/etc/profile.d/conda.sh
conda activate cfSNV


export dilutionseriesfolder=$1
export normal=$2
export config_file=/home/ubuntu/cfdna_snv_benchmark/callers/cfSNV/config_generate_fastq.yaml

for plasma in ${dilutionseriesfolder}/*/*[Tx].bam ; do

echo $plasma

export plasmafastqoutdir=$(dirname $plasma)
export plasmafastqoutdir=$(echo "${plasmafastqoutdir/mnt\/cfdnaseries/fastq}")  
echo $plasmafastqoutdir
if [ ! -d $plasmafastqoutdir ] ; then mkdir -p $plasmafastqoutdir ; fi
if [ ! -f $plasmafastqoutdir/$(basename $plasma .bam)_R1.fastq.gz ] ; then
	python /home/ubuntu/cfdna_snv_benchmark/callers/cfSNV/generate_fastqs_yaml.py  -i $plasma -t $plasmafastqoutdir/tmp -o $plasmafastqoutdir -c $config_file  
fi

done

export normalfastqoutdir=$(dirname $normal)
export normalfastqoutdir=$(echo "${normalfastqoutdir/mnt\/cfdnaseries/fastq}") 
echo $normalfastqoutdir
if [ ! -d $normalfastqoutdir ] ; then mkdir -p $normalfastqoutdir ; fi
if [ ! -f $normalfastqoutdir/$(basename $normal .bam)_R1.fastq.gz ] ; then
	python /home/ubuntu/cfdna_snv_benchmark/callers/cfSNV/generate_fastqs_yaml.py -i $normal -t $normalfastqoutdir/tmp -o $normalfastqoutdir -c $config_file 
fi
