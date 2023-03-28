#!/bin/bash

cd /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/pooledhealthy

export chr=$1
echo $chr

ls /mnt/projects/zhug/cfDNA/cohorts-validate/data-c2i-healthy/fastq-MUX11652/WH*/WH*.sort.bam

#for file in /mnt/projects/zhug/cfDNA/cohorts-validate/data-c2i-healthy/fastq-MUX11652/WH*/WH*.sort.bam ;
#do echo $file ; echo  $(dirname $file)_${chr}.sort.bam ;
#/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools view -b $file $chr > $(basename $file .sort.bam)_${chr}.sort.bam ; done
 
ls WH*.bam

/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools merge pooledhealthy_${chr}.bam WH*.bam 
/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools sort pooledhealthy_${chr}.bam > pooledhealthy_${chr}.sort.bam
/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools index pooledhealthy_${chr}.sort.bam
#rm  pooledhealthy_${chr}.bam
