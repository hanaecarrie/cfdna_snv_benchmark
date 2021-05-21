
#!/bin/bash
cd  /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/
/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/bwa mem -t 10 /mnt/projects/zhug/cfDNA/hg19-index/index /mnt/projects/skanderupamj/wgs/data/cfdna.crc/mix.samples/without.gatk.hap/qualfqs/healthy_chr22_merged_R1.fq.gz /mnt/projects/skanderupamj/wgs/data/cfdna.crc/mix.samples/without.gatk.hap/qualfqs/healthy_chr22_merged_R2.fq.gz  | /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools sort -@ 10 -o  healthy_chr22_merged.bam -
/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/bamsormadup level=-1 threads=16 indexfilename=healthy_chr22_merged.bam.bai tmpfile=healthy_chr22_merged_tmp  < healthy_chr22_merged.bam   > healthy_chr22_merged_rc.bam
