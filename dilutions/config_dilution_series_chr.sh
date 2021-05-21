##### Config file #####

samples:
  healthy: /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/healthy/healthy/healthy_chr22_merged-ready.bam 
  tumor: "/mnt/projects/skanderupamj/wgs/data/training/ready.bams/cfdna/NCC_CRC-809_110914-CW-T-ready.bam /mnt/projects/skanderupamj/wgs/data/training/ready.bams/cfdna/NCC_CRC-986_100215-CW-T-ready.bam"
  buffycoat: "/mnt/projects/skanderupamj/wgs/data/training/ready.bams/cfdna/NCC_CRC-809_110914-CW-N-ready.bam /mnt/projects/skanderupamj/wgs/data/training/ready.bams/cfdna/NCC_CRC-986_100215-CW-N-ready.bam"

samplenames:
  healthy: healthy
  tumor: "CRC-809_110914 CRC-986_100215"
  buffycoat: "CRC-809_110914-BC CRC-986_100215-BC"

dilutionfactors: "1-0 1-0.72 0.75-0.765 0.5-0.81 0.25-0.855 0.125-0.875"

chr: 22

healthydir: /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/dilutions_chr22/healthymerged
tumordirs: "/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/dilutions_chr22/NCC_CRC-809_110914-CW-T-ready /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/dilutions_chr22/NCC_CRC-986_100215-CW-T-ready"
buffycoatdirs: "/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/dilutions_chr22/NCC_CRC-809_110914-CW-N-ready /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/dilutions_chr22/NCC_CRC-986_100215-CW-N-ready"
outputfolder: /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/dilutions_series_chr22

vcffile: /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/PATIENTID-ensemble-annotated.vcf
