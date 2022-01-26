library(cfSNV)

# input files
outdir <- "/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/callers/cfSNV"
bamfile <- "/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/dilutions_series_chr22/dilutions_CRC-986_100215/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted.bam"
tmp.dir <- "/mnt/projects/carriehc/cfDNA/anaconda3/envs/cfSNV/lib/R/library/cfSNV/extdata/tmp/"
bam.file <- file.path(tmp.dir, basename(bamfile))
fastq1 <- file.path(tmp.dir, "fastq1.fq")
fastq2 <- file.path(tmp.dir, "fastq2.fq")
sample.id <- "dilution_chr22_CRC-986_100215_1_CRC-986_300316_0"
plasma.unmerged <-  file.path(tmp.dir, paste0(sample.id, '.recal.bam'))
plasma.merged.extendedFrags <- file.path(tmp.dir, paste0(sample.id, ".extendedFrags.recal.bam"))
plasma.merge.notCombined <- file.path(tmp.dir, paste0(sample.id, ".notCombined.recal.bam"))

bamnormal <- "/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/dilutions_series_chr22/NCC_CRC-986_100215-CW-N/CRC-986_100215-BC_chr22.bam"
normal.id <- "CRC-986_100215-BC_chr22"
bam.normal.file <- file.path(tmp.dir, basename(bamnormal))
fastq1.normal <- file.path(tmp.dir, "fastq1_normal.fq")
fastq2.normal <- file.path(tmp.dir, "fastq2_normal.fq")
normal <-  file.path(tmp.dir, paste0(normal.id, ".recal.bam"))

target.bed <- paste0(outdir, '/exome_chr22_hg19.bed')
reference <- "/mnt/projects/carriehc/cfDNA/data/refgenome/GRCh37/GRCh37.fa"
SNP.database <-  paste0(outdir, "/homo_sapiens-chr22.vcf")
samtools.dir <- "/mnt/projects/carriehc/cfDNA/anaconda3/envs/cfSNV/bin/samtools"
picard.dir <- "/mnt/projects/carriehc/cfDNA/utils/picard.jar"
bedtools.dir <- "/mnt/projects/carriehc/cfDNA/anaconda3/envs/cfSNV/bin/bedtools"
GATK.dir <- "/mnt/projects/carriehc/cfDNA/anaconda3/envs/cfSNV/opt/gatk-3.8/GenomeAnalysisTK.jar" #"/mnt/projects/carriehc/cfDNA/anaconda3/envs/cfSNV/share/gatk4-4.2.3.0-0/gatk-package-4.2.3.0-local.jar" #"/mnt/projects/carriehc/cfDNA/anaconda3/envs/cfSNV/bin/gatk3"
bwa.dir <- "/mnt/projects/carriehc/cfDNA/anaconda3/envs/cfSNV/bin/bwa"
flash.dir <- "/mnt/projects/carriehc/cfDNA/anaconda3/envs/cfSNV/bin/flash2"
output.dir <- "/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/callers/cfSNV"


# need to copy bam as issue in bwa mem cannot realign fastqs
#system2(command = "cp", args = paste(bamfile, tmp.dir))
#system2(command = "cp", args = paste(paste0(bamfile, '.bai'), tmp.dir))
#system2(command = bedtools.dir, args = paste("bamtofastq -i", bam.file, "-fq", fastq1, "-fq2", fastq2))
#system2(command = 'rm', args = bam.file)
#system2(command = 'rm', args = paste0(bam.file, '.bai'))
#getbam_align(fastq1, fastq2, reference, SNP.database, samtools.dir, picard.dir, bedtools.dir, GATK.dir, bwa.dir, sample.id, java.dir='java')

#getbam_align_after_merge(fastq1, fastq2, reference, SNP.database, samtools.dir, picard.dir, bedtools.dir, GATK.dir, bwa.dir, flash.dir, sample.id,  java.dir='java')

#system2(command = "cp", args = paste(bamnormal, tmp.dir))
#system2(command = "cp", args = paste(paste0(bamnormal, '.bai'), tmp.dir))
#system2(command = bedtools.dir, args = paste("bamtofastq -i", bam.normal.file, "-fq", fastq1.normal, "-fq2", fastq2.normal))
#getbam_align(fastq1.normal, fastq2.normal, reference, SNP.database, samtools.dir, picard.dir, bedtools.dir, GATK.dir, bwa.dir, normal.id, java.dir='java')

parameter_recommend(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, target.bed, reference, SNP.database, samtools.dir, sample.id, roughly_estimated_tf=TRUE)

#MIN_HOLD_SUPPORT_COUNT = 7
#MIN_PASS_SUPPORT_COUNT = 1

#results <- variant_calling(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, target.bed, reference, SNP.database, samtools.dir,picard.dir, bedtools.dir, sample.id, MIN_HOLD_SUPPORT_COUNT, MIN_PASS_SUPPORT_COUNT)

#print(results$variant.list)

#print(results$tumor.fraction)


