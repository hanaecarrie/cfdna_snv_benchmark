library(cfSNV)

# input files
tmp.dir <- "/home/ubuntu/anaconda3/envs/cfSNV/lib/R/library/cfSNV/extdata/tmp/"
output.dir <- "/home/ubuntu/cfSNV/results/"
fastq1 <- '/home/ubuntu/test_data/authors_data_cfSNV/P20_T4_cfDNA.chr22.paired.1.sorted.fastq.gz'
fastq2 <- '/home/ubuntu/test_data/authors_data_cfSNV/P20_T4_cfDNA.chr22.paired.2.sorted.fastq.gz'
sample.id <- 'plasma' #P20_T4_cfDNA.chr22.paired'
plasma.unmerged <-  file.path(output.dir, paste0(sample.id, '.recal.bam'))
print(sample.id)
plasma.merged.extendedFrags <- file.path(output.dir, paste0(sample.id, ".extendedFrags.recal.bam"))
plasma.merge.notCombined <- file.path(output.dir, paste0(sample.id, ".notCombined.recal.bam"))

normal.id <- 'normal' #'P20_WBC.chr22.paired'
fastq1.normal <- '/home/ubuntu/test_data/authors_data_cfSNV/P20_WBC.chr22.paired.1.sorted.fastq.gz' 
fastq2.normal <- '/home/ubuntu/test_data/authors_data_cfSNV/P20_WBC.chr22.paired.2.sorted.fastq.gz'
normal <-  file.path(output.dir, paste0(normal.id, ".recal.bam"))

target.bed <- '/home/ubuntu/test_data/authors_data_cfSNV/target_hg19_chr22.bed' #'/home/ubuntu/data/exome_chr22_hg19.bed' 
reference <- '/home/ubuntu/test_data/authors_data_cfSNV/demo/chr22.fa' #'/home/ubuntu/data/GRCh37/GRCh37.fa'
SNP.database <-  '/home/ubuntu/test_data/authors_data_cfSNV/demo/chr22_dbSNP.vcf' #'/home/ubuntu/data/homo_sapiens-chr22.vcf'
samtools.dir <- "/home/ubuntu/anaconda3/envs/cfSNV/bin/samtools"
picard.dir <- "/home/ubuntu/anaconda3/envs/cfSNV/share/picard-2.26.10-0/picard.jar"
bedtools.dir <- "/home/ubuntu/anaconda3/envs/cfSNV/bin/bedtools"
GATK.dir <- "/home/ubuntu/anaconda3/envs/cfSNV/opt/gatk-3.8/GenomeAnalysisTK.jar"
bwa.dir <- "/home/ubuntu/anaconda3/envs/cfSNV/bin/bwa"
flash.dir <- "/home/ubuntu/anaconda3/envs/cfSNV/bin/flash2"
results_file <- file.path(output.dir, "test.results.txt")

#getbam_align(fastq1, fastq2, reference, SNP.database, samtools.dir, picard.dir, bedtools.dir, GATK.dir, bwa.dir, sample.id, output.dir, java.dir='java')

#getbam_align_after_merge(fastq1, fastq2, reference, SNP.database, samtools.dir, picard.dir, bedtools.dir, GATK.dir, bwa.dir, flash.dir, sample.id, output.dir, java.dir='java')

#getbam_align(fastq1.normal, fastq2.normal, reference, SNP.database, samtools.dir, picard.dir, bedtools.dir, GATK.dir, bwa.dir, normal.id, output.dir, java.dir='java')

#parameter_recommend(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, target.bed, reference, SNP.database, samtools.dir, sample.id, roughly_estimated_tf=TRUE, python.dir='/home/ubuntu/anaconda3/envs/cfSNV/bin/python')

#The per base coverage of the plasma sample for each genomic region in the target bed file:
#average = 149.833, median = 139.046, 95th percentile = 280.356

#The roughly estimated tumor fraction in the plasma sample: 10.025%
#For a more accurate estimation, please run variant_calling().

#Lowest detectable VAF range under the default parameters: [1.783%, 4.28%]

#To detect different levels of lowest VAF,
#at 1% VAF: MIN_HOLD_SUPPORT_COUNT = 9, MIN_PASS_SUPPORT_COUNT = 3;
#at 5% VAF: MIN_HOLD_SUPPORT_COUNT = 20, MIN_PASS_SUPPORT_COUNT = 14
#Note: decreasing the parameters (i.e. MIN_HOLD_SUPPORT_COUNT and MIN_PASS_SUPPORT_COUNT)
#can lower the detection limit, but may also lower the variant quality.

MIN_HOLD_SUPPORT_COUNT = 20
MIN_PASS_SUPPORT_COUNT = 14

results <- variant_calling(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, target.bed, reference, SNP.database, samtools.dir,picard.dir, bedtools.dir, sample.id, MIN_HOLD_SUPPORT_COUNT, MIN_PASS_SUPPORT_COUNT, python.dir='/home/ubuntu/anaconda3/envs/cfSNV/bin/python')

write.table(results, file = results_file, sep='\t', row.names=F, quote=F)

print(results$variant.list)

print(results$tumor.fraction)


