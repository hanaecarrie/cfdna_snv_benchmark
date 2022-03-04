#!/usr/bin/env Rscript

library(cfSNV)
library(yaml)
library(optparse)
 
option_list = list(
  make_option(c("-c", "--config_file"), type='character', help="config yaml file path")); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

config = yaml.load_file(opt$config_file)

# input files
print(config$outputdir)
print(config$plasma$id)
print(config$plasma$fastq1)
print(config$plasma$fastq2)
plasma.unmerged <-  file.path(config$outputdir, paste0(config$plasma$id, '.recal.bam'))
plasma.merged.extendedFrags <- file.path(config$outputdir, paste0(config$plasma$id, ".extendedFrags.recal.bam"))
plasma.merge.notCombined <- file.path(config$outputdir, paste0(config$plasam$id, ".notCombined.recal.bam"))

print(config$normal$id)
print(config$normal$fastq1)
print(config$normal$fastq2)
normal <-  file.path(config$outputdir, paste0(config$normal$id, ".recal.bam"))

print(config$targetbeddir)
print(config$reference)
print(config$SNPdatabase)

print(config$dir$samtools)
print(config$dir$picard)
print(config$dir$bedtools)
print(config$dir$GATK)
print(config$dir$bwa)
print(config$dir$flash)
print(config$dir$python)
print(config$dir$java)

if (!file.exists(plasma.unmerged) ) {
        print('Get BAM align plasma')
        getbam_align(config$plasma$fastq1, config$plasma$fastq2, config$reference, config$SNPdatabase, config$dir$samtools, config$dir$picard, config$dir$bedtools, config$dir$GATK, config$dir$bwa, config$plasma$id, config$outputdir, java.dir=config$dir$java)
}

if (!file.exists(plasma.merged.extendedFrags) ) {
        print('Get BAM align after merge plasma')
        getbam_align_after_merge(config$plasma$fastq1, config$plasma$fastq2, config$reference, config$SNPdatabase, config$dir$samtools, config$dir$picard, config$dir$bedtools, config$dir$GATK, config$dir$bwa, config$dir$flash, config$plasma$id, config$outputdir, java.dir=config$dir$java)
}

if (!file.exists(normal)) {
        print('Get BAM align normal')
        getbam_align(config$normal$fastq1, config$normal$fastq2, config$reference, config$SNPdatabase, config$dir$samtools, config$dir$picard, config$dir$bedtools, config$dir$GATK, config$dir$bwa, config$plasma$id, config$outputdir, java.dir=config$dir$java)
}

print('Parameter recommend')
targetbedfull <- list.files(path = config$targetbeddir, pattern = 'chr[0-2]?[0-9].bed', all.files = TRUE, full.names = TRUE) # target bed directory with whole genome bed XXX_chr22.bed and split beds per chunks XXX_chr22_i.bed, i from 00 to 99 for instance
parameter_recommend(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, targetbedfull, config$reference, config$SNPdatabase, config$dir$samtools, config$plasma$id, roughly_estimated_tf=TRUE, python.dir=config$dir$python)
# read parameter recommended 
MIN_HOLD_SUPPORT_COUNT = as.integer(unlist(strsplit(unlist(strsplit(grep('at 1% VAF: ',readLines(file.path(config$outputdir, 'log.out'), value = TRUE), 'MIN_HOLD_SUPPORT_COUNT = '))[2], ','))[1]))
MIN_PASS_SUPPORT_COUNT = as.integer(unlist(strsplit(unlist(strsplit(grep('at 1% VAF: ',readLines(file.path(config$outputdir, 'log.out'), value = TRUE), 'MIN_PASS_SUPPORT_COUNT = '))[2], ';'))[1]))
print('MIN_HOLD_SUPPORT_COUNT', MIN_HOLD_SUPPORT_COUNT, 'MIN_PASS_SUPPORT_COUNT', MIN_PASS_SUPPORT_COUNT)

listtargetbed <- list.files(path = config$targetbeddir, pattern = 'chr[0-2]?[0-9]_[0-9][0-9].bed', all.files = TRUE, full.names = TRUE)
print(listtargetbed)

for targetbed in (listtargetbed) {

	print(targetbed)
	results_file <- file.path(config$outputdir, paste0(config$plasma$id, ".results_", str(i), ".txt"))
	print(results_file)
	results <- variant_calling(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, targetbed, reference, SNP.database, samtools.dir,picard.dir, bedtools.dir, sample.id, MIN_HOLD_SUPPORT_COUNT, MIN_PASS_SUPPORT_COUNT, python.dir=config$dir$python)
	write.table(results, file = results_file, sep='\t', row.names=F, quote=F)
	print(results$variant.list)
	print(results$tumor.fraction)

}


