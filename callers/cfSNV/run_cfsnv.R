#!/usr/bin/env Rscript

library(cfSNV)
library(yaml)
library(optparse)
library(callr)

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
plasma.merge.notCombined <- file.path(config$outputdir, paste0(config$plasma$id, ".notCombined.recal.bam"))

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
	plasma.unmerged.oldname <- list.files(path=config$outputdir, pattern=paste0(config$plasma$id, '[^a-zA-Z]*.recal.bam'), full.names=TRUE)
	print(plasma.unmerged.oldname)
	print(plasma.unmerged)
	file.rename(plasma.unmerged.oldname, plasma.unmerged)
	file.rename(paste0(substring(plasma.unmerged.oldname, 1, nchar(plasma.unmerged.oldname)-1), 'i'), paste0(substring(plasma.unmerged,1, nchar(plasma.unmerged)-1), 'i'))
}

if (!file.exists(plasma.merged.extendedFrags) ) {
        print('Get BAM align after merge plasma')
        getbam_align_after_merge(config$plasma$fastq1, config$plasma$fastq2, config$reference, config$SNPdatabase, config$dir$samtools, config$dir$picard, config$dir$bedtools, config$dir$GATK, config$dir$bwa, config$dir$flash, config$plasma$id, config$outputdir, java.dir=config$dir$java)
	plasma.merged.extendedFrags.oldname <- list.files(path=config$outputdir, pattern=paste0(config$plasma$id, '[^a-zA-Z]*.extendedFrags.recal.bam'), full.names=TRUE)
	print(plasma.merged.extendedFrags.oldname)
	print(plasma.merged.extendedFrags)
	file.rename(plasma.merged.extendedFrags.oldname, plasma.merged.extendedFrags)
	file.rename(paste0(substring(plasma.merged.extendedFrags.oldname, 1, nchar(plasma.merged.extendedFrags.oldname)-1), 'i'), paste0(substring(plasma.merged.extendedFrags,1, nchar(plasma.merged.extendedFrags)-1), 'i'))
	plasma.merge.notCombined.oldname <- list.files(path=config$outputdir, pattern=paste0(config$plasma$id, '[^a-zA-Z]*.notCombined.recal.bam'), full.names=TRUE)
	print(plasma.merge.notCombined.oldname)
	print(plasma.merge.notCombined)
	file.rename(plasma.merge.notCombined.oldname, plasma.merge.notCombined)
	file.rename(paste0(substring(plasma.merge.notCombined.oldname, 1, nchar(plasma.merge.notCombined.oldname)-1), 'i'), paste0(substring(plasma.merge.notCombined,1, nchar(plasma.merge.notCombined)-1), 'i'))
}

if (!file.exists(normal)) {
        print('Get BAM align normal')
        getbam_align(config$normal$fastq1, config$normal$fastq2, config$reference, config$SNPdatabase, config$dir$samtools, config$dir$picard, config$dir$bedtools, config$dir$GATK, config$dir$bwa, config$normal$id, config$outputdir, java.dir=config$dir$java)
	normal.oldname <- list.files(path=config$outputdir, pattern=paste0(config$normal$id, '[^a-zA-Z]*.recal.bam'), full.names=TRUE)
	print(normal.oldname)
	print(normal)
	file.rename(normal.oldname, normal)
	file.rename(paste0(substring(normal.oldname, 1, nchar(normal.oldname)-1), 'i'), paste0(substring(normal,1, nchar(normal)-1), 'i'))
}


print('Parameter recommend')
targetbedfull <- list.files(path = config$targetbeddir, pattern = 'chr[0-2]?[0-9].bed', all.files = TRUE, full.names = TRUE) # target bed directory with whole genome bed XXX_chr22.bed and split beds per chunks XXX_chr22_i.bed, i from 00 to 99 for instance
parameter_recommend(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, targetbedfull, config$reference, config$SNPdatabase, config$dir$samtools, config$plasma$id, roughly_estimated_tf=TRUE, python.dir=config$dir$python)
# read parameter recommended 
MIN_HOLD_SUPPORT_COUNT = as.integer(unlist(strsplit(unlist(strsplit(grep('at 1% VAF: ',readLines(file.path(config$outputdir, 'log.out'), warn=FALSE), value = TRUE), 'MIN_HOLD_SUPPORT_COUNT = '))[2], ","))[1])
MIN_PASS_SUPPORT_COUNT = as.integer(unlist(strsplit(unlist(strsplit(unlist(strsplit(grep('at 1% VAF: ',readLines(file.path(config$outputdir, 'log.out'), warn=FALSE), value = TRUE), 'MIN_HOLD_SUPPORT_COUNT = '))[2], 'MIN_PASS_SUPPORT_COUNT = '))[2], ';'))[1])
print(paste0('MIN_HOLD_SUPPORT_COUNT = ', MIN_HOLD_SUPPORT_COUNT, ', MIN_PASS_SUPPORT_COUNT = ', MIN_PASS_SUPPORT_COUNT))

"
listtargetbed <- list.files(path = config$targetbeddir, pattern = 'chr[0-2]?[0-9]_[0-9][0-9].bed', all.files = TRUE, full.names = TRUE)
print(listtargetbed)

for ( targetbed in listtargetbed ) {

        print(targetbed)
        i <- unlist(strsplit(substr(targetbed,1,nchar(targetbed)-4), '_'))
        i <- i[length(i)]
	#rscript('/home/ubuntu/cfSNV/run_variantcalling.R', cmdargs = c('--config_file', opt$config_file, '--targetbed', targetbed))
	options <- rscript_process_options(script = '/home/ubuntu/cfSNV/run_variantcalling.R', cmdargs = c('--config_file', opt$config_file, '--targetbed', targetbed))
	rp <- rscript_process$new(options)
	rp$wait()
	rp$read_output_lines()
}
"

"
targetbed <- paste0(config$targetbeddir, '/wholegenome_hg19_chr22_05.bed')
print(targetbed)
results <- variant_calling(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, targetbed, config$reference, config$SNPdatabase, config$dir$samtools, config$dir$picard, condif$dir$bedtools, config$plasma$id, MIN_HOLD_SUPPORT_COUNT, MIN_PASS_SUPPORT_COUNT, python.dir=config$dir$python)
write.table(results, file = results_file, sep='\t', row.names=F, quote=F)
print(results$variant.list)
print(results$tumor.fraction)
"
"
listtargetbed <- list.files(path = config$targetbeddir, pattern = 'chr[0-2]?[0-9]_[0-9][0-9].bed', all.files = TRUE, full.names = TRUE)
print(listtargetbed)

for ( targetbed in listtargetbed ) {

	print(targetbed)
	i <- unlist(strsplit(substr(targetbed,1,nchar(targetbed)-4), '_'))
	i <- i[length(i)]
	results_file <- file.path(config$outputdir, paste0(config$plasma$id, '.results_', i, '.txt'))
	print(results_file)
	results <- variant_calling(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, targetbed, config$reference, config$SNPdatabase, config$dir$samtools, config$dir$picard, condif$dir$bedtools, config$plasma$id, MIN_HOLD_SUPPORT_COUNT, MIN_PASS_SUPPORT_COUNT, python.dir=config$dir$python)
	write.table(results, file = results_file, sep='\t', row.names=F, quote=F)
	print(results$variant.list)
	print(results$tumor.fraction)

}
"
