#!/usr/bin/env Rscript

library(cfSNV)
library(yaml)
library(optparse)
 
option_list = list(
  make_option(c("-c", "--config_file"), type = 'character', help="config yaml file path"),
  make_option(c('=t', '--targetbed'), type = 'character', help='path to bedfile listing target positions')); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

config = yaml.load_file(opt$config_file)
targetbed = opt$targetbed

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

print(targetbed)

# read parameter recommended 
MIN_HOLD_SUPPORT_COUNT = as.integer(unlist(strsplit(unlist(strsplit(grep('at 1% VAF: ',readLines(file.path(config$outputdir, 'log.out'), warn=FALSE), value = TRUE), 'MIN_HOLD_SUPPORT_COUNT = '))[2], ","))[1])
MIN_PASS_SUPPORT_COUNT = as.integer(unlist(strsplit(unlist(strsplit(unlist(strsplit(grep('at 1% VAF: ',readLines(file.path(config$outputdir, 'log.out'), warn=FALSE), value = TRUE), 'MIN_HOLD_SUPPORT_COUNT = '))[2], 'MIN_PASS_SUPPORT_COUNT = '))[2], ';'))[1])
print(paste0('MIN_HOLD_SUPPORT_COUNT = ', MIN_HOLD_SUPPORT_COUNT, ', MIN_PASS_SUPPORT_COUNT = ', MIN_PASS_SUPPORT_COUNT))

print(targetbed)
i <- unlist(strsplit(substr(targetbed,1,nchar(targetbed)-4), '_'))
i <- i[length(i)]
results_file <- file.path(config$outputdir, paste0(config$plasma$id, '.results_', i, '.txt'))
print(results_file)

results <- variant_calling(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, targetbed, config$reference, config$SNPdatabase, config$dir$samtools, config$dir$picard, condif$dir$bedtools, config$plasma$id, MIN_HOLD_SUPPORT_COUNT, MIN_PASS_SUPPORT_COUNT, python.dir=config$dir$python)

write.table(results, file = results_file, sep='\t', row.names=F, quote=F)
print(results$variant.list)
print(results$tumor.fraction)
