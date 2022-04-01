#!/usr/bin/env Rscript

library(cfSNV, lib.loc='/outdir/Rlibs')
library(yaml)
library(optparse)
 
option_list = list(
  make_option(c("-c", "--config_file"), type = 'character', help="config yaml file path"),
  make_option(c("-p", "--plasmaid"), type='character', help="plasma sample id")); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

config = yaml.load_file(opt$config_file)
plasmaid = opt$plasmaid
print(plasmaid)

# input files
print(config$outdir)
outputdir <- file.path(config$outdir, plasmaid)
print(outputdir)
plasma.unmerged <-  file.path(outputdir, paste0(plasmaid, '.recal.bam'))
plasma.merged.extendedFrags <- file.path(outputdir, paste0(plasmaid, ".extendedFrags.recal.bam"))
plasma.merge.notCombined <- file.path(outputdir, paste0(plasmaid, ".notCombined.recal.bam"))
print(plasma.unmerged)
print(plasma.merged.extendedFrags)
print(plasma.merge.notCombined)

print(config$buffycoatbam)
normalid <- tools::file_path_sans_ext(basename(config$buffycoatbam))
normaldir <- dirname(config$buffycoatbam)
normaloutputdir <- file.path(config$outdir, normalid)
normal <-  file.path(normaloutputdir, paste0(normalid, ".recal.bam"))
print(normal)

targetbeddir <- file.path(config$extdata, 'wholegenome_bed')
reference <- file.path(config$extdat, 'GRCh37', 'GRCh37.fa')
SNPdatabase <- file.path(config$extdata, 'dbsnp_vcf', paste0('dbSNP_hg19_chr', config$chr, '.vcf'))
print(targetbeddir)
print(reference)
print(SNPdatabase)

print(config$dir$samtools)
print(config$dir$picard)
print(config$dir$bedtools)
print(config$dir$GATK)
print(config$dir$bwa)
print(config$dir$flash)
print(config$dir$python)
print(config$dir$java)

print('Parameter recommend')
targetbedfull <- file.path(targetbeddir, paste0('wholegenome_hg19_chr', config$chr, '.bed'))
parameter_recommend(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, targetbedfull, reference, SNPdatabase, config$dir$samtools, plasmaid, roughly_estimated_tf=TRUE, python.dir=config$dir$python)
# read parameter recommended 
MIN_HOLD_SUPPORT_COUNT = as.integer(unlist(strsplit(unlist(strsplit(grep('at 1% VAF: ',readLines(file.path(outputdir, 'log.out'), warn=FALSE), value = TRUE), 'MIN_HOLD_SUPPORT_COUNT = '))[2], ","))[1])
MIN_PASS_SUPPORT_COUNT = as.integer(unlist(strsplit(unlist(strsplit(unlist(strsplit(grep('at 1% VAF: ',readLines(file.path(outputdir, 'log.out'), warn=FALSE), value = TRUE), 'MIN_HOLD_SUPPORT_COUNT = '))[2], 'MIN_PASS_SUPPORT_COUNT = '))[2], ';'))[1])
print(paste0('MIN_HOLD_SUPPORT_COUNT = ', MIN_HOLD_SUPPORT_COUNT, ', MIN_PASS_SUPPORT_COUNT = ', MIN_PASS_SUPPORT_COUNT))


