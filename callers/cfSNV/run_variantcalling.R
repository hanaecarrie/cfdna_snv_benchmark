#!/usr/bin/env Rscript

library(cfSNV)
library(yaml)
library(optparse)
 
option_list = list(
  make_option(c("-c", "--config_file"), type = 'character', help="config yaml file path"),
  make_option(c('-t', '--targetbed'), type = 'character', help='path to bedfile listing target positions'),
  make_option(c("-p", "--plasmaid"), type='character', help="plasma sample id")); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

config = yaml.load_file(opt$config_file)
targetbed = opt$targetbed
plasmaid = opt$plasmaid

# input files
print(config$outdir)
outputdir <- file.path(config$outdir, plasmaid)
print(outputdir)
plasmafastq1 <- file.path(config$dilutionseriesfolder, plasmaid, paste0(plasmaid, '.sorted_R1.fastq.gz'))
plasmafastq2 <- file.path(config$dilutionseriesfolder, plasmaid, paste0(plasmaid, '.sorted_R2.fastq.gz'))
print(plasmafastq1)
print(plasmafastq2)
plasma.unmerged <-  file.path(outputdir, paste0(plasmaid, '.recal.bam'))
plasma.merged.extendedFrags <- file.path(outputdir, paste0(plasmaid, ".extendedFrags.recal.bam"))
plasma.merge.notCombined <- file.path(outputdir, paste0(plasmaid, ".notCombined.recal.bam"))
print(plasma.unmerged)
print(plasma.merged.extendedFrags)
print(plasma.merge.notCombined)

print(config$buffycoatbam)
normalid <- tools::file_path_sans_ext(config$buffycoatbam)
normaldir <- dirname(config$buffycoatbam)
normalfastq1 <- file.path(normaldir, paste0(normalid, '.sorted_R1.fastq.gz'))
normalfastq2 <- file.path(normaldir, paste0(normalid, '.sorted_R2.fastq.gz'))
print(normalfastq1)
print(normalfastq2)
normal <-  file.path(outputdir, paste0(normalid, ".recal.bam"))
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

# read parameter recommended 
MIN_HOLD_SUPPORT_COUNT = as.integer(unlist(strsplit(unlist(strsplit(grep('at 1% VAF: ',readLines(file.path(outputdir, 'log.out'), warn=FALSE), value = TRUE), 'MIN_HOLD_SUPPORT_COUNT = '))[2], ","))[1])
MIN_PASS_SUPPORT_COUNT = as.integer(unlist(strsplit(unlist(strsplit(unlist(strsplit(grep('at 1% VAF: ',readLines(file.path(outputdir, 'log.out'), warn=FALSE), value = TRUE), 'MIN_HOLD_SUPPORT_COUNT = '))[2], 'MIN_PASS_SUPPORT_COUNT = '))[2], ';'))[1])
print(paste0('MIN_HOLD_SUPPORT_COUNT = ', MIN_HOLD_SUPPORT_COUNT, ', MIN_PASS_SUPPORT_COUNT = ', MIN_PASS_SUPPORT_COUNT))

i <- unlist(strsplit(substr(targetbed,1,nchar(targetbed)-4), '_'))
i <- i[length(i)]
results_file <- file.path(outputdir, paste0(plasmaid, '.results_', i, '.txt'))
print(results_file)

results <- variant_calling(plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined, targetbed, reference, SNPdatabase, config$dir$samtools, config$dir$picard, config$dir$bedtools, plasmaid, MIN_HOLD_SUPPORT_COUNT, MIN_PASS_SUPPORT_COUNT, python.dir=config$dir$python)

write.table(results, file = results_file, sep='\t', row.names=F, quote=F)
print(results$variant.list)
print(results$tumor.fraction)
