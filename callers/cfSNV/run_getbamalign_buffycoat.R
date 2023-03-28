#!/usr/bin/env Rscript

library(cfSNV,  lib.loc='/rfs-storageservice/GIS/Projects/LOCCG/carriehc/Rlibs/')
library(yaml)
library(optparse)
library(callr)

option_list = list(
  make_option(c("-c", "--config_file"), type='character', help="config yaml file path")); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

config = yaml.load_file(opt$config_file)

# input files
print(config$outdir)

print(config$buffycoatbam)
normalid <- tools::file_path_sans_ext(basename(config$buffycoatbam))
normaldir <- dirname(config$buffycoatbam)
normalfastq1 <- file.path(normaldir, paste0(normalid, '_R1.fastq.gz'))
normalfastq2 <- file.path(normaldir, paste0(normalid, '_R2.fastq.gz'))
print(normalfastq1)
print(normalfastq2)
normaloutputdir <- file.path(config$outdir, normalid)
normal <-  file.path(normaloutputdir, paste0(normalid, ".recal.bam"))
print(normal)

targetbeddir <- file.path(config$extdata, 'exome_bed')
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

if (!file.exists(normal)) {
        print('Get BAM align normal')
        getbam_align(normalfastq1, normalfastq2, reference, SNPdatabase, config$dir$samtools, config$dir$picard, config$dir$bedtools, config$dir$GATK, config$dir$bwa, normalid, normaloutputdir, java.dir=config$dir$java)
	normal.oldname <- list.files(path=normaloutputdir, pattern=paste0(normalid, '[^a-zA-Z]*.recal.bam'), full.names=TRUE)
	print(normal.oldname)
	print(normal)
	file.rename(normal.oldname, normal)
	file.rename(paste0(substring(normal.oldname, 1, nchar(normal.oldname)-1), 'i'), paste0(substring(normal,1, nchar(normal)-1), 'i'))
}


