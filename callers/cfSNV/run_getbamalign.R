#!/usr/bin/env Rscript

library(cfSNV, lib.loc='/outdir/Rlibs')
library(yaml)
library(optparse)
library(callr)

option_list = list(
  make_option(c("-c", "--config_file"), type='character', help="config yaml file path"), 
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
plasmafastq1 <- file.path(config$dilutionseriesfolder, plasmaid, paste0(plasmaid, '_R1.fastq.gz'))
plasmafastq2 <- file.path(config$dilutionseriesfolder, plasmaid, paste0(plasmaid, '_R2.fastq.gz'))
print(plasmafastq1)
print(plasmafastq2)
plasma.unmerged <-  file.path(outputdir, paste0(plasmaid, '.recal.bam'))
plasma.merged.extendedFrags <- file.path(outputdir, paste0(plasmaid, ".extendedFrags.recal.bam"))
plasma.merge.notCombined <- file.path(outputdir, paste0(plasmaid, ".notCombined.recal.bam"))
print(plasma.unmerged)
print(plasma.merged.extendedFrags)
print(plasma.merge.notCombined)

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


if (!file.exists(plasma.unmerged) ) {
        print('Get BAM align plasma')
        getbam_align(plasmafastq1, plasmafastq2, reference, SNPdatabase, config$dir$samtools, config$dir$picard, config$dir$bedtools, config$dir$GATK, config$dir$bwa, plasmaid, outputdir, java.dir=config$dir$java)
	plasma.unmerged.oldname <- list.files(path=outputdir, pattern=paste0(plasmaid, '[^a-zA-Z]*.recal.bam'), full.names=TRUE)
	print(plasma.unmerged.oldname)
	print(plasma.unmerged)
	file.rename(plasma.unmerged.oldname, plasma.unmerged)
	file.rename(paste0(substring(plasma.unmerged.oldname, 1, nchar(plasma.unmerged.oldname)-1), 'i'), paste0(substring(plasma.unmerged,1, nchar(plasma.unmerged)-1), 'i'))

}

if (!file.exists(plasma.merged.extendedFrags) ) {
        print('Get BAM align after merge plasma')
        getbam_align_after_merge(plasmafastq1, plasmafastq2, reference, SNPdatabase, config$dir$samtools, config$dir$picard, config$dir$bedtools, config$dir$GATK, config$dir$bwa, config$dir$flash, plasmaid, outputdir, java.dir=config$dir$java)
	plasma.merge.notCombined.oldname <- list.files(path=outputdir, pattern=paste0(plasmaid, '[^a-zA-Z]*.notCombined.recal.bam'), full.names=TRUE)
        print(plasma.merge.notCombined.oldname)
        print(plasma.merge.notCombined)
        file.rename(plasma.merge.notCombined.oldname, plasma.merge.notCombined)
        file.rename(paste0(substring(plasma.merge.notCombined.oldname, 1, nchar(plasma.merge.notCombined.oldname)-1), 'i'), paste0(substring(plasma.merge.notCombined,1, nchar(plasma.merge.notCombined)-1), 'i'))
	plasma.merged.extendedFrags.oldname <- list.files(path=outputdir, pattern=paste0(plasmaid, '[^a-zA-Z]*.extendedFrags.recal.bam'), full.names=TRUE)
	print(plasma.merged.extendedFrags.oldname)
	print(plasma.merged.extendedFrags)
	file.rename(plasma.merged.extendedFrags.oldname, plasma.merged.extendedFrags)
	file.rename(paste0(substring(plasma.merged.extendedFrags.oldname, 1, nchar(plasma.merged.extendedFrags.oldname)-1), 'i'), paste0(substring(plasma.merged.extendedFrags,1, nchar(plasma.merged.extendedFrags)-1), 'i'))
}

if (!file.exists(normal)) {
        print('Get BAM align normal')
        getbam_align(normalfastq1, normalfastq2, reference, SNPdatabase, config$dir$samtools, config$dir$picard, config$dir$bedtools, config$dir$GATK, config$dir$bwa, normalid, normaloutputdir, java.dir=config$dir$java)
	normal.oldname <- list.files(path=normaloutputdir, pattern=paste0(normalid, '[^a-zA-Z]*.recal.bam'), full.names=TRUE)
	print(normal.oldname)
	print(normal)
	file.rename(normal.oldname, normal)
	file.rename(paste0(substring(normal.oldname, 1, nchar(normal.oldname)-1), 'i'), paste0(substring(normal,1, nchar(normal)-1), 'i'))
}


