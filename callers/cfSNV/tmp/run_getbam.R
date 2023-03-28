library(cfSNV)


library(yaml)
config = yaml.load_file("config.yml")

# input files
print(config$outputdir)
print(config$plasma$id)
print(config$plasma$fastq1, config$plasma$fastq2)
plasma.unmerged <-  file.path(config$outputdir, paste0(config$plasma$id, '.recal.bam'))
plasma.merged.extendedFrags <- file.path(config$outputdir, paste0(config$plasma$id, ".extendedFrags.recal.bam"))
plasma.merge.notCombined <- file.path(config$outputdir, paste0(config$plasam$id, ".notCombined.recal.bam"))

print(config$normal$id)
print(config$normal$fastq1, config$normal$fastq2)
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



