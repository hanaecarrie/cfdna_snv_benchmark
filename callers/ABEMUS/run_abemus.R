library(abemus)
library(data.table)
library(parallel)


args <- commandArgs(TRUE)
outdir <- args[1]   # "/data/abemus_outdir/"
sample.info.file <- args[2]   # "/home/ubuntu/ABEMUS/infofile_test.tsv"
targetbed <- args[3]  # "/data/chr22_data/wholegenome_chr22_hg19_00.bed"
pacbamfolder <- args[4]  # "/data/abemus_outdir/PaCBAM_outdir/"
pacbamfolder_bychrom <- args[5]  # "/data/abemus_outdir/pacbam_data_bychrom/"

print(outdir)
print(sample.info.file)
print(targetbed)
print(pacbamfolder)
print(pacbamfolder_bychrom)

setwd(outdir)

# split pacbam outputs by chrom

split_pacbam_bychrom(targetbed=targetbed, pacbamfolder=pacbamfolder, pacbamfolder_bychrom=pacbamfolder_bychrom)
targetbp_list <- bed2positions(targetbed=targetbed, get_only_chromosomes=TRUE)

# compute per-base error model

outpbem <- compute_pbem(sample.info.file=sample.info.file, targetbed=targetbed, outdir=outdir, pacbamfolder_bychrom=pacbamfolder_bychrom)

# outs

head(outpbem$pbem_tab)
outpbem$bperr_summary
outpbem$bgpbem
outpbem$mean_pbem

# compute coverage-based and not-coverage-based allelic fraction thresholds

outafth <- compute_afthreshold(outdir=outdir, pbem_dir=file.path(outdir,"BaseErrorModel"))

# outs

head( outafth$th_results )
head( outafth$th_results_bin )
head( outafth$datacount_bin )

# call snvs in case samples

calls <- callsnvs(sample.info.file=sample.info.file, controls_dir=file.path(outdir,"Controls"), pbem_dir=file.path(outdir,"BaseErrorModel"), detection.specificity=0.995, replicas=10, outdir=outdir, outdir.calls.name="Results", targetbed=targetbed, pacbamfolder_bychrom=pacbamfolder_bychrom)
head(calls$tabsnvs_index)

tabindex <- calls$tabsnvs_index

calls$tabsnvs_index_scalfact <- apply_scaling_factor(tabindex = tabindex, R = 1.2 )

# compute mean coverage

tabindex <- get_case_mean_coverage(tabindex = tabindex, pacbamfolder_bychrom = pacbamfolder_bychrom)
print(tabindex)
target_size <- get_target_size(targetbed=targetbed, Mbp = TRUE)
print(target_size)

# apply pbem scaling factor to calls

calls$tabsnvs_index_scalfact <- apply_scaling_factor(tabindex=tabindex, target_size=target_size, use.optimal.R=TRUE)

head(calls$tabsnvs_index_scalfact) 

