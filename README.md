# cfdna_snv_benchmark

## Requirements

XXX

## Install

XXX
add the repository to your Python path
export PYTHONPATH=$PYTHONPATH:"/Users/hanae/Repositories/cfdna_snv_benchmark/"

## Usage

XXX

### Create benchmark dataset

#### Create mixture series

$ qsub -pe OpenMP 1 -l mem_free=12G,h_rt=12:00:00 -o /mnt/projects/zhug/cfDNA/skandlab-public/carriehc/data/mixtures/mixtures_chr22/logs/ -e /mnt/projects/zhug/cfDNA/skandlab-public/carriehc/data/mixtures/mixtures_chr22/logs/ /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/mixtures/create_mixtures_series_chr.sh -c /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/mixtures/config/config_lowtfsample_1014_180816_chr22.yml
Save to AWS s3 bucket

#### Create spikein series

Prepare common cancer mutations to insert
$ bash create_spikeins_series_chr.sh -c /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/config/config_spikeins/config_lowtfsample_CRC-1014_180816_chr22_snv.yml
$ bash create_spikeins_series_chr.sh -c /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/config/config_spikeins/config_lowtfsample_CRC-1014_180816_chr22_indel.yml
Save to AWS s3 bucket

### Run callers

#### DNA-specific callers using bcbio pipeline

Prepare bcbio input file
Run bcbio pipeline in Aquila
Copy bcbio results to repository 
Save results

### Runs cfDNA-specific callers 

ABEMUS
SiNVICT
cfSNV

