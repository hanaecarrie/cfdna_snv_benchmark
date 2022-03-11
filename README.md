# cfdna_snv_benchmark

## Requirements

XXX

## Install

XXX
add the repository to your Python path
```
export PYTHONPATH=$PYTHONPATH:"/Users/hanae/Repositories/cfdna_snv_benchmark/"
```
## Usage

### 1. Create benchmark dataset

```
$ export repodir=/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark

$ export datadir=/mnt/projects/zhug/cfDNA/skandlab-public/carriehc
```

#### Create mixture series

Prepare dilution series ratios
XXX
Create series
```
$ qsub -pe OpenMP 1 -l mem_free=12G,h_rt=12:00:00 -o $datadir/data/mixtures/mixtures_chr22/logs/ -e $datadir/data/mixtures/mixtures_chr22/logs/ $repodir/mixtures/create_mixtures_series_chr.sh -c $repodir/mixtures/config/config_lowtfsample_1014_180816_chr22.yml
```
Save to AWS s3 bucket
XXX

#### Create spikein series

Prepare common cancer mutations to insert
XXX
Create series
```
$ bash $repodir/spikeins/create_spikeins_series_chr.sh -c $repodir/spikeins/config/config_lowtfsample_CRC-1014_180816_chr22.yml
```
Save to AWS s3 bucket
XXX

### 2. Run callers

#### DNA-specific callers using bcbio pipeline

Prepare bcbio input file
```
$ bash $repodir/bcbio/run_prepare_bcbio.sh
```
Run bcbio pipeline in Aquila
```
$ bash $repodir/bcbio/run_bcbio.sh
```
Copy bcbio results to repository
XXX
Save results
XXX

### Runs cfDNA-specific callers 

Transfer to Ronin machine
XXX
Run ABEMUS
```
$ screen -S abemus

$ bash ~/ABEMUS/run_abemus.sh -c ~/ABEMUS/config/config_mixtures_chr22_CRC-1014_180816-CW-T.yaml
```
Run SiNVICT
```
$ screen -S sinvict

$ bash ~/sinvict/run_sinvict.sh -c ~/sinvict/config/config_mixtures_chr22_CRC-1014_180816-CW-T.yaml
```
Run cfSNV
```
$ screen -S cfsnv

$ bash ~/cfSNV/run_cfsnv.sh -c ~/cfSNV/config/config_mixtures_chr22_CRC-1014_180816-CW-T.yaml
```


### 3. Analyse results
XXX
