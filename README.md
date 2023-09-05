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

Define the directory where the code repository is installed and the location of the dataset.
```
export repodir=/home/users/astar/gis/carriehc/cfdna_snv_benchmark
export datadir=/rfs-storageservice/GIS/Projects/LOCCG/carriehc
```

For deep Whole Genome Sequencing (WGS) mixture series, 
the samples are splitted by chromosomes due to storage space and memory limitations. \
For deep Whole Exome Sequencing (WES) mixture series, all the chromsomes are processed together.

#### 1. Prepare configuration files

Adapt the template '${repodir}/mixtures/config/config_template.yaml' to create your own.\
Here, the configuration file example is 'config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml' for a WGS mixture series. \
```
cd ${repodir}/mixtures
cp config/config_template.yml config/config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml
vi config/config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml
```
In the WGS case, write config file for chr1 and run the following to create the following config files for chr2-22: 
```
cd ${repodir}/mixtures
bash create_config_yaml.sh config/config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml
```

#### 2. Run mixture series creations in HPC as jobs (in Slurm)

Example for a WGS series splitted by chromosome. Here for chr1.
```
sbatch -p normal -J chr1_986_mixtures -t 24:00:00 --mem 48000 \
 --output=${datadir}/data/logs/chr1_986_mixture_wgs.o \
 --error=${datadir}/data/logs/chr1_986_mixtures_wgs.e \
 ${repodir}/cfdna_snv_benchmark/mixtures/create_mixtures_series_chr.sh \
  -c ${repodir}/mixtures/config/config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml
```

#### 3. Alternatively, download paper dataset from s3 bucket using the aws cli API.

Configure first the ~/.aws/credentials file to add the READ-ONLY access key of the database s3 bucket.
Here is displayed the database tree corresponding to the output of the mixture series creation detailed above.
A total of 8 samples were generated with different tumor fractions and restricted to chr1. 
For each mixture sample, the index bam file is present with the log files as well as the read group file, coverage, roughly estimated TF, and ichorCNA TF estimate.
```
aws s3 ls s3://cfdna-benchmark-dataset.store.genome.sg/mixtures/mixtures_chr1/mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T/ --profile [profilenickname] 
                           PRE mixture_chr1_CRC-986_100215-CW-T_10x_CRC-986_300316-CW-T_140x/
                           PRE mixture_chr1_CRC-986_100215-CW-T_20x_CRC-986_300316-CW-T_130x/
                           PRE mixture_chr1_CRC-986_100215-CW-T_30x_CRC-986_300316-CW-T_120x/
                           PRE mixture_chr1_CRC-986_100215-CW-T_50x_CRC-986_300316-CW-T_100x/
                           PRE mixture_chr1_CRC-986_100215-CW-T_5x_CRC-986_300316-CW-T_145x/
                           PRE mixture_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_0x/
                           PRE mixture_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_180x/
                           PRE mixture_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_80x/
2022-05-08 16:35:40       1301 config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml

aws s3 ls s3://cfdna-benchmark-dataset.store.genome.sg/mixtures/mixtures_chr1/mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T/mixture_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_80x/ --profile [profilenickname]  
                           PRE calls/
                           PRE ichorcna/
                           PRE logs/
2022-05-08 17:03:49          8 coverage_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_80x.txt
2022-05-08 17:03:51         22 estimated_tf_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_80x.txt
2022-05-08 17:04:03 23180309957 mixture_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_80x.bam
2022-05-08 17:12:59     825776 mixture_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_80x.bam.bai
2022-05-08 17:12:59         64 rg_chr1_CRC-986_100215-CW-T_70x_CRC-986_300316-CW-T_80x.txt

```


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
