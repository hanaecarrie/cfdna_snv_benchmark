# cfdna_snv_benchmark

## Requirements

XXX

## Install

XXX
add the repository to your Python path
```sh
export PYTHONPATH=$PYTHONPATH:"/Users/hanae/Repositories/cfdna_snv_benchmark/"
```

## Overview [figure 1]

![image](figures/figure1/Figures_1.png)

## Usage

#### 0. Create tree

Define the directory where the code repository is installed and the location of the dataset.
```sh
export repodir=/home/users/astar/gis/carriehc/cfdna_snv_benchmark
export datadir=/rfs-storageservice/GIS/Projects/LOCCG/carriehc/data
export extdatadir=/rfs-storageservice/GIS/Projects/LOCCG/carriehc/extdata
```
Create necessary folders and subfolders at the '$datadir' location below:
```sh
mkdir ${datadir} ${datadir}/initialsamples ${datadir}/logs ${datadir}/mixtures \
 ${datadir}/mixtures_ultradeep ${datadir}/mixtures_wholegenome \
 ${datadir}/PoNbuffycoat ${datadir}/mixtures_ultradeep

ls ${datadir}
```
```
├── initialsamples
├── logs
├── mixtures
├── mixtures_ultradeep
├── mixtures_wholegenome
├── PoNbuffycoat
└── PoNbuffycoat_ultradeep
```
```sh
mkdir ${extdatadir} ${extdatadir}/dbsnp_vcf ${datadir}/exome_bed ${extdatadir}/GRCh37 \
${datadir}/wholegenome_bed 

ls ${extdatadir}
```
```
├── dbsnp_vcf
├── exome_bed
├── GRCh37
└── wholegenome_bed
```
 
### 1. Design benchmark dataset [figure 2, table 1]

From a cohort of longitudinal cfDNA samples, with lpWGS and deep targeted sequencing done for most timepoints.

#### 1. Identify suitable candidates for deep WGS 

Generate patient timeline to screen for patients matching both criteria \
about ichorCNA tumor fraction estimated from lpWGS sample, \
and about VAF of called mutations in deep targeted sequencing samples.

```sh
mutationfolder=data/variant_calls/226\ PANEL\ VARIANTS\ CLASSIFICSATION\ EXCEL
listpatients=()
for f in ${mutationfolder}/*; do
    if [[ $(basename $f) == CCG* ]]; then
        patient=$(echo $(basename $f) | awk -F'_' '{print $3}')
        echo $patient
        listpatients+=("$patient")
        python $repodir/initialsamples/patient_timeline_analysis.py --patient $patient
    fi
done
echo $listpatients

```

Further check by generating pileups on raw bamfiles for the targeted samples matching criteria.
```sh
export filename="${repodir}/initialsamples/candidate_samples_pileup_targeted.txt"
export lines=$(cat $filename)
for $line in $lines
do
    IFS='\t'
    read -a strarr <<< $line 
    export bamfile=${strarr[0]}
    export bedfile=${strarr[1]}
    export outputfile=${strarr[2]}
    export refgenome=${strarr[3]}
    export condapath=${strarr[4]}
    echo $(basename $bamfile) $(basename $bedfile) $(basename $outputfile) $(basename $refgenome)
    bash ${repodir}/initialsamples/pileup.sh $bamfile $bedfile $outputfile $refgenome $condpath
done
```

#### 2. Check deep WGS samples match criteria 

Generate pileup for the selected samples which underwent deep WGS.
```sh
export filename="${repodir}/initialsamples/candidate_samples_pileup_deepWGS.txt"
export lines=$(cat $filename)
for $line in $lines
do
    IFS='\t'
    read -a strarr <<< $line 
    export bamfile=${strarr[0]}
    export bedfile=${strarr[1]}
    export outputfile=${strarr[2]}
    export refgenome=${strarr[3]}
    export condapath=${strarr[4]}
    echo $(basename $bamfile) $(basename $bedfile) $(basename $outputfile) $(basename $refgenome)
    bash ${repodir}/initialsamples/pileup.sh $bamfile $bedfile $outputfile $refgenome $condpath
done
```

Get paired plot.

```sh
mutationfolder=data/variant_calls/226\ PANEL\ VARIANTS\ CLASSIFICSATION\ EXCEL
listpatients=()
for f in ${mutationfolder}/*; do
    if [[ $(basename $f) == CCG* ]]; then
        patient=$(echo $(basename $f) | awk -F'_' '{print $3}')
        echo $patient
        listpatients+=("$patient")
        python $repodir/initialsamples/paireplots.py --patient $patient
    fi
done
```

To reproduce Figure 2 of manuscript, run:
```sh
cd $repodir/plots
python figure2.py
```

#### 3. Design mixtures series

Run
```sh
cd $repodir/intialsamples
python design_mixtures.py
```

### 2. Create benchmark dataset


For deep Whole Genome Sequencing (WGS) mixture series, 
the samples are splitted by chromosomes due to storage space and memory limitations. \
For deep Whole Exome Sequencing (WES) mixture series, all the chromsomes are processed together.

#### 1. Prepare configuration files

Adapt the template '${repodir}/mixtures/config/config_template.yaml' to create your own.\
Here, the configuration file example is 'config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml' for a WGS mixture series. \
```sh
cd ${repodir}/mixtures
cp config/config_template.yml config/config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml
vi config/config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml
```
In the WGS case, write config file for chr1. You need to specify chr1 in the config file name.
Then, run the following to create the following config files for chr2-22: 
```sh
cd ${repodir}/mixtures
bash create_config_yaml.sh config/config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml
```

#### 2. Run mixture series creations in HPC as jobs (in Slurm)

Example for a WGS series splitted by chromosome. Here for chr1.
```sh
mkdir -p ${datadir}/data/logs
sbatch -p normal -J chr1_986_mixtures -t 24:00:00 --mem 48000 \
 --output=${datadir}/data/logs/chr1_986_mixture_wgs.o \
 --error=${datadir}/data/logs/chr1_986_mixtures_wgs.e \
 ${repodir}/mixtures/create_mixtures_series_chr.sh \
  -c ${repodir}/mixtures/config/config_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml
```

#### 3. Alternatively, download paper dataset from s3 bucket using the aws cli API.

Configure first the ~/.aws/credentials file to add the READ-ONLY access key of the database s3 bucket.
Here is displayed the database tree corresponding to the output of the mixture series creation detailed above.
A total of 8 samples were generated with different tumor fractions and restricted to chr1. 
For each mixture sample, the index bam file is present with the log files as well as the read group file, coverage, roughly estimated TF, and ichorCNA TF estimate.
```sh
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

#### 2.1 Install callers

This benchmark study considers 9 different callers in the following versions:
- 5 callers included in the bcbio pipeline v.1.2.9
  - Freebayes: v.1.3.5 
  - Mutect2: v2.2
  - Strelka2: 2.9.10
  - Vardict: 1.8.2
  - Varscan: 2.4.4
- SMuRF: v2.0.12
- Varnet: v1.1.0
- ABEMUS: v.1.0.3
- SiNVICT: master branch commit id 8e87e8d25c19d287dd68c7daa7375095dc099fa5 (2020)
Those callers need to be installed before hand.
Please refer to the guidelines provided by each software.
You can also find some help in the installation notes located in '${repodir}/callers/installation_notes.md'.

#### 2.2 Prepare configuration file for each caller and each mixture series

For each caller named ${caller}, first adapt the template '${repodir}/callers/${caller}/config/config_template.yaml' to create your own.\
Here, the configuration file example is 'config_${caller}_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml' for a WGS mixture series. \
```sh
cd ${repodir}/callers/${caller}
cp config/config_template.yml config/config_${caller}_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml
vi config/config_${caller}_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml
```
In the WGS case, write config file for chr1. You need to specify chr1 in the config file name.
Then, run the following to create the following config files for chr2-22:
```sh
cd ${repodir}/callers/${caller}
bash create_config_yaml.sh config/config_${caller}_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml
```
Then, follow the following steps.

#### 2.3 DNA-specific callers using bcbio pipeline

* **Freebayes, Mutect2, Strelka2, Vardict, Varscan (via the bcbio-nextgen pipeline)**

*Install and Prepare bcbio input file* (see 2.1 and 2.2)

```
$ bash $repodir/bcbio/run_prepare_bcbio.sh
```
Run bcbio pipeline in Aquila
```
$ bash $repodir/bcbio/run_bcbio.sh
```

* **SMuRF**
  
*Install and Prepare SMuRF input file* (see 2.1 and 2.2)

*Run SMuRF*

* **VarNet**

*Install and Prepare VarNet input file* (see 2.1 and 2.2)

*Run VarNet* 
```sh
sbatch -p normal -J prepare_dbSNP -t 24:00:00 --mem 48000 \
 --output=${datadir}/data/logs/prepare_dbSNP.o \
 --error=${datadir}/data/logs/prepare_dbSNP.e \
  ${repodir}/callers/ABEMUS/prepare_dbSNP.sh $extdata $repopath
```


### 2.3 Runs cfDNA-specific callers 

* **ABEMUS**

*Install and Prepare ABEMUS input file* (see 2.1 and 2.2)

*Prepare dbSNP database*

```sh
sbatch -p normal -J prepare_dbSNP -t 24:00:00 --mem 48000 \
 --output=${datadir}/data/logs/prepare_dbSNP.o \
 --error=${datadir}/data/logs/prepare_dbSNP.e \
  ${repodir}/callers/ABEMUS/prepare_dbSNP.sh $extdata $repopath
```

*Prepare Panel of Normal buffycoats*

```sh
#TODO iterate over buffycoat bams and patient ids, declare condapath and mode
sbatch -p normal -J prepare_PoN -t 24:00:00 --mem 48000 \
 --output=${datadir}/data/logs/prepare_PoN.o \
 --error=${datadir}/data/logs/prepare_PoN.e \
  ${repodir}/callers/ABEMUS/prepare_PoN.sh $buffycoatbam $patientid $datadir $extdatadir $condapath $mode
```

*Run ABEMUS*

```sh
export chr='1'
sbatch -p normal -J abemus_${chr}_${mixtureid} -t 1-00:00:00 -N 1 --mem 64000 \
--output=${datadir}/logs/abemus_${mixtureid}_${chr}_${mode}.o \
--error=${datadir}/logs/abemus_${mixtureid}_${chr}_${mode}.o \
run_abemus.sh -c ${repodir}/callers/ABEMUS/config/config_abemus_mixtures_chr${chr}_${mixtureid}.yml ;

## WES calling on WGS data chrom 1: -p normal -t 1-00:00:00 -N 1 --mem 64000 
## WGS calling on WGS data chrom 1: -p largemem -t 3-00:00:00 -N 1 --mem 1000000 
## WES calling on WES data chrom all: -p largememlong  -t 14-00:00:00 --mem 1000000
```

* **SiNVICT**

*Install and Prepare SiNVICT input file* (see 2.1 and 2.2)

*Run SiNVICT*
```
$ screen -S sinvict

$ bash ~/sinvict/run_sinvict.sh -c ~/sinvict/config/config_mixtures_chr22_CRC-1014_180816-CW-T.yaml
```


### 3. Analyse results 

XXX calltables.py
XXX generate_groundtruth.py

#### 1. At intermediate depth (150x) [figure 3]

To reproduce Figure 3 of manuscript, run:
```sh
cd $repodir/plots
python figure3.py
```

#### 2. At high depth (2,000x) [figure 4]


To reproduce Figure 3 of manuscript, run:
```sh
cd $repodir/plots
python figure4.py
```


### 4. Feature analysis [figure 5]


To reproduce Figure 3 of manuscript, run:
```sh
cd $repodir/plots
python figure5.py
```






