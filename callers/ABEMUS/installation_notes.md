# Installation and Usage of ABEMUS version 1.0.3

## 1. Create conda environmnent

```sh
conda envs create -n ABEMUS
conda activate ABEMUS
conda install -c conda-forge r-base  # install R3.6 is ok
conda install -c confa-forge r-devtools # installing devtools from R install.packages('devtools') does not work
conda install -c conda-forge r-data.table  # rq: library parallel already installed by default in R
```

## 2. Install ABEMUS

```R
library( "devtools" )
devtools::install_github("cibiobcg/abemus") # does not work with build_vignettes = T 
library( "abemus" ) # check installation completed
```

## 3. Prepare data

Prepare inputs
```R
outdir <-  "/my_project/Abemus_analysis/"
sample.info.file <- "/my_project/info/sample_info_file.tsv"
targetbed <- "/my_project/info/regions.bed"
pacbamfolder_bychrom <- "/my_project/data/PaCBAM_outdir_bychrom"
```
Download dbSNP database
```sh
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_dbSNP_all.vcf.gz
tabix -p vcf GRCh37_latest_dbSNP_all.vcf.gz
```
Convert chromosome names
```sh
bcftools annotate --rename-chrs /data/extdata/GRCh37/convert_chrom_names.txt -o dbSNP.vcf GRCh37_latest_dbSNP_all.vcf.gz
```
Split per chromosome
```sh
conda install -c bioconda bcftools
conda install -c bioconda tabix
```
Need a SNP database file with single REF or ALT base per position, no duplicate positions and REF and ALT = 'A', 'T', 'G' or 'C'.
So there is a need to edit dbSNP database vcf file.
```sh
python3 edit_vcf.py
```

## 4. Install PaCBAM

```sh
git clone https://CibioBCG@bitbucket.org/CibioBCG/pacbam.git  
cd pacbam

sudo apt-get install gcc or conda install -c conda-forge gcc
sudo apt-get install zlib1g-dev or conda install -c anaconda zlib  # zlib library, zlib1g does not work

make -f Makefile.linux
```

## 5. Run PaCBAM to create pileup and pabs files

```sh
mkdir PaCBAM_outdir

~/bin/pacbam/pacbam \
    bam=/data/buffycoat_chr22/CRC-986_100215-BC_chr22.bam  \
    bed=/data/chr22_data/wholegenome_chr22_hg19.bed \
    vcf=homo_sapiens-chr22_edited.vcf \
    fasta=/data/GRCh37/GRCh37.fa \
    strandbias \
    mode=5 \
    out=/data/abemus_outdir/PaCBAM_outdir

for i in /data/dilutions_chr22/*.bam ; do 
    ~/bin/pacbam/pacbam \
    bam=$i \
    bed=/data/chr22_data/wholegenome_chr22_hg19.bed \
    vcf=homo_sapiens-chr22_edited.vcf \
    fasta=/data/GRCh37/GRCh37.fa \
    strandbias  mode=5  out=/data/abemus_outdir/PaCBAM_outdir
done

~/pacbam/pacbam \
    bam=../data/CRC-986_100215-BC_chr22.bam \
    bed=exome_chr22_hg19.bed \
    vcf=homo_sapiens-chr22_edited.vcf \
    fasta=../data/GRCh37/GRCh37.fa \
    strandbias mode=5 out=PaCBAM_outdir

~/pacbam/pacbam \
    bam=../data/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted.bam 
    bed=exome_chr22_hg19.bed \
    vcf=homo_sapiens-chr22_edited.vcf \
    fasta=../data/GRCh37/GRCh37.fa \
    strandbias mode=5 out=PaCBAM_outdir
```

## 6. Run ABEMUS

```sh
mkdir pacbam_data_bychrom 
bash run_abemus.sh &
```
OR run per chunck in parallel
```sh
for i in {00..10} ; do
    bash run_abemus_i.sh $i ;
done
```