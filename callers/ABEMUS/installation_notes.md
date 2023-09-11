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

In R
```R
library( "devtools" )
devtools::install_github("cibiobcg/abemus") # does not work with build_vignettes = T 
library( "abemus" ) # check installation completed
```

## 3. Install PaCBAM

In bash
```sh
cd $HOME/bin
git clone https://CibioBCG@bitbucket.org/CibioBCG/pacbam.git  
cd pacbam

sudo apt-get install gcc or conda install -c conda-forge gcc
sudo apt-get install zlib1g-dev or conda install -c anaconda zlib  # zlib library, zlib1g does not work

make -f Makefile.linux # if linux otherwise Makefile.macos or Makefile.mingw
```