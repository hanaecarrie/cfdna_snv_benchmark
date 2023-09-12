#!/bin/bash

#SBATCH -p express
#SBATCH -J XXX
#SBATCH -t 00:59:00
#SBATCH -N 1
#SBATCH --mem 8000
#SBATCH --output=bcbiodir/XXX/bcbio_work/smurf_workdir/smurf.XXX.o
#SBATCH --error=bcbiodir/XXX/bcbio_work/smurf_workdir/smurf.XXX.e

mkdir -p bcbiodir/XXX/bcbio_final/smurf

rscript smurfdir/run.smurf.R \
  -i bcbiodir/XXX/bcbio_final/2015-07-31_XXX \
  -o bcbiodir/XXX/bcbio_final/smurf \
  -r hg19 -p 50000 -a h2olib \
  -s smurflib