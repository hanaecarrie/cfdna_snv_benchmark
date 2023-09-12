#!/bin/bash

#SBATCH -p normal
#SBATCH -J jobid
#SBATCH -t 3-00:00:00
#SBATCH -N 1
#SBATCH --mem 60000
#SBATCH --output=bcbiodir/sampleid/bcbio_work/z.sampleid.o
#SBATCH --error=bcbiodir/sampleid/bcbio_work/z.sampleid.e

cd bcbiodir/sampleid/bcbio_work

bcbionextgenpath ../config.yaml -n 24 --timeout 60

