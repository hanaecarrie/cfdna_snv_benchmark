#!/bin/bash

#SBATCH -p normal
#SBATCH -J jobid
#SBATCH -t 2-00:00:00
#SBATCH -N 1
#SBATCH --mem 32000
#SBATCH --output=outputdir/z.sampleid.o
#SBATCH --error=outputdir/z.sampleid.e


