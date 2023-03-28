#!/bin/bash

#SBATCH -p normal
#SBATCH -J 5_986_sinvict
#SBATCH -t 3-00:00:00
#SBATCH -N 1
#SBATCH --mem 64000
#SBATCH --output=/rfs-storageservice/GIS/Projects/LOCCG/carriehc/sinvict_outdir/logs/z.5_986_cfsnv.o
#SBATCH --error=/rfs-storageservice/GIS/Projects/LOCCG/carriehc/sinvict_outdir/logs/z.5_986_cfsnv.e

export chr=5

cd /home/users/astar/gis/carriehc/cfdna_snv_benchmark/callers/sinvict

bash run_sinvict.sh -c config/config_sinvict_mixtures_chr${chr}_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml 

