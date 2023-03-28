#!/bin/bash

#SBATCH -p normal
#SBATCH -J chrall_986_mixtures_wes
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem 48000
#SBATCH --output=/rfs-storageservice/GIS/Projects/LOCCG/carriehc/data_wes/logs/z.chrall_986_mixture_wes.o
#SBATCH --error=/rfs-storageservice/GIS/Projects/LOCCG/carriehc/data_wes/logs/z.chrall_986_mixtures_wes.e

source /apps/anaconda3-individual-edition/2020.11/etc/profile.d/conda.sh

cd /home/users/astar/gis/carriehc/cfdna_snv_benchmark/mixtures

bash create_mixtures_series_chr.sh -c /home/users/astar/gis/carriehc/cfdna_snv_benchmark/mixtures/config/config_mixtures_chrall_CRC-986_100215-CW-T_CRC-986_300316-CW-T_WES.yml


