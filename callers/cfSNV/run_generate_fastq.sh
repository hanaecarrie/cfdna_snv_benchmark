#!/bin/bash

#SBATCH -p normal
#SBATCH -J fqchr15
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH --mem 48000
#SBATCH --output=/scratch/users/astar/gis/carriehc/logs/z.fastqchr15.o
#SBATCH --error=/scratch/users/astar/gis/carriehc/logs/z.fastqchr15.e

export chr=15

bash generate_fastq_bis.sh /rfs-storageservice/GIS/Projects/LOCCG/carriehc/data/mixtures/mixtures_chr${chr}/mixtures_chr${chr}_CRC-986_100215-CW-T_CRC-986_300316-CW-T /rfs-storageservice/GIS/Projects/LOCCG/carriehc/data/initialsamples/initialsamples_chr${chr}/initialsamples_chr${chr}_CRC-986-CW-N/CRC-986-CW-N_chr${chr}.bam

bash generate_fastq_bis.sh /rfs-storageservice/GIS/Projects/LOCCG/carriehc/data/mixtures/mixtures_chr${chr}/mixtures_chr${chr}_CRC-1014_180816-CW-T_CRC-1014_090516-CW-T /rfs-storageservice/GIS/Projects/LOCCG/carriehc/data/initialsamples/initialsamples_chr${chr}/initialsamples_chr${chr}_CRC-1014-CW-N/CRC-1014-CW-N_chr${chr}.bam

bash generate_fastq_bis.sh /rfs-storageservice/GIS/Projects/LOCCG/carriehc/data/mixtures/mixtures_chr${chr}/mixtures_chr${chr}_CRC-123_310715-CW-T_CRC-123_121115-CW-T /rfs-storageservice/GIS/Projects/LOCCG/carriehc/data/initialsamples/initialsamples_chr${chr}/initialsamples_chr${chr}_CRC-123-CW-N/CRC-123-CW-N_chr${chr}.bam

