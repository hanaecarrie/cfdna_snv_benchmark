#!/bin/bash

#SBATCH -p normal
#SBATCH -J 11_986_cfsnv
#SBATCH -t 3-00:00:00
#SBATCH -N 1
#SBATCH --mem 64000
#SBATCH --output=/rfs-storageservice/GIS/Projects/LOCCG/carriehc/cfsnv_outdir/logs/z.11_986_cfsnv.o
#SBATCH --error=/rfs-storageservice/GIS/Projects/LOCCG/carriehc/cfsnv_outdir/logs/z.11_986_cfsnv.e

export chr=11

cd /home/users/astar/gis/carriehc/cfdna_snv_benchmark/callers/cfSNV

export npid=0
for plasma in /rfs-storageservice/GIS/Projects/LOCCG/carriehc/data/mixtures/mixtures_chr${chr}/mixtures_chr${chr}_CRC-986_100215-CW-T_CRC-986_300316-CW-T/*/*[Tx].bam ; do 
	export npid=$((npid+1)) ;
	echo "n PID ${npid}" ;
	bash run_cfsnv_sample.sh -c config/config_cfsnv_mixtures_chr${chr}_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml -p $plasma  &  pids[${npid}]=$! ;
done

# wait for all jobs
for pid in ${pids[*]} ; do 
	wait $pid 
done

