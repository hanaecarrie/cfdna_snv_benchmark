#!/bin/bash

export patient=$1
export path_data=$2

# STEP 4: run bamsurgeon SNV

qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 00 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 01 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 02 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 03 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 04 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 05 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 06 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 07 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 08 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 09 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 10 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 11 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 12 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 13 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 14 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 15 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 16 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 17 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 18 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 19 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 20 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 21 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 22 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 23 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 24 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 25 32 $path_data
qsub -pe OpenMP 32 -l mem_free=128G,h_rt=144:00:00 ~/cfdna_snv_benchmark/analysis/preprocess_runbamsurgeonsnv.sh $patient 26 32 $path_data