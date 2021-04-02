
export patient_id=$1 # 986
export path_data=$2 # '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 ~/cfnda_snv_benchmark/preprocess_preparereadstoremove.sh $patient_id dbsnp 0 0 $oath_data
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 ~/cfnda_snv_benchmark/preprocess_preparereadstoremove.sh $patient_id genomead 0 10 $oath_data
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 ~/cfnda_snv_benchmark/preprocess_preparereadstoremove.sh $patient_id genomead 10 20 $oath_data
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 ~/cfnda_snv_benchmark/preprocess_preparereadstoremove.sh $patient_id genomead 20 30 $oath_data
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 ~/cfnda_snv_benchmark/preprocess_preparereadstoremove.sh $patient_id genomead 30 40 $oath_data
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 ~/cfnda_snv_benchmark/preprocess_preparereadstoremove.sh $patient_id genomead 40 50 $oath_data
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 ~/cfnda_snv_benchmark/preprocess_preparereadstoremove.sh $patient_id genomead 50 60 $oath_data
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 ~/cfnda_snv_benchmark/preprocess_preparereadstoremove.sh $patient_id genomead 60 70 $oath_data
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 ~/cfnda_snv_benchmark/preprocess_preparereadstoremove.sh $patient_id genomead 70 80 $oath_data
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 ~/cfnda_snv_benchmark/preprocess_preparereadstoremove.sh $patient_id genomead 80 90 $oath_data
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 ~/cfnda_snv_benchmark/preprocess_preparereadstoremove.sh $patient_id genomead 90 100 $oath_data
