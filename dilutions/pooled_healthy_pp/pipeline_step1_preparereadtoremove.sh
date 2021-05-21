
export patient_id=$1

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/preprocess_preparereadstoremove.sh $patient_id dbsnp 0 0 '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/' 
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/preprocess_preparereadstoremove.sh $patient_id genomead 0 10 '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/preprocess_preparereadstoremove.sh $patient_id genomead 10 20 '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/preprocess_preparereadstoremove.sh $patient_id genomead 20 30 '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/preprocess_preparereadstoremove.sh $patient_id genomead 30 40 '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/preprocess_preparereadstoremove.sh $patient_id genomead 40 50 '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/preprocess_preparereadstoremove.sh $patient_id genomead 50 60 '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/preprocess_preparereadstoremove.sh $patient_id genomead 60 70 '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/preprocess_preparereadstoremove.sh $patient_id genomead 70 80 '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/preprocess_preparereadstoremove.sh $patient_id genomead 80 90 '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'
qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/preprocess_preparereadstoremove.sh $patient_id genomead 90 100 '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'
