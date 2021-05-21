
touch qsub_pp_pooled_healthy.sh ;
for i in $(seq 0 45) ; do echo "qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/run_pp_pooled_healthy.sh "$(($i * 10))" "$((($i +1)* 10))" " >> qsub_pp_pooled_healthy.sh ; done
