
export patient=$1 # 986
export germline_vcf_name=$2 # buffycoat_CRC-986_chr22_nofilter
export path_data=$3 # '/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/'

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 -N dbsnp -m e -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg \
-o $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.o \
-e $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.e  \
~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.sh $patient $germline_vcf_name dbsnp 0 0 $oath_data

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 -N gnomad_0 -m e -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg \
-o $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.o \
-e $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.e  \
~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.sh $patient $germline_vcf_name gnomad 0 10 $oath_data

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 -N gnomad_1-m e -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg \
-o $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.o \
-e $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.e  \
~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.sh $patient $germline_vcf_name gnomad 10 20 $oath_data

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 -N gnomad_2 -m e -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg \
-o $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.o \
-e $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.e  \
~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.sh $patient $germline_vcf_name gnomad 20 30 $oath_data

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 -N gnomad_3 -m e -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg \
-o $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.o \
-e $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.e  \
~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.sh $patient $germline_vcf_name gnomad 30 40 $oath_data

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 -N gnomad_4 -m e -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg \
-o $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.o \
-e $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.e  \
~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.sh $patient $germline_vcf_name gnomad 40 50 $oath_data

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 -N gnomad_5 -m e -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg \
-o $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.o \
-e $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.e  \
~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.sh $patient $germline_vcf_name gnomad 50 60 $oath_data

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 -N gnomad_6 -m e -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg \
-o $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.o \
-e $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.e  \
~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.sh $patient $germline_vcf_name gnomad 60 70 $oath_data

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 -N gnomad_7 -m e -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg \
-o $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.o \
-e $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.e  \
~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.sh $patient $germline_vcf_name gnomad 70 80 $oath_data

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 -N gnomad_8 -m e -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg \
-o $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.o \
-e $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.e  \
~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.sh $patient $germline_vcf_name gnomad 80 90 $oath_data

qsub -pe OpenMP 2 -l mem_free=24G,h_rt=24:00:00 -N gnomad_9 -m e -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg \
-o $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.o \
-e $path_data/prepare_pooled_healthy/$patient/preprocess_preparereadstoremove_dbsnp.e  \
~/cfdna_snv_benchmark/preprocessing/preprocess_preparereadstoremove.sh $patient $germline_vcf_name gnomad 90 100 $oath_data
