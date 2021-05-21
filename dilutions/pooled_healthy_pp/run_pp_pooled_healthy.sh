
#!/bin/bash

# def path and environment
source /home/carriehc/miniconda2/etc/profile.d/conda.sh
conda activate default

python3 ~/cfdna_snv_benchmark/analysis/preprocess_pooled_healthy.py 986_100215 genomead $1 $2 /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/pooled_healthy_pp/
