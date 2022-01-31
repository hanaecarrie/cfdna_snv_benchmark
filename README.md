# cfdna_snv_benchmark

## Requirements


## Install


add the repository to your Python path
export PYTHONPATH=$PYTHONPATH:"/Users/hanae/Repositories/cfdna_snv_benchmark/"

## Demo



### Create dilution series

qsub -pe OpenMP 4 -l mem_free=32G,h_rt=24:00:00 -e /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/logs/ -o /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/logs/ -M hanae_camille_carrie_from.tp@gis.a-star.edu.sg /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/dilutions/create_dilution_series_chr.sh -c /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/config/config_dilutions/config_lowtfsample_1014_110116.yml


### Create spikein series



### Prepare bcbio input



### Run bcbio pipeline in Aquila



### Copy results on subfolder data of the repository


### Transfer results to Ronin s3

qsub -pe OpenMP 2 -l mem_free=4G,h_rt=24:00:00 -e /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/utils/logs/ -o /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/utils/logs/ /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/utils/cptos3.sh

### Copy to local machine from s3

