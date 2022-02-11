
for file in  /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/dilutions_series_chr22/dilutions_CRC-986_100215/dilution_chr22_CRC-986_100215_1*_CRC-986_300316*/*sorted.bam* ;
do
echo $file ;
scp -i ~/cfSNV.pem $file ubuntu@cfsnv.genome.sg:/home/ubuntu/dilutions_chr22/ ;
done
