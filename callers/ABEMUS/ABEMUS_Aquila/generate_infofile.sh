
touch infofile_test.tsv

for dil in /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/dilutions_series_chr22/dilutions_CRC-986_100215/dilution*_CRC-986_300316*/*.sorted.bam ;
do echo -e '986_100215_chr22\t'$(basename $dil .bam)'\t'$dil'\tCRC-986_100215-BC_chr22\`t/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/dilutions_series_chr22/NCC_CRC-986_100215-CW-N/CRC-986_100215-BC_chr22.bam'  >> infofile_test.tsv ;
done

