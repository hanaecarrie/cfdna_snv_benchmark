
#bash generate_fastq.sh /mnt/cfdnaseries/mixtures/mixtures_chr3/mixtures_chr3_CRC-986_100215-CW-T_CRC-986_300316-CW-T /mnt/cfdnaseries/initialsamples/initialsamples_chr3/initialsamples_chr3_CRC-986-CW-N/CRC-986-CW-N_chr3.bam
#bash generate_fastq.sh /mnt/cfdnaseries/mixtures/mixtures_chr4/mixtures_chr4_CRC-986_100215-CW-T_CRC-986_300316-CW-T /mnt/cfdnaseries/initialsamples/initialsamples_chr4/initialsamples_chr4_CRC-986-CW-N/CRC-986-CW-N_chr4.bam
#bash generate_fastq.sh /mnt/cfdnaseries/mixtures/mixtures_chr5/mixtures_chr5_CRC-986_100215-CW-T_CRC-986_300316-CW-T /mnt/cfdnaseries/initialsamples/initialsamples_chr5/initialsamples_chr5_CRC-986-CW-N/CRC-986-CW-N_chr5.bam
#bash generate_fastq.sh /mnt/cfdnaseries/mixtures/mixtures_chr6/mixtures_chr6_CRC-986_100215-CW-T_CRC-986_300316-CW-T /mnt/cfdnaseries/initialsamples/initialsamples_chr6/initialsamples_chr6_CRC-986-CW-N/CRC-986-CW-N_chr6.bam
for chr in {7..10} ; do 
	bash generate_fastq.sh /mnt/cfdnaseries/mixtures/mixtures_chr${chr}/mixtures_chr${chr}_CRC-986_100215-CW-T_CRC-986_300316-CW-T /mnt/cfdnaseries/initialsamples/initialsamples_chr${chr}/initialsamples_chr${chr}_CRC-986-CW-N/CRC-986-CW-N_chr${chr}.bam
done
#bash generate_fastq.sh /mnt/cfdnaseries/mixtures/mixtures_chr3/mixtures_chr3_CRC-1014_180816-CW-T_CRC-1014_090516-CW-T /mnt/cfdnaseries/initialsamples/initialsamples_chr3/initialsamples_chr3_CRC-1014-CW-N/CRC-1014-CW-N_chr3.bam
#bash generate_fastq.sh /mnt/cfdnaseries/mixtures/mixtures_chr3/mixtures_chr3_CRC-986_100215-CW-T_CRC-986_300316-CW-T /mnt/cfdnaseries/initialsamples/initialsamples_chr3/initialsamples_chr3_CRC-986-CW-N/CRC-986-CW-N_chr3.bam
#bash generate_fastq.sh /mnt/cfdnaseries/mixtures/mixtures_chr3/mixtures_chr3_CRC-123_310715-CW-T_CRC-123_121115-CW-T /mnt/cfdnaseries/initialsamples/initialsamples_chr3/initialsamples_chr3_CRC-123-CW-N/CRC-123-CW-N_chr3.bam

