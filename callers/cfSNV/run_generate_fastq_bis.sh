

for chr in {1..3} ; do 
	bash generate_fastq_bis.sh /mnt/cfdnaseries/mixtures/mixtures_chr${chr}/mixtures_chr${chr}_CRC-986_100215-CW-T_CRC-986_300316-CW-T /mnt/cfdnaseries/initialsamples/initialsamples_chr${chr}/initialsamples_chr${chr}_CRC-986-CW-N/CRC-986-CW-N_chr${chr}.bam &
done

