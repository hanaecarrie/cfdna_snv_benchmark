


java -Xmx4G -jar /home/ubuntu/sinvict/abra/temp/abra-0.94c.jar --in /home/ubuntu/data/CRC-986_100215-BC_chr22.bam,/home/ubuntu/data/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted.bam --out /home/ubuntu/sinvict/output/CRC-986_100215-BC_chr22.abra.bam,/home/ubuntu/sinvict/output/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted.abra.bam --ref /home/ubuntu/data/GRCh37/GRCh37.fa --threads 8 --targets /home/ubuntu/data/exome_chr22_hg19.bed  --working /home/ubuntu/sinvict/tmp/ > /home/ubuntu/sinvict/output/abra.log 2>&1

#/usr/bin/samtools index /home/ubuntu/sinvict/output/CRC-986_100215-BC_chr22.abra.bam
#/usr/bin/samtools index /home/ubuntu/sinvict/output/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted.abra.bam

#/home/ubuntu/bam-readcount/build/bin/bam-readcount -f /home/ubuntu/data/GRCh37/GRCh37.fa /home/ubuntu/sinvict/abra/CRC-986_100215-BC_chr22.abra.bam -l /home/ubuntu/data/exome_chr22_hg19.bed >  /home/ubuntu/sinvict/bam-readcount/CRC-986_100215-BC_chr22.tsv
#/home/ubuntu/bam-readcount/build/bin/bam-readcount -f /home/ubuntu/data/GRCh37/GRCh37.fa /home/ubuntu/sinvict/abra/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted.abra.bam -l /home/ubuntu/data/exome_chr22_hg19.bed >  /home/ubuntu/sinvict/bam-readcount//dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted.tsv

#/home/ubuntu/sinvict/sinvict -t /home/ubuntu/sinvict/bam-readcount -o /home/ubuntu/sinvict/results

