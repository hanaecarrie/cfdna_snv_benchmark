
export i=$1
echo $i

java -Xmx16G -jar /home/ubuntu/bin/abra2/target/abra2-2.24-jar-with-dependencies.jar --in /data/buffycoat_chr22/CRC-986_100215-BC_chr22.bam,/data/dilutions_chr22/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted.bam --out /data/sinvict_outdir/abra/CRC-986_100215-BC_chr22_${i}.abra.bam,/data/sinvict_outdir/abra/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted_${i}.abra.bam --ref /data/GRCh37/GRCh37.fa --threads 8 --targets /data/chr22_data/wholegenome_chr22_hg19_${i}.bed  --tmpdir /data/sinvict_outdir/tmp/ > /data/sinvict_outdir/abra/abra_${i}.log

#/usr/bin/samtools index /data/sinvict_outdir/abra/CRC-986_100215-BC_chr22_${i}.abra.bam
/usr/bin/samtools index /data/sinvict_outdir/abra/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted_${i}.abra.bam

#/home/ubuntu/bin/bam-readcount/build/bin/bam-readcount -f /data/GRCh37/GRCh37.fa /data/sinvict_outdir/abra/CRC-986_100215-BC_chr22.abra.bam -l /data/sinvict_outdir/exome_chr22_hg19.bed >  /data/sinvict_outdir/bam-readcount/CRC-986_100215-BC_chr22_${i}.tsv
/home/ubuntu/bin/bam-readcount/build/bin/bam-readcount -f /data/GRCh37/GRCh37.fa /data/sinvict_outdir/abra/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted_${i}.abra.bam -l /data/chr22_data/exome_chr22_hg19.bed >  /data/sinvict_outdir/bam-readcount/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted_${i}.tsv

#/home/ubuntu/sinvict/sinvict -t /data/sinvict_outdir/bam-readcount -o /data/sinvict_outdir/results

