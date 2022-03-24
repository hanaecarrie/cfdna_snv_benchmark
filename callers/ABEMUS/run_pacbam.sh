
#~/bin/pacbam/pacbam bam=../data/CRC-986_100215-BC_chr22.bam  bed=../chr22_data/exome_chr22_hg19.bed vcf=homo_sapiens-chr22_edited.vcf fasta=../data/GRCh37/GRCh37.fa strandbias mode=5 out=PaCBAM_outdir

for file in ~/dilutions_chr22/*.bam ; 
do echo $file ;
	~/bin/pacbam/pacbam bam=$file  bed=../chr22_data/exome_chr22_hg19.bed vcf=homo_sapiens-chr22_edited.vcf fasta=../data/GRCh37/GRCh37.fa strandbias mode=5 out=PaCBAM_outdir
done
