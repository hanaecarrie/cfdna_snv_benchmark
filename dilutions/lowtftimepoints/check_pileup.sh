

#export samples = ['2015-10-01', '2015-06-11', '2016-07-05', '2016-05-11', '2016-04-21', '2016-03-30'] 
export positionfile='/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/lowtftimepoints/positionfile.txt'
declare -a samples=('986.011015' '986.110615' '986.050716' '986.110516' '986.210416' '986.300316')
export lines=$(cat $positionfile)

for sample in ${samples[@]}; do

for line in $lines; do
export chr="$(cut -d',' -f1 <<<"$line")"
export pos_start="$(($(cut -d',' -f2 <<<"$line")))" 
export pos_end=$pos_start
export gene="$(cut -d',' -f3 <<<"$line")"
echo $sample
echo $gene
#  -f /mnt/projects/carriehc/cfDNA/data/refgenome/hg19/hg19.fa \
/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/samtools mpileup \
-f /mnt/projects/ngbhs/cirqseq/references/hg38.fa \
-r $chr:$pos_start-$pos_end \
 /mnt/projects/huangwt/wgs/hanae/sequenza_resort_bams/CCG_226_${sample}.reordered.bam >> /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/lowtftimepoints/output_pileup_${sample}.txt
#/mnt/projects/zwpoh/cfDNA/bulk/ccg/ccg_batch2/lpwgs/hg19_bam/patient/${sample}-CCG-sorted.bam >> /mnt/projects/carriehc/cfDNA/cfSNV/benchmark/lowtftimepoints/output_pileup_${sample}.txt ;
done
done
