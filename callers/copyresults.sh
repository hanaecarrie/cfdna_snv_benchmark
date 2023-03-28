
export dilutionseries='mixtures_chr22_CRC-1014_180816-CW-T'
export mixtureseriesoutdir='/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/callers_output/mixtures/mixtures_chr22/mixtures_chr22_CRC-1014_180816-CW-T/'

# cfdna callers on Ronin
echo "scp -i /home/carriehc/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/output/cfdna_callers/$dilutionseries/ $mixtureseriesoutdir"
#scp -i /home/carriehc/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/output/cfdna_callers/$dilutionseries/ $mixtureseriesoutdir

# bcbio pipeline on Aquila
for sample in $mixtureseriesoutdir/* ; do export sample=$(basename $sample .sorted) ; echo $sample ; 
if [ ! -d $mixtureseriesoutdir/${sample}.sorted/bcbio ] ; then mkdir  $mixtureseriesoutdir/${sample}.sorted/bcbio ; fi
scp -r /mnt/projects/skanderupamj/wgs/data/cfdna.crc/benchmark/$sample/bcbio_final/2015-07-31_${sample}/* $mixtureseriesoutdir/${sample}.sorted/bcbio/
done
