
export abemusoutdir='/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/callers_output/abemus/mixtures_chr22_CRC-1014_180816-CW-T'
if [ ! -d $abemusoutdir ] ; then mkdir $abemusoutdir ; fi

scp -i /home/carriehc/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/output/abemus_outdir/mixtures_chr22_CRC-1014_180816-CW-T/results/ $abemusoutdir

export sinvictoutdir='/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/callers_output/sinvict/mixtures_chr22_CRC-1014_180816-CW-T'
if [ ! -d $sinvictoutdir ] ; then mkdir $sinvictoutdir ; fi

scp -i /home/carriehc/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/output/sinvict_outdir/mixtures_chr22_CRC-1014_180816-CW-T/*/results/ $sinvictoutdir

export cfsnvoutdir='/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/callers_output/cfsnv//mixtures_chr22_CRC-1014_180816-CW-T'
if [ ! -d $cfsnvoutdir ] ; then mkdir $cfsnvoutdir ; fi

scp -i /home/carriehc/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/output/cfsnv_outdir/test/dilution_chr22_CRC-986_100215/testpicard/results/ $cfsnvoutdir


