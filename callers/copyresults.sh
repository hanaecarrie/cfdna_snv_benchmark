
export abemusoutdir='/Users/hanae/Repositories/cfdna_snv_benchmark/data/callers_output/abemus/mixtures_chr22_CRC-1014_180816-CW-T'
echo $abemusoutdir
if [ ! -d $abemusoutdir ] ; then mkdir $abemusoutdir ; fi

#scp -i  /Users/hanae/Documents/3_Monde\ professionnel/6_PhD/1_Research/Project\ 3\ -\ cfDNA\ SNV\ Calling/Ronin/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/output/abemus_outdir/mixtures_chr22_CRC-1014_180816-CW-T/results/ $abemusoutdir

export sinvictoutdir='/Users/hanae/Repositories/cfdna_snv_benchmark/data/callers_output/sinvict/mixtures_chr22_CRC-1014_180816-CW-T'
echo $sinvictoutdir
if [ ! -d $sinvictoutdir ] ; then mkdir $sinvictoutdir ; fi
find . -type d -wholename '/output/sinvict_outdir/mixtures_chr22_CRC-1014_180816-CW-T/*/results/' \
| xargs tar cf - \
| ssh -i  /Users/hanae/Documents/3_Monde\ professionnel/6_PhD/1_Research/Project\ 3\ -\ cfDNA\ SNV\ Calling/Ronin/ronin_cfdna_benchmark_hanae.pem tar xf - -C $sinvictoutdir 
#scp -i  /Users/hanae/Documents/3_Monde\ professionnel/6_PhD/1_Research/Project\ 3\ -\ cfDNA\ SNV\ Calling/Ronin/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/output/sinvict_outdir/mixtures_chr22_CRC-1014_180816-CW-T/*/results/ $sinvictoutdir

export cfsnvoutdir='/Users/hanae/Repositories/cfdna_snv_benchmark/data/callers_output/cfsnv//mixtures_chr22_CRC-1014_180816-CW-T'
echo $cfsnvoutdir
if [ ! -d $cfsnvoutdir ] ; then mkdir $cfsnvoutdir ; fi

#scp -i  /Users/hanae/Documents/3_Monde\ professionnel/6_PhD/1_Research/Project\ 3\ -\ cfDNA\ SNV\ Calling/Ronin/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/output/cfsnv_outdir/test/dilution_chr22_CRC-986_100215/testpicard/results/ $cfsnvoutdir


