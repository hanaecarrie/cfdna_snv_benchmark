
export datapath=$1 #'/mnt/projects/huangwt/wgs/hanae/'
export seriespath=$2 #'data/mixtures/mixtures_chr22/mixtures_chr22_CRC-1014_180816-CW-T/'
export dilutionseries=$(basename $seriespath) #'mixtures_chr22_CRC-1014_180816-CW-T'
echo $datapath
echo $seriespath
echo $dilutionseries

# cfdna callers on Ronin
for sample in $datapath/$seriespath/*[Tx] ; 
do export sample=$(basename $sample)
echo $sample

if [ ! -d $datapath/$seriespath/$sample/calls ] ; then mkdir $datapath/$seriespath/$sample/calls ; fi
echo 'copy ABEMUS from Ronin...'
echo "scp -i /home/carriehc/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/outdir/cfdna_callers/$dilutionseries/$sample/abemus $datapath/$seriespath/$sample/calls/"
scp -i /home/carriehc/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/outdir/cfdna_callers/$dilutionseries/$sample/abemus $datapath/$seriespath/$sample/calls/

echo 'copy SiNVICT from Ronin...'
echo "scp -i /home/carriehc/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/outdir/cfdna_callers/$dilutionseries/$sample/sinvict $datapath/$seriespath/$sample/calls/"
scp -i /home/carriehc/ronin_cfdna_benchmark_hanae.pem -r ubuntu@cfdna_callers.genome.sg:/outdir/cfdna_callers/$dilutionseries/$sample/sinvict $datapath/$seriespath/$sample/calls/

#TODO copy cfSNV

echo 'copy bcbio from Aquila...'
if [ ! -d  $datapath/$seriespath/$sample/calls/bcbio ] ; then mkdir $datapath/$seriespath/$sample/calls/bcbio ; fi
echo "scp -r /mnt/projects/skanderupamj/wgs/data/cfdna.crc/benchmark/$sample/bcbio_final/2015-07-31_${sample}/* $datapath/$seriespath/$sample/calls/bcbio/"
scp -r /mnt/projects/skanderupamj/wgs/data/cfdna.crc/benchmark/$sample/bcbio_final/2015-07-31_${sample}/* $datapath/$seriespath/$sample/calls/bcbio/

done
