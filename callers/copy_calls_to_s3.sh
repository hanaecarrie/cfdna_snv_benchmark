
export datapath=$1 #'/mnt/projects/huangwt/wgs/hanae/'
export seriespath=$2 #'/mnt/projects/huangwt/wgs/hanae/data/mixtures/mixtures_chr22/mixtures_chr22_CRC-1014_180816-CW-T/'

# cfdna callers on Ronin
for sample in $datapath/$seriespath/*[Tx] ; 
do export sample=$(basename $sample)
echo $sample
echo /mnt/projects/skanderupamj/wgs/bcbio_v107/bin/aws s3 cp $datapath/${seriespath}${sample}/calls s3://cfdna-mutation-calling.store.genome.sg/${seriespath}${sample}/calls --profile hanae --recursive
/mnt/projects/skanderupamj/wgs/bcbio_v107/bin/aws s3 cp $datapath/${seriespath}${sample}/calls s3://cfdna-mutation-calling.store.genome.sg/${seriespath}${sample}/calls --profile hanae --recursive
done
