
for i in /mnt/projects/skanderupamj/wgs/data/cfdna.crc/mix.samples/dilution_chr22_CRC-* ; do export fname=$(basename $i) ; echo $fname ;  if [ -d $i/bcbio_final/2015-07-31_${fname} ] ; then cp -r $i/bcbio_final/2015-07-31_${fname} /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/bcbio_output/${fname} ; fi ;  done
