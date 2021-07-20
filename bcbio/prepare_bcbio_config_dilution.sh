
cd /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/bcbio/

export dilfolder=$1
export bcpath=$2
export bcbioconfigpath=$3 #/mnt/projects/skanderupamj/wgs/data/cfdna.crc/mix.samples/mix.samples.nofilters.yaml
export bcbiooutput=$4 #.tsv

touch $bcbiooutput

for dildir in $dilfolder/*;

do echo "$dildir"
# dilution cfdna sample with pseudo matched healthy
export dilname=$(basename $dildir)
export dilfile=${dilname}.sorted.bam
# buffy coat
export bcdir=$(dirname $bcpath)
export bcfile=$(basename $bcpath)

echo "$dilname	symlink::$dildir	$dilfile	symlink::$bcdir	$bcfile	MIX	NA	NA	$bcbioconfigpath" >> $bcbiooutput
done
