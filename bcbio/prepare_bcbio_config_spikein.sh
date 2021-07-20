
cd /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/bcbio/

export typefolder=$1
export bcpath=$2
export bcbioconfigpath=$3 #/mnt/projects/skanderupamj/wgs/data/cfdna.crc/mix.samples/mix.samples.nofilters.yaml
export bcbiooutput=$4

touch $bcbiooutput

for vafdir in $typefolder/*;

do echo "$vafdir"
# dilution cfdna sample with pseudo matched healthy
export spikeinpaths=$(ls $vafdir/*sorted.bam)
for spikeinpath in $spikeinpaths 
do echo $spikeinpath
export spikename=$(basename $spikeinpath .sorted.bam)
# buffy coat
export bcdir=$(dirname $bcpath)
export bcfile=$(basename $bcpath)

echo "$spikename	symlink::$vafdir	$(basename $spikeinpath)	symlink::$bcdir	$bcfile	MIX	NA	NA	$bcbioconfigpath" >> $bcbiooutput
done
done
