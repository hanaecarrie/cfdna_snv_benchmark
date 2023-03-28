#!/bin/bash

# function to parse config file
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}


while getopts ":c:i:p:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
    i) i=${OPTARG} ;;
    p) plasma=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

source $condapath
conda activate default

cd ${repopath}/sinvict


echo $config_file
echo $dilutionseriesfolder
echo $buffycoatbam
echo $chr
echo $extdata
echo $outdir
echo $abra
echo $mode

echo $i

echo $plasma

export outdirplasma=$outdir/$(basename $plasma .bam)
echo $outdirplasma
if [ ! -d $outdirplasma ] ; then mkdir $outdirplasma ; fi
if [ ! -d $outdirplasma/log ] ; then mkdir $outdirplasma/log ; fi
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>$outdirplasma/log/log_${i}.out 2>&1
# Everything below will go to the file 'log.out':

if [ ! -f $outdirplasma/log/logtime_${i}.out ] ; then touch $outdirplasma/log/logtime_${i}.out ; fi

### ABRA ###

if [ $abra = "True" ] ; then

if [ ! -d $outdirplasma/abra ] ; then mkdir $outdirplasma/abra ; fi
if [ ! -d $outdirplasma/tmp ] ; then mkdir $outdirplasma/tmp ; fi
startabra=$(date +%s)
if [ ! -f ${outdirplasma}/abra/$(basename $plasma .bam)_${i}.abra.bam ] ;
	then java -Xmx16G -jar /home/users/astar/gis/carriehc/bin/abra2/target/abra2-2.24-jar-with-dependencies.jar \
	--in $buffycoatbam,$plasma \
	--out ${outdirplasma}/abra/$(basename $buffycoatbam .bam)_${i}.abra.bam,${outdirplasma}/abra/$(basename $plasma .bam)_${i}.abra.bam \
	--ref ${extdata}/GRCh37/GRCh37.fa --threads 8 --targets ${extdata}/${mode}_bed/${mode}_hg19_chr${chr}_${i}.bed  --tmpdir ${outdirplasma}/tmp/ > ${outdirplasma}/abra/abra_${i}.log
fi
endabra=$(date +%s)
timeabra=$(($endabra-$startabra))
hmsabra=$(printf '%02dh:%02dm:%02ds\n' $((timeabra/3600)) $((timeabra%3600/60)) $((timeabra%60)))
echo "Elapsed Time ABRA on ${plasma} for chunk ${i} of chr${chr}: ${hmsabra}"
echo "Elapsed Time ABRA on ${plasma} for chunk ${i} of chr${chr}: ${hmsabra}" >> $outdirplasma/log/logtime_${i}.out
if [ ! -f ${outdirplasma}/abra/$(basename $plasma .bam)_${i}.abra.bam.bai ] ; then samtools index ${outdirplasma}/abra/$(basename $plasma .bam)_${i}.abra.bam ; fi

fi

### BAM-READCOUNT ###
startreadcount=$(date +%s)
if [ $abra = 'True' ] ; 
then export bamfile=${outdirplasma}/abra/$(basename $plasma .bam)_${i}.abra.bam ;
else export bamfile=$plasma
fi
if [ ! -d $outdirplasma/bam-readcount ] ; then mkdir $outdirplasma/bam-readcount ; fi
if [ ! -f ${outdirplasma}/bam-readcount/$(basename $plasma .bam)_${i}.tsv ] ; then bam-readcount \
	-f ${extdata}/GRCh37/GRCh37.fa \
	$bamfile \
	-w 1 \
	-l ${extdata}/${mode}_bed/${mode}_hg19_chr${chr}_${i}.bed  >  ${outdirplasma}/bam-readcount/$(basename $plasma .bam)_${i}.tsv
fi
endreadcount=$(date +%s)
timereadcount=$(($endreadcount-$startreadcount))
hmsreadcount=$(printf '%02dh:%02dm:%02ds\n' $((timereadcount/3600)) $((timereadcount%3600/60)) $((timereadcount%60)))
echo "Elapsed Time READCOUNT on ${plasma} for chunk ${i} of chr${chr}: ${hmsreadcount}"
echo "Elapsed Time READCOUNT on ${plasma} for chunk ${i} of chr${chr}: ${hmsreadcount}" >> $outdirplasma/log/logtime_${i}.out

