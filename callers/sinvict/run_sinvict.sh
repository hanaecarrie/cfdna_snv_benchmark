#!/bin/bash

cd /home/ubuntu/sinvict

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

while getopts ":c:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

echo $config_file
echo $dilutionseriesfolder
echo $buffycoatbam
echo $chr
echo $extdata
echo $outdir

####### generate info file #######
if [ ! -d $outdir ] ; then mkdir $outdir ; fi
cp $config_file $outdir/

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>$outdir/log.out 2>&1
# Everything below will go to the file 'log.out'

###### prepare bedfile ######
# Download reference genome (here all with hg19)
# Index it
# Create dictionary
# Create bed file chr
export nbbasechr=$(cat $extdata/GRCh37/GRCh37.dict | grep SN:${chr} | awk 'BEGIN { FS = "LN:" } ; {print $2}' | awk 'BEGIN { FS = "\t" } ; {print $1}')
echo $nbbasechr
if [ ! -d $extdata/wholegenome_bed ] ; then mkdir $extdata/wholegenome_bed ; fi
if [ ! -f $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed ] ;
then touch $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed ;
for ((i = 1; i <= $nbbasechr; i+=1000)) ; 
do echo -e "$chr\t$i\t$(($i+999))" >> $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed ; done
fi
# Split bed file chr by chunks of 5,000 lines
if [ ! -f $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}_00.bed ] ; then split -l 5000 --numeric-suffixes --additional-suffix='.bed' $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}_ ; fi

###### run SINVICT per chunk in parallel ######
export nchunk=$(ls $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}_*.bed | wc -l)
echo $nchunk
echo "dilfolder ${dilutionseriesfolder}"
for plasma in ${dilutionseriesfolder}/*/*.bam ; do
	echo "plasma ${plasma}" ;
	export outdirplasma=$outdir/$(basename $plasma .bam)
        echo $outdirplasma
	if [ ! -d $outdirplasma ] ; then mkdir $outdirplasma ; fi
	
	#if [ ! -f $outdirplasma/log/logtime_${i}.out ] ; then touch $outdirplasma/log/logtime_${i}.out ; fi
	#### ABRA ###
	#if [ ! -d $outdirplasma/abra ] ; then mkdir $outdirplasma/abra ; fi
	#if [ ! -d $outdirplasma/tmp ] ; then mkdir $outdirplasma/tmp ; fi
	#startabra=$(date +%s)
	#if [ ! -f ${outdirplasma}/abra/$(basename $plasma .bam).abra.bam ] ;
        #then java -Xmx16G -jar /home/ubuntu/bin/abra2/target/abra2-2.24-jar-with-dependencies.jar \
        #--in $buffycoatbam,$plasma \
        #--out ${outdirplasma}/abra/$(basename $buffycoatbam .bam).abra.bam,${outdirplasma}/abra/$(basename $plasma .bam).abra.bam \
        #--ref ${extdata}/GRCh37/GRCh37.fa --threads 8 --targets $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed --tmpdir ${outdirplasma}/tmp/ > ${outdirplasma}/abra/abra.log
	#fi
	#endabra=$(date +%s)
	#timeabra=$(($endabra-$startabra))
	#hmsabra=$(printf '%02dh:%02dm:%02ds\n' $((timeabra/3600)) $((timeabra%3600/60)) $((timeabra%60)))
	#echo "Elapsed Time ABRA on ${plasma} for chr${chr}: ${hmsabra}"
	#echo "Elapsed Time ABRA on ${plasma} for chr${chr}: ${hmsabra}" >> $outdirplasma/log/logtime_${i}.out
	#if [ ! -f ${outdirplasma}/abra/$(basename $plasma .bam).abra.bam.bai ] ; then /usr/bin/samtools index ${outdirplasma}/abra/$(basename $plasma .bam).abra.bam ; fi

	for n in $(seq -f "%02g" 0 $(($nchunk - 1))) ; do 
		echo "run sinvict on chunk ${n}"
		export npid=$((${n#0} + 1))
	        bash /home/ubuntu/sinvict/run_sinvict_i.sh -c $config_file -i $n -p $plasma &  pids[${npid}]=$!
	done
	# wait for all pids
	for pid in ${pids[*]}; do
		wait $pid
	done
	### SINVICT ###
	if [ ! -d $outdirplasma/results ] ; then mkdir $outdirplasma/results ; fi
	/home/ubuntu/sinvict/sinvict -t ${outdirplasma}/bam-readcount -o ${outdirplasma}/results
done
