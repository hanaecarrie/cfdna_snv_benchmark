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
#export nbbasechr=$(cat $extdata/GRCh37/GRCh37.dict | grep SN:${chr} | awk 'BEGIN { FS = "LN:" } ; {print $2}' | awk 'BEGIN { FS = "\t" } ; {print $1}')
#echo $nbbasechr
if [ ! -d $extdata/wholegenome_bed ] ; then mkdir $extdata/wholegenome_bed ; fi
if [ ! -f $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed ] ;
then touch $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed ;
for ((i = 1; i <= $nbbasechr; i+=1000)) ; 
do echo -e "$chr\t$i\t$(($i+999))" >> $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed ; done
fi
# Split bed file chr by chunks of 5,000 lines
if [ ! -f $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}_00.bed ] ; then split -l 5000 --numeric-suffixes --additional-suffix='.bed' $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}_ ; fi

###### run SINVICT per chunk in parallel ######
#export nchunk=$(ls $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}_*.bed | wc -l)
#echo $nchunk
echo "dilfolder ${dilutionseriesfolder}"
for plasma in ${dilutionseriesfolder}/*/*.bam ; do
	echo "plasma ${plasma}" ;
	export outdirplasma=$outdir/$(basename $plasma .bam)
        echo $outdirplasma
	if [ ! -d $outdirplasma ] ; then mkdir $outdirplasma ; fi
	for n in $(seq -f "%02g" 0 1) ; do #$nchunk
		echo "run sinvict on chunk ${n}"
		if [ "$n" -eq "01" ] ; then echo $n ; echo "bash /home/ubuntu/sinvict/run_sinvict_i.sh -c $config_file -i $n -p $plasma & "
		else echo "last chunk" ; echo "bash /home/ubuntu/sinvict/run_sinvict_i.sh -c $config_file -i $n -p $plasma"
		fi
	done
	### SINVICT ###
	#if [ ! -d $outdirplasma/results ] ; then mkdir $outdirplasma/results ; fi
	#/home/ubuntu/sinvict/sinvict -t ${outdirplasma}/bam-readcount -o ${outdirplasma}/results
done
