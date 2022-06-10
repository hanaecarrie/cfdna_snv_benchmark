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

while getopts ":c:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

source $condapath
conda activate default

echo $sample_pseudohealthy
echo $sample_buffycoat
echo $samplename_pseudohealthy
echo $samplename_buffycoat
echo $chr
echo $outputfolder
echo $mutsbed_snv
echo $mutsbed_indel

export outputdir=$outputfolder/spikeins_chr${chr}_${samplename_knowntumormuts}_${samplename_pseudohealthy}
if [ ! -d $outputfolder ] ; then mkdir $outputfolder ; fi
echo $outputdir
echo $pseudohealthydir
echo $buffycoatdir

export sample_pseudohealthy_chr=$pseudohealthydir/${samplename_pseudohealthy}_chr${chr}.bam
export sample_buffycoat_chr=$buffycoatdir/${samplename_buffycoat}_chr${chr}.bam
echo $sample_pseudohealthy_chr

if [ ! -d $outputdir ] ; then mkdir $outputdir ; fi

cp $config_file $outputdir/

echo $dilutionfactor_pseudohealthy # 150x
export sample_pseudohealthy_chr_downsample=$pseudohealthydir/${samplename_pseudohealthy}_chr${chr}_${dilutionfactor_pseudohealthy}x.bam
echo $sample_pseudohealthy_chr_downsample

if [ ! -f $sample_pseudohealthy_chr_downsample ] ; then
# Select chr only
echo "Select chr ${chr} only for the pseudohealthy sample..."
if [ ! -d $pseudohealthydir ] ; then mkdir $pseudohealthydir ; fi
if [ ! -f $sample_pseudohealthy_chr ] ; then $samtools view -b $sample_pseudohealthy $chr > $sample_pseudohealthy_chr ; fi
if [ ! -f ${sample_pseudohealthy_chr}.bai ] ; then $samtools index $sample_pseudohealthy_chr ; fi

export pseudohealthy_chr_coverage=$pseudohealthydir/coverage_${samplename_pseudohealthy}_chr${chr}.txt
if  [ ! -f $pseudohealthy_chr_coverage ] ; then export pseudohealthy_cov=$($samtools depth -a $sample_pseudohealthy_chr | awk '{sum+=$3} END {print sum/NR}') ; echo $pseudohealthy_cov >  $pseudohealthy_chr_coverage ; else export pseudohealthy_cov=$(cat $pseudohealthy_chr_coverage) ; fi
echo $pseudohealthy_cov

fi

export dilutionfraction_pseudohealthy=$(echo "scale=4; $dilutionfactor_pseudohealthy / $pseudohealthy_cov" | bc -l)
echo $dilutionfraction_pseudohealthy
if [ $(echo "$dilutionfraction_pseudohealthy>1" | bc) == 1 ]; then dilutionfraction_pseudohealthy=1 ; fi
echo 'Dilution fraction pseudohealthy'
echo $dilutionfraction_pseudohealthy

echo "subsample pseudohealthy sample as desired..."

if [ $dilutionfraction_pseudohealthy != 1 ] ;
then if [ ! -f $sample_pseudohealthy_chr_downsample ] ; then $samtools view -b -s $dilutionfraction_pseudohealthy -o $sample_pseudohealthy_chr_downsample $sample_pseudohealthy_chr ; fi
else if [ ! -f $sample_pseudohealthy_chr_downsample ] ; then cp $sample_pseudohealthy_chr $sample_pseudohealthy_chr_downsample ; fi
fi
if [ ! -f ${sample_pseudohealthy_chr_downsample}.bai ] ; then $samtools index $sample_pseudohealthy_chr_downsample ; fi

for vaf in $vafs ;

do echo $vaf
if [ ! -d $outputfolder/spikeins_chr${chr}_${samplename_knowntumormuts}_${samplename_pseudohealthy}/logs ] ; then mkdir $outputfolder/spikeins_chr${chr}_${samplename_knowntumormuts}_${samplename_pseudohealthy}/logs ; fi
export outputdirnew=$outputfolder/spikeins_chr${chr}_${samplename_knowntumormuts}_${samplename_pseudohealthy}/spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}
if [ ! -d $outputdirnew ] ; then mkdir $outputdirnew ; fi
if [ ! -d $outputdirnew/logs ] ; then mkdir $outputdirnew/logs ; fi
if [ $machine == 'Aquila' ] ; then
	/opt/uge-8.1.7p3/bin/lx-amd64/qsub -pe OpenMP 4 -l mem_free=32G,h_rt=24:00:00 -o $outputdirnew/logs/ -e $outputdirnew/logs/ $repopath/cfdna_snv_benchmark/spikeins/create_spikeins_chr.sh -c $config_file -v $vaf
fi
if [ $machine == 'Ronin' ] ; then
        nohup $repopath/cfdna_snv_benchmark/spikeins/create_spikeins_chr.sh -c $config_file -v $vaf > $outputdir/logs/log_spikein_chr${chr}_vaf${vaf}.out 
fi
if [ $machine == 'Pluto' ] ; then
	if [ ! -d $outputdir/runs ] ; then mkdir $outputdir/runs ; fi
	export sampleid=spikein_chr${chr}_${samplename_knowntumormuts}_vaf${vaf}_${samplename_pseudohealthy}
	cp $defaultbash $outputdir/runs/run_${sampleid}.sh
	echo "bash $repopath/cfdna_snv_benchmark/spikeins/create_spikeins_chr.sh -c $config_file -v $vaf" >> $outputdir/runs/run_${sampleid}.sh  
	export jobid=s_${chr}_${vaf}
	export odir=${outputdir}/logs
	sed -i "s/jobid/${jobid}/g" $outputdir/runs/run_${sampleid}.sh
	sed -i "s/sampleid/${sampleid}/g" $outputdir/runs/run_${sampleid}.sh
      	sed -i "s!outputdir!${odir}!g" $outputdir/runs/run_${sampleid}.sh
	sbatch $outputdir/runs/run_${sampleid}.sh
fi

done

# check buffy coat select chr exists
#echo "buffy coat..."
#if [ ! -d $buffycoatdir ] ; then mkdir $buffycoatdir ; fi
#if [ ! -f $sample_buffycoat_chr ] ; then $samtools view -b $sample_buffycoat $chr > $sample_buffycoat_chr ; fi
#if [ ! -f ${sample_buffycoat_chr}.bai ] ; then $samtools index $sample_buffycoat_chr ; fi
