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

export bcbioconfigpath='/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/callers/bcbio/dilution.chr22.yaml'
export bcbiooutput='/mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/callers/bcbio/config/bcbio_config_chr22_mix_spik.tsv'

touch $bcbiooutput
export datadir='/mnt/projects/huangwt/wgs/hanae/data'

# mixtures
for mixturechromdir in $datadir/mixtures/* ;
	do echo $mixturechromdir 
	for mixtureseries in $mixturechromdir/mixtures* ;
		do echo $mixtureseries ; 
		if [ "$(ls $mixtureseries/*.yml | wc -l)" = "1" ] ; then export configfile=$(ls $mixtureseries/*.yml)
	        echo $configfile
        	# parse config file
        	eval $(parse_yaml $configfile)
        	echo $outputfolder
        	echo $sample_buffycoat
        	# buffy coat
        	export bcdir=$buffycoatdir
        	export bcfile=${samplename_buffycoat}_chr${chr}.bam
		for mixture in $mixtureseries/*/*x.bam ;
			do echo $mixture
			# dilution cfdna sample with pseudo matched healthy
			export mixdir=$(dirname $mixture)
			export mixfile=$(basename $mixture)
			echo -e "$(basename $mixfile .bam)\tsymlink::$mixdir\t$mixfile\tsymlink::$bcdir\t$bcfile\tMIX\tNA\tNA\t$bcbioconfigpath" >> $bcbiooutput
		
		done
		fi
	done
done


# spikeins
for spikeinchromdir in $datadir/spikeins/* ;
        do echo $spikeinchromdir
        for spikeinseries in $spikeinchromdir/spikeins* ;
                do echo $spikeinseries ;
		if [ "$(ls $spikeinseries/*.yml | wc -l)" = "1" ] ; then export configfile=$(ls $spikeinseries/*.yml)
	        echo $configfile
        	# parse config file
        	eval $(parse_yaml $configfile)
       		echo $outputfolder
        	echo $sample_buffycoat
        	# buffy coat
                export bcdir=$buffycoatdir
                export bcfile=${samplename_buffycoat}_chr${chr}.bam
		for spikein in $spikeinseries/*/*T.bam ;
                        do echo $spikein
                        # dilution cfdna sample with pseudo matched healthy
                        export spikedir=$(dirname $spikein)
                        export spikefile=$(basename $spikein)
                        echo -e "$(basename $spikefile .bam)\tsymlink::$spikedir\t$spikefile\tsymlink::$bcdir\t$bcfile\tMIX\tNA\tNA\t$bcbioconfigpath" >> $bcbiooutput
                done
        	fi
	done
done

