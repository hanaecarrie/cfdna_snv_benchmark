#!/bin/bash

cd config

for i in {4..22}
do
	export oldconfigfile="config_sinvict_spikeins_chr3_CRC-COSMIC-5p_CRC-986_300316-CW-T.yaml"
	export configfile=$(echo "${oldconfigfile/chr3/"chr${i}"}")    
	echo $configfile
	cp $oldconfigfile $configfile
	export search="chr3"
	export replace="chr${i}"
	sed -i "s/$search/$replace/g" $configfile
	export search="_chr3_"
        export replace="_chr${i}_"
        sed -i "s/$search/$replace/g" $configfile
	export search="chr: 3"
        export replace="chr: ${i}"
        sed -i "s/$search/$replace/g" $configfile
done
