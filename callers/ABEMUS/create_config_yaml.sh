#!/bin/bash

cd config

for i in {2..22}
do
	export oldconfigfile="config_abemus_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T_WGS.yml"
	export configfile=$(echo "${oldconfigfile/chr1/"chr${i}"}")    
	echo $configfile
	cp $oldconfigfile $configfile
	export search="chr1"
	export replace="chr${i}"
	sed -i "s/$search/$replace/g" $configfile
	export search="_chr1_"
        export replace="_chr${i}_"
        sed -i "s/$search/$replace/g" $configfile
	export search="chr: 1"
        export replace="chr: ${i}"
        sed -i "s/$search/$replace/g" $configfile
done
