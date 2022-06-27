#!/bin/bash

cd config

for i in {1..22}
do
	export oldconfigfile="config_mixtures_chr3_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml"
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

