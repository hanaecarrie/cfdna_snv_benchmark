#!/bin/bash

cd config

for i in {2..22}
do
        export oldconfigfile="config_varnet_mixtures_chr1_CRC-123_310715-CW-T_CRC-123_121115-CW-T.yml"
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


"""
for i in {2..22}
do
	export oldconfigfile="config_varnet_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml"
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
"""
"""
for i in {2..22}
do
        export oldconfigfile="config_varnet_mixtures_chr1_CRC-1014_180816-CW-T_CRC-1014_090516-CW-T.yml"
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
"""
