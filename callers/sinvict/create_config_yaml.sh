#!/bin/bash

cd config

for i in {2..22}
do
        export oldconfigfile="config_sinvict_mixtures_chr1_BRA-412_240820-CW-T_BRA-412_060220-CW-T.yml"
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
        export oldconfigfile="config_sinvict_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T_WGS.yml"
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

for i in {2..22}
do
	export oldconfigfile="config_sinvict_mixtures_chr1_CRC-1014_180816-CW-T_CRC-1014_090516-CW-T.yml"
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

for i in {4..21}
do
        export oldconfigfile="config_mixtures_chr3_CRC-1014_180816-CW-T_CRC-1014_090516-CW-T.yml"
        export configfile=$(echo "${oldconfigfile/chr3/"chr${i}"}")
        echo $configfile
        cp $oldconfigfile $configfile
        export search="chr3"
        export replace="chr${i}"
        sed -i "s/$search/$replace/" $configfile
	export search="_chr3_"
        export replace="_chr${i}_"
        sed -i "s/$search/$replace/" $configfile
        export search="chr: 3"
        export replace="chr: ${i}"
        sed -i "s/$search/$replace/" $configfile
done


for i in {4..21}
do
        export oldconfigfile="config_mixtures_chr3_CRC-123_310715-CW-T_CRC-123_121115-CW-T.yml"
        export configfile=$(echo "${oldconfigfile/chr3/"chr${i}"}")
        echo $configfile
        cp $oldconfigfile $configfile
        export search="chr3"
        export replace="chr${i}"
        sed -i "s/$search/$replace/" $configfile
	export search="_chr3_"
        export replace="_chr${i}_"
        sed -i "s/$search/$replace/" $configfile
        export search="chr: 3"
        export replace="chr: ${i}"
        sed -i "s/$search/$replace/" $configfile
done
"""
