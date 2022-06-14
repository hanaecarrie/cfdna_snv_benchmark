#!/bin/bash

cd config

"""
for i in {4..22}
do
	export oldconfigfile="config_spikeins_chr3_CRC-COSMIC-5p_CRC-986_300316-CW-T.yml"
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

for i in {1..22}
do
        export oldconfigfile="config_spikeins_chr5_CRC-COSMIC-5p_CRC-1014_090516-CW-T.yml"
        export configfile=$(echo "${oldconfigfile/chr5/"chr${i}"}")
        echo $configfile
        cp $oldconfigfile $configfile
        export search="chr5"
        export replace="chr${i}"
        sed -i "s/$search/$replace/" $configfile
	export search="_chr5_"
        export replace="_chr${i}_"
        sed -i "s/$search/$replace/" $configfile
        export search="chr: 5"
        export replace="chr: ${i}"
        sed -i "s/$search/$replace/" $configfile
done

"""
for i in {1..21}
do
        export oldconfigfile="config_mixtures_chr3_CRC-123_310715-CW-T_CRC-123_121115-CW-T.yml"
        export configfile=$(echo "${oldconfigfile/chr3/"chr${i}"}")
        echo $configfile
        cp $oldconfigfile $configfile
        export search="chr22"
        export replace="chr${i}"
        sed -i "s/$search/$replace/" $configfile
	export search="_chr22_"
        export replace="_chr${i}_"
        sed -i "s/$search/$replace/" $configfile
        export search="chr: 22"
        export replace="chr: ${i}"
        sed -i "s/$search/$replace/" $configfile
done
"""
