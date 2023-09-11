#!/bin/bash

# Date: 2023
# Author: Hanae Carrie
# This script takes as input a configuration file for the chr1 of a mixture series to generate the configuration files
# of the following autosomal chromosomes from 2 to 22.
# It needs to be in the cfdna_snv_benchmark/callers/sinvict/ subfolder directly
# and requires the configuration file of chr1 as input parameter.

if [ $# == 0 ]; then
    echo "Usage: $0 param1"
    echo "* param1: initial configuration file name for chr1 used as template. chr1 needs to be present in file name."
    echo "Example:"
    echo "$ cd cfdna_snv_benchmark/callers/sinvict"
    echo "$ bash $0 config_sinvict_mixtures_chr1_CRC-986_100215-CW-T_CRC-986_300316-CW-T.yml"
    exit 1
fi

# Parse input to get initial configuration file
export templateconfigfile=$1
echo 'initial config file' $1

# Move to the cfdna_snv_benchmark/mixtures/sinvict/config subfolder
cd config

# Generate corresponding configuration file for each chromosomes from 2 to 22
for i in {2..22}
do
        export oldconfigfile=$templateconfigfile
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

