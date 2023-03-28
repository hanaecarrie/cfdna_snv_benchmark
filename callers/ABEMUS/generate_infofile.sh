#!/bin/bash

export dilutionseries_folder=$1
export buffycoat_bam=$2
export outdir=$3
export panelofnormalbcdir=$4
export chr=$5

echo $dilutionseries_folder
echo $outdir

if [ ! -d $outdir ] ; then mkdir $outdir ; fi
if [ ! -f $outdir/infofile.tsv ] ; then touch $outdir/infofile.tsv ; fi

export patientid=$($(basename $dilutionseries_folder) | cut -d 'dilutions_' -f1)'chr'$(dirname $dilutionseries_folder | cut -d 'dilutions_series_')
echo $patientid

for dil in $dilutionseries_folder/*.bam ; 
do echo $dil ; 	
	echo -e $patientid'\t'$(basename $dil .bam)'\t'$dil'\t'$(basename $buffycoat_bam .bam)'\t'$buffycoat_bam >> $outdir/infofile.tsv ;
done

# add other buffy coat samples to build global error-sequencing (GSE) distribution.
for bc in $panelofnormalbcdir/*chr${chr}.bam ;
	export patientid=$(echo $bc | cut -d 'chr' -f1)
	echo $patientid
	echo -e $patientid'\tNA\tNA\t'$(basename bc .bam)'\t'$bc >> $outdir/infofile.tsv
	
