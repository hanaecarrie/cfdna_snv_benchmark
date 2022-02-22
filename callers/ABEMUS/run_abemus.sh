#!/bin/bash

source /home/ubuntu/anaconda3/etc/profile.d/conda.sh
conda activate ABEMUS

cd /home/ubuntu/ABEMUS

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

echo $config_file
echo $dilutionseriesfolder
echo $buffycoatbam
echo $chr
echo $extdata
echo $outdir


####### generate info file #######
if [ ! -d $outdir ] ; then mkdir $outdir ; fi
cp $config_file $outdir/

if [ ! -f $outdir/infofile.tsv ] ; then touch $outdir/infofile.tsv ; fi

export patientid=$($(basename $dilutionseriesfolder) | cut -d 'dilutions_' -f1)'chr'$(dirname $dilutionseriesfolder | cut -d 'dilutions_series_')
echo $patientid

for dil in $dilutionseriesfolder/*.bam ;
do echo $dil ;
        echo -e $patientid'\t'$(basename $dil .bam)'\t'$dil'\t'$(basename $buffycoatbam .bam)'\t'$buffycoatbam >> $outdir/infofile.tsv ;
done


###### prepare bedfile ######
# Download reference genome (here all with hg19)
# Index it
# Create dictionary
# Create bed file chr
#export nbbasechr=$(cat $extdata/GRCh37/GRCh37.dict | grep SN:${chr} | awk 'BEGIN { FS = "LN:" } ; {print $2}' | awk 'BEGIN { FS = "\t" } ; {print $1}')
#echo $nbbasechr
if [ ! -d $extdata/wholegenome_bed ] ; then mkdir $extdata/wholegenome_bed ; fi
if [ ! -f $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed ] ;
then touch $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed ;
for ((i = 1; i <= $nbbasechr; i+=1000)) ;
do echo -e "$chr\t$i\t$(($i+999))" >> $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed ; done
fi
# Split bed file chr by chunks of 5,000 lines
if [ ! -f $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}_00.bed ] ; then split -l 5000 --numeric-suffixes --additional-suffix='.bed' $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}.bed $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}_ ; fi

###### prepare dbSNP database ######
# 1. Download it in /data/extdata folder
# $ cd $extdata/dbsnp_vcf
# $ wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_dbSNP_all.vcf.gz
# $ tabix GRCh37_latest_dbSNP_all.vcf.gz
# # to get correspondance refs and chromosome: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/
# $ bcftools annotate --rename-chrs /data/extdata/GRCh37/convert_chrom_names.txt -o dbSNP.vcf GRCh37_latest_dbSNP_all.vcf.gz
# $ bcftools view dbSNP.vcf -Oz -o dbSNP.vcf.gz
# $ bcftools index dbSNP.vcf.gz
# 2. Split per chromosome
# $ bcftools filter -r chr${chr} $extdata/dbsnp_vcf/dbSNP.vcf.gz -o $extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}.vcf
# 3. Edit vcf to keep only SNPs with one Alt base
if [ -f $extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}_edited.vcf ] ; then python edit_vcf.py $extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}.vcf ; fi


###### run ABEMUS per chunk in parallel ######
export nchunk=$(ls $extdata/wholegenome_bed/wholegenome_hg19_chr${chr}_*.bed | wc -l)
echo $nchunk
for ((n=00; i<$nchunk; n++)) ; do echo "run abemus on chunk ${n}" ; bash run_abemus_i.sh $config_file $n & ; done
bash run_abemus_i.sh $config_file $nchunk 
for n in $(seq -f "%02g" 0 $(($nchunk - 1))) ; do
	echo "run abemus on chunk ${n}"
        export npid=$((${n#0} + 1))
        bash /home/ubuntu/ABEMUS/run_abemus_i.sh -c $config_file -i $n  &  pids[${npid}]=$!
done
# wait for all pids
for pid in ${pids[*]}; do
	wait $pid
done

##### concat results #######
for dil in $dilutionseriesfolder/*.bam ;
do echo $i ;
mkdir   ${outdir}/abemus_outdir_chr${chr}/$(basename $dil .bam)
for j in {1..3} ;
        do ls ${outdir}/abemus_outdir_chr${chr}_*/Results/$(basename $dil .bam)/pmtab_F${j}_*.tsv
awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' ${outdir}/abemus_outdir_chr${chr}_*/Results/$(basename $dil .bam)/pmtab_F${j}_*.tsv  > ${outdir}/abemus_outdir_chr${chr}/$(basename $dil .bam)/pmtab_F${j}_$(basename $dil .bam).tsv
echo ${outdir}/abemus_outdir_chr${chr}/$(basename $dil .bam)/pmtab_F${j}_$(basename $dil .bam).tsv
done
done


