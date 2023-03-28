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


while getopts ":c:i:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
    i) i=${OPTARG} ;;
  esac
done

# parse config file
eval $(parse_yaml $config_file)

source $condapath
conda activate ABEMUS

cd ${repopath}/ABEMUS

echo $config_file
echo $dilutionseriesfolder
echo $buffycoatbam
echo $chr
echo $extdata
echo $outdir
echo $mode

echo $i

if [ ! -d $outdir ] ; then mkdir $outdir ; fi ; 
if [ ! -d $outdir/chunk_${i} ] ; then mkdir $outdir/chunk_${i} ; fi ;
if [ ! -d $outdir/chunk_${i}/PaCBAM_outdir ] ; then mkdir $outdir/chunk_${i}/PaCBAM_outdir ; fi

# PaCBAM buffy coat file
if [ ! -f $extdata/${mode}_bed/${mode}_hg19_merge_chr${chr}_${i}.bed ] ; then bedtools merge -i $extdata/${mode}_bed/${mode}_hg19_chr${chr}_${i}.bed > $extdata/${mode}_bed/${mode}_hg19_merge_chr${chr}_${i}.bed ; fi
~/bin/pacbam/pacbam bam=$buffycoatbam  bed=$extdata/${mode}_bed/${mode}_hg19_merge_chr${chr}_${i}.bed vcf=$extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}_edited.vcf fasta=$extdata/GRCh37/GRCh37.fa strandbias mode=5 out=$outdir/chunk_${i}/PaCBAM_outdir ;
# PaCBAM plasma files
for bamfile in ${dilutionseriesfolder}/*/*.bam; do 
	echo $bamfile ; 
	~/bin/pacbam/pacbam bam=$bamfile bed=$extdata/${mode}_bed/${mode}_hg19_merge_chr${chr}_${i}.bed vcf=$extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}_edited.vcf fasta=$extdata/GRCh37/GRCh37.fa strandbias mode=5 out=$outdir/chunk_${i}/PaCBAM_outdir ;
done
# PaCBAM control buffycoat files
for bamfile in ${panelofnormalbcdir}/PoNbuffycoat_*/*.bam; do
        echo $bamfile ;
        ~/bin/pacbam/pacbam bam=$bamfile bed=$extdata/${mode}_bed/${mode}_hg19_merge_chr${chr}_${i}.bed vcf=$extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}_edited.vcf fasta=$extdata/GRCh37/GRCh37.fa strandbias mode=5 out=$outdir/chunk_${i}/PaCBAM_outdir ;
done

if [ ! -d $outdir/chunk_${i}/pacbam_data_bychrom ] ; then mkdir $outdir/chunk_${i}/pacbam_data_bychrom ; fi

Rscript run_abemus.R "${outdir}/chunk_${i}/" "${outdir}/infofile.tsv" "${extdata}/${mode}_bed/${mode}_hg19_merge_chr${chr}_${i}.bed" "${outdir}/chunk_${i}/PaCBAM_outdir" "${outdir}/chunk_${i}/pacbam_data_bychrom"

