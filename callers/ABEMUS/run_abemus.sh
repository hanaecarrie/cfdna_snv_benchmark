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

while getopts ":c:" opt; do
  case $opt in
    c) config_file=${OPTARG} ;;
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
echo $finaloutdir
echo $mode

####### generate info file #######
if [ ! -d $outdir ] ; then mkdir $outdir ; fi
cp $config_file $outdir/

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>$outdir/log.out 2>&1
# Everything below will go to the file 'log.out'

# Start logging the RAM usage and CPU usage
bash ${repopath}/log_mem_cpu.sh $outdir/log_mem_cpu.out  & export logpid=$!

if [ ! -f $outdir/infofile.tsv ] ; then touch $outdir/infofile.tsv ; fi

export patientid=$(echo $(basename $dilutionseriesfolder) | cut -d '_' -f2-3)
echo $patientid

echo $dilutionseriesfolder
for dil in ${dilutionseriesfolder}/*/*[Tx].bam ; do
	echo $dil ;
        echo -e $patientid'\t'$(basename $dil .bam)'\t'$dil'\t'$(basename $buffycoatbam .bam)'\t'$buffycoatbam >> $outdir/infofile.tsv ;
done
# add list of buffycoat to build GSE distribution
echo $panelofnormalbcdir
for bc in ${panelofnormalbcdir}/PoNbuffycoat_chr${chr}/*chr${chr}.bam ; do
        echo $bc ;	
	echo -e $(echo $(basename $bc .bam) | cut -d'-' -f2)'\t'NA'\t'NA'\t'$(basename $bc .bam)'\t'$bc >> $outdir/infofile.tsv
done

###### prepare bedfile ######
# Download reference genome (here all with hg19)
# Index it
# Create dictionary
# Create bed file chr
# Split bed file chr by chunks of 5,000 lines
if [ ! -f $extdata/${mode}_bed/${mode}_hg19_chr${chr}_00.bed ] ; then split -l 5000 --numeric-suffixes --additional-suffix='.bed' $extdata/${mode}_bed/${mode}_hg19_chr${chr}.bed $extdata/${mode}_bed/${mode}_hg19_chr${chr}_ ; fi

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
if [ ! -f $extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}.vcf ] ; then bcftools filter -r ${chr} $extdata/dbsnp_vcf/dbSNP.vcf.gz -o $extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}.vcf ; fi
# 3. Edit vcf to keep only SNPs with one Alt base
echo $extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}_edited.vcf
if [ ! -f $extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}_edited.vcf ] ; then python edit_vcf.py $extdata/dbsnp_vcf/dbSNP_hg19_chr${chr}.vcf ; fi

###### run ABEMUS per chunk in parallel ######
export nchunk=$(ls $extdata/exome_bed/exome_hg19_chr${chr}_*.bed | wc -l)
echo $nchunk
for n in $(seq -f "%02g" 0 $(($nchunk - 1))) ; do 
	echo "run abemus on chunk ${n}" ; 
	export npid=$((${n#0} + 1))
	bash ${repopath}/ABEMUS/run_abemus_i.sh -c $config_file -i $n  &  pids[${npid}]=$!
done

# wait for all pids
for pid in ${pids[*]}; do
	wait $pid
done

if [ ! -d ${outdir}/results ] ; then mkdir ${outdir}/results ; fi

###### concat results #####
python ${repopath}/ABEMUS/concat_results.py $dilutionseriesfolder $outdir

####### copy results to common folder ########
echo $finaloutdir
if [ ! -d $finaloutdir ] ; then mkdir $finaloutdir ; fi

for plasma in ${dilutionseriesfolder}/*/*[Tx].bam ; 
	do echo "plasma ${plasma}" ;
	if [ ! -d $finaloutdir/$(basename $plasma .bam) ] ; then mkdir $finaloutdir/$(basename $plasma .bam) ; fi
	if [ ! -d $finaloutdir/$(basename $plasma .bam)/abemus ] ; then mkdir $finaloutdir/$(basename $plasma .bam)/abemus ; fi
	scp $outdir/results/$(basename $plasma .bam)/* $finaloutdir/$(basename $plasma .bam)/abemus/
done

kill -9 $logpid
