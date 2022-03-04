#!/bin/bash

source /home/ubuntu/anaconda3/etc/profile.d/conda.sh
conda activate cfSNV

cd /home/ubuntu/cfSNV

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

#parse config file
eval $(parse_yaml $config_file)

echo $config_file
echo $outputdir
echo $tmpdir
echo $sampleid

if [ ! -d $outputdir ]; then mkdir $outputdir ; fi

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>$outputdir/log.out 2>&1
# Everything below will go to the file 'log.out':

echo 'Clean tmp folder'
rm -r $tmpdir

touch $outputdir/logtime.out

startcfsnv=$(date +%s)

Rscript run_cfsnv.R --config_file $config_file

endcfsnv=$(date +%s)
timecfsnv=$(($endcfsnv-$startcfsnv))
hmscfsnv=$(printf '%02dh:%02dm:%02ds\n' $((timecfsnv/3600)) $((timecfsnv%3600/60)) $((timecfsnv%60)))
        echo "Elapsed Time cfSNV on ${sampleid} for chr${chr}: ${hmscfsnv}"
        echo "Elapsed Time cfSNV on ${sampleid} for chr${chr}: ${hmscfsnv}" >> $outputdir/logtime.out
