

cd /mnt/data

for mix in mixtures/mix*/*/*/*.fastq.gz ; do
	echo $mix
	aws s3 cp $mix s3://cfdna-benchmark-dataset.store.genome.sg/${mix} --profile hanaecarrie
done

for init in intialsamples/ini*/*/*.fastq.gz ; do
        echo $init
        aws s3 cp $init s3://cfdna-benchmark-dataset.store.genome.sg/${init} --profile hanaecarrie
done

