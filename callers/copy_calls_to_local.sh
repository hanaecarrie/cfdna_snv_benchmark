
export dilutionseries=$1
echo aws s3 ls s3://cfdna-mutation-calling.store.genome.sg/$dilutionseries/ --profile hanae
aws s3 ls s3://cfdna-mutation-calling.store.genome.sg/$dilutionseries/ --profile hanae
echo aws s3 cp s3://cfdna-mutation-calling.store.genome.sg/$dilutionseries /Users/hanae/Repositories/cfdna_snv_benchmark/$dilutionseries --profile hanae --recursive --exclude '*.bam' 
aws s3 cp s3://cfdna-mutation-calling.store.genome.sg/$dilutionseries /Users/hanae/Repositories/cfdna_snv_benchmark/$dilutionseries --profile hanae --recursive --exclude '*.bam' 
