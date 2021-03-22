import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('patient', help='string describing the sample : patientid_date')
parser.add_argument('path_data', help='path to the data folder')

args = parser.parse_args()
print(args)

from supporting_reads import prepare_bamsurgeon_inputs

patient_date = args.patient  # '986_100215'
# load patient's SNPs detected on the deep WGS buffy coat sample with GATK Haplotype
patient_snps_df = pd.read_csv(args.path_data+'/data/patient_SNPs/patient_'+patient_date.split('_')[0]+'_snps.csv')
# get log file concatenate
if not os.path.exists(args.path_data+'/data/prepare_pooled_healthy/log_'+patient_date.split('_')[0]+'_total.csv'):
    log_df = pd.read_csv(args.path_data+'/data/prepare_pooled_healthy/log_'+patient_date.split('_')[0]+'_dbsnp.csv')
    for log_path_genomad in [f for f in os.listdir(args.path_data+'/data/prepare_pooled_healthy/') if f.startswith('log_')]:
        log_df_genomad = pd.read_csv(args.path_data+'/data/prepare_pooled_healthy/'+log_path_genomad)
        log_df = pd.concat([log_df, log_df_genomad])
    log_df.to_csv(args.path_data+'/data/prepare_pooled_healthy/log_'+patient_date.split('_')[0]+'_total.csv', index=False)
log_df = pd.read_csv(args.path_data+'/data/prepare_pooled_healthy/log_'+patient_date.split('_')[0]+'_total.csv')

# prepare bamsurgeon inputs
bamsurgeon_snv_pd, bamsurgeon_indel_pd = prepare_bamsurgeon_inputs(patient_snps_df, log_df, max_vaf=0.1)
# save file
bamsurgeon_snv_pd.to_csv(args.path_data+'/data/prepare_pooled_healthy/varfile_snv_'+patient_date.split('_')[0]+'_total.bed',
                         sep='\t', header=False, index=False)
bamsurgeon_indel_pd.to_csv(args.path_data+'/data/prepare_pooled_healthy/varfile_indel_'+patient_date.split('_')[0]+'_total.bed',
                           sep='\t', header=False, index=False)