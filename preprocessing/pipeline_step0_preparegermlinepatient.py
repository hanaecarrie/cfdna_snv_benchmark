import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('patient', help='string describing the sample : patientid_date')
parser.add_argument('germline_vcf_name', help='prefix name like buffycoat_CRC-986_chr22_nofilter')
parser.add_argument('path_data', help='path to the data folder')

args = parser.parse_args()
print(args)

from utils import read_vcf

patient = str(args.patient)  # 986

if not os.path.exists(args.path_data+'/data/patient_SNPs/patient_'+patient+'_snps.csv'):
    # Read SNPs detected in cancer patient
    patient_snps_df = read_vcf(args.path_data+'/data/2015-07-31_'+args.germline_vcf_name+'/'+args.germline_vcf_name+'-gatk-haplotype-annotated.vcf')
    print(patient_snps_df.shape)
    patient_snps_df = patient_snps_df[patient_snps_df['#CHROM'] == '22']
    print(patient_snps_df.shape)
    foo_vaf = lambda x: pd.Series(x.split(';AF=')[1].split(';')[0])
    patient_snps_df['VAF'] = patient_snps_df['INFO'].apply(foo_vaf)
    patient_snps_df = patient_snps_df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'VAF']]

    patient_snps_df.to_csv(args.path_data+'/data/patient_SNPs/patient_'+patient+'_snps.csv', index=False)