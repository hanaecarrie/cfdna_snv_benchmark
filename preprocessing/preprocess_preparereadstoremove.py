import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('patient', help='string describing the sample : patientid_date')
parser.add_argument('germline_vcf_name', help='prefix name like buffycoat_CRC-986_chr22_nofilter')
parser.add_argument('snp_database', help='either dbsnp or genomead')
parser.add_argument('chunk_start', help='int of genomead 500,000 size chunk to start processing')
parser.add_argument('chunk_end', help='int of genomead 500,000 size chunk to end processing')
parser.add_argument('path_data', help='path to the data folder')

args = parser.parse_args()
print(args)

from supporting_reads import list_reads_to_remove
from utils import read_vcf

patient = str(args.patient)  # 986
# load patient's SNPs detected on the deep WGS buffy coat sample with GATK Haplotype

if not os.path.exists(args.path_data+'/data/patient_SNPs/patient_'+patient+'_snps.csv'):
    # Read SNPs detected in cancer patient
    if os.path.exists(args.path_data+'/data/2015-07-31_'+args.germline_vcf_name+'/'+args.germline_vcf_name+'-gatk-haplotype-annotated.vcf'):
        patient_snps_df = read_vcf(args.path_data+'/data/2015-07-31_'+args.germline_vcf_name+'/'+args.germline_vcf_name+'-gatk-haplotype-annotated.vcf')
    else:
        patient_snps_df = read_vcf(args.path_data+'/data/2015-07-31_'+args.germline_vcf_name+'/'+args.germline_vcf_name+'-gatk-haplotype-annotated.vcf.gz')
    print(patient_snps_df.shape)
    patient_snps_df = patient_snps_df[patient_snps_df['#CHROM'] == '22']
    print(patient_snps_df.shape)
    foo_vaf = lambda x: pd.Series(x.split(';AF=')[1].split(';')[0])
    patient_snps_df['VAF'] = patient_snps_df['INFO'].apply(foo_vaf)
    patient_snps_df = patient_snps_df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'VAF']]

    patient_snps_df.to_csv(args.path_data+'/data/patient_SNPs/patient_'+patient+'_snps.csv', index=False)

patient_snps_df = pd.read_csv(args.path_data+'/data/patient_SNPs/patient_'+patient+'_snps.csv')

if not os.path.exists(args.path_data+'/data/prepare_pooled_healthy/'+patient):
    os.mkdir(args.path_data+'/data/prepare_pooled_healthy/'+patient)

if args.snp_database == 'dbsnp':
    # load known SNPs database
    dbsnp_df = pd.read_csv(args.path_data+'/data/common_SNPs/dbsnp_df.csv')
    reads2remove, log_df = list_reads_to_remove(args.path_data+"/data/healthy_chr22_merged-ready.bam",
                                                dbsnp_df, patient_snps_df, args.path_data+'/data/reference_genome/chr22.fa',
                                                verbose=-1)
    # save list of reads to remove and log dataframe
    log_df.to_csv(args.path_data+'/data/prepare_pooled_healthy/'+patient+'/log_'+patient+'_dbsnp.csv', index=False)
    with open(args.path_data+'/data/prepare_pooled_healthy/'+patient+'/readfile_'+patient+'_dbsnp.txt', "w") as output:
        for r in reads2remove:
            output.write(str(r) + "\n")

elif args.snp_database == 'gnomad':
    # load known SNPs database
    genomad_df_iterator = pd.read_csv(args.path_data+'/data/common_SNPs/gnomad_df.csv', iterator=True, chunksize=50000)
    ci = 0
    for genomad_df_chunk in genomad_df_iterator:
        ci += 1
        print(ci)
        if int(args.chunk_start) <= ci < int(args.chunk_end):
            print('chunk '+str(ci))
            genomad_df_chunk = genomad_df_chunk.drop('Unnamed: 0', axis=1)
            reads2remove, log_df = list_reads_to_remove(args.path_data+"/data/healthy_chr22_merged-ready.bam",
                                                        genomad_df_chunk, patient_snps_df, args.path_data+'/data/reference_genome/chr22.fa',
                                                        verbose=1)
            # save list of reads to remove and log dataframe
            log_df.to_csv(args.path_data+'/data/prepare_pooled_healthy/'+patient+'/log_'+patient+'_gnomad_'+str(ci)+'.csv', index=False)
            with open(args.path_data+'/data/prepare_pooled_healthy/'+patient+'/readfile_'+patient+'_gnomad_'+str(ci)+'.txt', "w") as output:
                for r in reads2remove:
                    output.write(str(r) + "\n")
