import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('patient', help='string describing the sample : patientid_date')
parser.add_argument('snp_database', help='either dbsnp or genomead')
parser.add_argument('chunk_start', help='int of genomead 500,000 size chunk to start processing')
parser.add_argument('chunk_end', help='int of genomead 500,000 size chunk to end processing')
parser.add_argument('path_data', help='path to the data folder')

args = parser.parse_args()
print(args)

from supporting_reads import list_reads_to_remove, prepare_bamsurgeon_inputs

patient_date = args.patient  # '986_100215'
# load patient's SNPs detected on the deep WGS buffy coat sample with GATK Haplotype
patient_snps_df = pd.read_csv(args.path_data+'/data/patient_SNPs/patient_'+patient_date.split('_')[0]+'_snps.csv')

if args.snp_database == 'dbsnp':
    # load known SNPs database
    dbsnp_df = pd.read_csv(args.path_data+'/data/common_SNPs/dbsnp_df.csv')
    reads2remove, log_df = list_reads_to_remove(args.path_data+"/data/healthy_chr22_merged-ready.bam",
                                                dbsnp_df, patient_snps_df, verbose=-1)
    # save list of reads to remove and log dataframe
    log_df.to_csv(args.path_data+'/data/prepare_pooled_healthy/log_'+patient_date.split('_')[0]+'_dbsnp.csv', index=False)
    with open(args.path_data+'/data/prepare_pooled_healthy/readfile_'+patient_date.split('_')[0]+'_dbsnp.txt', "w") as output:
        for r in reads2remove:
            output.write(str(r) + "\n")
    bamsurgeon_snv_pd, bamsurgeon_indel_pd = prepare_bamsurgeon_inputs(patient_snps_df, log_df, max_vaf=0.1)
    # save file
    bamsurgeon_snv_pd.to_csv(args.path_data+'/data/prepare_pooled_healthy/varfile_snv_'+patient_date.split('_')[0]+'_dbsnp.bed',
                             sep='\t', header=False, index=False)
    bamsurgeon_indel_pd.to_csv(args.path_data+'/data/prepare_pooled_healthy/varfile_indel_'+patient_date.split('_')[0]+'_dbsnp.bed',
                               sep='\t', header=False, index=False)

elif args.snp_database == 'genomead':
    # load known SNPs database
    genomad_df_iterator = pd.read_csv(args.path_data+'/data/common_SNPs/genomad_df.csv', iterator=True, chunksize=500000)
    ci = 0
    for genomad_df_chunk in genomad_df_iterator:
        ci += 1
        if int(args.chunk_start) <= ci < int(args.chunk_end):
            print('chunk '+str(ci))
            genomad_df_chunk = genomad_df_chunk.drop('Unnamed: 0', axis=1)
            reads2remove, log_df = list_reads_to_remove(args.path_data+"/data/healthy_chr22_merged-ready.bam",
                                                        genomad_df_chunk, patient_snps_df, verbose=-1)
            # save list of reads to remove and log dataframe
            log_df.to_csv(args.path_data+'/data/prepare_pooled_healthy/log_'+patient_date.split('_')[0]+'_genomad_'+str(ci)+'.csv', index=False)
            with open(args.path_data+'/data/prepare_pooled_healthy/readfile_'+patient_date.split('_')[0]+'_genomad_'+str(ci)+'.txt', "w") as output:
                output.write(str(reads2remove))

            bamsurgeon_snv_pd, bamsurgeon_indel_pd = prepare_bamsurgeon_inputs(patient_snps_df, log_df, max_vaf=0.1)
            # save file
            bamsurgeon_snv_pd.to_csv(args.path_data+'/data/prepare_pooled_healthy/varfile_snv_'+patient_date.split('_')[0]+'genomad_'+str(ci)+'.bed',
                                     sep='\t', header=False, index=False)
            bamsurgeon_indel_pd.to_csv(args.path_data+'/data/prepare_pooled_healthy/varfile_indel_'+patient_date.split('_')[0]+'genomad_'+str(ci)+'.bed',
                                       sep='\t', header=False, index=False)
