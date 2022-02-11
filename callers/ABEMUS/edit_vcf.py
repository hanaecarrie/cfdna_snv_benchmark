import pandas as pd
import numpy as np
from itertools import chain

vcf_df = pd.read_csv('homo_sapiens-chr22.vcf', skiprows=36, sep='\t', memory_map=True)
print(vcf_df.head())
print(vcf_df.shape)

# return list from series of comma-separated strings
def chainer(s):
    return list(chain.from_iterable(s.str.split(',')))

# calculate lengths of splits
lens = vcf_df['ALT'].str.split(',').map(len)

print(lens)

# create new dataframe, repeating or chaining as appropriate
res = pd.DataFrame({'#CHROM': np.repeat(vcf_df['#CHROM'], lens),
                    'POS': np.repeat(vcf_df['POS'], lens),
                    'ID': np.repeat(vcf_df['ID'], lens),
                    'REF': np.repeat(vcf_df['REF'], lens),
                    'ALT': chainer(vcf_df['ALT']),
                    'QUAL': np.repeat(vcf_df['QUAL'], lens),
                    'FILTER': np.repeat(vcf_df['FILTER'], lens),
                    'INFO': np.repeat(vcf_df['INFO'], lens)})

res.drop_duplicates(subset=['#CHROM', 'POS'], inplace=True)
res = res[(res['REF'] == 'A') | (res['REF'] == 'C') | (res['REF'] == 'G') | (res['REF'] == 'T')]
res = res[(res['ALT'] == 'A') | (res['ALT'] == 'C') | (res['ALT'] == 'G') | (res['ALT'] == 'T')]

print(res)
res.to_csv('homo_sapiens-chr22_edited.vcf', sep='\t', header=False, index=False)


