import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def check_pileup(pileup_path, genes='all'):
    # pileup_df = pd.read_csv('../data/raw_data/NCC_CRC-986_300316-CW-T_pileup.txt', sep='\t', header=None)
    # pileup_df = pd.read_csv('../data/raw_data/NCC_CRC-1014_090516-CW-T_pileup.txt', sep='\t', header=None)
    pileup_df = pd.read_csv(pileup_path, sep='\t', header=None)
    if genes != 'all':  # name of 1 gene
        pileup_df = pileup_df.iloc[pileup_df[(pileup_df[3].str.contains(genes))].index+1]
    else:
        pileup_df = pileup_df.iloc[1::2, :]
    nmut = 0
    cov_list = []
    vaf_list = []
    for prow in pileup_df.iterrows():
        pileup = prow[1][4]
        chrtocheck = ['a', 't', 'g', 'c', 'A', 'T', 'G', 'C', '-', '+']
        cov = len(pileup)
        if any(ext in pileup for ext in chrtocheck):
            nmut += 1
            for ext in chrtocheck:
                if ext in pileup:
                    alt = (pileup.count(ext.upper()) + pileup.count(ext.lower()))
        else:
            alt = 0
        vaf = alt/cov
        vaf_list.append(vaf)
        cov_list.append(cov)
        if genes != 'all':
            print(genes, 'chr_'+ prow[1][0], 'pos_'+str(prow[1][1]), 'n_alt='+str(alt), 'cov='+str(cov), 'VAF={:.3f}'.format(vaf))
    print(nmut)
    print(nmut/pileup_df.shape[0])
    samplename = os.path.basename(pileup_path).split('.txt')[0]

    plt.figure(figsize=(20, 15))
    sns.scatterplot(cov_list, vaf_list, s=100)
    plt.axhline(y=0, c='k')
    plt.axhline(y=0.01, c='k', ls='--', label='VAF=1%')
    plt.axhline(y=0.03, c='k', ls='--', label='VAF=3%')
    plt.axhline(y=0.05, c='k', ls='--', label='VAF=5%')
    plt.axvline(x=7, c='k')
    plt.xlabel('depth of coverage at given loci')
    plt.ylabel('VAF')
    plt.ylim([0, 0.25])
    plt.legend()
    plt.title(samplename+' - {}/{} ({}%) loci from 226 panel have variants'.format(nmut, pileup_df.shape[0], 100*nmut/pileup_df.shape[0]))
    plt.show()


if __name__ == "__main__":
    print('PATIENT 986')
    pileup_path = '../data/raw_data/NCC_CRC-986_300316-CW-T_pileup.txt'
    check_pileup(pileup_path, genes='APC')
    check_pileup(pileup_path, genes='KRAS')
    check_pileup(pileup_path, genes='all')
    print('PATIENT 1014')
    pileup_path = '../data/raw_data/NCC_CRC-1014_090516-CW-T_pileup.txt'
    check_pileup(pileup_path, genes='APC')
    check_pileup(pileup_path, genes='KRAS')
    check_pileup(pileup_path, genes='all')
