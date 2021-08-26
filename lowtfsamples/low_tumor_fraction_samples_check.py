


import pandas as pd
pileup_df = pd.read_csv('../data/raw_data/NCC_CRC-986_300316-CW-T_pileup.txt', sep='\t', header=None)
pileup_df = pileup_df.iloc[1::2, :]

nmut = 0
cov_list = []
vaf_list = []
for prow in pileup_df.iterrows():
    #print(prow)
    pileup = prow[1][4]
    chrtocheck = ['a', 't', 'g', 'c', 'A', 'T', 'G', 'C']
    cov = len(pileup)
    if any(ext in pileup for ext in chrtocheck):
        nmut += 1
        for ext in chrtocheck:
            if ext in pileup:
                vaf = (pileup.count(ext.upper()) + pileup.count(ext.lower()))/len(pileup)
        vaf_list.append(vaf)
        cov_list.append(cov)
    else:
        vaf = 0


print(nmut)
print(nmut/pileup_df.shape[0])

import matplotlib.pyplot as plt
plt.figure(figsize=(20,15))
sns.scatterplot(cov_list, vaf_list, s=100)
plt.axhline(y=0, c='k')
plt.axhline(y=0.01, c='k', ls='--', label='VAF=1%, expected tumor fraction')
plt.axvline(x=7, c='k')
plt.xlabel('depth of coverage at given loci')
plt.ylabel('VAF')
plt.ylim([0, 0.25])
plt.legend()
plt.title('986_300316 - {}/{} ({}%) loci from 226 panel have variants'.format(nmut, pileup_df.shape[0], 100*nmut/pileup_df.shape[0]))


import pandas as pd
pileup_df = pd.read_csv('../data/raw_data/NCC_CRC-1014_090516-CW-T_pileup.txt', sep='\t', header=None)
pileup_df = pileup_df.iloc[1::2, :]

nmut = 0
cov_list = []
vaf_list = []
for prow in pileup_df.iterrows():
    #print(prow)
    pileup = prow[1][4]
    chrtocheck = ['a', 't', 'g', 'c', 'A', 'T', 'G', 'C']
    cov = len(pileup)
    if any(ext in pileup for ext in chrtocheck):
        nmut += 1
        for ext in chrtocheck:
            if ext in pileup:
                vaf = (pileup.count(ext.upper()) + pileup.count(ext.lower()))/len(pileup)
        vaf_list.append(vaf)
        cov_list.append(cov)
    else:
        vaf = 0

print(nmut)
print(nmut/pileup_df.shape[0])

import matplotlib.pyplot as plt
plt.figure(figsize=(20,15))
sns.scatterplot(cov_list, vaf_list, s=100)
plt.axhline(y=0, c='k')
plt.axhline(y=0.03, c='k', ls='--', label='VAF=3%, expected tumor fraction')
plt.axvline(x=7, c='k')
plt.xlabel('depth of coverage at given loci')
plt.ylabel('VAF')
plt.ylim([0, 0.25])
plt.legend()
plt.title('1014_090516 - {}/{} ({:.2f}%) from 226 panel have variants'.format(nmut, pileup_df.shape[0], 100*nmut/pileup_df.shape[0]))

print('PATIENT 1014')
seq = '.$...,..,..,,,,,,....,..,,,,,,,,,,,,,,.....,,.....A..,...,.,...,..,.,,.........,.,,,,,,......,.............,............,.....,........,.........,..,'
print('AKT1: VAF={:.3f}, cov={}'.format((seq.count('A'.upper())+seq.count('A'.lower()))/len(seq), len(seq)))
seq = '.$,.........,..,,,..,..,,......,,,.,..,.,,...,.....,-2tg,.,.....,.,.,.,...,.......,..,,..,.......,..,...,,.,,......,,,,...-2TG.......'
print('TP53: VAF={:.3f}, cov={}'.format((seq.count('-2tg'.upper())+seq.count('-2tg'.lower()))/len(seq), len(seq)))
seq = ',,.,,.,..,.,,,,..,........,,.,...t.......,...,,.....,....,....,,..,...........,.,....,...,.....,..,.T.........,.....,..,...,^].'
print('SOX9: VAF={:.3f}, cov={}'.format((seq.count('T'.upper())+seq.count('T'.lower()))/len(seq), len(seq)))
seq = ',$.$.$....,.,,,,....,,,.........,.....,.,.....................,....,.,....'
print('PIK3CA: VAF={:.3f}, cov={}'.format((seq.count('G'.upper())+seq.count('G'.lower()))/len(seq), len(seq)))
seq = ',$,$.,..,$.,...,.................,t..,,,,,,......,.,...,..,......,T.........,........,..,....,..............,.'
print('APC: VAF={:.3f}, cov={}'.format((seq.count('t'.upper())+seq.count('t'.lower()))/len(seq), len(seq)))

print('PATIENT 986')
seq = '.$,,,,.....,..,,..,,,...,.,.,..,,...,...,,...,.......,,.c,,,..,.....,...................,............,....,,,,..................,.,...................,....,,...............................,......,.................'
print('EPHB2: VAF={:.3f}, cov={}'.format((seq.count('g'.upper())+seq.count('c'.lower()))/len(seq), len(seq)))
seq = '.$.$,,,,,,..,,,,,,,,.,....,,.....,,.....,,.,,.,......,....,..,....,..,..................,................,.......,...,.......,..,......,..,,,..,......,,.,..,.,....,,...,...,,.......,.......,,.,....'
print('TP53: VAF={:.3f}, cov={}'.format((seq.count('c'.upper())+seq.count('c'.lower()))/len(seq), len(seq)))
seq = ',$,$.,,,,,,.,,..,..,,,,,,,..,,,,,....,......,...,......,,..,,,,,,,,,,.,,,,,,.,,,,,,.....,.,.....,....,,,,,,,,,,.,,,,..,,..............,..,.........................^].^].'
print('SOX9: VAF={:.3f}, cov={}'.format((seq.count('+'.upper())+seq.count('+'.lower()))/len(seq), len(seq)))
seq = '.$.,........................,....,...............,......,.,...............,...,..........,,,..'
print('PIK3CA: VAF={:.3f}, cov={}'.format((seq.count('G'.upper())+seq.count('G'.lower()))/len(seq), len(seq)))
seq ='...................,.,....,.....,,.........,.............,......,'
print('APC: VAF={:.3f}, cov={}'.format((seq.count('G'.upper())+seq.count('G'.lower()))/len(seq), len(seq)))

pileup_df = pd.read_csv('../data/raw_data/NCC_CRC-1014_090516-CW-T_pileup.txt', sep='\t', header=None)
pileup_APC_df = pileup_df.iloc[pileup_df[(pileup_df[3].str.contains('APC'))].index+1]
chrtocheck = ['a', 't', 'g', 'c', 'A', 'T', 'G', 'C']
for mpileup in pileup_APC_df.iterrows():
    pileup = mpileup[1][4]
    cov = len(pileup)
    if any(ext in pileup for ext in chrtocheck):
        for ext in chrtocheck:
            if ext in pileup:
                alt = (pileup.count(ext.upper()) + pileup.count(ext.lower()))
    else:
        alt = 0
    vaf = alt/cov
    print('KRAS', 'chr_'+ mpileup[1][0], 'pos_'+str(mpileup[1][1]), 'n_alt='+str(alt), 'cov='+str(cov), 'VAF={:.3f}'.format(vaf))

pileup_df = pd.read_csv('../data/raw_data/NCC_CRC-986_300316-CW-T_pileup.txt', sep='\t', header=None)
pileup_APC_df = pileup_df.iloc[pileup_df[(pileup_df[3].str.contains('KRAS'))].index+1]
chrtocheck = ['a', 't', 'g', 'c', 'A', 'T', 'G', 'C']
for mpileup in pileup_APC_df.iterrows():
    pileup = mpileup[1][4]
    cov = len(pileup)
    if any(ext in pileup for ext in chrtocheck):
        for ext in chrtocheck:
            if ext in pileup:
                alt = (pileup.count(ext.upper()) + pileup.count(ext.lower()))
    else:
        alt = 0
    vaf = alt/cov
    print('APC', 'chr_'+ mpileup[1][0], 'pos_'+str(mpileup[1][1]), 'n_alt='+str(alt), 'cov='+str(cov), 'VAF={:.3f}'.format(vaf))