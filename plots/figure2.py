# Imports
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# set working directory
if not os.getcwd().endswith('cfdna_snv_benchmark'):
    os.chdir('../')
print('Current working directory: {}'.format(os.getcwd()))

from utils.config import Config
from utils.viz import set_display_params, function_to_split
from initialsamples.patient_timeline_analysis import plot_patient_timeline, get_mutations_stats
from initialsamples.pairedplots import paireplot

# Config and Display paramaters
config = Config("config/", "config_viz.yaml")
set_display_params(config)

# order of samples = 1) high tb sample 1, 2) high tb sample 2, 3) low tb sample
patientsample_dict = {
    '1014': ['NCC_CRC-1014_180816-CW-T', 'NCC_CRC-1014_110116-CW-T', 'NCC_CRC-1014_090516-CW-T'],
    '986': ['NCC_CRC-986_100215-CW-T', 'NCC_CRC-986_261016-CW-T', 'NCC_CRC-986_300316-CW-T'],
    '123': ['NCC_CRC-123_310715-CW-T', 'NCC_CRC-123_070116-CW-T', 'NCC_CRC-123_121115-CW-T']
}
patients = list(patientsample_dict.keys())
print(patients)

targetbedhg19 = pd.read_csv(os.path.join(*config.extdatafolder, 'Cancer226-targets_hg19.bed'), sep='\t')
targetbedhg38 = pd.read_csv(os.path.join(*config.extdatafolder, 'Cancer226-targets_hg38.bed'), sep='\t')

targetbedhg19list = []
for i, r in targetbedhg19.iterrows():
    ra = r[2] + 1 - r[1]
    rb = targetbedhg38.iloc[i, 2]+1-targetbedhg38.iloc[i, 1]
    if r[2] <= r[1]:
        raise ValueError('issue with {}'.format(r))
    for p in range(r[1], r[2]+1):
        targetbedhg19list.append('{}_{}'.format(r[0][3:], p))
    if ra != rb:
        print(ra, rb, i)
        for a in range(rb - ra):
            targetbedhg19list.append('{}_{}'.format(r[0][3:], p))
print("length of target bed file is {}".format(len(targetbedhg19list)))
print("example bed file locus for nomenclature {}".format(targetbedhg19list[0]))

targetbedhg38list = []
for _, r in targetbedhg38.iterrows():
    for p in range(r[1], r[2]+1):
        targetbedhg38list.append('{}_{}'.format(r[0], p))
print("length of target bed file is {}".format(len(targetbedhg38list)))
print("example bed file locus for nomenclature {}".format(targetbedhg38list[0]))

if not os.path.exists(os.path.join(*config.outputpath, 'figure2')):
    os.mkdir(os.path.join(*config.outputpath, 'figure2'))

#####################################################################################################################
# Figure 2 top: Identify elligible patients
#####################################################################################################################

for patient in patients:
    res = plot_patient_timeline(config, int(patient), figsize=(30, 8), mutations=True, highlight='discovery', treatment=False,
                                save=False, savepath=os.path.join(*config.outputpath, 'figure2'))
    lowtftimepoints_pd = get_mutations_stats(config, patient)
    lowtftimepoints_pd.dropna()

#####################################################################################################################
# Figure 2 bottom: Paired plots logscale and same scale
#####################################################################################################################

targetbedhg19 = pd.read_csv(os.path.join(*config.extdatafolder, 'Cancer226-targets_hg19.bed'), sep='\t')
targetbedhg38 = pd.read_csv(os.path.join(*config.extdatafolder, 'Cancer226-targets_hg38.bed'), sep='\t')

for patient in patients:
    for sc in ['samescale', 'logscale']:
        paireplot(config, patient, targetbedhg19, targetbedhg38, sc=sc, nbhightfsamples=1, save=False, savepath=None)

"""
# display params
markers = ['o', '^']
alphas = [1, 0.5]
fillstyles = ['full', 'none']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:grey']

for sc in ['logscale', 'samescale']:
    print(sc)
    for patient in patients:
        plt.figure(figsize=(10, 10))
        plt.grid()

        maxx = 0
        for si, seqtype in enumerate(['deepWGS', 'targeted']):

            mutationfile = pd.read_excel(os.path.join(*config.mutationfolder, [f for f in os.listdir(os.path.join(*config.mutationfolder)) if f.startswith('CCG_226_'+patient)][0]))
            mutationfile['chrom_pos'] = mutationfile['#CHROM'].astype('str').str.cat(mutationfile['POS'].astype('str'), sep="_")
            mutationfile.set_index('chrom_pos', inplace=True)
            mutationfile = mutationfile.loc[~mutationfile.index.duplicated()]
            for c in list(mutationfile.columns[13:-1]):
                if (sum([p.split('-')[1].split('_')[1] in c for p in patientsample_dict[patient]]) == 0) and (c != 'CCG_226_986.110215.P'):
                    mutationfile.drop(c, axis=1, inplace=True)
            ctokeep = [c for c in mutationfile.columns if c.startswith('CCG') or c.startswith(patient)]
            mutationfile[mutationfile['TIERS'] == 'Trusted'][['REF', 'ALT', 'GENE', *ctokeep]]

            print(seqtype)
            lowtfsample_vaf_df = pd.read_csv(os.path.join(*config.lowtfsamplesfolder, 'targeted', patientsample_dict[patient][2] + '_' + seqtype + '_vaf.txt'), index_col=0)
            lowtfsample_vaf_df['chrom_pos'] = lowtfsample_vaf_df['chrom'].astype('str').str.cat(lowtfsample_vaf_df['pos'].astype('str'), sep="_")
            lowtfsample_vaf_df.set_index('chrom_pos', inplace=True)
            lowtfsample_vaf_df = lowtfsample_vaf_df.loc[~lowtfsample_vaf_df.index.duplicated()]
            lowtfsample_vaf_df['annotation'] = 'No annotation'
            if seqtype == 'deepWGS':
                mutationfile.reset_index(inplace=True)
                for i in range(mutationfile.shape[0]):
                    if mutationfile['chrom_pos'].iloc[i] in targetbedhg38list:
                        val = targetbedhg19list[targetbedhg38list.index(mutationfile['chrom_pos'].iloc[i])]
                        mutationfile.iat[i, 0] = val
                    else:
                        mutationfile.iloc[i]['chrom_pos'] = np.nan
                mutationfile['chrom_pos'] = mutationfile['chrom_pos'].apply(lambda x: x.split('_')[0][:3] + '_' + x.split('_')[1]) ###
                mutationfile.set_index('chrom_pos', inplace=True)
                mutationfile = mutationfile.loc[~mutationfile.index.duplicated()]
            mutationfile = mutationfile.loc[~mutationfile.index.duplicated()]
            lowtfsample_vaf_df.loc[list(set(mutationfile.index) & set(lowtfsample_vaf_df.index)), 'annotation'] = mutationfile['TIERS']
            lowtfsample_vaf_df = lowtfsample_vaf_df[lowtfsample_vaf_df['annotation'] == 'Trusted'] #  != 'No annotation'
            lowtfsample_vaf_df['gene'] = mutationfile.loc[lowtfsample_vaf_df.index, 'GENE']
            lowtfsample_vaf_df['ref'] = mutationfile.loc[lowtfsample_vaf_df.index, 'REF']
            lowtfsample_vaf_df['alt'] = mutationfile.loc[lowtfsample_vaf_df.index, 'ALT']
            lowtfsample_vaf_df['pileup'].fillna('', inplace=True)
            lowtfsample_vaf_df['pileup'] = lowtfsample_vaf_df['pileup'].astype(str)
            lowtfsample_vaf_df['supporting altcov'] = lowtfsample_vaf_df[['pileup', 'alt', 'ref']].apply(lambda x: x[0].count(x[1]) if (len(x[1]) == 1) and (len(x[2]) == 1) else x[0].count('+'+str(len(x[1][1:]))+x[1][1:]) if len(x[1]) > 1 else x[0].count('-'+str(len(x[2][1:]))+x[2][1:]), axis=1)
            lowtfsample_vaf_df['supporting vaf'] = lowtfsample_vaf_df['supporting altcov']/lowtfsample_vaf_df['totcov']
            for pi in range(1):
                print(pi)
                hightfsample_vaf_df = pd.read_csv(os.path.join(*config.lowtfsamplesfolder, 'targeted', patientsample_dict[patient][pi] + '_' + seqtype + '_vaf.txt'), index_col=0)
                hightfsample_vaf_df['chrom_pos'] = hightfsample_vaf_df['chrom'].astype('str').str.cat(hightfsample_vaf_df['pos'].astype('str'), sep="_")
                hightfsample_vaf_df.set_index('chrom_pos', inplace=True)
                hightfsample_vaf_df = hightfsample_vaf_df.loc[~hightfsample_vaf_df.index.duplicated()]
                hightfsample_vaf_df['annotation'] = 'No annotation'
                if seqtype == 'deepWGS':
                    mutationfile.reset_index(inplace=True)
                    for i in range(mutationfile.shape[0]):
                        if mutationfile['chrom_pos'].iloc[i] in targetbedhg38list:
                            val = targetbedhg19list[targetbedhg38list.index(mutationfile['chrom_pos'].iloc[i])]
                            mutationfile.iat[i, 0] = val
                        else:
                            mutationfile.iloc[i]['chrom_pos'] = np.nan
                    mutationfile['chrom_pos'] = mutationfile['chrom_pos'].apply(lambda x: x.split('_')[0][:3] + '_' + x.split('_')[1]) ###
                    mutationfile.set_index('chrom_pos', inplace=True)
                    mutationfile = mutationfile.loc[~mutationfile.index.duplicated()]
                mutationfile = mutationfile.loc[~mutationfile.index.duplicated()]
                hightfsample_vaf_df.loc[list(set(mutationfile.index) & set(hightfsample_vaf_df.index)), 'annotation'] = mutationfile['TIERS']
                hightfsample_vaf_df = hightfsample_vaf_df[hightfsample_vaf_df['annotation'] == 'Trusted'] #  != 'No annotation'
                hightfsample_vaf_df['gene'] = mutationfile.loc[hightfsample_vaf_df.index, 'GENE']
                hightfsample_vaf_df['ref'] = mutationfile.loc[hightfsample_vaf_df.index, 'REF']
                hightfsample_vaf_df['alt'] = mutationfile.loc[hightfsample_vaf_df.index, 'ALT']
                hightfsample_vaf_df['pileup'].fillna('', inplace=True)
                hightfsample_vaf_df['pileup'] = hightfsample_vaf_df['pileup'].astype(str)
                hightfsample_vaf_df['supporting altcov'] = hightfsample_vaf_df[['pileup', 'alt', 'ref']].apply(lambda x: x[0].count(x[1]) if (len(x[1]) == 1) and (len(x[2]) == 1) else x[0].count('+'+str(len(x[1][1:]))+x[1][1:]) if len(x[1]) > 1 else x[0].count('-'+str(len(x[2][1:]))+x[2][1:]), axis=1)
                hightfsample_vaf_df['supporting vaf'] = hightfsample_vaf_df['supporting altcov']/hightfsample_vaf_df['totcov']

                paired_vaf = pd.concat([lowtfsample_vaf_df[lowtfsample_vaf_df['annotation'] == 'Trusted'][['supporting vaf']],
                                        hightfsample_vaf_df[hightfsample_vaf_df['annotation'] == 'Trusted'][['supporting vaf']]], axis=1)
                paired_vaf.columns = [patientsample_dict[patient][2].split('-')[1], patientsample_dict[patient][pi].split('-')[1]]
                maxx = paired_vaf[patientsample_dict[patient][2].split('-')[1]].max() if paired_vaf[patientsample_dict[patient][2].split('-')[1]].max() > maxx else maxx
                print(maxx)
                paired_vaf['gene'] = mutationfile.loc[paired_vaf.index, 'GENE']
                if patient == '986' and seqtype == 'targeted':
                    paired_vaf.iat[0, 0] = -0.00015  # (intead of 0 for APC, for visualisation)
                print(paired_vaf)
                print(patientsample_dict[patient][pi].split('-')[1])
                for gi, g in enumerate(list(paired_vaf['gene'])):
                    if g in config.genelist:
                        print(g, config.genelist.index(g))
                        c = sns.color_palette('tab10') + [sns.color_palette('Accent')[5]] + [sns.color_palette('tab20b')[14]] + [sns.color_palette('tab20b')[0]]
                        c = c[config.genelist.index(g)]
                    else:
                        c = colors[gi % len(colors)]
                    plt.plot(paired_vaf[paired_vaf['gene'] == g][patientsample_dict[patient][2].split('-')[1]].values[0],
                             paired_vaf[paired_vaf['gene'] == g][patientsample_dict[patient][pi].split('-')[1]].values[0], lw=5, markeredgewidth=1,
                             c=c, markersize=20, marker=markers[si], fillstyle=fillstyles[0], markeredgecolor='k', label=g + '*' + seqtype + ' sample')
        plt.ylim([-.02, 0.60])
        plt.axvline(x=0.000, ls='--', c='k')
        plt.title('Paired VAF plot for patient ' +patient)
        ax = plt.gca()
        hand, labl = ax.get_legend_handles_labels()
        hand, labl = function_to_split(hand, labl, '*')
        ax.legend(hand, labl, bbox_to_anchor=(1, 1), loc="upper left")
        plt.xlabel('VAF in ultra low TF sample')
        plt.ylabel('VAF in high TF sample')
        if sc == 'logscale':
            plt.xscale('symlog')
            plt.xlim([-0.002, min(0.05,round(maxx+0.01, 2))])
            if patient == '986':
                plt.xticks([i/200 for i in range(0, 2*int(100*(round(maxx+0.01, 2)+0.01)))], [i/200 for i in range(0, 2*int(100*(round(maxx+0.01, 2)+0.01)))])
            else:
                plt.xticks([i/100 for i in range(0, int(100*(round(maxx+0.01, 2)+0.01)))], [i/100 for i in range(0, int(100*(round(maxx+0.01, 2)+0.01)))])
            #plt.savefig(os.path.join(*config.outputpath, 'figure1c', 'pairedplot_vaf_trusted_'+patient+'_logscale.svg'), bbox_inches='tight')
            #plt.savefig(os.path.join(*config.outputpath, 'lowtbsamples', 'check', 'pairedplot_vaf_trusted_'+patient+'_logscale.svg'), bbox_inches='tight')
        elif sc == 'samescale':
            plt.xlim([-.02, 0.6])
            #plt.savefig(os.path.join(*config.outputpath, 'figure1c', 'pairedplot_vaf_trusted_'+patient+'_samescale.svg'), bbox_inches='tight')
        else:
            raise ValueError('Unknown sc parameter. should be samescale or logscale but here is {}'.format(sc))
        plt.show()
"""