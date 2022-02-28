import os
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt


def cosmictsv_to_bamsurgeonbed(extdatafolder, cancer_type, chrom, target='coding', threshold=5):
    if target == 'coding':
        cosmic_pd = pd.read_csv(os.path.join(extdatafolder, 'CosmicMutantExport.tsv'), encoding='ISO-8859-1', sep='\t', memory_map=True)
    else:  # target == 'noncoding'
        cosmic_pd = pd.read_csv(os.path.join(extdatafolder, 'CosmicNCV.tsv'), encoding='ISO-8859-1', sep='\t', memory_map=True)
    print(cosmic_pd.shape)
    # cancer type selection
    if cancer_type == 'CRC':
        if target == 'coding':
            cosmic_muts = cosmic_pd[(cosmic_pd['Primary site'] == 'large_intestine') & (cosmic_pd['Site subtype 1'] == 'colon') & (cosmic_pd['Tumour origin'] == 'primary')]
        else: # target == 'noncoding'
            cosmic_muts = cosmic_pd[(cosmic_pd['Primary site'] == 'large_intestine') & (cosmic_pd['Site subtype 1'] == 'colon')]
    else:
        raise ValueError('NEED TO INCLUDE THIS CANCER TYPE')
    print(cosmic_muts.shape)
    # remove nan positions in GRCh37
    cosmic_muts = cosmic_muts[~cosmic_muts["GRCh"].isna()]
    if len(cosmic_muts["GRCh"].unique()) > 1:
        raise ValueError('several reference human genome coordinates in database: GRCh{}'.format(cosmic_muts["GRCh"].unique()))
    print('reference human genome version used is GRCh{}'.format(str(cosmic_muts["GRCh"].unique()[0])))
    print(cosmic_muts.shape)
    # annotate mutation type
    cosmic_muts['type'] = np.nan
    if target == 'coding':
        print(cosmic_muts['Mutation Description'].unique())  # SNV or INDEL
        cosmic_muts['type'][cosmic_muts['Mutation Description'].str.contains('Substitution') == True] = 'SNV'
        cosmic_muts['type'][cosmic_muts['Mutation Description'].str.contains('Insertion') == True] = 'INS'
        cosmic_muts['type'][cosmic_muts['Mutation Description'].str.contains('Deletion') == True] = 'DEL'
        cosmic_muts.dropna(subset=['type'], inplace=True)
    else:  # target == 'noncoding'
        cosmic_muts['type'][cosmic_muts['MUT_SEQ'].str.len()-cosmic_muts['WT_SEQ'].str.len() == 0] = 'SNV'
        cosmic_muts['type'][cosmic_muts['MUT_SEQ'].str.len()-cosmic_muts['WT_SEQ'].str.len() > 0] = 'INS'
        cosmic_muts['type'][cosmic_muts['MUT_SEQ'].str.len()-cosmic_muts['WT_SEQ'].str.len() < 0] = 'DEL'
        cosmic_muts = cosmic_muts[~cosmic_muts['type'].isna()]
    print(cosmic_muts.shape)
    # encode mutation position
    if target == 'coding':
        mutposname = 'Mutation genome position'
    else:  # target == 'noncoding'
        mutposname = 'genome position'
    cosmic_muts['chrom'] = cosmic_muts[mutposname].str.split(':').str[0]
    cosmic_muts['chrom'][cosmic_muts['chrom'] == 23] = 'X'
    cosmic_muts['chrom'][cosmic_muts['chrom'] == 24] = 'Y'
    cosmic_muts['startpos'] = cosmic_muts[mutposname].str.split(':').str[1].str.split('-').str[0]
    cosmic_muts['endpos'] = cosmic_muts[mutposname].str.split(':').str[1].str.split('-').str[1]
    cosmic_muts['chrom_pos'] = cosmic_muts['chrom'] + '_' + cosmic_muts['startpos'] + '_' + cosmic_muts['endpos']
    cosmic_muts['startpos'] = cosmic_muts['startpos'].astype(int)
    cosmic_muts['endpos'] = cosmic_muts['endpos'].astype(int)
    cosmic_muts.dropna(subset=['chrom', 'startpos', 'endpos'], inplace=True)
    # get occurences
    occurences = cosmic_muts[mutposname].value_counts().values
    #plt.figure()
    #plt.hist(occurences[occurences <= 10])
    print(np.quantile(occurences, 0.25), np.mean(occurences), np.median(occurences), np.quantile(occurences, 0.75))
    occ = cosmic_muts['chrom_pos'].value_counts()
    common_muts_occ = occ[occ >= threshold]
    common_muts = list(common_muts_occ.index)
    print(len(common_muts))
    cosmic_muts = cosmic_muts[cosmic_muts['chrom_pos'].isin(common_muts)]
    print(cosmic_muts.shape)
    #plt.figure(figsize=(14, 10))
    #ax = sns.countplot(x='chrom', data=cosmic_muts[['chrom', 'startpos', 'endpos', 'chrom_pos']].drop_duplicates(), order=np.arange(1,25).astype(str))
    #plt.title('number of common {} mutations (SNV + INDEL) per chromosome found in at least {} patients'.format(cancer_type, threshold))
    #for p in ax.patches:
    #    ax.annotate('{}'.format(p.get_height()), (p.get_x(), p.get_height()))
    #plt.draw()
    cosmic_muts_chr = cosmic_muts[cosmic_muts['chrom'] == chrom]
    print(cosmic_muts_chr.shape)
    if target == 'coding':
        # info about nomenclature http://varnomen.hgvs.org/recommendations/DNA/
        cosmic_muts_chr['ref'] = np.nan
        cosmic_muts_chr['alt'] = np.nan
        cosmic_muts_chr['ref'][cosmic_muts_chr['type'] == 'SNV'] = cosmic_muts_chr[cosmic_muts_chr['type'] == 'SNV']['Mutation CDS'].str.split('>').str[0].str[-1]
        cosmic_muts_chr['alt'][cosmic_muts_chr['type'] == 'SNV'] = cosmic_muts_chr[cosmic_muts_chr['type'] == 'SNV']['Mutation CDS'].str.split('>').str[1]
        cosmic_muts_chr['alt'][cosmic_muts_chr['type'] == 'INS'] = cosmic_muts_chr[cosmic_muts_chr['type'] == 'INS']['Mutation CDS'].str.split('ins').str[1]
        for ri, row in cosmic_muts_chr.iterrows():
            if (row['type'] == 'INS') and (row['Mutation CDS'].endswith('dup')):
                chromrow, startpos, endpos = row['chrom_pos'].split('_')
                if chromrow.startswith('chr'):
                    chromrow = chromrow[3:]
                if chromrow == '23':
                    chromrow = 'X'
                elif chromrow == '24':
                    chromrow = 'Y'
                startpos, endpos = int(startpos), int(endpos)
                fasta = pysam.FastaFile(os.path.join('data', 'GRCh37', 'GRCh37.fa'))
                ref_seq = fasta.fetch(chromrow, startpos-1, startpos)
                fasta.close()
                cosmic_muts_chr.at[ri, 'ref'] = ref_seq
                cosmic_muts_chr.at[ri, 'alt'] = ref_seq
            # needs shift in deletions for bamsurgeon to indicate length of deletion
        cosmic_muts_chr.loc[cosmic_muts_chr['type'] == 'DEL', 'endpos'] += 1
    else:  # target noncoding
        cosmic_muts_chr['ref'] = cosmic_muts_chr['WT_SEQ']
        cosmic_muts_chr['alt'] = cosmic_muts_chr['MUT_SEQ']
    cosmic_bed_chr = cosmic_muts_chr[['chrom', 'startpos', 'endpos', 'alt', 'type']].reset_index(drop=True)
    cosmic_bed_chr.drop_duplicates(inplace=True, ignore_index=True)
    # check for errors
    # SNV
    print(sum(cosmic_bed_chr[cosmic_bed_chr['type'] == 'SNV']['endpos'] != cosmic_bed_chr[cosmic_bed_chr['type'] == 'SNV']['startpos']))
    # INS
    print(sum(cosmic_bed_chr[cosmic_bed_chr['type'] == 'INS']['endpos'] - cosmic_bed_chr[cosmic_bed_chr['type'] == 'INS']['startpos'] != cosmic_bed_chr[cosmic_bed_chr['type'] == 'INS']['alt'].str.len()))
    print(sum(cosmic_bed_chr[cosmic_bed_chr['type'] == 'INS']['endpos'] - cosmic_bed_chr[cosmic_bed_chr['type'] == 'INS']['startpos'] <= 0))
    # DEL
    print(sum(cosmic_bed_chr[cosmic_bed_chr['type'] == 'DEL']['startpos'] - cosmic_bed_chr[cosmic_bed_chr['type'] == 'DEL']['endpos'] >= 0))
    print(cosmic_bed_chr[cosmic_bed_chr['type'] == 'SNV'].shape[0])
    print(cosmic_bed_chr[(cosmic_bed_chr['type'] == 'INS') | (cosmic_bed_chr['type'] == 'DEL')].shape[0])
    # /!\ 1-based index for SNV bamsurgeon
    cosmic_bed_chr_snv = cosmic_bed_chr[cosmic_bed_chr['type'] == 'SNV'].drop('type', axis=1)
    cosmic_bed_chr_snv.insert(loc=3, column='vaf', value=1)  # Insert default VAF column
    # /!\ 0-based index for INDEL bamsurgeon
    cosmic_bed_chr_indel = cosmic_bed_chr[(cosmic_bed_chr['type'] == 'INS') | (cosmic_bed_chr['type'] == 'DEL')].drop('type', axis=1)
    cosmic_bed_chr_indel.insert(loc=3, column='vaf', value=1)  # Insert default VAF column
    cosmic_bed_chr_indel['startpos'] -= 1
    cosmic_bed_chr_indel['endpos'] -= 1
    return cosmic_bed_chr_snv, cosmic_bed_chr_indel


if __name__ == "__main__":
    from utils.config import Config

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")
    extdatafolder = os.path.join('data', 'extdata')
    cancer_type = 'CRC'

    threshold = 5

    for chrom in range(1, 25):
        chrom = str(chrom)
        print('#########')
        print(chrom)
        print('#########')
        if not os.path.exists(os.path.join('data', 'extdata', 'cosmic_mutations_atleast'+str(threshold)+'patients', cancer_type+'_chr'+chrom+'_SNV_tf1.bed')) \
                or not os.path.join('data', 'extdata', 'cosmic_mutations_atleast'+str(threshold)+'patients', cancer_type+'_chr'+chrom+'_INDEL_tf1.bed'):
            cosmic_bed_chr_snv_coding, cosmic_bed_chr_indel_coding = cosmictsv_to_bamsurgeonbed(extdatafolder, cancer_type, chrom, target='coding', threshold=threshold)
            cosmic_bed_chr_snv_noncoding, cosmic_bed_chr_indel_noncoding = cosmictsv_to_bamsurgeonbed(extdatafolder, cancer_type, chrom, target='noncoding', threshold=threshold)
            cosmic_bed_chr_snv = pd.concat([cosmic_bed_chr_snv_coding, cosmic_bed_chr_snv_noncoding], ignore_index=True).sort_values(by=['chrom', 'startpos', 'endpos'])
            cosmic_bed_chr_indel = pd.concat([cosmic_bed_chr_indel_coding, cosmic_bed_chr_indel_noncoding], ignore_index=True).sort_values(by=['chrom', 'startpos', 'endpos'])

            print(cosmic_bed_chr_snv.shape)
            print(cosmic_bed_chr_snv)
            print(cosmic_bed_chr_indel.shape)
            print(cosmic_bed_chr_indel)
            # save bed files
            if not os.path.exists(os.path.join('data', 'extdata')):
                os.mkdir(os.path.join('data', 'extdata'))
            if not os.path.exists(os.path.join('data', 'extdata', 'cosmic_mutations_atleast'+str(threshold)+'patients')):
                os.mkdir(os.path.join('data', 'extdata', 'cosmic_mutations_atleast'+str(threshold)+'patients'))
            cosmic_bed_chr_snv.to_csv(os.path.join('data', 'extdata', 'cosmic_mutations_atleast'+str(threshold)+'patients', cancer_type+'_chr'+chrom+'_SNV_tf1.bed'), sep='\t', header=False, index=False)
            cosmic_bed_chr_indel.to_csv(os.path.join('data', 'extdata', 'cosmic_mutations_atleast'+str(threshold)+'patients', cancer_type+'_chr'+chrom+'_INDEL_tf1.bed'), sep='\t', header=False, index=False)

    res_df = pd.DataFrame(columns=['chrom', 'startpos', 'endpos', 'vaf', 'alt',  'type'])
    for chrom in range(1, 25):
        chrom = str(chrom)
        snv_file = os.path.join('data', 'extdata', 'cosmic_mutations_atleast'+str(threshold)+'patients', cancer_type+'_chr'+chrom+'_SNV_tf1.bed')
        indel_file = os.path.join('data', 'extdata', 'cosmic_mutations_atleast'+str(threshold)+'patients', cancer_type+'_chr'+chrom+'_INDEL_tf1.bed')
        if not os.stat(snv_file).st_size == 0:
            aux_snv = pd.read_csv(snv_file, sep='\t', header=None)
            aux_snv.columns = ['chrom', 'startpos', 'endpos', 'vaf', 'alt']
            aux_snv['type'] = 'SNV'
        else:
            aux_snv = pd.DataFrame()
        if not os.stat(indel_file).st_size == 0:
            aux_indel = pd.read_csv(indel_file, sep='\t', header=None)
            aux_indel.columns = ['chrom', 'startpos', 'endpos', 'vaf', 'alt']
            aux_indel['type'] = 'INDEL'
        else:
            aux_indel = pd.DataFrame()
        res_df = pd.concat([res_df, aux_snv, aux_indel])
    print(res_df)
    df_plot = res_df.groupby(['chrom', 'type']).size().reset_index().pivot(columns='type', index='chrom', values=0)
    ax = df_plot.plot(kind='bar', stacked=True, figsize=(20, 7))
    for p in ax.patches:
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy()
        ax.text(x+width/2, y+height+25, int(height), horizontalalignment='center', verticalalignment='center')
    plt.title('Common {} mutations found in at least {} patients in the COSMIC database'.format(cancer_type, threshold))
    if not os.path.exists(os.path.join('figures', 'common_{}_mutations_COSMIC_{}patients.png'.format(cancer_type, threshold))):
        plt.savefig(os.path.join('figures', 'common_{}_mutations_COSMIC_{}patients.png'.format(cancer_type, threshold)))
    plt.show()



