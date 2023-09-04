import os
import numpy as np
import pandas as pd

from utils.calltable import read_vcf


def parse_caller_feature(calldir, method, save=False):
    # path to calls
    sampleid = os.path.basename(calldir)
    print(sampleid)
    print(method)
    if method in ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']:
        calltablemethod_path = os.path.join(calldir, 'calls', 'bcbio', sampleid+'-'+method+'-annotated.vcf.gz')
        if '.' in os.path.basename(calltablemethod_path[:-7]):
            calltablemethod_path = os.path.join(os.path.dirname(calltablemethod_path), os.path.basename(calltablemethod_path)[:-7].replace('.', '_') + '.vcf.gz')
        if not os.path.exists(calltablemethod_path):
            print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
            callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
        else:
            callmethod = read_vcf(calltablemethod_path)
            callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'FORMAT', sampleid.replace('.', '_')+'-T', 'ID', 'INFO']]
            callmethod.columns = ['chrom', 'pos', 'ref', 'alt', method, 'format', 'formatvalue', 'ID', 'INFO']
            #print(callmethod.shape)
            # callmethod[method].fillna('', inplace=True)
            callmethod[method] = callmethod[method].astype(str)
            if method == 'mutect2' or method == 'strelka2':
                # print(callmethod.loc[callmethod[method] == 'MinAF', method].value_counts())
                print('retrieving {} {} calls with MinAF tags out of {}'.format(method,
                                                                                callmethod[callmethod[method] =='MinAF'].shape[0], callmethod[callmethod[method] == "PASS"].shape[0]))
                callmethod.loc[callmethod[method] == 'MinAF', method] = True
            elif method == 'vardict':
                #print(callmethod.loc[(callmethod[method] == 'f0.01'), method].value_counts())
                print('retrieving {} {} calls with f0.01;REJECT;REJECT tags out of {}'.format(method,
                                                                                              callmethod[(callmethod[method] == 'f0.01;REJECT;REJECT')].shape[0], callmethod[callmethod[method] == "PASS"].shape[0]))
                callmethod.loc[(callmethod[method] == 'f0.01;REJECT;REJECT'), method] = True
        #print(callmethod.shape)
        #callmethod.loc[callmethod[method] == 'PASS', method] = True
        #print(callmethod.shape)
        #callmethod = callmethod[callmethod[method] == True]
        #print(callmethod.shape)
        info = callmethod['INFO']
        callmethod.drop('INFO', axis=1, inplace=True)
        if method == 'freebayes':  # logodds to probability score prob = exp(logODDS)/(1+exp(logODDS))
            callmethod[method + '_score'] = [np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0]))) / (1 + np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0]))))
                                             if 'ODDS' in i else np.nan for i in info.to_list()]
        elif method == 'mutect2':  # logodds to probability score prob = exp(TLOD)/(1+exp(TLOD))
            callmethod[method + '_score'] = [np.exp(float(i.split('TLOD=')[1].split(';')[0])) / (1 + np.exp(float(i.split('TLOD=')[1].split(';')[0])))
                                             if 'TLOD' in i and ',' not in i.split('TLOD=')[1].split(';')[0]
                                             else np.exp(float(i.split('TLOD=')[1].split(';')[0].split(',')[0])) / (1 + np.exp(float(i.split('TLOD=')[1].split(';')[0].split(',')[0])))
            if 'TLOD' in i and ',' in i.split('TLOD=')[1].split(';')[0]
            else np.nan for i in info.to_list()]
        if method == 'strelka2':  # phred score to probability, prob = 1-10^(-SomaticEVS/10)
            callmethod[method + '_score'] = [1-(10 ** (-float(i.split('SomaticEVS=')[1].split(';')[0]) / 10))
                                             if 'SomaticEVS' in i else np.nan for i in info.to_list()]
        if method == 'vardict':  # P-value /!\ opposite direction of score variation /!\
            callmethod[method + '_score'] = [1-float(i.split('SSF=')[1].split(';')[0])
                                             if 'SSF' in i else np.nan for i in info.tolist()]
        if method == 'varscan':  # P-value  /!\ opposite direction of score variation /!\  # [1-float(i.split('SSC=')[1].split(';')[0])/255
            callmethod[method + '_score'] = [1-float(i.split('SPV=')[1].split(';')[0])
                                             if 'SPV' in i else np.nan for i in info.tolist()]
        DPpos = [int(a.index('DP')) for a in callmethod['format'].str.split(':').values]
        callmethod['formatvalue'] = callmethod['formatvalue'].str.split(':')
        callmethod['totcov'] = [int(callmethod['formatvalue'].iloc[a][DPpos[a]]) if callmethod['formatvalue'].iloc[a][DPpos[a]] != '.' else 0 for a in range(callmethod.shape[0])]
        if method == 'freebayes' or method == 'mutect2' or method == 'vardict' or method == 'varscan':
            if method == 'freebayes':
                AOpos = [int(a.index('AO')) for a in callmethod['format'].str.split(':').values]
                callmethod['altcov'] = [callmethod['formatvalue'].iloc[a][AOpos[a]] if callmethod['formatvalue'].iloc[a][AOpos[a]] != '.' else 0 for a in range(callmethod.shape[0])]
            elif method == 'mutect2' or method == 'vardict':
                ADpos = [int(a.index('AD')) for a in callmethod['format'].str.split(':').values]
                callmethod['altcov'] = [','.join(callmethod['formatvalue'].iloc[a][ADpos[a]].split(',')[1:]) for a in range(callmethod.shape[0])]
            elif method == 'varscan':
                ADpos = [int(a.index('AD')) for a in callmethod['format'].str.split(':').values]
                callmethod['altcov'] = [callmethod['formatvalue'].iloc[a][ADpos[a]] for a in range(callmethod.shape[0])]
            callmethod.drop(['format', 'formatvalue'], axis=1, inplace=True)
        elif method == 'strelka2':
            callmethod['XUTIR'] = callmethod['alt'].copy()
            callmethod.loc[(callmethod.ref.str.len() - callmethod.alt.str.len() != 0), 'XUTIR'] = 'TIR'
            callmethod.loc[callmethod['XUTIR'] != 'TIR', 'XUTIR'] = callmethod.loc[callmethod['XUTIR'] != 'TIR', 'XUTIR'] + 'U'
            XUTIRpos = [int(callmethod.iloc[a]['format'].split(':').index(callmethod.iloc[a]['XUTIR'])) for a in range(callmethod.shape[0])]
            callmethod['altcov'] = [callmethod['formatvalue'].iloc[a][XUTIRpos[a]][0] for a in range(callmethod.shape[0])]
            callmethod.drop(['format', 'formatvalue', 'XUTIR'], axis=1, inplace=True)
        #print(info.shape)
        #print(callmethod.shape)
        colnames = []
        for i in info.to_list():
            # print([j.split('=')[0] for j in str(i).split(';')])
            for j in str(i).split(';'):
                colnames.append(j.split('=')[0])
        colnames = np.unique(colnames)
        #print(colnames)
        #print(len(colnames))
        for c in colnames:
            #print(c)
            #print(len([i.split(c+'=')[1].split(';')[0] if c+'=' in i else np.nan for i in info.to_list()]))
            #print(callmethod.shape)
            callmethod[c] = [i.split(c+'=')[1].split(';')[0] if c+'=' in i else np.nan for i in info.to_list()]
        if method == 'freebayes':
            # callmethod['altcov'] = callmethod['altcov'].astype(str)
            if not callmethod.empty:
                callmethod.loc[(callmethod.altcov.str.contains(',', regex=False)) & (~callmethod.alt.str.contains(',', regex=False)), 'altcov'] = callmethod.loc[(callmethod.altcov.str.contains(',', regex=False)) & (~callmethod.alt.str.contains(',', regex=False)), 'altcov'].apply(lambda x: str(max(x)))
        if not callmethod.empty:
            callmethod = callmethod.assign(alt=callmethod.alt.str.split(",")).assign(altcov=callmethod.altcov.str.split(","))
            for ci in callmethod['altcov'].isna().index:
                callmethod.at[ci, 'altcov'] = [0]
            callmethod.loc[(callmethod.alt.str.len() > callmethod.altcov.str.len()), 'altcov'] = \
                callmethod.loc[(callmethod.alt.str.len() > callmethod.altcov.str.len()), 'altcov'] * \
                callmethod.loc[(callmethod.alt.str.len() > callmethod.altcov.str.len()), 'alt'].str.len()
            checkindex = list(callmethod.loc[(callmethod.alt.str.len() < callmethod.altcov.str.len())].index)
            for ci in checkindex:
                callmethod.at[ci, 'altcov'] = list(map(int, callmethod.loc[ci, 'altcov']))[:int(len(callmethod.loc[ci, 'alt']))]
            callmethod = callmethod.set_index(callmethod.columns.difference(['alt', 'altcov']).tolist()).apply(pd.Series.explode).reset_index()
            callmethod['altcov'] = callmethod['altcov'].astype(int)
            callmethod['vaf'] = callmethod['altcov']/callmethod['totcov']
            callmethod.loc[callmethod['vaf'] > 1, 'vaf'] = 1  # bug for some complex indels
            callmethod['type'] = np.nan
            callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() == 0, 'type'] = 'SNV'
            callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() > 0, 'type'] = 'INS'
            callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() < 0, 'type'] = 'DEL'
            callmethod.loc[callmethod['ID'].str.contains('rs'), 'type'] = 'SNP'
            callmethod.drop('ID', axis=1, inplace=True)
        else:
            finalcolumns = ['chrom', method, method+'_score', 'pos', 'ref', 'totcov', 'alt', 'altcov', 'vaf', 'type']
            for c in colnames:
                finalcolumns.append(c)
            callmethod = pd.DataFrame(columns=finalcolumns)
        callmethod['chrom_pos_ref_alt'] = callmethod['chrom'].astype('str').str.cat(callmethod['pos'].astype('str'), sep="_").str.cat(callmethod['ref'].astype('str'), sep='_').str.cat(callmethod['alt'].astype('str'), sep='_')
        callmethod.set_index('chrom_pos_ref_alt', inplace=True)
        #callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score']] # need to rename
        callmethod.rename(columns={'totcov': method+'_totcov', 'altcov': method+'_altcov', 'vaf': method+'_vaf'}, inplace=True)
        callmethod_snv = callmethod[callmethod['type'] == 'SNV']
        callmethod_indel = callmethod[(callmethod['type'] == 'INS') | (callmethod['type'] == 'DEL')]
        callmethod_snp = callmethod[callmethod['type'] == 'SNP']
    elif method == 'smurf':
        calltablemethod_path = os.path.join(calldir, 'calls', method, 'snv-parse.txt')
        if not os.path.exists(calltablemethod_path):
            print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
            callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', method+'_totcov', method+'_altcov', method+'_vaf', method, method+'_score'])
        else:
            callmethod = pd.read_csv(calltablemethod_path, sep='\t')
            callmethod['T_refDepth'] = callmethod['T_refDepth'].astype(int) + callmethod['T_altDepth'].astype(int)
            callmethod.rename(columns={'Chr': 'chrom', 'START_POS_REF': 'pos', 'REF': 'ref', 'ALT': 'alt',
                                       'T_refDepth': method+'_totcov', 'T_altDepth':  method+'_altcov', 'Alt_Allele_Freq':  method+'_vaf',
                                       'predict': method, 'TRUE.': method+'_score'}, inplace=True)
            callmethod[method] = callmethod[method].astype(bool)
            callmethod['type'] = 'SNV'
            callmethod['chrom_pos_ref_alt'] = callmethod['chrom'].astype('str').str.cat(callmethod['pos'].astype('str'), sep="_").str.cat(callmethod['ref'].astype('str'), sep='_').str.cat(callmethod['alt'].astype('str'), sep='_')
            callmethod.set_index('chrom_pos_ref_alt', inplace=True)
        calltablemethod_path = os.path.join(calldir, 'calls', method, 'indel-parse.txt')
        if not os.path.exists(calltablemethod_path):
            print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
        else:
            callmethodindel = pd.read_csv(calltablemethod_path, sep='\t')
            callmethodindel['T_refDepth'] = callmethodindel['T_refDepth'].astype(int) + callmethodindel['T_altDepth'].astype(int)
            callmethodindel.rename(columns={'Chr': 'chrom', 'START_POS_REF': 'pos', 'REF': 'ref', 'ALT': 'alt',
                                       'T_refDepth': method+'_totcov', 'T_altDepth': method+'_altcov', 'Alt_Allele_Freq': method+'_vaf',
                                       'predict': method,'TRUE.': method+'_score'}, inplace=True)
            callmethodindel[method] = callmethodindel[method].astype(bool)
            callmethodindel.loc[callmethodindel['alt'].str.len() - callmethodindel['ref'].str.len() > 0, 'type'] = 'INS'
            callmethodindel.loc[callmethodindel['alt'].str.len() - callmethodindel['ref'].str.len() < 0, 'type'] = 'DEL'
            callmethodindel['chrom_pos_ref_alt'] = callmethodindel['chrom'].astype('str').str.cat(callmethodindel['pos'].astype('str'), sep="_").str.cat(callmethodindel['ref'].astype('str'), sep='_').str.cat(callmethodindel['alt'].astype('str'), sep='_')
            callmethodindel.set_index('chrom_pos_ref_alt', inplace=True)
            callmethod = pd.concat([callmethod, callmethodindel])
        callmethod_snv = callmethod[callmethod['type'] == 'SNV']
        callmethod_indel = callmethod[(callmethod['type'] == 'INS') | (callmethod['type'] == 'DEL')]
        callmethod_snp = callmethod[callmethod['type'] == 'SNP']
    # bcbio parsing INFO
    elif method == 'varnet':  # NB: no indel method
        calltablemethod_path = os.path.join(calldir, 'calls', method, sampleid+'_v1.vcf')
        if not os.path.exists(calltablemethod_path):
            print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
            callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
        else:
            callmethod = read_vcf(calltablemethod_path, varnet=True)
            if filter == 'PASS':
                callmethod = callmethod[callmethod['FILTER'] == "PASS"]
            elif filter == 'REJECT':
                callmethod = callmethod[callmethod['FILTER'] != "PASS"]
            #callmethod['#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  SAMPLE']
            #callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO', 'SAMPLE']]
            callmethod = callmethod.rename(columns={'CHROM': 'chrom', 'POS': 'pos', 'REF':'ref', 'ALT':'alt', 'SAMPLE':'sample'})
            callmethod[method] = callmethod['FILTER']
            callmethod[method+'_score'] = callmethod['INFO']
            callmethod['type'] = callmethod[method+'_score'].str.split('TYPE=').str[1].str.split(';').str[0].astype(str)
            callmethod[method+'_score'] = callmethod[method+'_score'].str.split('SCORE=').str[1].str.split(';').str[0].astype(float)
            callmethod['totcov'] = callmethod['sample'].str.split(':').str[1].astype(int)
            callmethod['altcov'] = callmethod['sample'].str.split(':').str[3].astype(int)
            callmethod['vaf'] = callmethod['sample'].str.split(':').str[4].astype(float)
            callmethod[method][callmethod[method] == 'PASS'] = True
            callmethod[method][callmethod[method] == 'REJECT'] = False
            callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() > 0, 'type'] = 'INS'
            callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() < 0, 'type'] = 'DEL'
            colnames = []
            info = callmethod['INFO']
            for i in info.to_list():
                # print([j.split('=')[0] for j in str(i).split(';')])
                for j in str(i).split(';'):
                    colnames.append(j.split('=')[0])
            colnames = np.unique(colnames)
            print(colnames)
            print(len(colnames))
            for c in colnames:
                #print(c)
                #print(len([i.split(c+'=')[1].split(';')[0] if c+'=' in i else np.nan for i in info.to_list()]))
                #print(callmethod.shape)
                callmethod[c] = [i.split(c+'=')[1].split(';')[0] if c+'=' in i else np.nan for i in info.to_list()]
            callmethod.drop('INFO', axis=1, inplace=True)
            callmethod_snv = callmethod[callmethod['type'] == 'SNV']
            callmethod_indel = callmethod[(callmethod['type'] == 'INS') | (callmethod['type'] == 'DEL')]
            callmethod_snp = callmethod[callmethod['type'] == 'SNP']
            #callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', 'varnet', 'varnet_score']]
    # cfdna callers parsing
    elif method == 'abemus':  # NB: no indel method
        calltablemethod_path = os.path.join(calldir, 'calls', method, 'pmtab_F3_optimalR_'+os.path.basename(calldir)+'.csv')
        if not os.path.exists(calltablemethod_path):
            print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
            callmethod_snv = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            callmethod_indel = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            callmethod_snp = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
        else:
            callmethod = pd.read_csv(calltablemethod_path, index_col=0)
            callmethod = callmethod.rename(columns={'chr': 'chrom', 'dbsnp': 'type', 'cov_case': 'abemus_totcov', 'cov.alt': 'abemus_altcov', 'af_case': 'abemus_vaf', 'filter.pbem_coverage': 'abemus_score'})
            callmethod['abemus_score'] = callmethod['abemus_vaf'] / (callmethod['abemus_vaf'] + callmethod['abemus_score'])  # af_case / (af_case + filter.pbem_coverage)
            callmethod.loc[~callmethod['type'].isna(), 'type'] = 'SNV'
            callmethod.loc[callmethod['type'].isna(), 'type'] = 'SNP'
            callmethod['abemus'] = True
            callmethod_snv = callmethod[callmethod['type'] == 'SNV']
            callmethod_indel = callmethod[(callmethod['type'] == 'INS') | (callmethod['type'] == 'DEL')]
            callmethod_snp = callmethod[callmethod['type'] == 'SNP']
    return callmethod_snv, callmethod_indel, callmethod_snp


if __name__ == "__main__":
    import os
    from utils.config import Config

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")

    callmethod_snv, callmethod_indel, callmethod_snp = parse_caller_feature(
        'data/mixtures/mixtures_chr22/mixtures_chr22_CRC-1014_180816-CW-T_CRC-1014_090516-CW-T/mixture_chr22_CRC-1014_180816-CW-T_50x_CRC-1014_090516-CW-T_100x',
        'freebayes', save=False)
    print(list(callmethod_snv.columns))
    print(callmethod_snv.head())