import io
import os
import gzip
import numpy as np
import pandas as pd


def read_vcf(path):
    if path.endswith('.gz'):
        path = path[:-3]
    if not os.path.exists(path) and os.path.exists(path + '.gz'):
        fp = open(path, "wb")
        with gzip.open(path + '.gz', "rb") as f:
            bindata = f.read()
        fp.write(bindata)
        fp.close()
    if os.path.exists(path):
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        res = pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})
        return res


def get_calltable(calldir, methods, save=False, filter='PASS'):
    callmethods_snv, callmethods_indel, callmethods_snp = {}, {}, {}
    sampleid = os.path.basename(calldir)
    print(sampleid)
    for method in methods:
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
                if filter == 'PASS':
                    callmethod = callmethod[callmethod['FILTER'] == "PASS"]  # TODO: get also low vaf
                elif filter == 'REJECT':
                    callmethod = callmethod[callmethod['FILTER'] != "PASS"]
                info = callmethod['INFO']
                callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'FORMAT', sampleid.replace('.', '_')+'-T', 'ID']]
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', method, 'format', 'formatvalue', 'ID']
                callmethod[method] = True
                if method == 'freebayes':  # logodds to probability score prob = exp(logODDS)/(1+exp(logODDS))
                    callmethod[method + '_score'] = [np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0]))) / (1 + np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0]))))
                                                     if 'ODDS' in i else np.nan for i in info.to_list()]
                elif method == 'mutect2':  # logodds to probability score prob = exp(TLOD)/(1+exp(TLOD))
                    callmethod[method + '_score'] = [np.exp(float(i.split('TLOD=')[1].split(';')[0])) / (1 + np.exp(float(i.split('TLOD=')[1].split(';')[0])))
                                                     if 'TLOD' in i and ',' not in i.split('TLOD=')[1].split(';')[0]
                                                     else np.exp(float(i.split('TLOD=')[1].split(';')[0].split(',')[0])) / (1 + np.exp(float(i.split('TLOD=')[1].split(';')[0].split(',')[0])))
                                                     if 'TLOD' in i and ',' in i.split('TLOD=')[1].split(';')[0]
                                                     else np.nan for i in info.to_list()]
                if method == 'strelka2':  # phred score to probability, prob = 1 - 10^(-SomaticEVS/10)
                    callmethod[method + '_score'] = [1 - (10 ** (-float(i.split('SomaticEVS=')[1].split(';')[0]) / 10))
                                                     if 'SomaticEVS' in i else np.nan for i in info.to_list()]
                if method == 'vardict':  # P-value /!\ opposite direction of score variation /!\
                    callmethod[method + '_score'] = [1-float(i.split('SSF=')[1].split(';')[0])
                                                     if 'SSF' in i else np.nan for i in info.tolist()]
                if method == 'varscan':  # P-value
                    callmethod[method + '_score'] = [float(i.split('SSC=')[1].split(';')[0])
                                                     if 'SSC' in i else np.nan for i in info.tolist()]
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
                if method == 'freebayes':
                    callmethod.loc[(callmethod.altcov.str.contains(',', regex=False)) & (~callmethod.alt.str.contains(',', regex=False)), 'altcov'] = callmethod.loc[(callmethod.altcov.str.contains(',', regex=False)) & (~callmethod.alt.str.contains(',', regex=False)), 'altcov'].apply(lambda x: str(max(x)))
                callmethod = callmethod.assign(alt=callmethod.alt.str.split(",")).assign(altcov=callmethod.altcov.str.split(","))
                checkindex = callmethod.loc[(callmethod.alt.str.len() > callmethod.altcov.str.len())].index
                callmethod.loc[(callmethod.alt.str.len() > callmethod.altcov.str.len()), 'altcov'] = callmethod.loc[(callmethod.alt.str.len() > callmethod.altcov.str.len()), 'altcov']*callmethod.loc[(callmethod.alt.str.len() > callmethod.altcov.str.len()), 'alt'].str.len()
                #print(callmethod.loc[checkindex][['pos', 'ref', 'alt', 'totcov', 'altcov']])
                #print(callmethod.loc[checkindex, 'altcov'].dtype)
                checkindex2 = list(callmethod.loc[(callmethod.alt.str.len() < callmethod.altcov.str.len())].index)
                for ci in checkindex2:
                    print(list(map(int, callmethod.loc[ci, 'altcov'])))
                    print(len(callmethod.loc[ci, 'alt']))
                    print(list(map(int, callmethod.loc[ci, 'altcov']))[:len(callmethod.loc[ci, 'alt'])])
                    callmethod.loc[ci, 'altcov'] = list(map(int, callmethod.loc[ci, 'altcov']))[:int(len(callmethod.loc[ci, 'alt']))]
                #print(callmethod.loc[checkindex2][['pos', 'ref', 'alt', 'totcov', 'altcov']])
                #print(callmethod.loc[checkindex2, 'altcov'].dtype)
                #print(callmethod.loc[(callmethod.alt.str.len() - callmethod.altcov.str.len()) != 0][['pos', 'ref', 'alt', 'totcov', 'altcov']])
                #print(callmethod[['pos', 'ref', 'alt', 'totcov', 'altcov']])
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
        elif method == 'abemus':  # NB: no indel method
            calltablemethod_path = os.path.join(calldir, 'calls', method, 'pmtab_F3_optimalR_'+os.path.basename(calldir)+'.csv')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                callmethod = pd.read_csv(calltablemethod_path)
                callmethod = callmethod[['chr', 'pos', 'ref', 'alt', 'dbsnp', 'cov_case', 'cov.alt', 'af_case', 'filter.pbem_coverage']]
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', 'abemus_score']
                callmethod['abemus_score'] = callmethod['abemus_score']
                callmethod.loc[~callmethod['type'].isna(), 'type'] = 'SNV'
                callmethod.loc[callmethod['type'].isna(), 'type'] = 'SNP'
                callmethod['abemus'] = True
        elif method == 'cfsnv':  # NB: no indel method
            calltablemethod_paths = [os.path.join(os.getcwd(), calldir, 'calls', method,  l) for l in os.listdir(os.path.join(calldir, 'calls', method)) if l.endswith('.txt')]
            if not os.path.exists(os.path.join(calldir, 'calls', method)) or calltablemethod_paths == []:
                print('calls for caller {} do not exist. paths in folder {} not found.'.format(method, os.path.join(calldir, 'calls', method)))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                li = []
                for calltablemethod_path in calltablemethod_paths:
                    li.append(pd.read_csv(calltablemethod_path, index_col=None, header=0, sep='\t'))
                callmethod = pd.concat(li, axis=0, ignore_index=True)
                callmethod.columns = ['chrom', 'pos', 'type', 'ref', 'alt', 'cfsnv_score', 'filter', 'vaf', 'tumor fraction']
                callmethod = callmethod[callmethod['filter'] == 'PASS']
                callmethod.loc[callmethod['type'] == '.', 'type'] = 'SNV'
                callmethod.loc[callmethod['type'] != 'SNV', 'type'] = 'SNP'
                callmethod['totcov'] = np.nan
                callmethod['altcov'] = np.nan
                # phred score to probability, prob = 1 - 10^(-variant.list.QUAL/10)
                callmethod['cfsnv_score'] = 1 - (10 ** (-callmethod['cfsnv_score'].astype(float) / 10))
                callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', 'cfsnv_score']]
                callmethod[method] = True
        elif method.startswith('sinvict'):
            calltablemethod_path = os.path.join(calldir, 'calls', method.split('_')[0], 'calls_'+method.split('_')[1]+'.sinvict')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                callmethod = pd.read_csv(calltablemethod_path, sep='\t', header=None)
                callmethod.columns = ['chrom', 'pos', 'samplename', 'ref', 'totcov', 'alt', 'altcov', 'vaf', '+strand', '-strand', 'avgreadpos', 'type']
                callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf']]
                callmethod['chrom'] = callmethod['chrom'].astype(str)
                callmethod['pos'] = callmethod['pos'].astype(str)
                callmethod.loc[callmethod['type'] == 'Somatic', 'type'] = 'SNV'  # NB: based on vaf not on dbsnp
                callmethod.loc[callmethod['type'] == 'Germline', 'type'] = 'SNP'  # NB: based on vaf not on dbsnp
                callmethod.loc[callmethod['alt'].str.startswith('+'), 'type'] = 'INS'  # generates warning
                callmethod.loc[callmethod['alt'].str.startswith('+'), 'alt'] = callmethod.loc[callmethod['alt'].str.startswith('+'), 'ref'] + callmethod.loc[callmethod['alt'].str.startswith('+'), 'alt'].str.replace('+', '', regex=True)
                callmethod.loc[callmethod['alt'].str.startswith('-'), 'type'] = 'DEL'  # generates warning
                callmethod.loc[callmethod['alt'].str.startswith('-'), 'alt'] = callmethod.loc[callmethod['alt'].str.startswith('-'), 'ref']
                callmethod.loc[callmethod['alt'].str.startswith('-'), 'ref'] = callmethod.loc[callmethod['alt'].str.startswith('-'), 'ref'] + callmethod.loc[callmethod['alt'].str.startswith('-'), 'alt'].str.replace('-', '', regex=True)
                callmethod['vaf'] = callmethod['vaf'] / 100  # % to fraction
                callmethod[method+'_score'] = np.nan
                callmethod[method] = True
        else:
            raise ValueError('caller {} is unknown'.format(method))
        callmethod['chrom_pos_ref_alt'] = callmethod['chrom'].astype('str').str.cat(callmethod['pos'].astype('str'), sep="_").str.cat(callmethod['ref'].astype('str'), sep='_').str.cat(callmethod['alt'].astype('str'), sep='_')
        callmethod.set_index('chrom_pos_ref_alt', inplace=True)
        callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score']]
        callmethod.rename(columns={'totcov': method+'_totcov', 'altcov': method+'_altcov', 'vaf': method+'_vaf'}, inplace=True)
        callmethod_snv = callmethod[callmethod['type'] == 'SNV']
        callmethod_indel = callmethod[(callmethod['type'] == 'INS') | (callmethod['type'] == 'DEL')]
        callmethod_snp = callmethod[callmethod['type'] == 'SNP']
        callmethods_snv[method] = callmethod_snv
        callmethods_indel[method] = callmethod_indel
        callmethods_snp[method] = callmethod_snp
    calltable_snv = pd.concat([cm.drop(['chrom', 'pos', 'ref', 'alt', 'type'], axis=1) for cm in list(callmethods_snv.values())], axis=1, sort=True)
    calltable_indel = pd.concat([cm.drop(['chrom', 'pos', 'ref', 'alt', 'type'], axis=1) for cm in list(callmethods_indel.values())], axis=1, sort=True)
    calltable_snp = pd.concat([cm.drop(['chrom', 'pos', 'ref', 'alt', 'type'], axis=1) for cm in list(callmethods_snp.values())], axis=1, sort=True)
    calltable_snv['chrom_pos_ref_alt'] = list(calltable_snv.index)
    calltable_indel['chrom_pos_ref_alt'] = list(calltable_indel.index)
    calltable_snp['chrom_pos_ref_alt'] = list(calltable_snp.index)
    calltable_snv[['chrom', 'pos', 'ref', 'alt']] = calltable_snv['chrom_pos_ref_alt'].str.split('_', expand=True)
    calltable_indel[['chrom', 'pos', 'ref', 'alt']] = calltable_indel['chrom_pos_ref_alt'].str.split('_', expand=True)
    calltable_snp[['chrom', 'pos', 'ref', 'alt']] = calltable_snp['chrom_pos_ref_alt'].str.split('_', expand=True)
    calltable_snv['type'] = 'SNV'
    calltable_indel.loc[calltable_indel['alt'].str.len() - calltable_indel['ref'].str.len() > 0, 'type'] = 'INS'
    calltable_indel.loc[calltable_indel['alt'].str.len() - calltable_indel['ref'].str.len() < 0, 'type'] = 'DEL'
    calltable_snp['type'] = 'SNP'
    for m in methods:
        calltable_snv[m] = calltable_snv[m].fillna(False)
        calltable_snv[m] = calltable_snv[m].astype(bool)
        calltable_indel[m] = calltable_indel[m].fillna(False)
        calltable_indel[m] = calltable_indel[m].astype(bool)
        calltable_snp[m] = calltable_snp[m].fillna(False)
        calltable_snp[m] = calltable_snp[m].astype(bool)
    calltable_snv = calltable_snv[['chrom', 'pos', 'ref', 'alt', 'type'] + [m+suffix for m in methods for suffix in ['', '_score']] + [m+suffix for m in methods for suffix in ['_altcov', '_totcov', '_vaf']]]
    calltable_indel = calltable_indel[['chrom', 'pos', 'ref', 'alt', 'type'] + [m+suffix for m in methods for suffix in ['', '_score']] + [m+suffix for m in methods for suffix in ['_altcov', '_totcov', '_vaf']]]
    calltable_snp = calltable_snp[['chrom', 'pos', 'ref', 'alt', 'type'] + [m+suffix for m in methods  for suffix in ['', '_score']] + [m+suffix for m in methods for suffix in ['_altcov', '_totcov', '_vaf']]]
    # correct for germline calls with GATK Haplotype
    calltablemethod_path = os.path.join(calldir, 'calls', 'bcbio', sampleid+'-N-gatk-haplotype-annotated.vcf.gz')
    if '.' in os.path.basename(calltablemethod_path[:-7]):
        calltablemethod_path = os.path.join(os.path.dirname(calltablemethod_path), os.path.basename(calltablemethod_path)[:-7].replace('.', '_') + '.vcf.gz')
    if not os.path.exists(calltablemethod_path):
        print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
        print('cannot use GATK Haplotype calls to filter germline calls')
    else:
        callmethod = read_vcf(calltablemethod_path)
        callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER']]
        callmethod = callmethod[callmethod['FILTER'] == 'PASS']
        callmethod.columns = ['chrom', 'pos', 'ref', 'alt', 'gatk']
        callmethod['chrom_pos_ref_alt'] = callmethod['chrom'].astype('str').str.cat(callmethod['pos'].astype('str'), sep="_").str.cat(callmethod['ref'].astype('str'), sep='_').str.cat(callmethod['alt'].astype('str'), sep='_')
        callmethod.set_index('chrom_pos_ref_alt', inplace=True)
        print("# calls before using germline calls from GATK Haplotype: {} SNV, {} INDEL, {} SNP".format(
            calltable_snv.shape[0], calltable_indel.shape[0], calltable_snp.shape[0]))
        idx_snv = calltable_snv.index.difference(callmethod.index, sort=False)
        idx_snp_snv = calltable_snv.index.intersection(callmethod.index, sort=False)
        calltable_snp_snv = calltable_snv.loc[idx_snp_snv]
        calltable_snv = calltable_snv.loc[idx_snv]
        idx_indel = calltable_indel.index.difference(callmethod.index, sort=False)
        idx_snp_indel = calltable_indel.index.intersection(callmethod.index, sort=False)
        calltable_snp_indel = calltable_indel.loc[idx_snp_indel]
        calltable_indel = calltable_indel.loc[idx_indel]
        calltable_snp = pd.concat([calltable_snp, calltable_snp_snv, calltable_snp_indel])
        print("# calls after using germline calls from GATK Haplotype: {} SNV, {} INDEL, {} SNP".format(
            calltable_snv.shape[0], calltable_indel.shape[0], calltable_snp.shape[0]))
        print(sum(calltable_snp.duplicated()))
    print('final shape SNV: {}'.format(calltable_snv.shape))
    print('final shape INDEL: {}'.format(calltable_indel.shape))
    print('final shape SNP: {}'.format(calltable_snp.shape))
    if save:
        if not os.path.exists(os.path.join(calldir, 'calls')):
            os.mkdir(os.path.join(calldir, 'calls'))
        #if not os.path.exists(os.path.join(calldir, 'calls', sampleid+'_snv_calls_'+filter+'.csv')):
        calltable_snv.to_csv(os.path.join(calldir, 'calls', sampleid+'_snv_calls_'+filter+'.csv'))
        #else:
        #    print('snv_calls already exists')
        #if not os.path.exists(os.path.join(calldir, 'calls', sampleid+'_indel_calls_'+filter+'.csv')):
        calltable_indel.to_csv(os.path.join(calldir, 'calls', sampleid+'_indel_calls_'+filter+'.csv'))
        #else:
        #    print('indel_calls already exists')
        #if not os.path.exists(os.path.join(calldir, 'calls', sampleid+'_snp_calls_'+filter+'.csv')):
        calltable_snp.to_csv(os.path.join(calldir, 'calls', sampleid+'_snp_calls_'+filter+'.csv'))
        #else:
        #    print('snp_calls already exists')

    return calltable_snv, calltable_indel, calltable_snp


if __name__ == "__main__":
    import os
    from utils.config import Config

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")
    #calltable_snv_all, calltable_indel_all, calltable_snp_all = get_calltable(
    #    'data/mixtures/mixtures_chr22/mixtures_chr22_CRC-1014_180816-CW-T_CRC-1014_090516-CW-T/mixture_chr22_CRC-1014_180816-CW-T_50x_CRC-1014_090516-CW-T_100x',
    #    config.methods, save=False, filter='all')
    #print(calltable_snv_all)
    #print(calltable_snv_all.columns)
    #calltable_snv_pass, calltable_indel_pass, calltable_snp_pass = get_calltable(
    #    'data/mixtures/mixtures_chr22/mixtures_chr22_CRC-1014_180816-CW-T_CRC-1014_090516-CW-T/mixture_chr22_CRC-1014_180816-CW-T_50x_CRC-1014_090516-CW-T_100x',
    #    config.methods, save=False, filter='PASS')
    #print(calltable_snv_pass)
    #print(calltable_snv_pass.columns)

    #diff_index = calltable_snv_all.index.difference(calltable_snv_pass.index)
    #calltable_snv_diff = calltable_snv_all.loc[diff_index]
    #print(calltable_snv_diff.shape)
    #res = calltable_snv_diff[[m+'_vaf' for m in config.methods]].min(axis=1)
    #print(res)
    #print(res.describe())

    calltable_snv_reject, calltable_indel_reject, calltable_snp_reject = get_calltable(
        'data/mixtures/mixtures_chr22/mixtures_chr22_CRC-1014_180816-CW-T_CRC-1014_090516-CW-T/mixture_chr22_CRC-1014_180816-CW-T_50x_CRC-1014_090516-CW-T_100x',
        config.methods, save=False, filter='PASS')
    print(calltable_snv_reject)
    res = calltable_snv_reject[[m+'_vaf' for m in config.methods]].min(axis=1)
    print(res[res > 1])
    print(res)
    print(res.describe())
    """
    for irow, row in calltable_snv.iterrows():  # sanity check on vaf
        if (',' in list(row.values)[-1]) and int(sum(row[config.methods].values)) != (int(list(row.values)[-1].count(',')) + 1):
            print(calltable_snv.columns.tolist())
            print(list(row.values))
    print(calltable_indel)
    print(calltable_indel.columns)
    for irow, row in calltable_indel.iterrows():  # sanity check on vaf
        if (',' in list(row.values)[-1]) and int(sum(row[config.methods].values)) != (int(list(row.values)[-1].count(',')) + 1):
            print(calltable_indel.columns.tolist())
            print(list(row.values))
    print(calltable_snp)
    print(calltable_snp.columns)
    for irow, row in calltable_snp.iterrows():  # sanity check on vaf
        if (',' in list(row.values)[-1]) and int(sum(row[config.methods].values)) != (int(list(row.values)[-1].count(',')) + 1):
            print(calltable_snp.columns.tolist())
            print(list(row.values))
    """