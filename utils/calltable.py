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


def get_calltable(calldir, methods):
    callmethods_snv, callmethods_indel, callmethods_snp = {}, {}, {}
    sampleid = os.path.basename(calldir)[:-7]
    print(sampleid)
    for method in methods:
        print(method)
        if method in ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']:
            calltablemethod_path = os.path.join(calldir, 'bcbio', sampleid+'-'+method+'-annotated.vcf.gz')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
            else:
                callmethod = read_vcf(calltablemethod_path)
                callmethod = callmethod[callmethod['FILTER'] == "PASS"]
                callmethod = callmethod[callmethod['CHROM'] == '22']  # TODO: add bedfile to bcbio
                info = callmethod['INFO']
                sampleidbis = sampleid.replace('mixture', 'dilution')  # just for test sample with mismatch in IDs
                callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'FORMAT', sampleidbis+'-T', 'ID']]
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', method, 'format', 'formatvalue', 'ID']
                callmethod[method] = True
                if method == 'freebayes':  # logodds to probability score prob = exp(logODDS)/(1+exp(logODDS))
                    callmethod[method + '_score'] = [np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0]))) / (1 + np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0]))))
                                                     if 'ODDS' in i else np.nan for i in info.to_list()]
                elif method == 'mutect2':  # logodds to probability score prob = exp(logTLOD)/(1+exp(logTLOD))
                    callmethod[method + '_score'] = [np.exp(np.log(float(i.split('TLOD=')[1].split(';')[0]))) / (1 + np.exp(np.log(float(i.split('TLOD=')[1].split(';')[0]))))
                                                     if 'TLOD' in i else np.nan for i in info.to_list()]
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
                callmethod['totcov'] = [int(callmethod['formatvalue'].iloc[a][DPpos[a]]) for a in range(callmethod.shape[0])]
                if method == 'freebayes' or method == 'mutect2' or method == 'vardict' or method == 'varscan':
                    if method == 'freebayes':
                        AOpos = [int(a.index('AO')) for a in callmethod['format'].str.split(':').values]
                        callmethod['altcov'] = [callmethod['formatvalue'].iloc[a][AOpos[a]] for a in range(callmethod.shape[0])]
                    elif method == 'mutect2' or method == 'vardict':
                        ADpos = [int(a.index('AD')) for a in callmethod['format'].str.split(':').values]
                        callmethod['altcov'] = [callmethod['formatvalue'].iloc[a][ADpos[a]][1:] for a in range(callmethod.shape[0])]
                    elif method == 'varscan':
                        ADpos = [int(a.index('AD')) for a in callmethod['format'].str.split(':').values]
                        callmethod['altcov'] = [callmethod['formatvalue'].iloc[a][ADpos[a]] for a in range(callmethod.shape[0])]
                    callmethod.drop(['format', 'formatvalue'], axis=1, inplace=True)
                elif method == 'strelka2':
                    callmethod['XUTIR'] = callmethod['alt'].copy()
                    callmethod.loc[(callmethod.ref.str.len() != 1) | (callmethod.alt.str.len() != 1), 'XUTIR'] = 'TIR'
                    callmethod.loc[callmethod['XUTIR'] != 'TIR', 'XUTIR'] = callmethod.loc[callmethod['XUTIR'] != 'TIR', 'XUTIR'] + 'U'
                    XUTIRpos = [int(callmethod.iloc[a]['format'].split(':').index(callmethod.iloc[a]['XUTIR'])) for a in range(callmethod.shape[0])]
                    callmethod['altcov'] = [callmethod['formatvalue'].iloc[a][XUTIRpos[a]] for a in range(callmethod.shape[0])]
                    callmethod.drop(['format', 'formatvalue', 'XUTIR'], axis=1, inplace=True)
                callmethod = callmethod.assign(alt=callmethod.alt.str.split(",")).assign(altcov=callmethod.altcov.str.split(","))
                callmethod.loc[callmethod.alt.str.len() != callmethod.altcov.str.len(), 'altcov'] = callmethod.loc[callmethod.alt.str.len() != callmethod.altcov.str.len(), 'altcov'].apply(lambda x: max(x))
                callmethod = callmethod.set_index(callmethod.columns.difference(['alt', 'altcov']).tolist()).apply(pd.Series.explode).reset_index()
                callmethod['altcov'] = callmethod['altcov'].astype(int)
                callmethod['vaf'] = callmethod['altcov']/callmethod['totcov']
                callmethod['type'] = np.nan
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() == 0, 'type'] = 'SNV'
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() > 0, 'type'] = 'INS'
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() < 0, 'type'] = 'DEL'
                callmethod.loc[callmethod['ID'].str.contains('rs'), 'type'] = 'SNP'
                callmethod.drop('ID', axis=1, inplace=True)
        elif method == 'abemus':  # NB: no indel method
            calltablemethod_path = os.path.join(calldir, method, 'pmtab_F3_optimalR_'+os.path.basename(calldir)+'.tsv')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
            else:
                callmethod = pd.read_csv(calltablemethod_path, sep='\t')
                # 'pbem.allele', 'filter.pbem_coverage' ??? scores
                callmethod = callmethod[['chr', 'pos', 'ref', 'alt', 'dbsnp', 'cov_case', 'cov.alt', 'af_case', 'filter.pbem_coverage']]
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', 'abemus_score']
                callmethod.loc[~callmethod['type'].isna(), 'type'] = 'SNV'
                callmethod.loc[callmethod['type'].isna(), 'type'] = 'SNP'
                callmethod['abemus'] = True
                callmethod['chrom_pos_ref_alt'] = callmethod['chrom'].astype('str').str.cat(callmethod['pos'].astype('str'), sep="_").str.cat(callmethod['ref'].astype('str'), sep='_').str.cat(callmethod['alt'].astype('str'), sep='_')
                callmethod.set_index('chrom_pos_ref_alt', inplace=True)
        elif method == 'sinvict':
            calltablemethod_path = os.path.join(calldir, method, 'calls_level4.sinvict')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
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
                callmethod['sinvict_score'] = np.nan
                callmethod['sinvict'] = True
        else:
            raise ValueError('caller {} is unknown'.format(method))
        callmethod['chrom_pos_ref_alt'] = callmethod['chrom'].astype('str').str.cat(callmethod['pos'].astype('str'), sep="_").str.cat(callmethod['ref'].astype('str'), sep='_').str.cat(callmethod['alt'].astype('str'), sep='_')
        callmethod.set_index('chrom_pos_ref_alt', inplace=True)
        if method == 'abemus':
            callmethod.drop('chr_pos_ref_alt', inplace=True)  # TODO: fix abemus concatenation
        callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score']]
        # callmethod.columns = ['chrom', 'pos', 'ref', 'alt', 'type', 'totcov_'+method, 'altcov_'+method, 'vaf_'+method, method, method+'_score']
        callmethod_snv = callmethod[callmethod['type'] == 'SNV']
        callmethod_indel = callmethod[(callmethod['type'] == 'INS') | (callmethod['type'] == 'DEL')]
        callmethod_snp = callmethod[callmethod['type'] == 'SNP']
        callmethods_snv[method] = callmethod_snv
        callmethods_indel[method] = callmethod_indel
        callmethods_snp[method] = callmethod_snp
    calltable_snv = pd.concat([cm.set_index(['chrom', 'pos', 'ref', 'alt', 'type']) for cm in list(callmethods_snv.values())], axis=1)
    calltable_indel = pd.concat([cm.set_index(['chrom', 'pos', 'ref', 'alt', 'type']) for cm in list(callmethods_indel.values())], axis=1)
    calltable_snp = pd.concat([cm.set_index(['chrom', 'pos', 'ref', 'alt', 'type']) for cm in list(callmethods_snp.values())], axis=1)
    calltable_snv = calltable_snv.groupby(calltable_snv.columns, axis=1).agg(lambda x: x.apply(lambda y: ','.join([str(l) for l in y if str(l) != "nan"]), axis=1))
    calltable_indel = calltable_indel.groupby(calltable_indel.columns, axis=1).agg(lambda x: x.apply(lambda y: ','.join([str(l) for l in y if str(l) != "nan"]), axis=1))
    calltable_snp = calltable_snp.groupby(calltable_snp.columns, axis=1).agg(lambda x: x.apply(lambda y: ','.join([str(l) for l in y if str(l) != "nan"]), axis=1))
    for m in methods:
        calltable_snv[m] = calltable_snv[m].fillna(False)
        calltable_snv[m] = calltable_snv[m].astype(bool)
        calltable_indel[m] = calltable_indel[m].fillna(False)
        calltable_indel[m] = calltable_indel[m].astype(bool)
        calltable_snp[m] = calltable_snp[m].fillna(False)
        calltable_snp[m] = calltable_snp[m].astype(bool)
    print('final shape SNV: {}'.format(calltable_snv.shape))
    print('final shape INDEL: {}'.format(calltable_indel.shape))
    print('final shape SNP: {}'.format(calltable_snp.shape))
    return calltable_snv, calltable_indel, calltable_snp


if __name__ == "__main__":
    import os
    from utils.config import Config

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")
    calltable_snv, calltable_indel, calltable_snp = get_calltable(
        'data/callers_output/mixtures/mixtures_chr22/mixture_chr22_CRC-1014_180816-CW-T_10x_CRC-1014_090516-CW-T_140x.sorted',
        config.methods)
    print(calltable_snv)
    print(calltable_snv.columns)
    for irow, row in calltable_snv.iterrows():  # sanity check
        if (',' in list(row.values)[2]) and int(sum(row[config.methods].values)) != (int(list(row.values)[2].count(',')) + 1):
            print(calltable_snv.columns.tolist())
            print(list(row.values))
    print(calltable_indel)
    print(calltable_indel.columns)
    for irow, row in calltable_indel.iterrows():  # sanity check
        if (',' in list(row.values)[2]) and int(sum(row[config.methods].values)) != (int(list(row.values)[2].count(',')) + 1):
            print(calltable_indel.columns.tolist())
            print(list(row.values))
    print(calltable_snp)
    print(calltable_snp.columns)
    for irow, row in calltable_snp.iterrows():  # sanity check
        if (',' in list(row.values)[2]) and int(sum(row[config.methods].values)) != (int(list(row.values)[2].count(',')) + 1):
            print(calltable_snp.columns.tolist())
            print(list(row.values))