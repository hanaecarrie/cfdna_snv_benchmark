import io
import os
import gzip
import numpy as np
import pandas as pd


def read_vcf(path, varnet=False, memory_map=True):
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
        if varnet:
            l0 = [l for l in lines[0].split(' ') if l != '']
            ll = '\t'.join(l0) + ''.join(lines[1:])
            index_col=False
        else:
            ll = ''.join(lines)
            index_col=None
        res = pd.read_csv(
            io.StringIO(ll),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t',
            index_col=index_col,
            memory_map=memory_map,
        ).rename(columns={'#CHROM': 'CHROM'})
        return res


def get_calltable(calldir, methods, save=False, filter='PASS', bcbiovaf=1, gatkcorr=True, verbose=0):
    callmethods_snv, callmethods_indel, callmethods_snp = {}, {}, {}
    sampleid = os.path.basename(calldir)
    if verbose >= 1:
        print(sampleid)
    for method in methods:
        if verbose >= 1:
            print(method)
        # the 5 bcbio callers
        if method in ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']:
            if bcbiovaf == 1:
                calltablemethod_path = os.path.join(calldir, 'calls', 'bcbio', sampleid+'-'+method+'-annotated.vcf.gz')
            else:
                if "BRP" in sampleid or "IDT" in sampleid or "ILM" in sampleid or "TFS" in sampleid:
                    calltablemethod_path = os.path.join(calldir, 'calls', 'bcbio'+str(bcbiovaf)+'vaf', '-'.join(sampleid.split('-')[:-1]) + '_vaf-' + sampleid.split('-')[-1] +'-'+method+'-annotated.vcf.gz')
                else:
                    calltablemethod_path = os.path.join(calldir, 'calls', 'bcbio'+str(bcbiovaf)+'vaf', sampleid +'-'+method+'-annotated.vcf.gz')
            if verbose >= 1:
                print(calltablemethod_path)
            if '.' in os.path.basename(calltablemethod_path[:-7]):
                calltablemethod_path = os.path.join(os.path.dirname(calltablemethod_path), os.path.basename(calltablemethod_path)[:-7].replace('.', '_') + '.vcf.gz')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found. Returning an empty DataFrame.'.format(method, calltablemethod_path))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                callmethod = read_vcf(calltablemethod_path)
                if callmethod.empty:
                    callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
                else:
                    if filter == 'PASS':
                        callmethod = callmethod[callmethod['FILTER'] == "PASS"]
                    elif filter == 'REJECT':
                        callmethod = callmethod[callmethod['FILTER'] != "PASS"]
                    if (bcbiovaf != 1) and ("BRP" in sampleid or "IDT" in sampleid or "ILM" in sampleid or "TFS" in sampleid):
                        sid = '-'.join(sampleid.split('-')[:-1]) + '_vaf-' + sampleid.split('-')[-1]
                        sid = sid.replace('.', '_')+'-T'
                    else:
                        sid = sampleid.replace('.', '_')+'-T'
                    if verbose >= 1:
                        print(sid)
                    callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'FORMAT', sid, 'ID', 'INFO']]
                    callmethod.columns = ['chrom', 'pos', 'ref', 'alt', method, 'format', 'formatvalue', 'ID', 'INFO']
                    if filter == 'all':
                        callmethod[method] = callmethod[method].astype(str)
                        if method == 'mutect2' or method == 'strelka2':
                            # print(callmethod.loc[callmethod[method] == 'MinAF', method].value_counts())
                            print('retrieving {} {} calls with MinAF tags on top of {} PASS calls'.format(
                                callmethod[callmethod[method] == 'MinAF'].shape[0], method, callmethod[callmethod[method] == "PASS"].shape[0]))
                            callmethod.loc[callmethod[method] == 'MinAF', method] = True  # considering variant exluded just for low VAF true calls
                        elif method == 'vardict':
                            if bcbiovaf == 1:
                                print('retrieving {} {} calls with f0.01;REJECT;REJECT tags out of {}'.format(
                                    callmethod[(callmethod[method] == 'f0.01;REJECT;REJECT')].shape[0], method, callmethod[callmethod[method] == "PASS"].shape[0]))
                                callmethod.loc[(callmethod[method] == 'f0.01;REJECT;REJECT'), method] = True  # considering variant exluded just for low VAF true calls
                            elif bcbiovaf == 0.01:
                                print('retrieving {} {} calls with v3;f0.0001000000000000000020816681712;REJECT;REJECT tags out of {}'.format(method,
                                    callmethod[(callmethod[method] == 'v3;f0.0001000000000000000020816681712;REJECT;REJECT')].shape[0], callmethod[callmethod[method] == "PASS"].shape[0]))
                                callmethod.loc[(callmethod[method] == 'v3;f0.0001000000000000000020816681712;REJECT;REJECT'), method] = True  # considering variant exluded just for low VAF true calls
                        callmethod.loc[callmethod[method] == 'PASS', method] = True
                    # Extracting caller score
                    info = callmethod['INFO']
                    callmethod.drop('INFO', axis=1, inplace=True)
                    if method == 'freebayes':  # logodds to probability score prob = exp(logODDS)/(1+exp(logODDS))
                        callmethod[method + '_score'] = [np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0]))) / (1 + np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0]))))
                                                         if 'ODDS' in i else np.nan for i in info.to_list()]
                    elif method == 'mutect2':  # logodds to probability score prob = exp(TLOD)/(1+exp(TLOD))
                        callmethod[method + '_score'] = [min(1.0, np.exp(float(i.split('TLOD=')[1].split(';')[0])) / (1 + np.exp(float(i.split('TLOD=')[1].split(';')[0]))))
                                                         if 'TLOD' in i and ',' not in i.split('TLOD=')[1].split(';')[0]
                                                         else min(1.0, np.exp(float(i.split('TLOD=')[1].split(';')[0].split(',')[0])) / (1 + np.exp(float(i.split('TLOD=')[1].split(';')[0].split(',')[0]))))
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
                    # Getting totcov and altcov
                    callmethod.dropna(subset=['format', 'formatvalue'], inplace=True)
                    DPpos = [int(a.index('DP')) for a in callmethod['format'].str.split(':').values]  # DPpos=2 for Freebayes, DPpos=3 for Mutect2
                    callmethod['formatvalue'] = callmethod['formatvalue'].str.split(':')
                    # cases where no tot cov indicated
                    callmethod['totcov'] = [int(callmethod['formatvalue'].iloc[a][DPpos[a]]) if callmethod['formatvalue'].iloc[a][DPpos[a]] != '.' else 0 for a in range(callmethod.shape[0])]
                    if method == 'freebayes':
                        callmethod['altcov'] = callmethod['formatvalue'].str[-3]  # AO
                        if not callmethod[(callmethod['altcov'].str.count(',').subtract(callmethod['alt'].str.count(',')) != 0)].empty:
                            # dropping cases with issues in Freebayes output vcf where number of alt base is not consistent with number of alt coverage values indicated
                            callmethod.drop(callmethod[(callmethod['altcov'].str.count(',').subtract(callmethod['alt'].str.count(',')) != 0)].index, inplace=True)
                        # cases where no alt cov indicated
                        callmethod.loc[((callmethod['altcov'] == '.') | (callmethod['altcov'].isna())), 'altcov'] = '0'
                        # repeat alt alleles in different rows
                        callmethod['alt'] = callmethod['alt'].str.split(',')
                        callmethod['altcov'] = callmethod['altcov'].str.split(',')
                        callmethod = callmethod.explode(['alt', 'altcov'])  # requires pandas ≥ 1.3.0
                        callmethod['altcov'] = callmethod['altcov'].astype(int)
                        callmethod.drop(['format', 'formatvalue'], axis=1, inplace=True)
                    elif method == 'mutect2' or method == 'vardict':
                        ADpos = [int(a.index('AD')) for a in callmethod['format'].str.split(':').values]  # ADpos=1
                        callmethod['altcov'] = [','.join(callmethod['formatvalue'].iloc[a][ADpos[a]].split(',')[1:]) for a in range(callmethod.shape[0])]
                        callmethod['alt'] = callmethod['alt'].str.split(',')
                        callmethod['altcov'] = callmethod['altcov'].str.split(',')
                        callmethod = callmethod.explode(['alt', 'altcov'])  # requires pandas ≥ 1.3.0
                        callmethod['altcov'] = callmethod['altcov'].astype(int)
                        callmethod.drop(['format', 'formatvalue'], axis=1, inplace=True)
                    elif method == 'varscan':
                        ADpos = [int(a.index('AD')) for a in callmethod['format'].str.split(':').values]  # ADpos=4
                        callmethod['altcov'] = [int(callmethod['formatvalue'].iloc[a][ADpos[a]]) for a in range(callmethod.shape[0])]
                        # callmethod['vaf'] = callmethod['formatvalue'].str[-2].astype(float)  # FREQ
                        callmethod.drop(['format', 'formatvalue'], axis=1, inplace=True)
                    elif method == 'strelka2':
                        callmethod['XUTIR'] = callmethod['alt'].copy()
                        callmethod.loc[(callmethod.ref.str.len() - callmethod.alt.str.len() != 0), 'XUTIR'] = 'TIR' # indels TIR
                        callmethod.loc[callmethod['XUTIR'] != 'TIR', 'XUTIR'] = callmethod.loc[callmethod['XUTIR'] != 'TIR', 'XUTIR'] + 'U' # snvs <ALT>U
                        XUTIRpos = [int(callmethod.iloc[a]['format'].split(':').index(callmethod.iloc[a]['XUTIR'])) for a in range(callmethod.shape[0])]
                        callmethod['altcov'] = [callmethod['formatvalue'].iloc[a][XUTIRpos[a]][0] for a in range(callmethod.shape[0])]
                        callmethod['totcov'] = [callmethod['formatvalue'].iloc[a][XUTIRpos[a]][0] for a in range(callmethod.shape[0])]
                        # callmethod['vaf'] = callmethod['formatvalue'].str[-1].astype(float)
                        callmethod.drop(['format', 'formatvalue', 'XUTIR'], axis=1, inplace=True)
                    callmethod['altcov'] = callmethod['altcov'].astype(int)
                    callmethod['totcov'] = callmethod['totcov'].astype(int)
                    callmethod['vaf'] = callmethod['altcov']/callmethod['totcov'] # VAF = altcov/totcov
                    if not callmethod.loc[callmethod['vaf'] > 1].empty:
                        print("{} call(s) with VAF > 1. Usually some bugs in VCF for indels. Setting VAF to 1 instead.".format(callmethod.loc[callmethod['vaf'] > 1].shape[0]))
                        print(callmethod.loc[callmethod['vaf'] > 1])
                    callmethod.loc[callmethod['vaf'] > 1, 'vaf'] = 1
                    # Determine mutation type: SNV, INSertion or DELetion, or SNP
                    callmethod['type'] = np.nan
                    callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() == 0, 'type'] = 'SNV'
                    callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() > 0, 'type'] = 'INS'
                    callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() < 0, 'type'] = 'DEL'
                    callmethod.loc[callmethod['ID'].str.contains('rs'), 'type'] = 'SNP'
                    callmethod.drop('ID', axis=1, inplace=True)
        # SMuRF
        elif method == 'smurf':
            calltablemethod_path = os.path.join(calldir, 'calls', method, 'snv-parse.txt')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                callmethod = pd.read_csv(calltablemethod_path, sep='\t')
                callmethod = callmethod[['Chr', 'START_POS_REF', 'REF', 'ALT', 'T_refDepth', 'T_altDepth', 'Alt_Allele_Freq', 'predict','TRUE.']]
                callmethod['T_refDepth'] = callmethod['T_refDepth'].astype(int) + callmethod['T_altDepth'].astype(int)
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', 'totcov', 'altcov', 'vaf', method, method+'_score']
                callmethod[method] = callmethod[method].astype(bool)
                callmethod['type'] = 'SNV'
            calltablemethod_path = os.path.join(calldir, 'calls', method, 'indel-parse.txt')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
            else:
                callmethodindel = pd.read_csv(calltablemethod_path, sep='\t')
                callmethodindel = callmethodindel[['Chr', 'START_POS_REF', 'REF', 'ALT', 'T_refDepth', 'T_altDepth', 'Alt_Allele_Freq', 'predict','TRUE.']]
                callmethodindel['T_refDepth'] = callmethodindel['T_refDepth'].astype(int) + callmethodindel['T_altDepth'].astype(int)
                callmethodindel.columns = ['chrom', 'pos', 'ref', 'alt', 'totcov', 'altcov', 'vaf', method, method+'_score']
                callmethodindel[method] = callmethodindel[method].astype(bool)
                callmethodindel.loc[callmethodindel['alt'].str.len() - callmethodindel['ref'].str.len() > 0, 'type'] = 'INS'
                callmethodindel.loc[callmethodindel['alt'].str.len() - callmethodindel['ref'].str.len() < 0, 'type'] = 'DEL'
                callmethod = pd.concat([callmethod, callmethodindel])
        # VarNet
        elif method == 'varnet':
            calltablemethod_path = os.path.join(calldir, 'calls', method, sampleid+'.vcf')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                callmethod = read_vcf(calltablemethod_path, varnet=True)
                if filter == 'PASS':
                    callmethod = callmethod[callmethod['FILTER'] == "PASS"]
                elif filter == 'REJECT':
                    callmethod = callmethod[callmethod['FILTER'] != "PASS"]
                callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO', 'SAMPLE']]
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', method, method+'_score', 'sample']
                callmethod['type'] = callmethod[method+'_score'].str.split('TYPE=').str[1].str.split(';').str[0].astype(str)
                callmethod[method+'_score'] = callmethod[method+'_score'].str.split('SCORE=').str[1].str.split(';').str[0].astype(float)
                callmethod['totcov'] = callmethod['sample'].str.split(':').str[1].astype(int)
                callmethod['altcov'] = callmethod['sample'].str.split(':').str[3].astype(int)
                callmethod['vaf'] = callmethod['sample'].str.split(':').str[4].astype(float)
                callmethod.loc[callmethod[method] == 'PASS', method] = True
                callmethod.loc[callmethod[method] == 'REJECT', method] = False
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() > 0, 'type'] = 'INS'
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() < 0, 'type'] = 'DEL'
                callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', 'varnet', 'varnet_score']]
        # ABEMUS
        elif method == 'abemus':  # NB: no indel method
            calltablemethod_path = os.path.join(calldir, 'calls', method, 'pmtab_F3_optimalR_'+os.path.basename(calldir)+'.csv')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                callmethod = pd.read_csv(calltablemethod_path)
                callmethod = callmethod[['chr', 'pos', 'ref', 'alt', 'dbsnp', 'cov_case', 'cov.alt', 'af_case', 'filter.pbem_coverage']]
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', 'abemus_score']
                callmethod['abemus_score'] = callmethod['vaf'] / (callmethod['vaf'] + callmethod['abemus_score'])  # af_case / (af_case + filter.pbem_coverage)
                callmethod.loc[~callmethod['type'].isna(), 'type'] = 'SNV'
                callmethod.loc[callmethod['type'].isna(), 'type'] = 'SNP'
                callmethod['abemus'] = True
        # cfSNV
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
                # phred score to probability, prob = 10^(-SomaticEVS/10)
                callmethod['cfsnv_score'] = [1-10**(-float(cm)/10) for cm in callmethod['cfsnv_score']]
                callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', 'cfsnv_score']]
                callmethod[method] = True
        # SiNVICT
        elif method == 'sinvict':
            callmethod_dict = {}
            for lev in range(1, 7):
                calltablemethod_path = os.path.join(calldir, 'calls', method.split('_')[0], 'calls_level'+str(lev)+'.sinvict')
                if not os.path.exists(calltablemethod_path):
                    print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
                    callmethod_sinvict = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method+'_level'+str(lev)])
                elif os.path.getsize(calltablemethod_path) == 0:  # empty file like level 5 and 6
                    callmethod_sinvict = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method+'_level'+str(lev)])
                else:
                    callmethod_sinvict = pd.read_csv(calltablemethod_path, sep='\t', header=None)
                    callmethod_sinvict.columns = ['chrom', 'pos', 'samplename', 'ref', 'totcov', 'alt', 'altcov', 'vaf', '+strand', '-strand', 'avgreadpos', 'type']
                    callmethod_sinvict = callmethod_sinvict[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf']]
                    callmethod_sinvict['chrom'] = callmethod_sinvict['chrom'].astype(str)
                    callmethod_sinvict['pos'] = callmethod_sinvict['pos'].astype(str)
                    callmethod_sinvict.loc[callmethod_sinvict['type'] == 'Somatic', 'type'] = 'SNV'  # NB: based on vaf not on dbsnp
                    callmethod_sinvict.loc[callmethod_sinvict['type'] == 'Germline', 'type'] = 'SNP'  # NB: based on vaf not on dbsnp
                    callmethod_sinvict.loc[callmethod_sinvict['alt'].str.startswith('+'), 'type'] = 'INS'  # generates warning
                    callmethod_sinvict.loc[callmethod_sinvict['alt'].str.startswith('+'), 'alt'] = callmethod_sinvict.loc[callmethod_sinvict['alt'].str.startswith('+'), 'ref'] + callmethod_sinvict.loc[callmethod_sinvict['alt'].str.startswith('+'), 'alt'].str.replace('+', '', regex=True)
                    callmethod_sinvict.loc[callmethod_sinvict['alt'].str.startswith('-'), 'type'] = 'DEL'  # generates warning
                    callmethod_sinvict.loc[callmethod_sinvict['alt'].str.startswith('-'), 'alt'] = callmethod_sinvict.loc[callmethod_sinvict['alt'].str.startswith('-'), 'ref']
                    callmethod_sinvict.loc[callmethod_sinvict['alt'].str.startswith('-'), 'ref'] = callmethod_sinvict.loc[callmethod_sinvict['alt'].str.startswith('-'), 'ref'] + callmethod_sinvict.loc[callmethod_sinvict['alt'].str.startswith('-'), 'alt'].str.replace('-', '', regex=True)
                    callmethod_sinvict['vaf'] = callmethod_sinvict['vaf'] / 100  # % to fraction
                    callmethod_sinvict[method+'_level'+str(lev)] = True
                callmethod_dict[method+'_level'+str(lev)] = callmethod_sinvict
            callmethod = pd.concat(list(callmethod_dict.values()), axis=1)
            callmethod = callmethod.loc[:, ~callmethod.columns.duplicated()]
            callmethod[method+'_score'] = callmethod[[method+'_level'+str(le) for le in range(1, 7)]].sum(axis=1) / 6
            callmethod[method] = False
            callmethod.loc[callmethod[method+'_score'] > 0, method] = True
            callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', 'sinvict', 'sinvict_score']]
        elif method == 'BRP':
            calltablemethod_path = os.path.join(calldir, 'calls', method, sampleid[:-2]+'.vcf.gz')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                callmethod = read_vcf(calltablemethod_path)
                if filter == 'PASS':
                    callmethod = callmethod[callmethod['FILTER'] == "PASS"]
                elif filter == 'REJECT':
                    callmethod = callmethod[callmethod['FILTER'] != "PASS"]
                callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'Sample1']]
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', method, 'sample']
                callmethod['type'] = 'SNV'
                callmethod[method+'_score'] = callmethod['sample'].str.split(':').str[7].astype(float)
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() > 0, 'type'] = 'INS'
                callmethod['totcov'] = callmethod['sample'].str.split(':').str[3].astype(int)
                callmethod['altcov'] = callmethod['sample'].str.split(':').str[5].astype(int)
                callmethod['vaf'] = callmethod['sample'].str.split(':').str[6].str.split('%').str[0].astype(float) / 100 # % percentage
                callmethod[method][callmethod[method] == 'PASS'] = True
                callmethod[method][callmethod[method] == 'REJECT'] = False
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() > 0, 'type'] = 'INS'
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() < 0, 'type'] = 'DEL'
                callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score']]
        elif method == 'IDT':
            calltablemethod_path = os.path.join(calldir, 'calls', method, sampleid[:-2]+'.vcf.gz')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                callmethod = read_vcf(calltablemethod_path)
                if filter == 'PASS':
                    callmethod = callmethod[callmethod['FILTER'] == "PASS"]
                elif filter == 'REJECT':
                    callmethod = callmethod[callmethod['FILTER'] != "PASS"]
                callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'tumor']]
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', method, 'tumor']
                callmethod['type'] = 'SNV'
                callmethod['totcov'] = callmethod['tumor'].str.split(':').str[1].astype(int)
                callmethod['altcov'] = callmethod['tumor'].str.split(':').str[2].astype(int)
                callmethod['vaf'] = callmethod['altcov']/callmethod['totcov']
                callmethod[method+'_score'] = callmethod['vaf'] # vardict not present
                callmethod[method][callmethod[method] == 'PASS'] = True
                callmethod[method][callmethod[method] == 'REJECT'] = False
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() > 0, 'type'] = 'INS'
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() < 0, 'type'] = 'DEL'
                callmethod['chrom'] = callmethod['chrom'].str[3:]
                callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score']]
        elif method == 'ILM':
            calltablemethod_path = os.path.join(calldir, 'calls', method, sampleid[:-2]+'.vcf.gz')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                callmethod = read_vcf(calltablemethod_path)
                if filter == 'PASS':
                    callmethod = callmethod[callmethod['FILTER'] == "LowSupport"]
                elif filter == 'REJECT':
                    callmethod = callmethod[callmethod['FILTER'] != "LowSupport"]
                callmethod = callmethod[callmethod['FILTER'] == "LowSupport"]
                samplename, met, st, ng, libp = sampleid.split('_')
                print('WG2_'+met+'_'+ st + '_' + samplename[-2:]+'_' + ng[:-2]+'_' + libp[3]+'.sorted.stitched.bam')
                samplen = 'WG2_'+met+'_'+ st + '_' + samplename[-2:]+'_' + ng[:-2]+'_' + libp[3]+'.sorted.stitched.bam'
                callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', samplen]]
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', method, 'sample']
                callmethod['type'] = 'SNV'
                callmethod['totcov'] = callmethod['sample'].str.split(':').str[3].astype(int)
                callmethod['altcov'] = callmethod['sample'].str.split(':').str[2].str.split(',').str[1].astype(int)
                callmethod['vaf'] = callmethod['sample'].str.split(':').str[4].astype(float)
                callmethod[method+'_score'] = callmethod['sample'].str.split(':').str[-1].astype(float)
                callmethod = callmethod[callmethod[method+'_score'] > 0]
                callmethod[method][callmethod[method] == 'LowSupport'] = True
                callmethod[method][callmethod[method] != True] = False
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() > 0, 'type'] = 'INS'
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() < 0, 'type'] = 'DEL'
                callmethod['chrom'] = callmethod['chrom'].str[3:]
                callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score']]
        elif method == 'ROC':
            calltablemethod_path = os.path.join(calldir, 'calls', method, sampleid[:-3]+'.vcf.gz')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                callmethod = read_vcf(calltablemethod_path)
                if filter == 'PASS':
                    callmethod = callmethod[callmethod['FILTER'] == "PASS"]
                elif filter == 'REJECT':
                    callmethod = callmethod[callmethod['FILTER'] != "PASS"]
                samplename, met, st, ng, libp = sampleid.split('_')
                samplename = samplename[-2:] + '_' + ng[:-2] + '_Rep' + libp[3]
                callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO']]
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', method, 'INFO']
                callmethod['type'] = 'SNV'
                callmethod['totcov'] = callmethod['INFO'].str.split(';').str[0].str.split('=').str[1].astype(int)
                callmethod['altcov'] = callmethod['INFO'].str.split(';').str[1].str.split('=').str[1].astype(int)
                callmethod['vaf'] = callmethod['INFO'].str.split(';').str[5].str.split('=').str[1].astype(float)
                callmethod[method+'_score'] = callmethod['INFO'].str.split('SCORE').str[1].str.split(';').str[0].str.split('=').str[1].astype(float)
                callmethod[method][callmethod[method] == 'PASS'] = True
                callmethod[method][callmethod[method] == 'REJECT'] = False
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() > 0, 'type'] = 'INS'
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() < 0, 'type'] = 'DEL'
                callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score']]
        elif method == 'TFS':
            calltablemethod_path = os.path.join(calldir, 'calls', method, sampleid[:-2]+'.vcf.gz')
            if not os.path.exists(calltablemethod_path):
                print('calls for caller {} do not exist. path {} not found.'.format(method, calltablemethod_path))
                callmethod = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score'])
            else:
                callmethod = read_vcf(calltablemethod_path)
                if filter == 'PASS':
                    callmethod = callmethod[callmethod['FILTER'] == "PASS"]
                elif filter == 'REJECT':
                    callmethod = callmethod[callmethod['FILTER'] != "PASS"]
                samplen = callmethod.columns[-1]
                callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', samplen]]
                callmethod.columns = ['chrom', 'pos', 'ref', 'alt', method, 'INFO']
                callmethod['type'] = 'SNV'
                callmethod = callmethod[~callmethod['alt'].str.contains('[', regex=False)]
                print(callmethod)
                callmethod['totcov'] = callmethod['INFO'].str.split(':').str[6]
                callmethod['totcov'][callmethod['totcov'] == '.'] = 0
                callmethod['totcov'].fillna(0, inplace=True)
                callmethod['totcov'] = callmethod['totcov'].astype(int)
                callmethod['altcov'] = callmethod['INFO'].str.split(':').str[10]
                callmethod['altcov'][callmethod['altcov'] == '.'] = 0
                callmethod['altcov'].fillna(0, inplace=True)
                callmethod['altcov'] = callmethod['altcov'].astype(int)
                callmethod['vaf'] = callmethod['INFO'].str.split(':').str[12]
                callmethod['vaf'][callmethod['vaf'] == '.'] = 0
                callmethod['vaf'].fillna(0, inplace=True)
                callmethod['vaf'] = callmethod['vaf'].astype(float)
                callmethod[method+'_score'] = callmethod['vaf']
                callmethod[method][callmethod[method] == 'PASS'] = True
                callmethod[method][callmethod[method] == 'REJECT'] = False
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() > 0, 'type'] = 'INS'
                callmethod.loc[callmethod['alt'].str.len() - callmethod['ref'].str.len() < 0, 'type'] = 'DEL'
                callmethod['chrom'] = callmethod['chrom'].str[3:]
                callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score']]
        else:
            raise ValueError('caller {} is unknown'.format(method))
        # add calls of the given method to results dictonary
        callmethod['chrom_pos_ref_alt'] = callmethod['chrom'].astype('str').str.cat(callmethod['pos'].astype('str'), sep="_").str.cat(callmethod['ref'].astype('str'), sep='_').str.cat(callmethod['alt'].astype('str'), sep='_')
        callmethod.set_index('chrom_pos_ref_alt', inplace=True)
        callmethod = callmethod[['chrom', 'pos', 'ref', 'alt', 'type', 'totcov', 'altcov', 'vaf', method, method+'_score']]
        callmethod.rename(columns={'totcov': method+'_totcov', 'altcov': method+'_altcov', 'vaf': method+'_vaf'}, inplace=True)
        callmethod_snv = callmethod[callmethod['type'] == 'SNV']
        callmethod_indel = callmethod[(callmethod['type'] == 'INS') | (callmethod['type'] == 'DEL')]
        callmethod_snp = callmethod[callmethod['type'] == 'SNP']
        callmethods_snv[method] = callmethod_snv.drop_duplicates()
        callmethods_indel[method] = callmethod_indel.drop_duplicates()
        callmethods_snp[method] = callmethod_snp.drop_duplicates()
        callmethods_snv[method] = callmethods_snv[method].loc[~callmethods_snv[method].index.duplicated(keep='last')]
        callmethods_indel[method] = callmethods_indel[method].loc[~callmethods_indel[method].index.duplicated(keep='last')]
        if not callmethod_snp[method].empty:
           callmethods_snp[method] = callmethods_snp[method].loc[~callmethods_snp[method].index.duplicated(keep='last')]
        print("{}: {} SNV calls, {} INDEL calls, {} SNP calls".format(method, callmethods_snv[method].shape[0], callmethods_indel[method].shape[0], callmethods_snp[method].shape[0]))
    # concat results methods
    calltable_snv = pd.concat([cm.drop(['chrom', 'pos', 'ref', 'alt', 'type'], axis=1) for cm in list(callmethods_snv.values())], axis=1, sort=True)
    calltable_indel = pd.concat([cm.drop(['chrom', 'pos', 'ref', 'alt', 'type'], axis=1) for cm in list(callmethods_indel.values())], axis=1, sort=True)
    calltable_snp = pd.concat([cm.drop(['chrom', 'pos', 'ref', 'alt', 'type'], axis=1) for cm in list(callmethods_snp.values())], axis=1, sort=True)
    calltable_snv['chrom_pos_ref_alt'] = list(calltable_snv.index)
    calltable_indel['chrom_pos_ref_alt'] = list(calltable_indel.index)
    calltable_snp['chrom_pos_ref_alt'] = list(calltable_snp.index)
    calltable_snv[['chrom', 'pos', 'ref', 'alt']] = calltable_snv['chrom_pos_ref_alt'].str.split('_', expand=True)
    if not calltable_indel.empty:
        calltable_indel[['chrom', 'pos', 'ref', 'alt']] = calltable_indel['chrom_pos_ref_alt'].str.split('_', expand=True)
    else:
        calltable_indel = calltable_indel.reindex(columns=['chrom', 'pos', 'ref', 'alt', 'type'] + list(calltable_indel.columns)[:-1])
    if not calltable_snp.empty:
        calltable_snp[['chrom', 'pos', 'ref', 'alt']] = calltable_snp['chrom_pos_ref_alt'].str.split('_', expand=True)
    else:
        calltable_snp = calltable_snp.reindex(columns = ['chrom', 'pos', 'ref', 'alt', 'type'] + list(calltable_snp.columns)[:-1])
    calltable_snv['type'] = 'SNV'
    if not calltable_indel.empty:
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
    if not calltable_indel.empty:
        calltable_indel = calltable_indel[['chrom', 'pos', 'ref', 'alt', 'type'] + [m+suffix for m in methods for suffix in ['', '_score']] + [m+suffix for m in methods for suffix in ['_altcov', '_totcov', '_vaf']]]
    if not calltable_snp.empty:
        calltable_snp = calltable_snp[['chrom', 'pos', 'ref', 'alt', 'type'] + [m+suffix for m in methods  for suffix in ['', '_score']] + [m+suffix for m in methods for suffix in ['_altcov', '_totcov', '_vaf']]]
    else:
        calltable_snp = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'type'] + [m+suffix for m in methods  for suffix in ['', '_score']] + [m+suffix for m in methods for suffix in ['_altcov', '_totcov', '_vaf']])
    # correct for germline calls with GATK Haplotype
    calltablemethod_path = os.path.join(calldir, 'calls', 'bcbio', sampleid+'-N-germline-gatk-haplotype-annotated.vcf.gz')
    if '.' in os.path.basename(calltablemethod_path[:-7]):
        calltablemethod_path = os.path.join(os.path.dirname(calltablemethod_path), os.path.basename(calltablemethod_path)[:-7].replace('.', '_') + '.vcf.gz')
    if not os.path.exists(calltablemethod_path) or not gatkcorr:
        calltablemethod_path = os.path.join(calldir, 'calls', 'bcbio', sampleid+'-N-gatk-haplotype-annotated.vcf.gz')
        if '.' in os.path.basename(calltablemethod_path[:-7]):
            calltablemethod_path = os.path.join(os.path.dirname(calltablemethod_path), os.path.basename(calltablemethod_path)[:-7].replace('.', '_') + '.vcf.gz')
        if not os.path.exists(calltablemethod_path):
            print('calls for caller GATK Haplotype do not exist. path {} not found.'.format(calltablemethod_path))
            print('cannot use GATK Haplotype calls to filter germline calls')
    else:
        callmethod = read_vcf(calltablemethod_path)
        callmethod = callmethod[['CHROM', 'POS', 'REF', 'ALT', 'FILTER']]
        callmethod = callmethod[callmethod['FILTER'] == 'PASS']
        callmethod.columns = ['chrom', 'pos', 'ref', 'alt', 'gatk']
        callmethod['chrom_pos_ref_alt'] = callmethod['chrom'].astype('str').str.cat(callmethod['pos'].astype('str'), sep="_").str.cat(callmethod['ref'].astype('str'), sep='_').str.cat(callmethod['alt'].astype('str'), sep='_')
        callmethod.set_index('chrom_pos_ref_alt', inplace=True)
        beforegatk_snv,  beforegatk_indel, beforegatk_snp = calltable_snv.shape[0], calltable_indel.shape[0], calltable_snp.shape[0]
        print("# calls before using germline calls from GATK Haplotype: {} SNV, {} INDEL, {} SNP".format(beforegatk_snv, beforegatk_indel, beforegatk_snp))
        idx_snv = calltable_snv.index.difference(callmethod.index, sort=False)
        idx_snp_snv = calltable_snv.index.intersection(callmethod.index, sort=False)
        calltable_snp_snv = calltable_snv.loc[idx_snp_snv]
        calltable_snv = calltable_snv.loc[idx_snv]
        idx_indel = calltable_indel.index.difference(callmethod.index, sort=False)
        idx_snp_indel = calltable_indel.index.intersection(callmethod.index, sort=False)
        calltable_snp_indel = calltable_indel.loc[idx_snp_indel]
        calltable_indel = calltable_indel.loc[idx_indel]
        calltable_snp = pd.concat([calltable_snp, calltable_snp_snv, calltable_snp_indel])
        print("# calls after using germline calls from GATK Haplotype: {} SNV, {} INDEL, +{} SNP".format(
            calltable_snv.shape[0]-beforegatk_snv, calltable_indel.shape[0]-beforegatk_indel, calltable_snp.shape[0]-beforegatk_snp))
        if beforegatk_snv + beforegatk_indel + beforegatk_snp != calltable_snv.shape[0] + calltable_indel.shape[0] +  calltable_snp.shape[0]:
            raise ValueError('some calls are missed after using GATK Haplotype.')
    print('Final # calls SNV: {}, INDEL {}, SNP {} in dataframes of {} columns.'.format(calltable_snv.shape[0], calltable_indel.shape[0], calltable_snv.shape[0], calltable_snv.shape[1]))
    if save:
        if not os.path.exists(os.path.join(calldir, 'calls')):
            os.mkdir(os.path.join(calldir, 'calls'))
        calltable_snv.to_csv(os.path.join(calldir, 'calls', sampleid+'_snv_calls_'+filter+'.csv'))
        calltable_indel.to_csv(os.path.join(calldir, 'calls', sampleid+'_indel_calls_'+filter+'.csv'))
        calltable_snp.to_csv(os.path.join(calldir, 'calls', sampleid+'_snp_calls_'+filter+'.csv'))
    return calltable_snv, calltable_indel, calltable_snp


if __name__ == "__main__":
    import os
    from utils.config import Config

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")
    print(config.methods)

    calltable_snv, calltable_indel, calltable_snp = get_calltable(
        'data/mixtures/mixtures_chr3/mixtures_chr3_CRC-1014_180816-CW-T_CRC-1014_090516-CW-T/mixture_chr3_CRC-1014_180816-CW-T_50x_CRC-1014_090516-CW-T_100x',
        config.methods, save=False, filter='all')
    print(calltable_snv)
    res = calltable_snv[[m+'_vaf' for m in config.methods]].min(axis=1)
    print(res[res > 1])
    print(res)
    print(res.describe())
