
import warnings
warnings.filterwarnings('ignore')

from benchmark.metrics import *
from benchmark.calltable import *


def generate_groundtruth(config, calltablesseries, calltablestf, ground_truth_method=5, muttype='snv', matchedtissuepath=None, methods=None):
    if methods is not None:
        refmethods = methods
    else:
        refmethods = list(np.copy(config.methods))
    print(refmethods)
    # Approach 1: CONSENSUS
    # pseudo ground truth = mutations found by at least k callers
    if type(ground_truth_method) == int:
        calltablesseries['truth'] = False
        #refmethods = list(np.copy(config.methods))
        print(refmethods)
        if '{:.2f}_{}'.format(max(calltablestf), refmethods[0]) in list(calltablesseries.columns):
            print('test')
            ncallsinundiluted = calltablesseries[['{:.2f}_{}'.format(max(calltablestf), m) for m in refmethods]].sum(axis=0)
            callsinundiluted = calltablesseries[['{:.2f}_{}_score'.format(max(calltablestf), m) for m in refmethods]]
        else:
            ncallsinundiluted = calltablesseries[calltablesseries['sampletf'] == max(calltablestf)][refmethods].sum(axis=0)
            callsinundiluted = calltablesseries[calltablesseries['sampletf'] == max(calltablestf)][[m+'_score' for m in refmethods]]
        print(ncallsinundiluted)
        callsinundiluted.columns = refmethods
        callsinundiluted = callsinundiluted.stack().reset_index(level=0, drop=False).reset_index()
        callsinundiluted.set_index('chrom_pos_ref_alt', inplace=True)
        callsinundiluted.columns = ['method', 'score']
        for mi, m in enumerate(refmethods):
            plt.figure()
            sns.histplot(callsinundiluted[callsinundiluted['method'] == m], x='score', stat="probability", color=config.colors[config.methods.index(m)], binwidth=0.01)
            plt.xlim([0, 1])
            plt.ylim([0, 1])
            plt.title(m)
        if '{:.2f}_{}'.format(max(calltablestf), refmethods[0]) in list(calltablesseries.columns):
            truthpos = list(calltablesseries[calltablesseries[['{:.2f}_{}'.format(max(calltablestf), m) for m in refmethods]].sum(axis=1) >= ground_truth_method].index)
        else:
            aux = calltablesseries[calltablesseries['sampletf'] == max(calltablestf)]
            print(aux[refmethods])
            truthpos = list(aux[aux[refmethods].sum(axis=1) >= ground_truth_method].index)
        calltablesseries.loc[truthpos, 'truth'] = True
    elif ground_truth_method == 'spikein':
        chroms = list(calltablesseries.chrom.unique().astype(str))
        if chroms in [str(c) for c in range(1, 23)]:
            gt = pd.read_csv(os.path.join(*config.extdatafolder, 'cosmic_mutations_atleast5patients', 'CRC_chr'+str(chroms)+'_'+muttype.upper()+'_tf1.bed'), sep='\t', header=None)
        else:
            gt_list = [pd.read_csv(os.path.join(*config.extdatafolder, 'cosmic_mutations_atleast5patients', 'CRC_chr'+str(chrom)+'_'+muttype.upper()+'_tf1.bed'), sep='\t', header=None) for chrom in chroms]
            gt = pd.concat(gt_list)
            print(gt)
            calltablesseries['truth'] = False
        if muttype == 'snv':
            gt.columns = ['chrom', 'startpos', 'endpos', 'vaf', 'alt']
        else:  # indel
            gt.columns = ['chrom', 'startpos', 'endpos', 'vaf', 'type', 'alt']
        truthpos = []
        posvalues = calltablesseries['pos'].values
        c = 0
        for igt, gtv in enumerate(gt['startpos'].values):
            if gtv in posvalues:
                if calltablesseries[calltablesseries['pos'] == gtv].shape[0] > 1:
                    print('ISSUE: cannot retrieve reference easily')
                    print(calltablesseries[calltablesseries['pos'] == gtv].shape[0])
                # print(gt)
                truthpos.append(str(gt.iloc[igt]['chrom']) +'_'+ str(gt.iloc[igt]['startpos']) +'_'+calltablesseries[calltablesseries['pos'] == gtv]['ref'].values[0] +'_'+str(gt.iloc[igt]['alt']))
            c +=1
        print(c)
        calltablesseries['truth'] = False
        calltablesseries = calltablesseries.reindex(list(set(list(calltablesseries.index) + truthpos)))
        calltablesseries['truth'].loc[truthpos] = True
    # Approach 2: rank mutations
    # pseudo ground truth = best K mutations found by each caller
    # number mutations found by each method
    elif ground_truth_method == 'ranked':
        # calltablesseries['truth'] = False
        if muttype == 'indel':
            refmethods.remove('abemus')  # not an indel method
            refmethods.remove('cfsnv')  # not an indel method
        ncallsinundiluted = calltablesseries[['{:.2f}_{}'.format(max(calltablestf), m) for m in refmethods]].sum(axis=0)
        callsinundiluted = calltablesseries[['{:.2f}_{}_score'.format(max(calltablestf), m) for m in refmethods]]
        print(ncallsinundiluted)
        callsinundiluted.columns = refmethods
        callsinundiluted = callsinundiluted.stack().reset_index(level=0, drop=False).reset_index()
        callsinundiluted.set_index('chrom_pos_ref_alt', inplace=True)
        callsinundiluted.columns = ['method', 'score']
        for mi, m in enumerate(refmethods):
            plt.figure()
            sns.histplot(callsinundiluted[callsinundiluted['method'] == m], x='score', stat="probability", color=config.colors[config.methods.index(m)], binwidth=0.01)
            plt.xlim([0, 1])
            plt.ylim([0, 1])
            plt.title(m)
        ncallsinundiluted = calltablesseries[['{:.2f}_{}'.format(max(calltablestf), m) for m in refmethods]].sum(axis=0)
        print(ncallsinundiluted)
        ncallsinundiluted = calltablesseries[['{:.2f}_{}'.format(max(calltablestf), m) for m in refmethods]].sum(axis=0)
        ncallsinundiluted = ncallsinundiluted.max()/ncallsinundiluted
        ncallsinundiluted = ncallsinundiluted/ncallsinundiluted.max()
        print(ncallsinundiluted.shape)
        print(ncallsinundiluted)
        callsinundiluted = calltablesseries[['{:.2f}_{}_score'.format(max(calltablestf), m) for m in refmethods]]
        callsinundiluted.sort_values(by=list(callsinundiluted.columns), ascending=False).head()
        #callsinundiluted.head()
        for c in callsinundiluted.columns:
            callsinundiluted[c] *= callsinundiluted[c] * ncallsinundiluted.loc[c[:-6]]
        callsinundiluted = callsinundiluted.fillna(0)
        callsinundiluted['score'] = callsinundiluted[['{:.2f}_{}_score'.format(max(calltablestf), m) for m in refmethods]].sum(axis=1) / len(refmethods)
        print(callsinundiluted['score'].describe())
        print(callsinundiluted[callsinundiluted['score'] > 1/(len(refmethods))].shape)
        plt.figure()
        sns.histplot(callsinundiluted[callsinundiluted['score'] > 1/(len(refmethods))], x='score', stat='probability', bins=10)
        plt.xlim([0, 1])
        plt.title('score')
        callsinundiluted.head()
        calltablesseries['truth'] = False
        truthpos = list(callsinundiluted[callsinundiluted['score'] > 1/(len(config.methods))].index)
        calltablesseries.loc[truthpos, 'truth'] = True

        print(calltablesseries[calltablesseries['truth'] == True][['{:.2f}_{}'.format(max(calltablestf), m) for m in config.methods]].sum(axis=1).value_counts())

    # Approach 3: tissue as matched ground truth
    elif ground_truth_method == 'tissue':
        print('tissue')
        calltablesseries['truth'] = False
        tissuetable = pd.read_csv(matchedtissuepath, index_col=0)  # 10% VAF filter fine for tissue
        # which chroms as represented in cfDNA
        print(tissuetable['chrom'])
        print(calltablesseries['chrom'].values)
        calltablesseries['chrom'] = calltablesseries['chrom'].astype(str)
        tissuetable['chrom'] = tissuetable['chrom'].astype(str)
        chromlist = list(set(calltablesseries['chrom'].values))
        #try:
        #    chromlist = [str(c) for c in chromlist]
        #except:
        #    print(chromlist)
        tissuetable = tissuetable[tissuetable['chrom'].isin(chromlist)]
        refmethods = list(np.copy(config.methods_tissue))
        print(refmethods)
        ncallsinundiluted = tissuetable[['{}'.format(m) for m in refmethods]].sum(axis=0)
        callsinundiluted = tissuetable[['{}_score'.format(m) for m in refmethods]]
        print(ncallsinundiluted)
        callsinundiluted.columns = refmethods
        callsinundiluted = callsinundiluted.stack().reset_index(level=0, drop=False).reset_index()
        callsinundiluted.set_index('chrom_pos_ref_alt', inplace=True)
        callsinundiluted.columns = ['method', 'score']
        for mi, m in enumerate(refmethods):
            plt.figure()
            sns.histplot(callsinundiluted[callsinundiluted['method'] == m], x='score', stat="probability", color=config.colors[config.methods.index(m)], binwidth=0.01)
            plt.xlim([0, 1])
            plt.ylim([0, 1])
            plt.title(m)
        truthpos = list(tissuetable[tissuetable[['{}'.format(m) for m in refmethods]].sum(axis=1) >= 3].index)  # 3/5 callers on tissue
        print(len(truthpos))
        print(calltablesseries.shape[0])
        calltablesseries.reindex[truthpos, 'truth'] = True
        print(calltablesseries.shape[0])
        print(calltablesseries['truth'].sum())
        #if calltablesseries.index.duplicated().sum() > 0:
        #    truthpos = [t for t in truthpos if t in list(calltablesseries.index)] # take intersection
        #print(len(truthpos))
        #calltablesseries = calltablesseries.reindex(truthpos)
        #calltablesseries.loc[truthpos, 'truth'] = True
        #truthposout = [t for t in truthpos if t not in list(calltablesseries.index)] # take missed
        #for tp in truthposout:
        #    calltablesseries.loc[tp, 'truth'] = True
    print(calltablesseries['truth'].value_counts())
    if '{:.2f}_{}'.format(max(calltablestf), refmethods[0]) in list(calltablesseries.columns):
        print(calltablesseries[calltablesseries['truth'] == True][['{:.2f}_{}'.format(max(calltablestf), m) for m in config.methods]].sum(axis=0))
    else:
        calltablesseries[(calltablesseries['truth'] == True) & (calltablesseries['sampletf'] == max(calltablestf))][refmethods].sum(axis=0)
    return calltablesseries


def compare_groundtruth(calltabledict):
    res = {}
    for j in range(1, 8): # 7 methods
        print(j)
        c = 0
        for gttype, calltable in calltabledict.items():
            print(gttype)
            refmethods = np.unique([r.split('_')[0] for r in calltable.columns[5:]])
            if c == 0:
                refcalltable = calltable[calltable[['{}'.format(m) for m in refmethods]].sum(axis=1) >= j]
            print(refmethods)
            print(calltable.shape[0])
            for i in range(1, len(refmethods)+1):
                aux = list(calltable[calltable[['{}'.format(m) for m in refmethods]].sum(axis=1) >= i].index)
                auxf = [t for t in aux if t in list(refcalltable.index)]
                print(i, len(aux), len(auxf))
                res[gttype + '_' + str(j) + '_' + str(i)] = aux
            c += 1
    return res


def compare_groundtruth_venn(calltabledict):
    res = {}
    for j in range(0, 8): # 7 methods
        print(j)
        c = 0
        for gttype, calltable in calltabledict.items():
            print(gttype)
            refmethods = np.unique([r.split('_')[0] for r in calltable.columns[5:]])
            if c == 0:
                refcalltable = calltable[calltable[['{}'.format(m) for m in refmethods]].sum(axis=1) >= j]
            print(refmethods)
            print(calltable.shape[0])
            for i in range(1, len(refmethods)+1):
                aux = list(calltable[calltable[['{}'.format(m) for m in refmethods]].sum(axis=1) >= i].index)
                auxf = [t for t in aux if t in list(refcalltable.index)]
                print(i, len(aux), len(auxf))
                res[gttype + '_' + str(j) + '_' + str(i)] = aux
            c += 1
    return res







if __name__ == "__main__":
    import os
    # set working directory
    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    # Config and Display paramaters
    from utils.config import Config
    config = Config("config/", "config_viz.yaml")
    set_display_params(config)
    print(config.methods)

    from benchmark.calltableseries import get_calltableseries

    chrom = '22'
    reload = False
    save = True
    filterparam = 'all'
    muttypes = ['snv', 'indel']
    mixtureids = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T']

    for muttype in muttypes:
        if muttype == 'snv':
            nref = 5
        else:  # elif muttype == 'indel':
            nref = 3
        for mixtureid in mixtureids:
            print('############# {} {} ############'.format(mixtureid, muttype))
            calltablesseries, calltables = get_calltableseries(config, mixtureid, chrom, muttype, filterparam, reload, save)
            print(calltablesseries.head())
            print(calltables['tf'])
            calltablesseries = generate_groundtruth(config, calltablesseries, calltables['tf'], ground_truth_method=nref, muttype=muttype)
            calltablesseries = generate_groundtruth(config, calltablesseries, calltables['tf'], ground_truth_method='ranked', muttype=muttype)
