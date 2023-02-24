

if __name__ == "__main__":
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import warnings
    from sklearn.metrics import precision_recall_curve, f1_score, average_precision_score
    warnings.filterwarnings('ignore')

    # set working directory
    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    from utils.config import Config
    from utils.viz import *
    from utils.table import *
    from utils.metrics import *
    from utils.calltable import *
    from utils.calltableseries import *
    from utils.groundtruth import *
    from utils.metricsseries import *

    # Config and Display paramaters

    config = Config("config/", "config_viz.yaml")
    set_display_params(config)
    print(config.methods)

    # Chomosome
    mixtureids = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T']
    mixtureid = 'CRC-986_100215-CW-T_CRC-986_300316-CW-T'
    reload = False
    save = True
    fixedvars=['coverage', 'ctdna']
    filterparam = 'all'

    markers = ['o', '^', 'X']
    linestyles = ['-', '-', '-']
    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}

    muttypes = ['snv', 'indel']
    metrics = ['auprc', 'precision', 'recall']

    chrom = 'all'

    ####### 150x WGS exome calling #######

    for fixedvar in fixedvars:
        if fixedvar == 'coverage':
            xaxis = 'tumor burden'
        elif fixedvar == 'ctdna':
            xaxis = 'coverage'
        mt = 'snv'
        if mt == 'snv':
            gtm = 5
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        else:  # elif mt == 'indel':
            gtm = 3
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        print(refname)
        # for metric in metrics:
        metric = 'auprc'
        # load results tables
        restables = {'snv': [], 'indel': []}
        mixtureid = 'CRC-986_100215-CW-T_CRC-986_300316-CW-T'
        #if mixtureid == 'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
        #    gtm = 3
        #    refname = 'intissuesamplebyatleast'+str(gtm)+'callers'
        #else:
        #    gtm = 4
        #    refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        plasmasample = '_'.join(mixtureid.split('_')[:2])
        print(mixtureid, plasmasample)
        xa = xaxis if xaxis != 'tumor burden' else 'tb'
        print(xa)
        restable = pd.read_csv(os.path.join(*config.mixturefolder, 'mixtures_allchr', 'results', mixtureid+'_'+mt+'_'+metric+'_'+refname+'_fixed'+fixedvar+'_'+ xa +'.csv'), index_col=0)
        restable['plasma sample'] = plasmasample
        restables[mt].append(restable)
        restables[mt] = pd.concat(restables[mt])
        res1 = plot_metricsseries(config, restables, mixtureids, 'all', metric=metric, muttype=mt,
                                  ground_truth_method='mixture', fixedvar=fixedvar, refname=refname, allpatients=True, logscale=False, save=False)
        if fixedvar == 'coverage':
            plt.xlim([22, 0])
            ax = plt.gca()
            # Major ticks every 20, minor ticks every 5
            major_ticks = np.arange(20, -1, -5)
            minor_ticks = np.arange(22, -1, -1.)
            ax.set_xticks(major_ticks)
            ax.set_xticks(minor_ticks, minor=True)
            # Or if you want different settings for the grids:
            ax.grid(which='minor', alpha=0.2)
            ax.grid(which='major', alpha=1.)
        else:
            plt.grid(linewidth=1)
        plt.ylim([0, .5])
        if not os.path.exists(os.path.join(*config.outputpath, 'figure2b')):
            os.mkdir(os.path.join(*config.outputpath, 'figure2b'))
        plt.savefig(os.path.join(*config.outputpath, 'figure2b', 'perf_auprc_986_150x_exomecalling_'+fixedvar+'_'+mt+'_'+str(gtm)+'callers.svg'), bbox_inches='tight')
        plt.show()

   ###### Ultra deep test: 2,000x WES exome calling #########

    mixtureids =  ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T']

    for fixedvar in fixedvars:
        if fixedvar == 'coverage':
            xaxis = 'tumor burden'
        elif fixedvar == 'ctdna':
            xaxis = 'coverage'
        mt = 'snv'
        if mt == 'snv':
            gtm = 5
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
            #if mixtureid == 'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
            #    gtm = 3
            #    refname = 'intissuesamplebyatleast'+str(gtm)+'callers'
            #else:
            #    gtm = 3
            #    refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        else:  # elif mt == 'indel':
            gtm = 3
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        print(refname)
        # for metric in metrics:
        metric = 'auprc'
        # load results tables
        restables = {'snv': [], 'indel': []}
        # for mixtureid in mixtureids:
        plasmasample = '_'.join(mixtureid.split('_')[:2])
        print(mixtureid, plasmasample)
        xa = xaxis if xaxis != 'tumor burden' else 'tb'
        print(xa)
        restable = pd.read_csv(os.path.join(*config.mixturefolderultradeep, 'mixtures_allchr', 'results', mixtureid+'_'+mt+'_'+metric+'_'+refname+'_fixed'+fixedvar+'_'+ xa +'.csv'), index_col=0)
        restable['plasma sample'] = plasmasample
        restables[mt].append(restable)
        restables[mt] = pd.concat(restables[mt])
        res1 = plot_metricsseries(config, restables, mixtureids, 'all', metric=metric, muttype=mt,
                                  ground_truth_method='mixture', fixedvar=fixedvar, refname=refname, allpatients=True, logscale=False, save=False)

        if fixedvar == 'coverage':
            plt.xlim([22, 0])
            ax = plt.gca()
            # Major ticks every 20, minor ticks every 5
            major_ticks = np.arange(20, -1, -5)
            minor_ticks = np.arange(22, -1, -1.)
            ax.set_xticks(major_ticks)
            ax.set_xticks(minor_ticks, minor=True)
            # Or if you want different settings for the grids:
            ax.grid(which='minor', alpha=0.2)
            ax.grid(which='major', alpha=1.)
        else:
            plt.grid(linewidth=1)
        plt.ylim([0, 0.85])
        if not os.path.exists(os.path.join(*config.outputpath, 'figure2b')):
            os.mkdir(os.path.join(*config.outputpath, 'figure2b'))
        #plt.savefig(os.path.join(*config.outputpath, 'figure2b', 'perf_auprc_986_2000x_exomecalling_'+fixedvar+'_'+mt+'.svg'), bbox_inches='tight')
        plt.show()

    ############ 150x WGS whole genome calling #############

    for fixedvar in fixedvars:
        print(fixedvar)
        if fixedvar == 'coverage':
            xaxis = 'tumor burden'
        elif fixedvar == 'ctdna':
            xaxis = 'coverage'
        mt = 'snv'
        if mt == 'snv':
            gtm = 5
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        else:  # elif mt == 'indel':
            gtm = 3
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        print(refname)
        # for metric in metrics:
        metric = 'auprc'
        # load results tables
        restables = {'snv': [], 'indel': []}
        mixtureid = 'CRC-986_100215-CW-T_CRC-986_300316-CW-T'
        #if mixtureid == 'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
        #    gtm = 3
        #    refname = 'intissuesamplebyatleast'+str(gtm)+'callers'
        #else:
        #    gtm = 4
        #    refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        plasmasample = '_'.join(mixtureid.split('_')[:2])
        print(mixtureid, plasmasample)
        xa = xaxis if xaxis != 'tumor burden' else 'tb'
        print(xa)
        restable = pd.read_csv(os.path.join(*config.mixturefolderwholegenome, 'mixtures_allchr', 'results', mixtureid+'_'+mt+'_'+metric+'_'+refname+'_fixed'+fixedvar+'_'+ xa +'.csv'), index_col=0)
        restable['plasma sample'] = plasmasample
        restables[mt].append(restable)
        restables[mt] = pd.concat(restables[mt])
        res1 = plot_metricsseries(config, restables, mixtureids, 'all', metric=metric, muttype=mt,
                                  ground_truth_method='mixture', fixedvar=fixedvar, refname=refname, allpatients=True, logscale=False, save=False, diltype='mixtures_wgs')
        plt.grid(linewidth=1)
        plt.ylim([0, 0.85])
        if not os.path.exists(os.path.join(*config.outputpath, 'figure2b')):
            os.mkdir(os.path.join(*config.outputpath, 'figure2b'))
        #plt.savefig(os.path.join(*config.outputpath, 'figure2b', 'perf_auprc_986_150x_chr22calling_'+fixedvar+'_'+mt+'.svg'), bbox_inches='tight')
        plt.show()


