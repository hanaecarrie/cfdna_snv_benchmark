
if __name__ == "__main__":

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import warnings
    from sklearn.metrics import precision_recall_curve, f1_score, average_precision_score
    warnings.filterwarnings('ignore')
    from utils.config import Config
    from utils.viz import *
    from utils.table import *
    from utils.metrics import *
    from utils.calltable import *
    from utils.calltableseries import *
    from utils.groundtruth import *
    from utils.metricsseries import *

    # set working directory
    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    # Config and Display paramaters
    config = Config("config/", "config_viz.yaml")
    set_display_params(config)
    print(config.methods)

    # Chomosome
    mixtureids = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T',
                  'CRC-986_100215-CW-T_CRC-986_300316-CW-T',
                  'CRC-123_310715-CW-T_CRC-123_121115-CW-T']
    reload = False
    save = True
    fixedvars = ['coverage', 'ctdna']
    filterparam = 'all'

    markers = ['o', '^', 'X']
    linestyles = ['-', '-', '-']
    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}

    muttypes = ['snv', 'indel']
    metrics = ['auprc', 'precision', 'recall']

    chrom = 'all'

    for fixedvar in fixedvars:
        if fixedvar == 'coverage':
            xaxis = 'tumor burden'
        elif fixedvar == 'ctdna':
            xaxis = 'coverage'
        #for mt in muttypes:
        mt = 'indel'
        print('####### ' + mt + ' #######')
        if mt == 'snv':
            gtm = 4
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        else:  # elif mt == 'indel':
            gtm = 2
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        print(refname)
        # for metric in metrics:
        metric = 'auprc'
        # load results tables
        restables = {'snv': [], 'indel': []}
        for mixtureid in mixtureids:
            if mixtureid == 'CRC-986_100215-CW-T_CRC-986_300316-CW-T' and mt == 'snv':
                gtm = 3
                refname = 'intissuesamplebyatleast'+str(gtm)+'callers'
            else:
                if mt == 'snv':
                    gtm = 4
                else:
                    gtm = 2
                refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
            plasmasample = '_'.join(mixtureid.split('_')[:2])
            print(mixtureid, plasmasample)
            xa = xaxis if xaxis != 'tumor burden' else 'tb'
            restable = pd.read_csv(os.path.join(*config.mixturefolder, 'mixtures_allchr', 'results', mixtureid + '_'
                                                + mt + '_'+metric+'_' + refname+'_fixed' + fixedvar + '_' + xa + '.csv'),
                                   index_col=0)
            restable['plasma sample'] = plasmasample
            restables[mt].append(restable)
        restables[mt] = pd.concat(restables[mt])
        res1 = plot_metricsseries(config, restables, mixtureids, 'all', metric=metric, muttype=mt,
                                  ground_truth_method='mixture', fixedvar=fixedvar, refname=refname,
                                  allpatients=True, logscale=False, save=False)
        resx = np.array([rx.values for rx in res1['x']])

        resx.mean(axis=0)
        resy = np.array([ry.values for ry in res1['y']])
        lmaux = list(set(config.methods) & set(restable['caller'].unique()))
        lm = config.methods[:]
        for lmi in lm:
            if lmi not in lmaux:
                lm.remove(lmi)
        if mt == 'indel':
            lm.remove('abemus')
        res = {m: [] for m in lm}
        for mi, m in enumerate(lm):
            resmean = np.mean([resy[mi], resy[mi+len(lm)], resy[mi+2*len(lm)]], axis=0)
            resstd = np.std([resy[mi], resy[mi+len(lm)], resy[mi+2*len(lm)]], axis=0)
            res[m] = [resmean, resstd, resx.mean(axis=0)]

        color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}

        plt.figure(figsize=(15, 10))
        plt.grid(linewidth=1)
        plt.grid()
        for m in lm:
            plt.errorbar(res[m][2], res[m][0], xerr = resx.std(axis=0), yerr=res[m][1], marker='s', c=color_dict[m],
                         label=m,  markersize=15, lw=2, fmt='-s')
        ax = plt.gca()
        if fixedvar == 'coverage':
            plt.gca().invert_xaxis()
            xlab='tumor burden (%)'
            plt.xlim([27, 0])
            # Major ticks every 20, minor ticks every 5
            major_ticks = np.arange(25, -1, -5)
            minor_ticks = np.arange(27, -1, -1.)
            ax.set_xticks(major_ticks)
            ax.set_xticks(minor_ticks, minor=True)
            # Or if you want different settings for the grids:
            ax.grid(which='minor', alpha=0.2)
            ax.grid(which='major', alpha=1.)
        else:
            xlab = 'effective coverage (x) increased by added noise'
            plt.xlim([60,260])
            plt.grid()
        hand, labl = ax.get_legend_handles_labels()
        ax.legend(hand, labl, bbox_to_anchor=(1, 1), loc="upper left")
        plt.xlabel(xlab)
        plt.ylabel(metric.upper()+' score')
        plt.title(metric.upper() + " score for {} calling in chr{} with ref {}".format(mt.upper(), chrom, refname),
                  pad=50)
        if mt == 'snv':
            plt.ylim([0, .9])
        if not os.path.exists(os.path.join(*config.outputpath, 'figure2a')):
            os.mkdir(os.path.join(*config.outputpath, 'figure2a'))
        plt.savefig(os.path.join(*config.outputpath, 'figure2a',
                                 'perf_auprc_3patients_150x_'+fixedvar+'_'+mt+'.svg'), bbox_inches='tight')
