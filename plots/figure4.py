

if __name__ == "__main__":
    import os
    import matplotlib.pyplot as plt
    import warnings
    from sklearn.metrics import precision_recall_curve, f1_score, average_precision_score
    warnings.filterwarnings('ignore')

    # set working directory
    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    from utils.config import Config
    from utils.viz import *
    from benchmark.table import *
    from benchmark.metrics import *
    from benchmark.calltable import *
    from benchmark.calltableseries import *
    from benchmark.groundtruth import *
    from benchmark.metricsseries import *

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

    """

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
            gtm = 4
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        print(refname)
        # for metric in metrics:
        metric = 'auprc'
        #metric = 'recall'
        #metric = 'maxrecallatleast0_03precision'
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
        if metric == 'auprc':
            #plt.ylim([0, .5])
            plt.ylim([0, .10])
            plt.xlim([5, 1])
            if mt == 'indel':
                plt.ylim([0, .4])
        else:
            plt.ylim([0, 1.])
        if not os.path.exists(os.path.join(*config.outputpath, 'figure2b')):
            os.mkdir(os.path.join(*config.outputpath, 'figure2b'))
        plt.savefig(os.path.join(*config.outputpath, 'figure2b', 'perf_'+metric+'_986_150x_exomecalling_'+fixedvar+'_'+mt+'_'+str(gtm)+'callers_zoom.svg'), bbox_inches='tight')
        plt.show()

    ###### Ultra deep test: 2,000x WES exome calling #########
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
            gtm = 4
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        print(refname)
        # for metric in metrics:
        metric = 'auprc'
        #metric = 'recall'
        #metric = 'maxrecallatleast0_03precision'
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
        if metric == 'auprc':
            #plt.ylim([0, .25])
            plt.ylim([0, .05])
            plt.xlim([5, 1])
        else:
            plt.ylim([0, 1.])
        if not os.path.exists(os.path.join(*config.outputpath, 'figure2b')):
            os.mkdir(os.path.join(*config.outputpath, 'figure2b'))
        plt.savefig(os.path.join(*config.outputpath, 'figure2b', 'perf_'+metric+'_986_2000x_exomecalling_'+fixedvar+'_'+mt+'_zoom.svg'), bbox_inches='tight')
        plt.show()
    """

    ####### 150x WGS whole genome calling #######
    for fixedvar in fixedvars:
        if fixedvar == 'coverage':
            xaxis = 'tumor burden'
        elif fixedvar == 'ctdna':
            xaxis = 'coverage'
        mt = 'snv'
        if mt == 'snv':
            gtm = 5
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        else:
            gtm = 4
            refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
        print(refname)
        #for metric in metrics:
        metric = 'auprc'
        #metric = 'recall'
        #metric = 'maxrecallatleast0_03precision'
        # load results tables
        restables = {'snv': [], 'indel': []}
        # for mixtureid in mixtureids:
        plasmasample = '_'.join(mixtureid.split('_')[:2])
        print(mixtureid, plasmasample)
        xa = xaxis if xaxis != 'tumor burden' else 'tb'
        print(xa)
        restable = pd.read_csv(os.path.join(*config.mixturefolderwholegenome, 'mixtures_allchr', 'results', mixtureid+'_'+mt+'_'+metric+'_'+refname+'_fixed'+fixedvar+'_'+ xa +'.csv'), index_col=0,  memory_map=True)
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
        if metric == 'auprc':
            plt.ylim([0, .25])
        else:
            plt.ylim([0, 1.])
        if not os.path.exists(os.path.join(*config.outputpath, 'figure2b')):
            os.mkdir(os.path.join(*config.outputpath, 'figure2b'))
        plt.savefig(os.path.join(*config.outputpath, 'figure2b', 'perf_'+metric+'_986_150x_wholegenomecalling_'+fixedvar+'_'+mt+'.svg'), bbox_inches='tight')
        plt.show()


    """
    ############# comp ground truths 150x and 2000x ##############

    A = pd.read_csv('figures/figure2b/gt_986_exome_150x_atleast5callersinundilutedsample_snv.csv', index_col=0)
    B = pd.read_csv('figures/figure2b/gt_986_exome_2000x_atleast5callersinundilutedsample_snv.csv', index_col=0)
    print(A.shape[0], B.shape[0])
    ab = list(set(set(list(A.index)) & set(list(B.index))))
    len(ab)
    Awithoutab = A.loc[list(set(A.index) - set(ab))]
    Bwithoutab = B.loc[list(set(B.index) - set(ab))]

    plt.figure(figsize=(10,10))

    calltablesseries = pd.read_csv(os.path.join('data', 'mixtures', 'mixtures_allchr', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T_snv_calls_all.csv'), index_col=0, memory_map=True)
    aux = pd.read_csv(os.path.join('data', 'mixtures_ultradeep', 'mixtures_chrall', 'mixtures_chrall_CRC-986_100215-CW-T_CRC-986_300316-CW-T' , 'calls', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T_tf_cov.csv'), index_col=0)
    ac = pd.DataFrame(calltablesseries.loc[[a for a in list(Bwithoutab.index) if a in calltablesseries.index]][['{:.2f}_{}_vaf'.format(aux['tf'].max(), m) for m in config.methods if ('{:.2f}_{}'.format(aux['tf'].max(), m) in calltablesseries.columns) and (m != 'smurf')]].median(skipna=True, axis=1))
    ac = ac[~ac.index.duplicated()]
    Bwithoutab = Bwithoutab[~Bwithoutab.index.duplicated()]
    comp = pd.concat([ac, Bwithoutab], axis=1)
    comp.columns = ['vaf 150x', 'vaf 2000x']
    comp.fillna(0, inplace=True)
    list_only2000x = list(comp.index)
    print(Bwithoutab.shape[0], ac.shape[0])
    sns.histplot(x='vaf 150x', y='vaf 2000x', data=comp, binwidth=0.01, binrange=[0,.5], alpha=1, color='tab:blue', label='2000x only')

    comp = pd.concat([A.loc[ab], B.loc[ab]], axis=1)
    comp.columns = ['vaf 150x', 'vaf 2000x']
    comp.fillna(0, inplace=True)
    list_both = list(comp.index)
    sns.histplot(x='vaf 150x', y='vaf 2000x', data=comp, binwidth=0.01, binrange=[0,.5], alpha=1, color='tab:red', label='both')

    table_ultradeep = pd.read_csv(os.path.join('data', 'mixtures_ultradeep', 'mixtures_chrall', 'mixtures_chrall_CRC-986_100215-CW-T_CRC-986_300316-CW-T' , 'calls', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T_snv_calls_all.csv'), index_col=0, memory_map=True)
    aux_ultradeep = pd.read_csv(os.path.join('data', 'mixtures_ultradeep', 'mixtures_chrall', 'mixtures_chrall_CRC-986_100215-CW-T_CRC-986_300316-CW-T' , 'calls', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T_tf_cov.csv'), index_col=0)
    ac = calltablesseries.loc[[a for a in list(Awithoutab.index) if a in table_ultradeep.index]][['{:.2f}_{}_vaf'.format(aux_ultradeep['tf'].max(), m) for m in config.methods if ('{:.2f}_{}'.format(aux_ultradeep['tf'].max(), m) in table_ultradeep.columns) and (m != 'smurf')]].median(skipna=True, axis=1)
    comp = pd.concat([Awithoutab, ac], axis=1)
    comp.columns = ['vaf 150x', 'vaf 2000x']
    comp.fillna(0, inplace=True)
    list_only150x = list(comp.index)
    print(Awithoutab.shape[0], ac.shape[0])
    sns.histplot(x='vaf 150x', y='vaf 2000x', data=comp, binwidth=0.01, binrange=[0,.5], alpha=1, color='tab:olive', label='150x only')

    from matplotlib.lines import Line2D
    a = Line2D([0], [0], color='tab:red', lw=4)
    b = Line2D([0], [0], color='tab:blue', lw=4)
    c = Line2D([0], [0], color='tab:olive', lw=4)

    plt.legend([a, b,  c], ['both', '2000x only', '150x only'], bbox_to_anchor=(1,1), loc="upper left")
    plt.savefig(os.path.join(*config.outputpath, 'figure2b', 'gt_vaf_150x_vs_2000x.svg'), bbox_inches='tight')
    plt.show()
    
    """


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



