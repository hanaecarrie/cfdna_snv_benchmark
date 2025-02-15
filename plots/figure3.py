
if __name__ == "__main__":

    import matplotlib.pyplot as plt
    import warnings
    from sklearn.metrics import precision_recall_curve, f1_score, average_precision_score
    warnings.filterwarnings('ignore')
    from utils.config import Config
    from utils.viz import *
    from benchmark.table import *
    from benchmark.metrics import *
    from benchmark.calltable import *
    from benchmark.calltableseries import *
    from benchmark.groundtruth import *
    from benchmark.metricsseries import *

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

    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}

    for fixedvar in fixedvars:
        if fixedvar == 'coverage':
            xaxis = 'tumor burden'
        elif fixedvar == 'ctdna':
            xaxis = 'coverage'
        for mt in muttypes:
            #mt = 'snv'
            print('####### ' + mt + ' #######')
            if mt == 'snv':
                gtm = 5
                refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
            else:  # elif mt == 'indel':
                gtm = 4
                refname = 'inundilutedsamplebyatleast'+str(gtm)+'callers'
            print(refname)
            # for metric in metrics:
            #metric = 'auprc'
            #metric = 'recall'
            metric = 'maxrecallatleast0_03precision'
            # load results tables
            restables = {'snv': [], 'indel': []}
            for mixtureid in mixtureids:
                if mt == 'snv':
                    gtm = 5
                else:
                    gtm = 4
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
            #tfdict = {'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T': [40.3, 3.2],
            #          'CRC-123_310715-CW-T_CRC-123_121115-CW-T': [62.1, 2.2],
            #          'CRC-986_100215-CW-T_CRC-986_300316-CW-T': [32.6, 0]}
            #for mi, mixtureid in enumerate(mixtureids):
            #    print(mixtureid)
            #    if fixedvar == 'coverage':
            #        seriesorder = [(70, 0), (70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]
            #        seriesorder = seriesorder[1:]
            #        #elif fixedvar == 'ctdna':
            #        #    seriesorder = [(70, 0), (70, 80), (70, 180)]
            #        print(seriesorder)
            #        esttflist = []
            #        for s in seriesorder:
            #            sh, sl = s
            #            esttf = round((tfdict[mixtureid][0]*sh + tfdict[mixtureid][1]*sl)/(sl+sh), 2)
            #            print(esttf)
            #            esttflist.append(esttf)
            #        if mt == 'snv':
            #            esttflist =[esttflist] * len(config.methods)
            #            print(esttflist)
            #            resx[mi*len(config.methods):(mi+1)*len(config.methods), :] = esttflist
            #        else:
            #            esttflist =[esttflist] * (len(config.methods)-1)
            #            print(esttflist)
            #            resx[mi*(len(config.methods)-1):(mi+1)*(len(config.methods)-1), :] = esttflist
            print(resx)
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


            plt.figure(figsize=(15, 10))
            plt.grid(linewidth=1)
            plt.grid()
            for m in lm:
                plt.errorbar(res[m][2], res[m][0],
                             xerr=resx.std(axis=0)/np.sqrt(int(resx.shape[0]/(len(config.methods)))),
                             yerr=res[m][1]/int(resx.shape[0]/(len(config.methods))),
                             marker='s', c=color_dict[m],
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
                plt.xlim([60, 260])
                plt.grid()
            hand, labl = ax.get_legend_handles_labels()
            ax.legend(hand, labl, bbox_to_anchor=(1, 1), loc="upper left")
            plt.xlabel(xlab)
            plt.ylabel(metric.upper()+' score')
            plt.title(metric.upper() + " score for {} calling in chr{} with ref {}".format(mt.upper(), chrom, refname),
                      pad=50)
            if metric == 'auprc':
                plt.ylim([0, .5])
                #plt.ylim([0, .1])
                #plt.xlim([6, 2.])
            if mt == 'indel' and metric == 'auprc':
                if fixedvar == 'coverage':
                    plt.ylim([0, .3])
                else:
                    plt.ylim([0, .3])
            if metric != 'auprc':
                plt.ylim([0, 1.])
                #plt.ylim([0, .5])
                #plt.xlim([6, 2.])
            if not os.path.exists(os.path.join(*config.outputpath, 'figure2a')):
                os.mkdir(os.path.join(*config.outputpath, 'figure2a'))
            plt.savefig(os.path.join(*config.outputpath, 'figure2a',
                                     'perf_'+metric+'_3patients_150x_'+fixedvar+'_'+mt+'_'+refname+'.svg'), bbox_inches='tight')
            plt.show()

    """

    # ground truth analysis
    fixedvar = 'coverage'
    #for fixedvar in fixedvars:
    if fixedvar == 'coverage':
        seriesorder = [(70, 0), (70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]
        xaxis = 'tumor burden'
    elif fixedvar == 'ctdna':
        seriesorder = [(70, 0), (70, 80), (70, 180)]
        xaxis = 'coverage'
    for mixtureid in mixtureids:
        # mixtureid = 'CRC-986_100215-CW-T_CRC-986_300316-CW-T'
        print('############# {} ############'.format(mixtureid))
        if mixtureid == 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T':
            chroms = [str(c) for c in range(1,23) if c != 2 and c!=6 and c !=17 and c!=19 and c!=20 and c!=21]
            #chroms = [str(c) for c in range(1,9) if c != 2 and c!=6]
        elif mixtureid ==  'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
            chroms = [str(c) for c in range(1,23) if c !=1 and c!= 2 and c !=8 and c!=20 and c!=21 and c!=22]
        else:
            chroms = [str(c) for c in range(1,23) if c !=6 and c!=19 and c!=20]  # c !=1 and c!= 2 and
        calltables = {'sampleid':[], 'tf':[], 'cov':[], 'snv':[], 'indel':[], 'snp':[]}
        aux_all = []
        calltable_snv, aux = get_calltableseries(config, mixtureid, chroms, muttype='snv', filterparam=filterparam, reload=reload, save=save)
        calltable_indel, aux = get_calltableseries(config, mixtureid, chroms, muttype='indel', filterparam=filterparam, reload=reload, save=save)
        calltable_snp, aux = get_calltableseries(config, mixtureid, chroms, muttype='snp', filterparam=filterparam, reload=reload, save=save)
        print(calltable_snv.shape, calltable_indel.shape, calltable_snp.shape)
        print(aux)
        plasmasample = '_'.join(mixtureid.split('_')[:2])
        print(plasmasample)
        healthysample = '_'.join(mixtureid.split('_')[2:])
        print(healthysample)
        calltables['snv'].append(calltable_snv)
        calltables['indel'].append(calltable_indel)
        calltables['snp'].append(calltable_snp)
        calltables['sampleid'] = mixtureid
        calltables['tf'] = np.unique([cn.split('_')[0] for cn in list(calltable_snv.columns)])[:-5].astype(float)
        calltables['snv'] = pd.concat(calltables['snv'])
        calltables['indel'] = pd.concat(calltables['indel'])
        calltables['snp'] = pd.concat(calltables['snp'])
        dilutionseries = aux.T[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]].T
        #for muttype in muttypes:
        muttype = 'snv'
        refsample = 'undiluted'
        if muttype == 'snv':
            gtm = 5
        else:  # elif muttype == 'indel':
            gtm = 4
        print(max(aux['tf']))
        #if mixtureid != 'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
        #    calltablesseries = generate_groundtruth(config, calltables[muttype], aux['tf'], ground_truth_method=gtm, muttype=muttype,
        #                                            matchedtissuepath=None, methods=['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan', 'abemus', 'sinvict'])
        #else:
        calltablesseries = generate_groundtruth(config, calltables[muttype], aux['tf'], ground_truth_method=gtm, muttype=muttype,
                                                matchedtissuepath=None, methods=['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan', 'varnet', 'abemus', 'sinvict'])

        import matplotlib.transforms as transforms

        patient = mixtureid.split('-')[1].split('_')[0]
        print(patient)
        gtanalysis = calltablesseries[calltablesseries['truth'] == True][['{:.2f}_{}'.format(aux['tf'].max(), m) for m in config.methods]]
        initialanalysis = calltablesseries[calltablesseries['truth'] == False][['{:.2f}_{}'.format(aux['tf'].max(), m) for m in config.methods]]
        ngt = gtanalysis.shape[0]
        print(ngt)
        gtanalysis = pd.DataFrame(gtanalysis.sum()).T
        gtanalysis.columns = [l.split('_')[1] for l in list(gtanalysis.columns)]
        gtanalysis.index = ['calls in GT']
        initialanalysis = pd.DataFrame(initialanalysis.sum()).T
        initialanalysis.columns = [l.split('_')[1] for l in list(initialanalysis.columns)]
        initialanalysis.index = ['calls not in GT']

        fig, ax = plt.subplots(figsize=(10,5))
        plt.bar(gtanalysis.columns, gtanalysis.values.flatten(), bottom=1, color=[config.colors[config.methods.index(m)] for m in gtanalysis.columns], width=1)
        plt.axhline(y=ngt, c='blue', lw='3')
        trans = transforms.blended_transform_factory(
            ax.get_yticklabels()[0].get_transform(), ax.transData)
        ax.text(0, ngt, "{:.0f}".format(ngt), color="blue", transform=trans, ha="right", va="center")
        for col in initialanalysis.columns:
            print(col, initialanalysis[col].values[0], gtanalysis[col].values[0])
            plt.bar(col, initialanalysis[col].values[0], bottom=gtanalysis[col].values[0], label=col, color='grey', alpha=0.5, width=1)

        for pi, p in enumerate(ax.patches):
            print(pi, p)
            if pi < (len(ax.patches)/2):
                width, height = p.get_width(), p.get_height()
                x, y = p.get_xy()
                ax.text(x+width/2,
                        20,
                        '{:.0f}'.format(height),
                        horizontalalignment='center')
            else:
                width, height = p.get_width(), p.get_height()
                x, y = p.get_xy()
                ax.text(x+width/2,
                        y+height/10,
                        '{:.0f}'.format(height),
                        horizontalalignment='center')
        ax.set_yscale('log')
        plt.ylabel('# SNV')
        ax.grid(axis='y')
        plt.ylim([1, 5e5])
        ax = plt.gca()
        ax.set_xticklabels(labels=gtanalysis.columns,rotation=90)
        plt.savefig(os.path.join(*config.outputpath, 'figure2a', 'gtanalysis_barplot_'+patient+'_SNV.svg'), bbox_inches='tight')
        plt.show()
        """
