# Imports
import warnings
warnings.filterwarnings('ignore')

from utils.metrics import *
from utils.calltable import *
from utils.viz import function_to_split


def plot_metricsseries(config, restables, mixtureids, chrom, metric, muttype='snv', ground_truth_method='mixture',
                       fixedvar='coverage', refname='inundilutedsamplebyatleast5callers',
                       allpatients=True, logscale=False,  save=False, diltype='mixtures'):
            xlab = 'tumor burden' if fixedvar == 'coverage' else 'coverage'
            print(xlab)
            res = {'x': [], 'y': [], 'label': []}
            color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
            lc = 'logscale' if logscale else 'linearscale'
            if allpatients:
                plt.figure(figsize=(15, 10))
                plt.grid(linewidth=1)
            for mi, mixtureid in enumerate(mixtureids):
                if mixtureid == '':
                    break
                if not allpatients:  # one patient per plot
                    plt.figure(figsize=(15, 10))
                    plt.grid(linewidth=1)
                plasmasample = '_'.join(mixtureid.split('_')[:2])
                print(plasmasample)
                restablesample = restables[muttype][restables[muttype]['plasma sample'] == plasmasample]
                for method in config.methods:
                    if muttype == 'indel' and method == 'abemus':
                        pass
                    else:
                        plt.plot(restablesample[restablesample['caller'] == method][xlab], restablesample[restablesample['caller'] == method][metric.upper()+' score'],
                                 c=color_dict[method], marker=config.markers[mi], markersize=15, lw=2, label = method + '*' + 'patient '+plasmasample.split('_')[0].split('-')[1])
                        plt.grid()
                        res['x'].append(restablesample[restablesample['caller'] == method][xlab])
                        res['y'].append(restablesample[restablesample['caller'] == method][metric.upper()+' score'])
                        res['label'].append(method + '*' + 'patient '+plasmasample.split('_')[0].split('-')[1])

                if not allpatients:
                    if fixedvar == 'coverage':
                        plt.gca().invert_xaxis()
                    ax = plt.gca()
                    hand, labl = ax.get_legend_handles_labels()
                    hand, labl = function_to_split(hand, labl, '*')
                    ax.legend(hand, labl, bbox_to_anchor=(1, 1), loc="upper left")
                    plt.xlabel(xlab)
                    plt.ylabel(metric.upper()+' score')
                    plt.title(metric.upper() + " score for {} calling in chr{} with ref {}".format(muttype.upper(), chrom, refname))
                    if logscale:
                        plt.semilogy()
                    if ground_truth_method == 'spikein':
                        dilfolder = config.spikeinfolder
                    elif diltype == 'mixture':
                        dilfolder = config.mixturefolder
                    elif diltype == 'mixture_wes':
                        dilfolder = config.mixturefolderultradeep
                    elif diltype == 'mixture_wgs':
                        dilfolder = config.mixturefolderwholegenome
                    if save:
                        if not os.path.exists(os.path.join(*dilfolder, 'figures')):
                            os.mkdir(os.path.join(*dilfolder, 'figures'))
                        plt.savefig(os.path.join(*dilfolder, 'figures', mixtureid + '_' + metric + '_' + muttype + '_chr' + chrom + '_' + refname + '_' + lc + '_' + config.context), bbox_inches='tight')
            if allpatients:
                if fixedvar == 'coverage':
                    plt.gca().invert_xaxis()
                ax = plt.gca()
                hand, labl = ax.get_legend_handles_labels()
                hand, labl = function_to_split(hand, labl, '*')
                ax.legend(hand, labl, bbox_to_anchor=(1, 1), loc="upper left")
                plt.xlabel(xlab)
                plt.ylabel(metric.upper()+' score')
                plt.title(metric.upper() + " score for {} calling in chr{} with ref {}".format(muttype.upper(), chrom, refname), pad=50)
                if logscale:
                    plt.semilogy()
                dilfolder = config.spikeinfolder if ground_truth_method == 'spikein' else config.mixturefolder
                if save:
                    if not os.path.exists(os.path.join(*dilfolder, 'figures')):
                        os.mkdir(os.path.join(*dilfolder, 'figures'))
                    plt.savefig(os.path.join(*dilfolder, 'figures', metric + '_' + muttype + '_chr' + chrom + '_' + refname + '_'  + fixedvar + '_' + lc + '_' + config.context), bbox_inches='tight')
            # plt.show()
            return res


if __name__ == '__main__':
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

    # Chomosome
    chrom = '22'
    mixtureids = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T',
                  'CRC-986_100215-CW-T_CRC-986_300316-CW-T',
                  'CRC-123_310715-CW-T_CRC-123_121115-CW-T']
    # refnamesnv = 'inundilutedsamplebyatleast5callers'
    # refnameindel = 'inundilutedsamplebyatleast3callers'
    refnamesnv = 'inundilutedsampleranked'
    refnameindel = 'inundilutedsampleranked'
    fixedvar='coverage'

    markers = ['o', '^', 'X']
    linestyles = ['-', '-', '-']
    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}

    muttypes = ['snv', 'indel']
    metrics = ['auprc', 'precision', 'recall']

    for muttype in muttypes:
        if muttype == 'snv':
            refname = refnamesnv
        else:  # muttype == 'indel':
            refname = refnameindel
        print(refname)
        for metric in metrics:
            # Save table if do not exist and load tables
            restables = {'snv': [], 'indel': []}
            for mt in ['snv', 'indel']:
                for mixtureid in mixtureids:
                    plasmasample = '_'.join(mixtureid.split('_')[:2])
                    restable = pd.read_csv(os.path.join(
                        *config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+mixtureid, 'results',
                        mixtureid+'_'+mt+'_'+metric+'_'+refname+'_fixed'+fixedvar+'.csv'), index_col=0)
                    restables[mt].append(restable)
                restables[mt] = pd.concat(restables[mt])
            plot_metricsseries(config, restables, mixtureids, chrom, metric=metric, muttype=muttype,
                               ground_truth_method='mixture', fixedvar='coverage', allpatients=True, logscale=False, save=True)
            plot_metricsseries(config, restables, mixtureids, chrom, metric=metric, muttype=muttype,
                               ground_truth_method='mixture', fixedvar='coverage', allpatients=True, logscale=True, save=True)
            plot_metricsseries(config, restables, mixtureids, chrom, metric=metric, muttype=muttype,
                               ground_truth_method='mixture', fixedvar='coverage', allpatients=False, logscale=False, save=True)
            plot_metricsseries(config, restables, mixtureids, chrom, metric=metric, muttype=muttype,
                               ground_truth_method='mixture', fixedvar='coverage', allpatients=False, logscale=True, save=True)

