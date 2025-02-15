import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import warnings
import os
from sklearn.metrics import precision_recall_curve, roc_curve, f1_score, precision_score, recall_score, average_precision_score, roc_auc_score

warnings.filterwarnings("ignore")


def plot_pr_curve(precision, recall, estimator_name=None, f1_score=None, figax=None, kwargs={}, plot='all'):
    if f1_score is not None and estimator_name is not None:
        kwargs["label"] = f"{estimator_name} (F1 = {f1_score:0.2f})"
    elif f1_score is not None:
        kwargs["label"] = f"AP = {f1_score:0.2f}"
    elif estimator_name is not None:
        kwargs["label"] = estimator_name
    fig, ax = figax if figax is not None else plt.subplots()
    print('before', recall, precision)
    if plot == 'all':
        kwargs["drawstyle"] = "steps-post"
        ax.plot(recall, precision, **kwargs)
    elif plot == 'partial':
        listtoremove = []  # index to remove
        invprecision = precision[::-1]
        for i in range(len(invprecision)-1):
            if (invprecision[i] == 1) and (invprecision[i+1] == 1):
                listtoremove.append(len(invprecision)-1-i)
        recall = [recall[i] for i in range(len(recall)) if i not in listtoremove]
        precision = [precision[i] for i in range(len(precision)) if i not in listtoremove]
        print(listtoremove)
        print('after', recall, precision)
        if (len(recall) == 2) and (recall[0] == 1) and (recall[1] != 0):
            ax.scatter(recall[1], precision[1],  **kwargs)
        else:
            kwargs["drawstyle"] = "steps-post"
            ax.plot(recall, precision, **kwargs)
    elif plot == 'average':
        kwargs["drawstyle"] = "steps-post"
        sns.lineplot(x=recall, y=precision, ax=ax, **kwargs)
    else:
        raise ValueError('plot argument should be either all, partial or average but here is {}'.format(plot))
    xlabel = "Recall"
    ylabel = "Precision"
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.legend()  # loc="lower left")
    f_scores = np.linspace(0.1, 0.9, num=9)
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.1, lw=1)
        plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.6, y[45] + 0.02), color='grey')
    return recall, precision


def figure_curve_allchr(config, df_table, dilutionseries, mixtureid, xy='pr', ground_truth_method=3, refsample='undiluted',
                        muttype='SNV', methods=None, fixedvar='coverage', save=True, diltype='mixture', splitby='method', plot='all', figax=None):
    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
    print(list(dilutionseries.index))
    alpha_dict = {d: 1-0.1*i for i, d in enumerate(list(dilutionseries.index))}
    baseline_dict = {}
    dilutionseries_present = []
    if methods is None:
        methods = config.methods
    if splitby == 'method':
        prall = {}
        for d in list(dilutionseries.index):
            prall[d] = {}
        for method in methods:
            if figax is None:
                fig, ax = plt.subplots(figsize=(10, 10))
            else:
                fig, ax = figax[0], figax[1]
            for i, d in enumerate(list(dilutionseries.index)):
                if ground_truth_method != 'spikein' and ground_truth_method != 'SEQC2':
                    factorprefix = '{:.2f}'.format(round(dilutionseries.loc[d, 'tf'], 2))
                elif ground_truth_method == 'SEQC2':
                    factorprefix = '{}'.format(dilutionseries.loc[d, 'samplename'])
                else:
                    factorprefix = '{:.2f}'.format(round(dilutionseries.loc[d, 'vaf'], 2))
                print(factorprefix + '_' + method + '_score')
                if factorprefix + '_' + method + '_score' in list(df_table.columns):
                    if type(ground_truth_method) == int or ground_truth_method == 'spikein' or ground_truth_method == 'ranked' or ground_truth_method == 'SEQC2':
                        truth_name = 'truth'
                    elif ground_truth_method == 'method':
                        truth_name = method +'_truth'
                    else:
                        raise ValueError('unknown ground truth {}'.format(ground_truth_method))
                    df_method = df_table[[truth_name, factorprefix + '_' + method + '_score']]
                    print(df_method.index[df_method.index.duplicated()])
                    df_method[truth_name].fillna(False, inplace=True)
                    df_method[factorprefix + '_' + method + '_score'].fillna(0, inplace=True)
                    baseline_dict[str(d)] = len(df_method[truth_name][df_method[truth_name]])/len(df_method[truth_name])
                    if str(d) not in dilutionseries_present:
                        dilutionseries_present.append(str(d))
                    if xy == 'pr':
                        precision, recall, thresholds = precision_recall_curve(df_method[truth_name], df_method[factorprefix + '_' + method + '_score'])
                        if i == 0:
                            p, r = plot_pr_curve(precision, recall, estimator_name=method, f1_score=None, figax=(fig, ax),
                                          kwargs={'color': color_dict[method], 'alpha': alpha_dict[str(d)], 'lw': 4-i/3}, plot=plot)
                        else:
                            p, r = plot_pr_curve(precision, recall, estimator_name='', f1_score=None, figax=(fig, ax),
                                          kwargs={'color': color_dict[method], 'alpha': alpha_dict[str(d)], 'lw': 4-i/3}, plot=plot)
                        prall[d][method] = [p, r]
                    elif xy == 'roc':
                        fpr, tpr, thresholds = roc_curve(df_method[truth_name], df_method[factorprefix + '_' + method + '_score'])
                        if i == 0:
                            plot_roc_curve(fpr, tpr, estimator_name=method, auc_score=None, figax=(fig, ax),
                                           kwargs={'color': color_dict[method], 'alpha': alpha_dict[str(d)], 'lw': 3.5-int(i/2)})
                        else:
                            plot_roc_curve(fpr, tpr, estimator_name='', auc_score=None, figax=(fig, ax),
                                           kwargs={'color': color_dict[method], 'alpha': alpha_dict[str(d)], 'lw': 3.5-int(i/2)})
            # print(baseline_dict)
            # print(tb_dict)
            print(dilutionseries_present)
            list_lines_baseline = []
            if xy == 'pr':
                if len(np.unique(baseline_dict.values())) == 1:
                    plt.axhline(y=baseline_dict[str(dilutionseries_present[0])], ls='--', c='k')
                    list_lines_baseline.append(Line2D([0], [0], color='black', ls='--', label="baseline = {:.2f}".format(baseline_dict[str(dilutionseries_present[0])])))
                else:
                    for d in dilutionseries_present:
                        plt.axhline(y=baseline_dict[str(d)], alpha=alpha_dict[str(d)], ls='--', c='k')
                        if ground_truth_method != 'spikein' and ground_truth_method != 'SEQC2':
                            list_lines_baseline.append(Line2D([0], [0], color='black', ls='--', alpha=alpha_dict[str(d)], label="baseline tf {:.2f}% = {:.2f}".format(dilutionseries.loc[d, 'tf'], baseline_dict[str(d)])))
                        elif ground_truth_method == 'SEQC2':
                            list_lines_baseline.append(Line2D([0], [0], color='black', ls='--', alpha=alpha_dict[str(d)], label="baseline sample {}% = {:.2f}".format(dilutionseries.loc[d, 'samplename'], baseline_dict[str(d)])))
                        else:
                            list_lines_baseline.append(Line2D([0], [0], color='black', ls='--', alpha=alpha_dict[str(d)], label="baseline vaf {:.2f}% = {:.2f}".format(dilutionseries.loc[d, 'cov'], baseline_dict[str(d)])))
            handles, labels = plt.gca().get_legend_handles_labels()
            if fixedvar == 'coverage':
                if ground_truth_method != 'spikein' and ground_truth_method != 'SEQC2':
                    list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[str(i)], label='tumor burden = {:.2f}%'.format(dilutionseries.loc[i, 'tf'])) for i in dilutionseries_present]
                elif ground_truth_method == 'SEQC2':
                    list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[str(i)], label='sample = {}'.format(dilutionseries.loc[i, 'samplename'])) for i in dilutionseries_present]
                else:
                    list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[str(i)], label='VAF = {:.2f}%'.format(float(i))) for i in dilutionseries_present]
            elif fixedvar == 'ctdna':
                list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[str(i)], label='cov = {}x'.format(int(round(dilutionseries.loc[i, 'cov']/10)*10))) for i in dilutionseries_present]
            if xy == 'pr':
                legend_list = handles + list_lines + list_lines_baseline
            else:
                legend_list = handles + list_lines
            # Creating legend with color box
            plt.legend(bbox_to_anchor=(1, 1), loc="upper left", handles=legend_list)
            if xy == 'pr':
                plt.title("Precision Recall curve for SNV calling in sample {}".format(mixtureid), pad=50)
            elif xy == 'roc':
                plt.title("Receiver Operation Characteristics curve for SNV calling in sample {}".format(mixtureid), pad=50)
            if xy == 'pr':
                plt.semilogx()
                plt.xlim([0.01, 1.01])
            else:
                plt.xlim([-0.01, 1.01])
            plt.ylim([-0.01, 1.01])
            if ground_truth_method == 'spikein':
                dilfolder = config.spikeinfolder
                dilution = 'spikeins'
            elif diltype == 'mixture':
                dilfolder = config.mixturefolder
                dilution = 'mixtures'
            elif diltype == 'mixture_wes':
                dilfolder = config.mixturefolderultradeep
                dilution = 'mixtures'
            elif diltype == 'mixture_wgs':
                dilfolder = config.mixturefolderwholegenome
                dilution = 'mixtures'
            elif diltype == 'SEQC2':
                dilfolder = config.mixturefolderSEQC2
                dilution = 'SEQC2'
            if save != False:
                if type(save) == list:
                    savepath = os.path.join(*save)
                    ext = '.svg'
                else:
                    savepath = os.path.join(*dilfolder, dilution+'_allchr', 'figures')
                    ext = '.png'
                if not os.path.exists(savepath):
                    os.mkdir(savepath)
                if type(ground_truth_method) == int:
                    refname = 'in'+refsample + 'samplebyatleast' + str(ground_truth_method) +'callers'
                elif ground_truth_method == 'ranked':
                    refname = 'in'+refsample + 'sampleranked'
                elif ground_truth_method == 'SEQC2':
                    refname = 'inKnownVariants'
                else:
                    refname = 'in'+refsample + 'samplebythesamecaller'
                print(os.path.join(savepath, mixtureid + '_' + muttype + '_' + xy.upper() + 'curve_' + refname + '_' + method + '_' + config.context + ext))
                plt.savefig(os.path.join(savepath, mixtureid + '_' + muttype + '_' + xy.upper() + 'curve_' + refname + '_' + method + '_' + config.context + ext), bbox_inches='tight')
                plt.savefig(os.path.join(savepath, mixtureid + '_' + muttype + '_' + xy.upper() + 'curve_' + refname + '_' + method + '_' + config.context + ext + '.svg'), bbox_inches='tight')
    elif splitby == 'dilution':
        prall = {}
        for d in list(dilutionseries.index):
            prall[d] = {}
        for i, d in enumerate(list(dilutionseries.index)):
            if figax is None:
                fig, ax = plt.subplots(figsize=(10, 10))
            else:
                fig, ax = figax[0], figax[1]
            for method in methods:
                if ground_truth_method != 'spikein' and ground_truth_method != 'SEQC2':
                    factorprefix = '{:.2f}'.format(round(dilutionseries.loc[d, 'tf'], 2))
                elif ground_truth_method == 'SEQC2':
                    factorprefix = '{}'.format(dilutionseries.loc[d, 'samplename'])
                else:
                    factorprefix = '{:.2f}'.format(round(dilutionseries.loc[d, 'vaf'], 2))
                print(factorprefix + '_' + method + '_score')
                if factorprefix + '_' + method + '_score' in list(df_table.columns):
                    if type(ground_truth_method) == int or ground_truth_method == 'spikein' or ground_truth_method == 'ranked' or ground_truth_method == 'SEQC2':
                        truth_name = 'truth'
                    elif ground_truth_method == 'method':
                        truth_name = method +'_truth'
                    else:
                        raise ValueError('unknown ground truth {}'.format(ground_truth_method))
                    df_method = df_table[[truth_name, factorprefix + '_' + method + '_score']]
                    print(df_method.index[df_method.index.duplicated()])
                    df_method[truth_name].fillna(False, inplace=True)
                    df_method[factorprefix + '_' + method + '_score'].fillna(0, inplace=True)
                    baseline_dict[str(d)] = len(df_method[truth_name][df_method[truth_name]])/len(df_method[truth_name])
                    if str(d) not in dilutionseries_present:
                        dilutionseries_present.append(str(d))
                    if xy == 'pr':
                        precision, recall, thresholds = precision_recall_curve(df_method[truth_name], df_method[factorprefix + '_' + method + '_score'])
                        if i == 0:
                            p, r = plot_pr_curve(precision, recall, estimator_name=method, f1_score=None, figax=(fig, ax),
                                          kwargs={'color': color_dict[method], 'alpha': alpha_dict[str(d)], 'lw': 4-i/3}, plot=plot)
                        else:
                            p, r = plot_pr_curve(precision, recall, estimator_name='', f1_score=None, figax=(fig, ax),
                                          kwargs={'color': color_dict[method], 'alpha': alpha_dict[str(d)], 'lw': 4-i/3}, plot=plot)
                        prall[d][method] = [p, r]
                    elif xy == 'roc':
                        fpr, tpr, thresholds = roc_curve(df_method[truth_name], df_method[factorprefix + '_' + method + '_score'])
                        if i == 0:
                            plot_roc_curve(fpr, tpr, estimator_name=method, auc_score=None, figax=(fig, ax),
                                           kwargs={'color': color_dict[method], 'alpha': alpha_dict[str(d)], 'lw': 3.5-int(i/2)})
                        else:
                            plot_roc_curve(fpr, tpr, estimator_name='', auc_score=None, figax=(fig, ax),
                                           kwargs={'color': color_dict[method], 'alpha': alpha_dict[str(d)], 'lw': 3.5-int(i/2)})
            # print(baseline_dict)
            # print(tb_dict)
            print(dilutionseries_present)
            list_lines_baseline = []
            if xy == 'pr':
                if len(np.unique(baseline_dict.values())) == 1:
                    plt.axhline(y=baseline_dict[str(dilutionseries_present[0])], ls='--', c='k')
                    list_lines_baseline.append(Line2D([0], [0], color='black', ls='--', label="baseline = {:.2f}".format(baseline_dict[str(dilutionseries_present[0])])))
                else:
                    for d in dilutionseries_present:
                        plt.axhline(y=baseline_dict[str(d)], alpha=alpha_dict[str(d)], ls='--', c='k')
                        if ground_truth_method != 'spikein' and ground_truth_method != 'SEQC2':
                            list_lines_baseline.append(Line2D([0], [0], color='black', ls='--', alpha=alpha_dict[str(d)], label="baseline tf {:.2f}% = {:.2f}".format(dilutionseries.loc[d, 'tf'], baseline_dict[str(d)])))
                        elif ground_truth_method == 'SEQC2':
                            list_lines_baseline.append(Line2D([0], [0], color='black', ls='--', alpha=alpha_dict[str(d)], label="baseline sample {}% = {:.2f}".format(dilutionseries.loc[d, 'samplename'], baseline_dict[str(d)])))
                        else:
                            list_lines_baseline.append(Line2D([0], [0], color='black', ls='--', alpha=alpha_dict[str(d)], label="baseline vaf {:.2f}% = {:.2f}".format(dilutionseries.loc[d, 'cov'], baseline_dict[str(d)])))
            handles, labels = plt.gca().get_legend_handles_labels()
            if fixedvar == 'coverage':
                if ground_truth_method != 'spikein' and ground_truth_method != 'SEQC2':
                    list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[str(i)], label='tumor burden = {:.2f}%'.format(dilutionseries.loc[i, 'tf'])) for i in dilutionseries_present]
                elif ground_truth_method == 'SEQC2':
                    list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[str(i)], label='sample = {}'.format(dilutionseries.loc[i, 'samplename'])) for i in dilutionseries_present]
                else:
                    list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[str(i)], label='VAF = {:.2f}%'.format(float(i))) for i in dilutionseries_present]
            elif fixedvar == 'ctdna':
                list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[str(i)], label='cov = {}x'.format(int(round(dilutionseries.loc[i, 'cov']/10)*10))) for i in dilutionseries_present]
            if xy == 'pr':
                legend_list = handles + list_lines + list_lines_baseline
            else:
                legend_list = handles + list_lines
            # Creating legend with color box
            plt.legend(bbox_to_anchor=(1, 1), loc="upper left", handles=legend_list)
            if xy == 'pr':
                plt.title("Precision Recall curve for SNV calling in sample {}".format(mixtureid), pad=50)
            elif xy == 'roc':
                plt.title("Receiver Operation Characteristics curve for SNV calling in sample {}".format(mixtureid), pad=50)
            if xy == 'pr':
                plt.semilogx()
                plt.xlim([0.01, 1.01])
            else:
                plt.xlim([-0.01, 1.01])
            plt.ylim([-0.01, 1.01])
            if ground_truth_method == 'spikein':
                dilfolder = config.spikeinfolder
                dilution = 'spikeins'
            elif diltype == 'mixture':
                dilfolder = config.mixturefolder
                dilution = 'mixtures'
            elif diltype == 'mixture_wes':
                dilfolder = config.mixturefolderultradeep
                dilution = 'mixtures'
            elif diltype == 'mixture_wgs':
                dilfolder = config.mixturefolderwholegenome
                dilution = 'mixtures'
            elif diltype == 'SEQC2':
                dilfolder = config.mixturefolderSEQC2
                dilution = 'SEQC2'
            if save != False:
                if type(save) == list:
                    savepath = os.path.join(*save)
                    ext = '.svg'
                else:
                    savepath = os.path.join(*dilfolder, dilution+'_allchr', 'figures')
                    ext = '.png'
                if not os.path.exists(savepath):
                    os.mkdir(savepath)
                if type(ground_truth_method) == int:
                    refname = 'in'+refsample + 'samplebyatleast' + str(ground_truth_method) +'callers'
                elif ground_truth_method == 'ranked':
                    refname = 'in'+refsample + 'sampleranked'
                elif ground_truth_method == 'SEQC2':
                    refname = 'inKnownVariants'
                else:
                    refname = 'in'+refsample + 'samplebythesamecaller'
                print(os.path.join(savepath, mixtureid + '_' + muttype + '_' + xy.upper() + 'curve_' + refname + '_' + method + '_' + config.context + '_splitby' + splitby + ext))
                plt.savefig(os.path.join(savepath, mixtureid + '_' + muttype + '_' + xy.upper() + 'curve_' + refname + '_' + method + '_' + config.context + '_splitby' + splitby + ext), bbox_inches='tight')
                plt.savefig(os.path.join(savepath, mixtureid + '_' + muttype + '_' + xy.upper() + 'curve_' + refname + '_' + method + '_' + config.context + '_splitby' + splitby + ext + '.svg'), bbox_inches='tight')
    if xy == 'pr':
        return prall




def plot_roc_curve(fpr, tpr, estimator_name=None, auc_score=None, figax=None, kwargs={}):
    kwargs["drawstyle"] = "steps-post"
    if auc_score is not None and estimator_name is not None:
        kwargs["label"] = f"{estimator_name} (AUC = {auc_score:0.2f})"
    elif auc_score is not None:
        kwargs["label"] = f"AUC = {auc_score:0.2f}"
    elif estimator_name is not None:
        kwargs["label"] = estimator_name
    fig, ax = figax if figax is not None else plt.subplots()
    ax.plot(fpr, tpr, **kwargs)
    xlabel = "False Positive Rate"
    ylabel = "True Positive Rate"
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.legend() #(loc="upper right")
    plt.plot([0, 1], [0, 1], ls='--', color='grey')
    plt.annotate('baseline', xy=(0.9, 0.9), color='grey')


def metric_curve(config, df_table, plasmasample, healthysample, dilutionseries, metric='auprc', ground_truth_method=3,
                 refsample='undiluted', muttype='SNV', chrom='22', methods=None, fixedvar='coverage', xaxis='tumor burden', save=True, diltype='mixture'):
    tb_dict = {}
    cov_dict = {}
    dilutionseries_present = []
    if methods is None:
        methods = config.methods
    color_dict = {}
    for i, m in enumerate(config.methods):
        if m in methods:
            color_dict[m] = config.colors[i]
    if ground_truth_method != 'spikein':
        for i, d in enumerate(dilutionseries):
            if chrom in [str(c) for c in range(1, 23)]:
                mixturepath = 'mixture_chr'+chrom+'_'+plasmasample +"_" + str(d[0]) +"x_" + healthysample + "_" + str(d[1]) + 'x'
                tb_dict[str(d)] = float(pd.read_csv(os.path.join(
                    *config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+plasmasample+'_'+healthysample,
                    mixturepath, 'estimated_tf_chr'+chrom+mixturepath[len(('mixture_chr'+chrom)):]+'.txt')).columns[0])
                cov_dict[str(d)] = float(pd.read_csv(os.path.join(
                    *config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+plasmasample+'_'+healthysample,
                    mixturepath, 'coverage_chr'+chrom+mixturepath[len(('mixture_chr'+chrom)):]+'.txt')).columns[0])
            #else:  # chrom == 'all'
            #    calltables_aux = pd.read_csv(os.path.join(mixturefolder, 'calls', mixtureid+'_tf_cov.csv'), index_col=0)
    results_df = pd.DataFrame()
    aux_metric = []
    aux_metricrelative = []
    aux_method = []
    aux_tb = []
    aux_cov = []
    baseline = {}
    if xaxis == 'tumor burden' or ground_truth_method == 'spikein':
        for i, d in enumerate(dilutionseries):
            print(d)
            if ground_truth_method != 'spikein':
                factorprefix = '{:.2f}'.format(100*round(tb_dict[str(d)], 4))
            else:
                factorprefix = '{:.2f}'.format(d)
            for method in methods:
                if factorprefix + '_' + method + '_score' in list(df_table.columns):
                    print(factorprefix + '_' + method + '_score')
                    if i != 0 or ground_truth_method == 'spikein' or xaxis == 'coverage':
                        if type(ground_truth_method) == int or ground_truth_method == 'spikein' or ground_truth_method == 'ranked':
                            truth_name = 'truth'
                        elif ground_truth_method == 'caller':
                            truth_name = method + '_truth'
                        else:
                            raise ValueError('unknown ground truth {}'.format(ground_truth_method))
                        if metric == 'auprc':
                            df_method = df_table[[truth_name, factorprefix + '_' + method + '_score', factorprefix + '_' + method]]
                            # df_method.dropna(how='all', inplace=True)
                            df_method[truth_name].fillna(False, inplace=True)
                            df_method[factorprefix + '_' + method + '_score'].fillna(0, inplace=True)
                            df_method.drop(factorprefix + '_' + method, axis=1, inplace=True)
                        else:
                            df_method = df_table[[truth_name, factorprefix + '_' + method]]
                            # df_method.dropna(how='all', inplace=True)
                            df_method[truth_name].fillna(False, inplace=True)
                            df_method[factorprefix + '_' + method].fillna(False, inplace=True)
                        df_method[truth_name] = df_method[truth_name].fillna(False)
                        baseline[method] = len(df_method[truth_name][df_method[truth_name]])/len(df_method[truth_name])
                        if str(d) not in dilutionseries_present:
                            dilutionseries_present.append(str(d))
                        if metric == 'auprc':
                            aux_metric.append(average_precision_score(df_method[truth_name], df_method[factorprefix + '_' + method + '_score']))
                            aux_metricrelative.append(average_precision_score(
                                df_method[truth_name], df_method[factorprefix + '_' + method + '_score']) - baseline[method])
                        elif metric == 'precision':
                            aux_metric.append(precision_score(df_method[truth_name], df_method[factorprefix + '_' + method]))
                        elif metric == 'precisionmaxf1':
                            aux_metric.append(precision_score(df_method[truth_name], df_method[factorprefix + '_' + method]))
                            # print(df_method[factorprefix + '_' + method][df_method[truth_name] == True])
                        elif metric == 'recall':
                            aux_metric.append(recall_score(df_method[truth_name], df_method[factorprefix + '_' + method]))
                        elif metric == 'f1':
                            aux_metric.append(f1_score(df_method[truth_name], df_method[factorprefix + '_' + method]))
                        elif metric == 'auc':
                            aux_metric.append(roc_auc_score(df_method[truth_name], df_method[factorprefix + '_' + method + '_score']))
                        aux_method.append(method)
                        if ground_truth_method != 'spikein':
                            aux_tb.append(round(100*tb_dict[str(d)], 3))
                            aux_cov.append(int(cov_dict[str(d)]))
                        else:
                            aux_tb.append(d)
    elif xaxis == 'vaf':
        np.sort(pd.qcut(calltable_snv[calltable_snv['truth']==True]['median_vaf'], q=10).unique())
        #vafranges = [1., 0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005, 0.0025, 0.001]
        #vafranges = [1., 0.016, 0.012, 0.0093, 0.0073, 0.0063, 0.0056, 0.0051, 0.0045, 0.0038, 0]
        vafranges = [.195, 0.0894, 0.0758, 0.0597, 0.0495, 0.0429, 0.0342, 0.0281, 0.024, 0.0172, 0.0052]
        vafranges = vafranges[::-1]
        res = {}
        calltable_snv = calltables['snv'][calltables['snv']['sampletf'] != calltables['snv']['sampletf'].unique().max()]
        for method in config.methods:
            print(method)
            x, y  = [], []
            for vi, vafrange in enumerate(vafranges):
                #print(vafranges[vi-1], vafranges[vi])
                if vi > 0:
                    aux = calltable_snv[(calltable_snv['median_vaf'] >= vafranges[vi-1]) & (calltable_snv['median_vaf'] < vafranges[vi])]
                    #if not aux.empty:
                    print(vafranges[vi-1], vafranges[vi], aux.shape[0])
                    precision, recall, thresholds = precision_recall_curve(aux['truth'], aux[method + '_score'].fillna(0))
                    f1list = 2*(precision * recall)/(precision + recall)
                    #print(len(f1list), len(precision))
                    #max(f1list)
                    #print(precision, recall, thresholds)
                    #print(average_precision_score(aux['truth'], aux[method + '_score'].fillna(0)))
                    #if np.nanmax(f1list) > 0:
                    #    x.append((vafranges[vi-1]+ vafranges[vi])/2)
                    #    y.append(recall[list(f1list).index(np.nanmax(f1list))])
                    #x.append((vafranges[vi-1]+ vafranges[vi])/2)
                    #y.append(average_precision_score(aux['truth'], aux[method + '_score'].fillna(0)))
                    x.append((vafranges[vi-1]+ vafranges[vi])/2)
                    y.append(np.nanmax(f1list))
            res[method] = {'x': x, 'y': y}
        color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
        plt.figure(figsize=(15,8))
        for k,v in res.items():
            plt.plot(v['x'], v['y'], marker='o', label=k, color=color_dict[k])
        plt.gca().invert_xaxis()
        #plt.xscale("log")
        plt.legend()
        plt.show()
        df_table['median VAF'] = df_table[[c for c in list(df_table.columns) if c.endswith('vaf')]].median()
        df_table['VAF'] = pd.cut(df_table['median VAF'],
                                 bins=[0, .01, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, 1],
                                 labels=['≤1%', '≤5%', '≤10%', '≤15%', '≤20%', '≤25%', '≤30%', '≤35%', '≤40%', '≤45%', '≤50%', '>50%'],
                                 include_lowest=True)
    # print(baseline)
    results_df[metric.upper() + ' score'] = aux_metric
    if metric == 'auprc':
        results_df[metric.upper() + ' score - baseline ' + metric.upper() + ' score'] = aux_metricrelative
    # print(aux_tb)
    if ground_truth_method != 'spikein':
        if fixedvar == 'coverage':
            results_df[xaxis] = aux_tb
            xvar = xaxis
            reverse = True
        else:  # fixedvar = ctdna
            results_df['coverage'] = aux_cov
            xvar = 'coverage'
            reverse = False
    else:
        results_df['VAF'] = aux_tb
        xvar = 'VAF'
        reverse = True
    results_df['caller'] = aux_method
    if config.context == 'paper':
        sns.set_style("whitegrid")
    else:
        sns.set_style("darkgrid")

    sns.catplot(x=xvar, y=metric.upper() + " score", hue="caller",
                capsize=.1, height=6, aspect=1.5, kind="point",
                order=sorted(results_df[xvar].unique(), reverse=reverse),
                palette=sns.color_palette(list(color_dict.values())), data=results_df)
    plt.title(metric.upper() + " score for {} calling in {} - chr{}".format(muttype, plasmasample, chrom))

    plt.ylim([0, max(0.5, max(aux_metric)+0.05)])
    if metric == 'auprc':
        if len(np.unique(baseline.values())) == 1:
            print(baseline[config.methods[0]])
            plt.axhline(y=baseline[config.methods[0]], color='k', linestyle='--')
    if type(ground_truth_method) == int:
        refname = 'in'+refsample + 'samplebyatleast' + str(ground_truth_method) +'callers'
    elif ground_truth_method == 'ranked':
        refname = 'in'+refsample + 'sampleranked'
    else:
        refname = 'in'+refsample + 'samplebythesamecaller'
    dilution = 'spikeins' if ground_truth_method == 'spikein' else 'mixtures'
    if ground_truth_method == 'spikein':
        dilfolder = config.spikeinfolder
        dilution = 'spikeins'
    elif diltype == 'mixture':
        dilfolder = config.mixturefolder
        dilution = 'mixtures'
    elif diltype == 'mixture_wes':
        dilfolder = config.mixturefolderultradeep
        dilution = 'mixtures'
    elif diltype == 'mixture_wgs':
        dilfolder = config.mixturefolderwholegenome
        dilution = 'mixtures'
    if save:
        if not os.path.exists(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'figures')):
            os.mkdir(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'figures'))
        if methods != config.methods:
            plt.savefig(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'figures',  plasmasample + '_' + healthysample + '_' + muttype + '_' + metric + '_' + refname + '_' + '_'.join(methods) + '_' + xaxis + '_' + config.context), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'figures', plasmasample + '_' + healthysample + '_' + muttype + '_' + metric + '_' + refname + '_' + xaxis + '_' + config.context), bbox_inches='tight')
    #plt.show()
    summary_df = results_df.copy()
    summary_df['plasma sample'] = plasmasample
    summary_df['healthy sample'] = healthysample
    summary_df['mutation type'] = muttype
    summary_df['metric'] = metric
    summary_df['mutation type'] = refname
    if save:
        if not os.path.exists(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'results')):
            os.mkdir(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'results'))
        if methods != config.methods:
            summary_df.to_csv(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'results',  plasmasample + '_' + healthysample + '_' + muttype + '_' + metric + '_' + refname + '_'.join(methods) + '_fixed' + fixedvar + '_' + xaxis + '.csv'))
        else:
            summary_df.to_csv(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'results', plasmasample + '_' + healthysample + '_' + muttype + '_' + metric + '_' + refname + '_fixed'+ fixedvar + '_' + xaxis + '.csv'))
    return summary_df


def metric_curve_allchr(config, df_table, dilutionseries, mixtureid, metric='auprc', ground_truth_method=4,
                 refsample='undiluted', muttype='SNV', methods=None, fixedvar='coverage', xaxis='tumor burden', save=True, diltype='mixture'):
    dilutionseries_present = []
    if methods is None:
        methods = config.methods
    color_dict = {}
    for i, m in enumerate(config.methods):
        if m in methods:
            color_dict[m] = config.colors[i]
    results_df = pd.DataFrame()
    aux_metric = []
    aux_metricrelative = []
    aux_method = []
    aux_tb = []
    aux_cov = []
    aux_samplename = []
    baseline = {}
    if xaxis == 'tumor burden' or xaxis == 'coverage' or ground_truth_method == 'spikein' or xaxis == 'vaf':
        for i, d in enumerate(list(dilutionseries.index)):
            print(d)
            if xaxis == 'vaf':
                factorprefix = ''
            elif ground_truth_method != 'spikein' and ground_truth_method != 'SEQC2':
                factorprefix = '{:.2f}'.format(round(dilutionseries.loc[d, 'tf'], 2))
            elif ground_truth_method == 'SEQC2':
                factorprefix = '{}'.format(dilutionseries.loc[d, 'samplename'])
            else:
                factorprefix = '{:.2f}'.format(dilutionseries.loc[d, 'vaf']) \
                    if d != 'spikein_CRC-COSMIC-5p_vaf0.025_CRC-986_300316-CW-T' \
                    else '0.03' #'{:.3f}'.format(dilutionseries.loc[d, 'vaf']) #TODO fix
                print(factorprefix)
            print(factorprefix)
            for method in methods:
                colname = factorprefix + '_' + method + '_score' if xaxis != 'vaf' else method + '_score'
                colnamebis = factorprefix + '_' + method if xaxis != 'vaf' else method
                print(colname)
                if colname in list(df_table.columns):
                    print('is present')
                    if i != 0 or ground_truth_method == 'spikein' or xaxis == 'coverage' or ground_truth_method == 'SEQC2':
                        if type(ground_truth_method) == int or ground_truth_method == 'spikein' or ground_truth_method == 'ranked' or ground_truth_method =='SEQC2':
                            truth_name = 'truth'
                        elif ground_truth_method == 'caller':
                            truth_name = method + '_truth'
                        else:
                            raise ValueError('unknown ground truth {}'.format(ground_truth_method))
                        if metric == 'auprc' or metric.startswith('max'):
                            df_method = df_table[[truth_name, colname, colnamebis]]
                            # df_method.dropna(how='all', inplace=True)
                            df_method[truth_name].fillna(False, inplace=True)
                            df_method[colname].fillna(0, inplace=True)
                            df_method.drop(colnamebis, axis=1, inplace=True)
                        else:
                            df_method = df_table[[truth_name, colnamebis]]
                            # df_method.dropna(how='all', inplace=True)
                            df_method[truth_name].fillna(False, inplace=True)
                            df_method[colnamebis].fillna(False, inplace=True)
                        df_method[truth_name] = df_method[truth_name].fillna(False)
                        baseline[method] = len(df_method[truth_name][df_method[truth_name]])/len(df_method[truth_name])
                        if str(d) not in dilutionseries_present:
                            dilutionseries_present.append(str(d))
                        if metric == 'auprc':
                            aux_metric.append(average_precision_score(df_method[truth_name], df_method[colname]))
                            aux_metricrelative.append(average_precision_score(
                                df_method[truth_name], df_method[colname]) - baseline[method])
                        elif metric.startswith('max'):
                            precision, recall, _ = precision_recall_curve(df_method[truth_name], df_method[colname])
                            #print(precision, recall)
                            f1list = 2*(precision * recall)/(precision + recall)
                            if metric == 'maxf1':
                                aux_metric.append(np.nanmax(f1list))
                            elif metric == 'maxf1precision':
                                #print(precision)
                                #print(baseline[method])
                                p = np.where(precision > baseline[method], precision, 0)
                                f1list = 2*(p * recall)/(p + recall)
                                #print(p)
                                if np.nanmax(f1list) >= 0.01:
                                    #print('f1max', np.nanmax(f1list))
                                    #print(list(f1list).index(np.nanmax(f1list)))
                                    aux_metric.append(precision[list(f1list).index(np.nanmax(f1list))])
                                else:
                                    aux_metric.append(0)
                            elif metric == 'maxf1recall':
                                #print(precision)
                                p = np.where(precision > baseline[method], precision, 0)
                                f1list = 2*(p * recall)/(p + recall)
                                #print(p)
                                if np.nanmax(f1list) >= 0.01:
                                    #print('f1max', np.nanmax(f1list))
                                    #print(list(f1list).index(np.nanmax(f1list)))
                                    aux_metric.append(recall[list(f1list).index(np.nanmax(f1list))])
                            elif 'maxrecallatleast' in metric:
                                threshold = float(metric.split('_')[1].split('precision')[0])/100
                                print(threshold)
                                p = np.where(precision >= threshold, 1, 0)
                                rlist = recall * p
                                #print(precision)
                                #print(rlist)
                                aux_metric.append(np.max(rlist))
                        elif metric == 'precision':
                            aux_metric.append(precision_score(df_method[truth_name], df_method[colnamebis]))
                            # print(df_method[factorprefix + '_' + method][df_method[truth_name] == True])
                        elif metric == 'recall':
                            aux_metric.append(recall_score(df_method[truth_name], df_method[colnamebis]))
                        elif metric == 'f1':
                            aux_metric.append(f1_score(df_method[truth_name], df_method[colnamebis]))
                        aux_method.append(method)
                        if ground_truth_method != 'spikein' and ground_truth_method != 'SEQC2':
                            aux_tb.append(round(dilutionseries.loc[d, 'tf'], 2))
                            aux_cov.append(int(dilutionseries.loc[d, 'cov']))
                        elif ground_truth_method == 'SEQC2':
                            aux_samplename.append(dilutionseries.loc[d, 'samplename'])
                        else:
                            aux_tb.append(round(dilutionseries.loc[d, 'vaf'], 3))
                else:
                    print('is not present')
                    # print(np.unique([cn.split('_')[0] for cn in list(df_table.columns)])[:-5].astype(float))
                    print(np.unique([cn.split('_')[0] for cn in list(df_table.columns[5:-1])]).astype(float))
    #elif xaxis == 'vaf':
        #df_table['median VAF'] = df_table[[c for c in list(df_table.columns) if c.endswith('vaf')]].median()
        #df_table['VAF'] = pd.cut(df_table['median VAF'],
        #                         bins=[0, .01, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, 1],
        #                         labels=['≤1%', '≤5%', '≤10%', '≤15%', '≤20%', '≤25%', '≤30%', '≤35%', '≤40%', '≤45%', '≤50%', '>50%'],
        #                         include_lowest=True)
        #df_table['VAF'] = df_table[[c for c in list(df_table.columns) if c.endswith('vaf')]]
    # print(baseline)
    results_df[metric.upper() + ' score'] = aux_metric
    if metric == 'auprc':
        results_df[metric.upper() + ' score - baseline ' + metric.upper() + ' score'] = aux_metricrelative
    # print(aux_tb)
    if ground_truth_method != 'spikein' and ground_truth_method != 'SEQC2':
        if fixedvar == 'coverage':
            results_df[xaxis] = aux_tb
            xvar = xaxis
            reverse = True
        elif fixedvar == 'ctdna':
            results_df['coverage'] = aux_cov
            xvar = 'coverage'
            reverse = False
        else: # fixedvar = 'vaf'
            xvar = 'median_vaf'
    elif ground_truth_method == 'SEQC2':
        results_df['samplename'] = aux_samplename
        xvar = 'samplename'
        reverse = False
    else:
        results_df['VAF'] = aux_tb
        xvar = 'VAF'
        reverse = True
    results_df['caller'] = aux_method
    if config.context == 'paper':
        sns.set_style("whitegrid")
    else:
        sns.set_style("darkgrid")
    print(results_df)
    sns.catplot(x=xvar, y=metric.upper() + " score", hue="caller",
                capsize=.1, height=6, aspect=1.5, kind="point",
                order=sorted(results_df[xvar].unique(), reverse=reverse),
                palette=sns.color_palette(list(color_dict.values())), data=results_df)
    plt.title(metric.upper() + " score for {} calling in {} - all chroms".format(muttype, mixtureid))

    # plt.ylim([0, max(0.5, max(aux_metric)+0.05)])
    if metric == 'auprc':
        if len(np.unique(baseline.values())) == 1:
            print(baseline[config.methods[0]])
            plt.axhline(y=baseline[config.methods[0]], color='k', linestyle='--')
    if type(ground_truth_method) == int:
        refname = 'in'+refsample + 'samplebyatleast' + str(ground_truth_method) +'callers'
    elif ground_truth_method == 'ranked':
        refname = 'in'+refsample + 'sampleranked'
    elif ground_truth_method == 'SEQC2':
        refname = 'inKnownVariants'
    else:
        refname = 'in'+refsample + 'samplebythesamecaller'
    if ground_truth_method == 'spikein':
        dilfolder = config.spikeinfolder
        dilution = 'spikeins'
    elif diltype == 'mixture':
        dilfolder = config.mixturefolder
        dilution = 'mixtures'
    elif diltype == 'mixture_wes':
        dilfolder = config.mixturefolderultradeep
        dilution = 'mixtures'
    elif diltype == 'mixture_wgs':
        dilfolder = config.mixturefolderwholegenome
        dilution = 'mixtures'
    elif diltype == 'SEQC2':
        dilfolder = config.mixturefolderSEQC2
        dilution = 'SEQC2'
    xa = xaxis if xaxis != 'tumor burden' else 'tb'
    if save:
        if not os.path.exists(os.path.join(*dilfolder, dilution+'_allchr')):
            os.mkdir(os.path.join(*dilfolder, dilution+'_allchr'))
        if not os.path.exists(os.path.join(*dilfolder, dilution+'_allchr', 'figures')):
            os.mkdir(os.path.join(*dilfolder, dilution+'_allchr', 'figures'))
        if methods != config.methods:
            plt.savefig(os.path.join(*dilfolder, dilution+'_allchr','figures', mixtureid + '_' + muttype + '_' + metric + '_' + refname + '_' + '_'.join(methods) + '_' + xa + '_' + config.context), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(*dilfolder, dilution+'_allchr', 'figures', mixtureid + '_' + muttype + '_' + metric + '_' + refname + '_' + xa + '_' + config.context), bbox_inches='tight')
    plt.show()
    summary_df = results_df.copy()
    summary_df['mutation type'] = muttype
    summary_df['metric'] = metric
    summary_df['mutation type'] = refname
    if save:
        if not os.path.exists(os.path.join(*dilfolder, dilution+'_allchr', 'results')):
            os.mkdir(os.path.join(*dilfolder, dilution+'_allchr', 'results'))
        if methods != config.methods:
            summary_df.to_csv(os.path.join(*dilfolder, dilution+'_allchr', 'results',  mixtureid + '_' + muttype + '_' + metric + '_' + refname + '_'.join(methods) + '_fixed' + fixedvar + '_' + xa + '.csv'))
        else:
            summary_df.to_csv(os.path.join(*dilfolder, dilution+'_allchr', 'results', mixtureid + '_' + muttype + '_' + metric + '_' + refname + '_fixed' + fixedvar + '_' + xa + '.csv'))
    return summary_df


if __name__ == "__main__":
    import os
    # set working directory
    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    # Config and Display paramaters
    from utils.config import Config
    from utils.viz import set_display_params
    config = Config("config/", "config_viz.yaml")
    set_display_params(config)
    print(config.methods)

    from benchmark.calltableseries import get_calltableseries
    from benchmark.groundtruth import generate_groundtruth

    chrom = '22'
    reload = False
    save = True
    filterparam = 'all'
    muttypes = ['snv', 'indel']
    mixtureids = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T']
    fixedvars = ['coverage', 'ctdna']

    for fixedvar in fixedvars:
        if fixedvar == 'coverage':
            dilutionseries = [(70, 0), (70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]
        elif fixedvar == 'ctdna':
            dilutionseries = [(70, 0), (70, 30), (70, 80), (70, 130), (70, 180)]
        for muttype in muttypes:
            if muttype == 'snv':
                nref = 5
            else:  # elif muttype == 'indel':
                nref = 3
            for mixtureid in mixtureids:
                print('############# {} {} ############'.format(mixtureid, muttype))
                plasmasample = '_'.join(mixtureid.split('_')[:2])
                print(plasmasample)
                healthysample = '_'.join(mixtureid.split('_')[2:])
                print(healthysample)
                calltablesseries, calltables = get_calltableseries(config, mixtureid, chrom, muttype, filterparam, reload, save)
                print(calltablesseries.head())
                print(calltables['tf'])
                for gtm in [nref, 'ranked']:
                    calltablesseries = generate_groundtruth(config, calltablesseries, calltables['tf'], ground_truth_method=gtm, muttype=muttype)
                    results_auprc_df = metric_curve(config, calltablesseries, plasmasample, healthysample, dilutionseries,
                                                    metric='auprc', ground_truth_method=gtm, refsample='undiluted',
                                                    muttype=muttype, chrom=chrom, methods=config.methods,
                                                    fixedvar=fixedvar, save=save)
                    results_recall_df = metric_curve(config, calltablesseries, plasmasample, healthysample, dilutionseries,
                                                     metric='recall', ground_truth_method=gtm, refsample='undiluted',
                                                     muttype=muttype, chrom=chrom, methods=config.methods,
                                                     fixedvar=fixedvar, save=save)
                    results_precision_df = metric_curve(config, calltablesseries, plasmasample, healthysample, dilutionseries,
                                                        metric='precision', ground_truth_method=gtm, refsample='undiluted',
                                                        muttype=muttype, chrom=chrom,
                                                        methods=config.methods,
                                                        fixedvar=fixedvar, save=save)
