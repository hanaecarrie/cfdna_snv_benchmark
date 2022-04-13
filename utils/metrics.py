import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import warnings
from sklearn.metrics import precision_recall_curve, roc_curve, f1_score, precision_score, recall_score, average_precision_score

warnings.filterwarnings("ignore")


def plot_pr_curve(precision, recall, estimator_name=None, f1_score=None, figax=None, kwargs={}):
    kwargs["drawstyle"] = "steps-post"
    if f1_score is not None and estimator_name is not None:
        kwargs["label"] = f"{estimator_name} (F1 = {f1_score:0.2f})"
    elif f1_score is not None:
        kwargs["label"] = f"AP = {f1_score:0.2f}"
    elif estimator_name is not None:
        kwargs["label"] = estimator_name
    fig, ax = figax if figax is not None else plt.subplots()
    ax.plot(recall, precision, **kwargs)
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


def figure_curve(config, df_table, plasmasample, healthysample, dilutionseries, xy='pr', ground_truth_method=3, refsample='undiluted', muttype='SNV', chrom='22', methods=None, fixedvar='coverage', save=True):
    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
    alpha_dict = {str(dilutionseries[i]): 1-0.1*i for i in range(len(dilutionseries))}
    tb_dict, cov_dict = {}, {}
    baseline_dict = {}
    dilutionseries_present = []
    if methods is None:
        methods = config.methods
    for i, d in enumerate(dilutionseries):
        if ground_truth_method != 'spikein':
            mixturepath = 'mixture_chr'+chrom+'_'+plasmasample +"_" + str(d[0]) +"x_" + healthysample + "_" + str(d[1]) + 'x'
            tb_dict[str(d)] = float(pd.read_csv(os.path.join(
                *config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+plasmasample+'_'+healthysample,
                mixturepath, 'estimated_tf_chr'+chrom+mixturepath[len(('mixture_chr'+chrom)):]+'.txt')).columns[0])
            cov_dict[str(d)] = float(pd.read_csv(os.path.join(
                *config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+plasmasample+'_'+healthysample,
                mixturepath, 'coverage_chr'+chrom+mixturepath[len(('mixture_chr'+chrom)):]+'.txt')).columns[0])
            print(tb_dict)
    for method in methods:
        fig, ax = plt.subplots(figsize=(10, 10))
        for i, d in enumerate(dilutionseries):
            if ground_truth_method != 'spikein':
                factorprefix = '{:.2f}'.format(100*round(tb_dict[str(d)], 4))
            else:
                factorprefix = '{:.2f}'.format(d)
            print(factorprefix + '_' + method + '_score')
            if factorprefix + '_' + method + '_score' in list(df_table.columns):
                if type(ground_truth_method) == int or ground_truth_method == 'spikein':
                    truth_name = 'truth'
                elif ground_truth_method == 'method':
                    truth_name = method +'_truth'
                else:
                    raise ValueError('unknown ground truth {}'.format(ground_truth_method))
                df_method = df_table[[truth_name, factorprefix + '_' + method + '_score']]
                df_method[truth_name].fillna(False, inplace=True)
                df_method[factorprefix + '_' + method + '_score'].fillna(0, inplace=True)
                baseline_dict[str(d)] = len(df_method[truth_name][df_method[truth_name]])/len(df_method[truth_name])
                if str(d) not in dilutionseries_present:
                    dilutionseries_present.append(str(d))
                if xy == 'pr':
                    precision, recall, thresholds = precision_recall_curve(df_method[truth_name], df_method[factorprefix + '_' + method + '_score'])
                    if i == 0:
                        plot_pr_curve(precision, recall, estimator_name=method, f1_score=None, figax=(fig, ax),
                                      kwargs={'color': color_dict[method], 'alpha': alpha_dict[str(d)], 'lw': 4-i/3})
                    else:
                        plot_pr_curve(precision, recall, estimator_name='', f1_score=None, figax=(fig, ax),
                                      kwargs={'color': color_dict[method], 'alpha': alpha_dict[str(d)], 'lw': 4-i/3})
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
                    if ground_truth_method != 'spikein':
                        list_lines_baseline.append(Line2D([0], [0], color='black', ls='--', alpha=alpha_dict[str(d)], label="baseline tf {:.2f}% = {:.2f}".format(100*tb_dict[str(d)], baseline_dict[str(d)])))
                    else:
                        list_lines_baseline.append(Line2D([0], [0], color='black', ls='--', alpha=alpha_dict[str(d)], label="baseline vaf {:.2f}% = {:.2f}".format(100*tb_dict[str(d)], baseline_dict[str(d)])))
        handles, labels = plt.gca().get_legend_handles_labels()
        if fixedvar == 'coverage':
            if ground_truth_method != 'spikein':
                list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[str(i)], label='tumor burden = {:.2f}%'.format(100*tb_dict[str(i)])) for i in dilutionseries_present]
            else:
                list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[str(i)], label='VAF = {:.2f}%'.format(float(i))) for i in dilutionseries_present]
        elif fixedvar == 'ctdna':
            list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[str(i)], label='cov = {}x'.format(int(round(cov_dict[str(i)]/10)*10))) for i in dilutionseries_present]
        if xy == 'pr':
            legend_list = handles + list_lines + list_lines_baseline
        else:
            legend_list = handles + list_lines
        # Creating legend with color box
        plt.legend(bbox_to_anchor=(1, 1), loc="upper left", handles=legend_list)
        if xy == 'pr':
            plt.title("Precision Recall curve for SNV calling in sample {}".format(plasmasample + '_' + healthysample))
        elif xy == 'roc':
            plt.title("Receiver Operation Characteristics curve for SNV calling in sample {}".format(plasmasample + '_' + healthysample))
        if xy == 'pr':
            plt.semilogx()
            plt.xlim([0.01, 1.01])
        else:
            plt.xlim([-0.01, 1.01])
        plt.ylim([-0.01, 1.01])
        dilution = 'spikeins' if ground_truth_method == 'spikein' else 'mixtures'
        dilfolder = config.spikeinfolder if ground_truth_method == 'spikein' else config.mixturefolder
        if save:
            if not os.path.exists(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'figures')):
                os.mkdir(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'figures'))
            if type(ground_truth_method) == int:
                refname = 'in'+refsample + 'samplebyatleast' + str(ground_truth_method) +'callers'
            else:
                refname = 'in'+refsample + 'samplebythesamecaller'
            plt.savefig(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'figures', plasmasample + '_' + healthysample + '_' + muttype + '_' + xy.upper() + 'curve_' + refname + '_' + method + '_' + config.context), bbox_inches='tight')
        plt.show()


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
    ax.legend(loc="upper right")
    plt.plot([0, 1], [0, 1], ls='--', color='grey')
    plt.annotate('baseline', xy=(0.9, 0.9), color='grey')


def metric_curve(config, df_table, plasmasample, healthysample, dilutionseries, metric='auprc', ground_truth_method=3,
                 refsample='undiluted', muttype='SNV', chrom='22', methods=None, fixedvar='coverage', save=True):
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
            mixturepath = 'mixture_chr'+chrom+'_'+plasmasample +"_" + str(d[0]) +"x_" + healthysample + "_" + str(d[1]) + 'x'
            tb_dict[str(d)] = float(pd.read_csv(os.path.join(
                *config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+plasmasample+'_'+healthysample,
                mixturepath, 'estimated_tf_chr'+chrom+mixturepath[len(('mixture_chr'+chrom)):]+'.txt')).columns[0])
            cov_dict[str(d)] = float(pd.read_csv(os.path.join(
                *config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+plasmasample+'_'+healthysample,
                mixturepath, 'coverage_chr'+chrom+mixturepath[len(('mixture_chr'+chrom)):]+'.txt')).columns[0])
    results_df = pd.DataFrame()
    aux_metric = []
    aux_metricrelative = []
    aux_method = []
    aux_tb = []
    aux_cov = []
    baseline = {}
    for i, d in enumerate(dilutionseries):
        if ground_truth_method != 'spikein':
            factorprefix = '{:.2f}'.format(100*round(tb_dict[str(d)], 4))
        else:
            factorprefix = '{:.2f}'.format(d)
        for method in methods:
            if factorprefix + '_' + method + '_score' in list(df_table.columns):
                if i != 0 or ground_truth_method == 'spikein':
                    if type(ground_truth_method) == int or ground_truth_method == 'spikein':
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
                        # print(df_method[factorprefix + '_' + method][df_method[truth_name] == True])
                    elif metric == 'recall':
                        aux_metric.append(recall_score(df_method[truth_name], df_method[factorprefix + '_' + method]))
                    elif metric == 'f1':
                        aux_metric.append(f1_score(df_method[truth_name], df_method[factorprefix + '_' + method]))
                    aux_method.append(method)
                    if ground_truth_method != 'spikein':
                        aux_tb.append(round(100*tb_dict[str(d)], 3))
                        aux_cov.append(int(cov_dict[str(d)]))
                    else:
                        aux_tb.append(d)
    # print(baseline)
    results_df[metric.upper() + ' score'] = aux_metric
    if metric == 'auprc':
        results_df[metric.upper() + ' score - baseline ' + metric.upper() + ' score'] = aux_metricrelative
    # print(aux_tb)
    if ground_truth_method != 'spikein':
        if fixedvar == 'coverage':
            results_df['tumor burden'] = aux_tb
            xvar = 'tumor burden'
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
    else:
        refname = 'in'+refsample + 'samplebythesamecaller'
    dilution = 'spikeins' if ground_truth_method == 'spikein' else 'mixtures'
    dilfolder = config.spikeinfolder if ground_truth_method == 'spikein' else config.mixturefolder
    if save:
        if not os.path.exists(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'figures')):
            os.mkdir(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'figures'))
        if methods != config.methods:
            plt.savefig(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'figures',  plasmasample + '_' + healthysample + '_' + muttype + '_' + metric + '_' + refname + '_' + '_'.join(methods) + '_' + config.context), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'figures', plasmasample + '_' + healthysample + '_' + muttype + '_' + metric + '_' + refname + '_' + config.context), bbox_inches='tight')
    plt.show()
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
            summary_df.to_csv(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'results',  plasmasample + '_' + healthysample + '_' + muttype + '_' + metric + '_' + refname + '_'.join(methods) + '_fixed'+ fixedvar + '.csv'))
        else:
            summary_df.to_csv(os.path.join(*dilfolder, dilution+'_chr'+chrom, dilution+'_chr'+chrom+'_'+plasmasample+'_'+healthysample, 'results', plasmasample + '_' + healthysample + '_' + muttype + '_' + metric + '_' + refname + '_fixed'+ fixedvar + '.csv'))
    return summary_df


if __name__ == "__main__":
    print('TODO')