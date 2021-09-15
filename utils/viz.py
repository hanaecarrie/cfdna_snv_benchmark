import gzip
import io
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings("ignore")


def load_files(filenames):
    for filename in filenames:
        yield pd.read_csv(filename, names=['sample_id', 'tumor_burden'])


def set_display_params(config):
    import seaborn as sns
    import matplotlib.pyplot as plt
    if config.context == 'talk':
        sns.set(style=config.context_talk.style, context=config.context_talk.context,
                rc={"lines.linewidth": config.context_talk.llw, "legend.fontsize": config.context_talk.llw})
        plt.style.use(config.context_talk.styleuse)
        plt.rcParams.update({"grid.linewidth": config.context_talk.glw, "grid.alpha": config.context_talk.galpha,
                             'font.size': config.context_talk.fs})
        sns.set_palette(config.context_talk.palette)
    elif config.context == 'paper':
        print(config.context_paper)
        sns.set(style=config.context_paper.style, context=config.context_paper.context,
                rc={"lines.linewidth": config.context_paper.llw, "legend.fontsize": config.context_paper.llw})
        plt.style.use(config.context_paper.styleuse)
        plt.rcParams.update({"grid.linewidth": config.context_paper.glw, "grid.alpha": config.context_paper.galpha,
                             'font.size': config.context_paper.fs})
        sns.set_palette(config.context_paper.palette)
    else:
        raise ValueError("unknown context {}. Should be either 'talk' or 'paper'".format(config.context))


def read_vcf(path):
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


def load_calls_from_vcf(vcf_path, methods, chrom='all'):
    sample = None
    for m in methods:
        vcf_path_m = vcf_path.split('ensemble')[0] + m + vcf_path.split('ensemble')[1]
        if os.path.exists(vcf_path_m) or os.path.exists(vcf_path_m + '.gz'):
            #print(vcf_path_m)
            res = read_vcf(vcf_path_m)
            res = res[res['FILTER'] == "PASS"]
            # res['callers'] = res['INFO'].apply(lambda x: pd.Series(x.split('CALLERS=')[1].split(';')[0]))
            res['type'] = np.nan
            res['type'][res['ALT'].str.len() - res['REF'].str.len() == 0] = 'SNV'
            res['type'][res['ALT'].str.len() - res['REF'].str.len() > 0] = 'INS'
            res['type'][res['ALT'].str.len() - res['REF'].str.len() < 0] = 'DEL'
            res['type'][res['ID'].str.contains('rs')] = 'SNP'
            # for m in methods:
            #    res[m] = res['INFO'].str.contains(m)
            info = res['INFO']
            res = res[['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'type']]
            res.rename(columns={'FILTER': m}, inplace=True)
            res[m] = True
            # sample = res[['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'type', *methods]]
            # for m in methods:
            if m == 'vardict':  # P-value
                res[m + '_score'] = [float(i.split('SSF=')[1].split(';')[0]) if 'SSF' in i else np.nan for i in info.tolist()]
            elif m == 'varscan':  # P-value
                res[m + '_score'] = [float(i.split('SSC=')[1].split(';')[0]) if 'SSC' in i else np.nan for i in info.tolist()]
            elif m == 'mutect2':  # logodds to probability score prob = exp(logTLOD)/(1+exp(logTLOD))
                res[m + '_score'] = [np.exp(np.log(float(i.split('TLOD=')[1].split(';')[0]))) / (1 + np.exp(np.log(float(i.split('TLOD=')[1].split(';')[0])))) if 'TLOD' in i else np.nan for i in info.to_list()]
                # [float(i.split('TLOD=')[1].split(';')[0]) if 'TLOD' in i else np.nan for i in info.tolist()]
                # [np.exp(float(i.split('TLOD=')[1].split(';')[0])) / (1 + np.exp(float(i.split('TLOD=')[1].split(';')[0]))) if 'TLOD' in i else np.nan for i in res['INFO']]
            elif m == 'freebayes':  # logodds to probability score prob = exp(logODDS)/(1+exp(logODDS))
                res[m + '_score'] = [np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0]))) / (1 + np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0])))) if 'ODDS' in i else np.nan for i in info.to_list()]
                    # [float(i.split('ODDS=')[1].split(';')[0]) if 'ODDS' in i else np.nan for i in info.tolist()]
                # [float(i.split('ODDS=')[1].split(';')[0]) if 'ODDS' in i else np.nan for i in res['INFO']]
                # [np.exp(float(i.split('ODDS=')[1].split(';')[0])) / (1 + np.exp(float(i.split('ODDS=')[1].split(';')[0]))) if 'ODDS' in i else np.nan for i in res['INFO']]
            elif m == 'strelka2':  # phred score to probability, prob = 1 - 10^(-SomaticEVS/10)
                res[m + '_score'] = [1 - (10 ** (-float(i.split('SomaticEVS=')[1].split(';')[0]) / 10)) if 'SomaticEVS' in i else np.nan for i in info.to_list()]
                # [float(i.split('SomaticEVS=')[1].split(';')[0]) if 'SomaticEVS' in i else np.nan for i in info.tolist()]
                # [1 - (10 ** (-float(i.split('SomaticEVS=')[1].split(';')[0]) / 10)) if 'SomaticEVS' in i else np.nan for i in res['INFO']]
            res['CHROM_POS'] = res['CHROM'].astype('str').str.cat(res['POS'].astype('str'), sep="_")
            res.set_index('CHROM_POS', inplace=True)
            if sample is None:
                sample = res.copy()
            else:
                res.drop(['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'type'], axis=1, inplace=True)
                sample = sample.join(res, how='outer')  # pd.concat([sample, res], axis=1)
    if chrom != 'all' and sample is not None:
        print('select a single chrom = {} for analysis'.format(chrom))
        sample = sample[sample['CHROM'] == chrom]
        return sample
    else:
        print("sample is not present with path {}".format(vcf_path))
        return None


def load_calls_from_vcf_dilutionseries(dirpath, plasmasample, reference, dilutionseries, methods,
                                       prefix='dilution_chr22', chrom='all'):
    vcf_samples_dict = {}
    samples_dict = {}
    ci = 0
    dilutionseries_new = []
    for i, d in enumerate(dilutionseries):
        print('vcf_pd_' + str(ci), d)
        d0 = str(d[0]).replace('.', '_')
        d1 = str(d[1]).replace('.', '_')
        if d == (1, 0):
            ref = 'pooledhealthy'
        vcf_path = os.path.join(*dirpath, prefix + plasmasample + "_" + str(d[0]) + "_" + ref + "_" + str(d[1]),
                                prefix + plasmasample + "_" + d0 + "_" + ref + "_" + d1 + "-ensemble-annotated.vcf")
        if d == (1, 0):
            ref = reference

        vcf_sample = load_calls_from_vcf(vcf_path, methods, chrom=chrom)
        if vcf_sample is None:
            print("dilution {} is not present with path {}".format(d, vcf_path))
        else:
            vcf_samples_dict['sample_' + str(ci)] = vcf_sample
            ci += 1
            dilutionseries_new.append(d)
    return vcf_samples_dict, dilutionseries_new


def count_mutations(samples, methods, samples_tf, mutationtypes=['all', 'INDEL', 'SNV', 'SNP'], threshold=0.5):
    for mutationtype in mutationtypes:
        if mutationtype not in ['all', 'INDEL', 'SNV', 'SNP']:
            raise ValueError(
                'mutation type {} is unknown. It should be in {}'.format(mutationtype, ['all', 'INDEL', 'SNV', 'SNP']))
        numbersnvs_pd = pd.DataFrame()
        for si, s in enumerate(samples):
            nb_snv = []
            for method in methods:
                if mutationtype == 'all':
                    nb_snv.append(s[s[method] > 0].shape[0])
                elif mutationtype == 'INDEL':
                    nb_snv.append(s[(s[method] > 0) & ((s['type'] == 'INS') | (s['type'] == 'DEL'))].shape[0])
                elif mutationtype == 'SNV':
                    nb_snv.append(s[(s[method] > 0) & (s['type'] == 'SNV')].shape[0])
                elif mutationtype == 'SNP':
                    nb_snv.append(s[(s[method] > 0) & (s['type'] == 'SNP')].shape[0])
                else:
                    raise ValueError('mutation type {} unknown'.format(mutationtype))
            if si == 0:
                numbersnvs_pd = pd.DataFrame.from_dict({'sample_' + str(si): nb_snv}).T
                numbersnvs_pd.columns = methods
            else:
                numbersnvs_pd.loc['sample_' + str(si)] = nb_snv
            numbersnvs_pd = numbersnvs_pd.rename(index=samples_tf)
        return numbersnvs_pd


def pr_table(results_df, value_colname, truth_colname, verbose=False, decreasing=False):
    # get PR table of values
    res_df = results_df[[value_colname, truth_colname]]
    res_df.dropna(inplace=True)
    res_df.columns = ['value', 'truth']
    thresholds = np.sort(np.unique(res_df['value'].values))
    # print(len(thresholds))
    thresholds = thresholds[:len(thresholds)]
    if len(thresholds) > 10000:  # TODO
        seq_idx = np.linspace(thresholds[0], len(thresholds), 10000)  # length=1000
        thresholds = thresholds[seq_idx]

    table = pd.DataFrame(columns=['Threshold', 'TP', 'FP', 'FN', 'Precision', 'Recall', 'F1'])
    table['Threshold'] = thresholds
    print('Total thresholds: {}'.format(len(thresholds)))

    for i in range(len(thresholds)):
        if verbose:
            print('THRESHOLD ', thresholds[i])
        res_df['predict'] = False
        if decreasing:
            res_df.loc[res_df['value'] >= thresholds[i], 'predict'] = True
        else:
            res_df.loc[res_df['value'] <= thresholds[i], 'predict'] = True

        # print(res_df[(res_df['predict'] == True)])
        # print(res_df[(res_df['predict'] == True) & (res_df['truth'] == True)])

        table['TP'].iloc[i] = res_df.loc[(res_df['predict'] == True) & (res_df['truth'] == True)].shape[0]
        table['FP'].iloc[i] = res_df.loc[(res_df['predict'] == True) & (res_df['truth'] == False)].shape[0]
        table['FN'].iloc[i] = res_df.loc[(res_df['predict'] == False) & (res_df['truth'] == True)].shape[0]
        print(table.iloc[i])

    table['Precision'] = table['TP'] / (table['TP'] + table['FP'])
    table['Recall'] = table['TP'] / (table['TP'] + table['FN'])
    table['F1'] = (2 * table['Precision'] * table['Recall']) / (table['Precision'] + table['Recall'])

    return table


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
        plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02), color='grey')


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


if __name__ == "__main__":
    """
    dilutiondirpath = ["..", "data", "dilutions_chr22"]
    bcbiooutputdirpath = ["..", "data", "bcbio_output"]
    prefix = 'dilution_chr22_'
    chrom = '22'
    plasmasample1 = 'CRC-986_100215'
    methods = ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']
    vcf_path = os.path.join(*bcbiooutputdirpath, prefix + plasmasample1 + "_1_pooledhealthy_0",
                            prefix + plasmasample1 + "_1_pooledhealthy_0-ensemble-annotated.vcf")
    vcf_ref = load_calls_from_vcf(vcf_path, methods, chrom=chrom)
    vcf_ref = vcf_ref[vcf_ref['type'] == 'SNV']
    y_true = pd.DataFrame(vcf_ref[methods].T.sum())
    y_true.columns = ['truth']
    y_true['truth'][y_true['truth'] <= 1] = 0
    y_true = y_true.astype(bool)
    #y_true['CHROM'] = y_true.index.str.split('_').str[0]
    #y_true['POS'] = y_true.index.str.split('_').str[1]
    #print(y_true)
    #vcf_path = os.path.join(*bcbiooutputdirpath, prefix + plasmasample1 + "_1_pooledhealthy_0",
    #                        prefix + plasmasample1 + "_1_pooledhealthy_0-ensemble-annotated.vcf")
    vcf_path = os.path.join(*bcbiooutputdirpath, prefix + plasmasample1 + "_1_CRC-986_300316_0.72",
                            prefix + plasmasample1 + "_1_CRC-986_300316_0_72-ensemble-annotated.vcf")
    vcf_sample = load_calls_from_vcf(vcf_path, methods, chrom=chrom)
    vcf_sample = vcf_sample[vcf_sample['type'] == 'SNV']
    print(vcf_sample.shape, y_true.shape)
    #vcf_sample = pd.merge(vcf_sample, y_true, how='outer', on=["CHROM", "POS"])
    #df_sample = pd.concat([vcf_sample, y_true], join='outer', axis=1)
    df_sample = vcf_sample.join(y_true, how='outer') # not concat because duplicated axis
    # prt = pr_table(vcf_sample, 'freebayes_score', 'truth', verbose=True, decreasing=False)
    # print(prt)

    from sklearn.metrics import precision_recall_curve

    fig, ax = plt.subplots()
    for method in methods:
        df_sample_method = df_sample[['truth', method + '_score']]
        df_sample_method['truth'].fillna(False, inplace=True)
        df_sample_method[method + '_score'].fillna(0, inplace=True)
        print(df_sample_method)
        precision, recall, thresholds = precision_recall_curve(df_sample_method['truth'], df_sample_method[method + '_score'])
        plot_pr_curve(precision, recall, estimator_name=method, f1_score=None, figax=(fig, ax))
    plt.show()
    

    dilutiondirpath = ["..", "data", "dilutions_chr22"]
    bcbiooutputdirpath = ["..", "data", "bcbio_output"]
    prefix = 'dilution_chr22_'
    chrom = '22'
    plasmasample1 = 'CRC-986_100215'
    methods = ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']
    #vcf_path = os.path.join(*bcbiooutputdirpath, prefix + plasmasample1 + "_1_CRC-986_300316_0.72",
    #                        prefix + plasmasample1 + "_1_CRC-986_300316_0_72-ensemble-annotated.vcf")
    vcf_path = os.path.join(*bcbiooutputdirpath, prefix + plasmasample1 + "_0.125_CRC-986_300316_0.875",
                            prefix + plasmasample1 + "_0_125_CRC-986_300316_0_875-ensemble-annotated.vcf")
    vcf_sample = load_calls_from_vcf(vcf_path, methods, chrom=chrom)
    vcf_sample = vcf_sample[vcf_sample['type'] == 'SNV']
    vcf_path = os.path.join(*bcbiooutputdirpath, prefix + plasmasample1 + "_1_pooledhealthy_0",
                            prefix + plasmasample1 + "_1_pooledhealthy_0-ensemble-annotated.vcf")
    vcf_ref = load_calls_from_vcf(vcf_path, methods, chrom=chrom)
    vcf_ref = vcf_ref[vcf_ref['type'] == 'SNV']

    from sklearn.metrics import precision_recall_curve

    fig, ax = plt.subplots()
    for method in methods:
        y_true = vcf_ref[[method]]
        y_true.columns = ['truth']
        y_true = y_true.astype(bool)
        df_sample = vcf_sample.join(y_true, how='outer') # not concat because duplicated axis
        df_sample_method = df_sample[['truth', method + '_score']]
        df_sample_method['truth'].fillna(False, inplace=True)
        df_sample_method[method + '_score'].fillna(0, inplace=True)
        print(method)
        print(df_sample_method[method + '_score'].describe())
        precision, recall, thresholds = precision_recall_curve(df_sample_method['truth'], df_sample_method[method + '_score'])
        plot_pr_curve(precision, recall, estimator_name=method, f1_score=None, figax=(fig, ax))
    plt.show()
    """
    bcbiooutputdirpath = ["..", "data"]
    methods = ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']
    samples = ['icgc_cll_tumour']# , 'icgc_cll_T40_tumour', 'icgc_cll_T20_tumour']
    vcf_path = os.path.join(*bcbiooutputdirpath, "SMURF benchmark", "snv-smurf-test20-v7-new5-maxdepth-2.tsv")
    vcf_samples = pd.read_csv(vcf_path, sep='\t')
    alpha_dict = {'icgc_cll_tumour': 1, "icgc_cll_T40_tumour": 0.8, 'icgc_cll_T20_tumour': 0.5}
    color_dict = {'freebayes': 'tab:blue', 'mutect2': 'tab:orange', 'strelka2': 'tab:green', 'vardict': 'tab:red', 'varscan': 'tab:purple' }
    fig, ax = plt.subplots()
    for sample in samples:
        vcf_sample = vcf_samples[vcf_samples['Sample_Name'] == sample]
        vcf_sample.reset_index(inplace=True)
        vcf_sample['CHROM_POS'] = vcf_sample['X.CHROM'].astype('str').str.cat(vcf_sample['POS'].astype('str'), sep="_")
        vcf_sample.set_index('CHROM_POS', inplace=True)
        vcf_sample['mutect2_score'] = vcf_sample['m2_TLOD']
        vcf_sample['freebayes_score'] = vcf_sample['f_ODDS']
        vcf_sample['strelka2_score'] = vcf_sample['s2_SomaticEVS']
        vcf_sample['varscan_score'] = vcf_sample['vs_SSC']
        vcf_sample['vardict_score'] = vcf_sample['vd_SSF']
        from sklearn.metrics import precision_recall_curve

        for method in methods:
            df_sample_method = vcf_sample[[method+'_score', 'TRUTH']]
            # df_sample = vcf_sample.join(y_true, how='outer') # not concat because duplicated axis
            # df_sample_method = df_sample[['truth', method + '_score']]
            # df_sample_method['truth'].fillna(False, inplace=True)
            df_sample_method[method + '_score'].fillna(0, inplace=True)
            print(method)
            print(df_sample_method[method + '_score'].describe())
            precision, recall, thresholds = precision_recall_curve(df_sample_method['TRUTH'], df_sample_method[method + '_score'])
            plot_pr_curve(precision, recall, estimator_name=method, f1_score=None, figax=(fig, ax), color=color_dict[method], alpha=alpha_dict[sample])
    plt.show()
