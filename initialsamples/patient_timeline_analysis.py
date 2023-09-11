# patient_timeline_analysis.py

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


def load_files(filenames):
    """Loads csv file indicated as input and yields dataframe with sample id and matched ichorCNA tumor fraction.

    Parameters
    ----------
    filenames : list
        list of str.

    Yields
    ------
    pandas DataFrame instance
        a dataframe containing sample ID = patient ID + date, ichorCNA tumor fraction, patient ID and date.
    """

    for filename in filenames:
        yield pd.read_csv(filename, names=['sample_id', 'tumor_burden'])


def get_tf(config, batch='all'):
    """Gets the ichorCNA tumor fractions of one or all batches of longitudinal cfDNA samples.

    Parameters
    ----------
    config : object
        configuration class instance containing paths and batch names
    batch : str, optional
        indicate a single batch or all batches of cfDNA samples. Can be 'batch1', 'batch2', 'deepwgs or 'all'. (default is 'all')

    Returns
    -------
    pandas DataFrame instance
        a dataframe containing sample ID = patient ID + date, ichorCNA tumor fraction, patient ID and date.
    """

    if batch == 'deepwgs':
        tffile = os.path.join(os.getcwd(), *config.tumorfractionfiles.deepwgs)
        tf_df = pd.read_csv(tffile, names=['sample_id', 'tumor_burden'])
    elif batch == 'batch1':
        tffile = os.path.join(os.getcwd(), *config.tumorfractionfiles.batch1)
        tf_df = pd.read_csv(tffile, names=['sample_id', 'tumor_burden'])
    elif batch == 'batch2':
        tffile = os.path.join(os.getcwd(), *config.tumorfractionfiles.batch2)
        tf_df = pd.read_csv(tffile, names=['sample_id', 'tumor_burden'])
    elif batch == 'all':
        tffiles = [os.path.join(os.getcwd(), *config.tumorfractionfiles[batch]) for batch in config.batches]
        tf_df = pd.concat(load_files(tffiles))
    else:
        raise ValueError("unknown batch {}. batch should be either 'batch1', 'batch2', 'deepwgs or 'all'".format(batch))
    tf_df['patient'] = [int(tf_df['sample_id'].iloc[i].split('-')[1].split('_')[0])
                        if tf_df['sample_id'].iloc[i].startswith('sortpanel')
                        else int(tf_df['sample_id'].iloc[i].split('_')[0]) for i in range(tf_df.shape[0])]
    tf_df['date'] = pd.to_datetime(tf_df['sample_id'].str.split('_').str[1], format='%d%m%y')
    return tf_df


def get_patient_tf_treatment(config, patient):
    """Gets the ichorCNA tumor fractions of one or all batches of longitudinal cfDNA samples.

    Parameters
    ----------
    config : object
        configuration class instance containing paths and batch names
    patient : str or int
        patient number

    Returns
    -------
    pandas DataFrame instance
        a dataframe containing sample ID = patient ID + date, ichorCNA tumor fraction, patient ID and date.
    """

    patient = int(patient)
    # Tumor burden
    tf_df = get_tf(config, batch='all')
    tf_df['patient'] = [int(tf_df['sample_id'].iloc[i].split('-')[1].split('_')[0])
                        if tf_df['sample_id'].iloc[i].startswith('sortpanel')
                        else int(tf_df['sample_id'].iloc[i].split('_')[0]) for i in range(tf_df.shape[0])]
    tf_df['date'] = pd.to_datetime(tf_df['sample_id'].str.split('_').str[1], format='%d%m%y')
    tf_patient = tf_df[tf_df['patient'] == patient].sort_values('date')
    tf_patient['date'] = tf_patient['date'].astype(str)
    # Treatment file
    treatment_df = pd.read_csv(os.path.join(os.getcwd(), *config.treatmentfile), sep='\t')
    treatment_df['patient'] = treatment_df['patient'].astype(int)
    treatment_df['date'] = pd.to_datetime(treatment_df['date'], format='%Y-%m-%d')
    treatment_df.rename(columns={'value': 'treatment'}, inplace=True)
    treatment_file = treatment_df[['patient', 'date', 'treatment']]
    treatment_patient = treatment_file[treatment_file['patient'] == patient].sort_values('date')
    treatment_patient['date'] = treatment_patient['date'].astype(str)
    # Patient treatment + tumor fraction infos
    df_patient = pd.concat([treatment_patient, tf_patient])
    df_patient = df_patient.sort_values('date')
    return df_patient, tf_patient, tf_df


def plot_patient_timeline(config, patient, treatment=True, mutations=False, highlight='discovery', figsize=(40, 10), save=False, savepath=None):
    """Plot timeline of patient with ichorCNA tumor fraction estimate evolution from lpWGS, and eventually,
    treatment information and/or mutation VAF evolution from deep targeted samples.
    Highlight in orange the high TF cfDNA samples and in blue the low low cfDNA sample.

    Parameters
    ----------
    config : object
        configuration class instance containing paths and batch names
    patient : str or int
        patient number
    treatment : bool
        whether to display treatment information or not annotated along the left y-axis. (default is True)
    mutations : bool
        whether to display mutation VAF information or not annotated along the right y-axis. (default is True)
    highlight : str
        highlights the high and low TF cfDNA samples that match criteria below if 'discovery' mode
        or that are already listed in configuration file if 'validation' mode
        If 'discovery', the high TF cfDNA samples have ichorCNA TF = 0% (or if None has, criteria is relaxed to TF ≤ 10%)
        the low TF cfDNA samples have
    figsize : tuple
         width and height of the timeline plot. (default is (40, 10))
    save : bool
        whether to save the plot or not. (default is False)
    savepath : None
        path where plot is saved if save=True. If None, plot is saved to config outputpath folder. (default is None)
    """

    patient = int(patient)
    # Display parameter
    lc = 'w' if config.context == 'talk' else 'k'
    if config.context == 'talk':
        fs = config.context_talk.fs
        lwd = config.context_talk.llw
    elif config.context == 'paper':
        fs = config.context_paper.fs
        lwd = config.context_paper.llw
    # Patient treatment + tumor fraction information
    df_patient, tf_patient, _ = get_patient_tf_treatment(config, patient)
    alldates = sorted(list(set(list(df_patient['date'].values) + list(tf_patient['date'].values))))
    tumorburden_dates = tf_patient[tf_patient['patient'] == patient]['date'].sort_values().astype(str).unique()
    # Start plot
    plt.figure()
    fig, ax2 = plt.subplots(figsize=figsize)
    df_patient['date'] = df_patient['date'].astype(str)
    df_patient['treatment'] = df_patient['treatment'].astype(str)
    df_patient['treatment'][df_patient['treatment'] == 'nan'] = 'cfDNA'
    # Display treatment information on left y-axis when treatment = True
    if treatment:
        sns.stripplot(y='treatment', x='date', data=df_patient, color='grey', marker='X', size=10, ax=ax2)
    ax2.grid(False)
    plt.legend((), ())
    # Display ichorCNA tumor fraction estimate in range [0,1] from lpWGS sample
    # on left y-axis if available, on right y-axis otherwise
    if treatment:
        ax = ax2.twinx()
    else:
        ax = ax2
    ele0 = ax.plot(df_patient['date'], df_patient['tumor_burden'], lc + '.', marker='s', markersize=10,
                   label='ichorCNA tumor burden')
    ax.plot(df_patient['date'][~df_patient['tumor_burden'].isna()],
            df_patient['tumor_burden'][~df_patient['tumor_burden'].isna()], lc + '-', linewidth=4)
    ax.plot(df_patient['date'][~df_patient['tumor_burden'].isna()],
            df_patient['tumor_burden'][~df_patient['tumor_burden'].isna()], lc + '.', marker='s', markersize=10)
    ax.set_ylabel('ctDNA or VAF fraction', fontsize=30)
    if treatment:
        ax.set_ylim(-0.01, 1)
    else:
        ax.set_ylim(-0.01, 0.8)
    fig.legend([ele0], ['ichorCNA tumor burden'], loc='upper left')
    labels = [ad if ad in tumorburden_dates else '' for ad in alldates]
    ax2.set_xticklabels(labels, rotation=90, fontsize=fs)
    ax.set_xticklabels(labels, rotation=90, fontsize=fs)
    plt.title('Patient {}: VAF evolution across timepoints'.format(patient))
    # Display mutations VAF in range [0,1] detected from targeted sequencing sample
    # on left y-axis if available, on right y-axis otherwise
    if mutations:
        # get timepoints with low ichorCNA tumor fraction estimate from lpWGS sample
        date_lowtftimepoints = list(tf_patient[tf_patient['tumor_burden'] == 0]['date'].unique())
        if not date_lowtftimepoints:  # if no samples with TF = 0%
            print('no zero ichorCNA estimate tumor burden')
            print('min tumor burden is {}'.format(min(tf_patient['tumor_burden'])))
            if min(tf_patient['tumor_burden']) < 0.1:  # relax criteria to samples with TF ≤ 10% instead
                date_lowtftimepoints = list(
                    tf_patient[tf_patient['tumor_burden'] == min(tf_patient['tumor_burden'])]['date'].unique())
        # get mutations using Saranya's variant calling pipeline (VarDict based) from deep 226-gene panel targeted seq
        # annotated as 'Trusted'
        mutations_targeted_df = pd.read_excel(os.path.join(*config.mutationfolder, 'CCG_226_' + str(patient) + '_reGeno.VEP.readable_tiers.xls'))
        col = ['#CHROM', 'POS', 'REF', 'ALT', 'GENE', 'TIERS']
        if mutations_targeted_df[mutations_targeted_df["TIERS"] == 'Trusted'].shape[0] != 0:  # Trusted mutations
            mutations_targeted_df = mutations_targeted_df[mutations_targeted_df["TIERS"] == 'Trusted']
            mutations_targeted_df.drop_duplicates(['#CHROM', 'POS'], inplace=True)
        else:  # if not Trusted Mutation, relax criteria to Low Evidence mutations
            mutations_targeted_df = mutations_targeted_df[mutations_targeted_df["TIERS"] == 'LowEvidence']
        if mutations_targeted_df.shape[0] > 0:  # if Trusted or Low Evidence mutations are detected in targeted seq
            for c in list(mutations_targeted_df.columns[6:]):
                if c.startswith('CCG_226_' + str(patient)):
                    col.append(c)
            mutations_targeted_df = mutations_targeted_df[col]
            mutations_targeted_df.insert(loc=6,
                                   column='helper',
                                   value='hello')
            try:
                mutations_targeted_acrosstime = (mutations_targeted_df.set_index(col[:6] + ['helper'])
                                            .stack()
                                            .unstack(-2)
                                            .ffill(axis=1)
                                            .bfill(axis=1, downcast='infer')
                                            .add_prefix('new_')
                                            .reset_index()
                                            .rename({'level_6': 'date'}, axis=1))
            except:
                mutations_targeted_acrosstime = (mutations_targeted_df.set_index(col[:6] + ['helper'])
                                            .stack()
                                            .drop_duplicates()
                                            .unstack(-2)
                                            .ffill(axis=1)
                                            .bfill(axis=1, downcast='infer')
                                            .add_prefix('new_')
                                            .reset_index()
                                            .rename({'level_6': 'date'}, axis=1))
            mutations_targeted_acrosstime['date'] = mutations_targeted_acrosstime['date'].str.split('.').str[1]
            mutations_targeted_acrosstime['date'] = pd.to_datetime(mutations_targeted_acrosstime['date'], format='%d%m%y').astype(str)
            foo2 = lambda x: pd.Series(float(x.split(' / ')[0]) / float(x.split(' / ')[1].split(' = ')[0]))
            mutations_targeted_acrosstime['VAF'] = mutations_targeted_acrosstime['new_hello'].apply(foo2)
            mutations_targeted_acrosstime.drop('new_hello', axis=1, inplace=True)
            mutations_targeted_acrosstime = mutations_targeted_acrosstime.pivot_table(values='VAF', index='GENE', columns='date', aggfunc='first')
            mutations_targeted_acrosstime = mutations_targeted_acrosstime.T
            mutations_acrosstime = mutations_targeted_acrosstime
            print(mutations_acrosstime)
            print(list(mutations_acrosstime.columns))
            if patient == 1014:  # add 101 panel mutation calls by Sarah timepoint 2016-08-18 for patient 1014
                add = pd.read_excel(os.path.join(*config.mutationfolder, 'CHB1461_CHC1126_variantcalls.xlsx'))
                add = add[['Symbol', 'freebayes_VAFTumour', 'mutect_VAFTumour', 'varscan_VAFTumour']].iloc[:83,:]
                add = add[add['Symbol'].isin(mutations_acrosstime.columns)]
                add['VAF'] = add[['freebayes_VAFTumour', 'mutect_VAFTumour', 'varscan_VAFTumour']].median(axis=1, skipna=True)
                add.set_index('Symbol', inplace=True)
                add = add.reindex(list(mutations_acrosstime.columns))
                print(add)
                mutations_acrosstime.loc['2016-08-18'] = add['VAF'].values
            xacrosstime = [i for i in mutations_acrosstime.index if i in df_patient['date'].values]
            eles = [ele0]
            collist = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
            if mutations_targeted_df["TIERS"].iloc[0] == 'Trusted':
                lstype = '-'
            elif mutations_targeted_df["TIERS"].iloc[0] == 'LowEvidence':
                lstype = '--'
            for gi, gene in enumerate(mutations_acrosstime.columns):
                yacrosstime = [m for i, m in enumerate(mutations_acrosstime[gene].values) if mutations_acrosstime.index[i] in df_patient['date'].values]
                if gene in config.genelist:
                    c = sns.color_palette('tab10') + [sns.color_palette('Accent')[5]] + [sns.color_palette('tab20b')[14]] + [sns.color_palette('tab20b')[0]]
                    c = c[config.genelist.index(gene)]
                else:
                    c = collist[gi % len(collist)]
                elei = ax.plot(xacrosstime, yacrosstime, color=c, ls=lstype, linewidth=lwd, label=gene)
                ax.plot(xacrosstime, yacrosstime, color=c, ls=lstype, marker='^', markersize=10)
                eles.append(elei)
            # indicate dates when tumor burden is available, that is to say timepoints with lpWGS
            labels = [ad if ad in tumorburden_dates else '' for ad in alldates]
            ax2.set_xticklabels(labels, rotation=90, fontsize=fs)
            ax.set_xticklabels(labels, rotation=90, fontsize=fs)
            # timepoints with high tumor burden deepWGS samples in red
            tf_deepwgs_df = get_tf(config, batch='deepwgs')
            listdeepwgs = list(pd.to_datetime(tf_deepwgs_df[tf_deepwgs_df['patient'] == patient]['date'],
                                              format='%d%m%y').astype(str).values)
            print(listdeepwgs)
            # timepoints with low tumor burden lpWGS samples in blue (ichorCNA tumor fraction estimate)
            if highlight == 'discovery':
                for dltbt in date_lowtftimepoints:
                    ax.get_xticklabels()[labels.index(dltbt)].set_color('blue')
                    ax2.get_xticklabels()[labels.index(dltbt)].set_color('blue')
                for ldw in listdeepwgs:
                    ax.get_xticklabels()[labels.index(ldw)].set_color('orange')
                    ax2.get_xticklabels()[labels.index(ldw)].set_color('orange')
            elif highlight == 'validation':
                ax.get_xticklabels()[labels.index(config.lowtfsamples[patient])].set_color('blue')
                ax2.get_xticklabels()[labels.index(config.lowtfsamples[patient])].set_color('blue')
                for htfsample in config.hightfsamples[patient]:
                    ax.get_xticklabels()[labels.index(htfsample)].set_color('orange')
                    ax2.get_xticklabels()[labels.index(htfsample)].set_color('orange')
            else:
                raise ValueError("highlight should be 'discovery' or 'validation' but here is {}".format(highlight))
            if treatment:
                ax.legend(loc='upper left')
            else:
                ax.legend(bbox_to_anchor=(1, 1),  loc='upper left')
                plt.ylim([-0.01, 0.8])
    options = ''
    if treatment:
        options += '_treatment'
    if mutations:
        options += '_mutations'
    if save and savepath is None:
        plt.savefig(os.path.join(os.getcwd(), *config.outputpath, 'timeline_patient', 'timeline_patient_' + str(patient) + options + '.png'), bbox_inches='tight')
        plt.savefig(os.path.join(os.getcwd(), *config.outputpath, 'timeline_patient', 'timeline_patient_' + str(patient) + options + '.svg'), bbox_inches='tight')
    if save and savepath is not None:
        plt.savefig(os.path.join(savepath, 'timeline_patient_' + str(patient) + options + '.png'), bbox_inches='tight')
        plt.savefig(os.path.join(savepath, 'timeline_patient_' + str(patient) + options + '.svg'), bbox_inches='tight')
    plt.show()


def get_mutations_stats(config, patient, mutation_folder=None, prefix='CCG_226_', suffix='_reGeno.VEP.readable_tiers'):
    """Gets the ichorCNA tumor fractions of one or all batches of longitudinal cfDNA samples.

    Parameters
    ----------
    config : object
        configuration class instance containing paths and batch names
    patient : str or int
        patient number ID
    mutation_folder : str or None
        folder containing excel files of mutation calls. If None, uses location indicated in config. (default is None).
    prefix : str, optional
        prefix in mutation file excel nomenclature. (default is 'CCG_226')
    suffix : str, optional
        suffix in mutation file excel nomenclature without 'xls' extension. (default is '_reGeno.VEP.readable_tiers').

    Returns
    -------
    pandas DataFrame instance
        a dataframe containing for each cfDNA timepoint samples of a patient
        the nb of mutations (annotated as Trusted or not), the median VAF over the 226-panel of genes and
        of mutations (annotated as Trusted or not).
    """

    if mutation_folder is None:
        mutation_folder = os.path.join(os.getcwd(), *config.mutationfolder)
    mutation_df = pd.read_excel(os.path.join(mutation_folder, prefix + str(patient) + suffix + '.xls'))
    col = ['#CHROM', 'POS', 'REF', 'ALT', 'GENE', 'TIERS']
    targeted_lowtf = []
    # Patient treatment + tumor fraction infos
    df_patient, tf_patient, _ = get_patient_tf_treatment(config, patient)
    # get timepoints with low tumor burden estimate (ichorCNA)
    date_lowtftimepoints = list(tf_patient[tf_patient['tumor_burden'] == 0]['date'].unique())
    if not date_lowtftimepoints:
        print('no zero ichorCNA estimate tumor burden')
        print('min tumor burden is {}'.format(min(tf_patient['tumor_burden'])))
        if min(tf_patient['tumor_burden']) < 0.1:
            date_lowtftimepoints = list(
                tf_patient[tf_patient['tumor_burden'] == min(tf_patient['tumor_burden'])]['date'].unique())
    for it in list(date_lowtftimepoints):
        aux1 = 'CCG_226_' + str(patient) + '.' + it.split('-')[-1] + it.split('-')[1] + it.split('-')[0][-2:] + '.P'
        aux2 = 'CCG_226_' + str(patient) + '.' + it.split('-')[-1] + it.split('-')[1] + it.split('-')[0][-2:]
        if sum(mutation_df.columns.str.contains(aux1)) == 1 or sum(mutation_df.columns.str.contains(aux2)) == 1:
            if sum(mutation_df.columns.str.contains(aux1)) == 1:
                idx = mutation_df.columns.str.contains(aux1).tolist().index(True)
            else:  # sum(mutation_df.columns.str.contains(aux2)) == 1:
                idx = mutation_df.columns.str.contains(aux2).tolist().index(True)
            col.append(mutation_df.columns[idx])
            targeted_lowtf.append(
                str(pd.to_datetime(mutation_df.columns[idx].split('.')[1], format='%d%m%y')).split(' ')[0])
    mutation_df = mutation_df[col]
    mutation_df.insert(loc=6, column='helper', value='hello')
    mutation_lowtftimepoints_226 = (mutation_df.set_index(col[:6] + ['helper'])
                                    .stack()
                                    .unstack(-2)
                                    .ffill(axis=1)
                                    .bfill(axis=1, downcast='infer')
                                    .add_prefix('new_')
                                    .reset_index()
                                    .rename({'level_6': 'date'}, axis=1))
    mutation_lowtftimepoints_226['date'] = mutation_lowtftimepoints_226['date'].str.split('.').str[1]
    mutation_lowtftimepoints_226['date'] = pd.to_datetime(mutation_lowtftimepoints_226['date'],
                                                          format='%d%m%y').astype(str)
    helfperfunction = lambda x: pd.Series(float(x.split(' / ')[0]) / float(x.split(' / ')[1].split(' = ')[0]))
    mutation_lowtftimepoints_226['VAF'] = mutation_lowtftimepoints_226['new_hello'].apply(helfperfunction)
    mutation_lowtftimepoints_226.drop('new_hello', axis=1, inplace=True)
    mutation_lowtftimepoints_226 = mutation_lowtftimepoints_226[mutation_lowtftimepoints_226['TIERS'] == 'Trusted']
    mutation_lowtftimepoints = mutation_lowtftimepoints_226
    lowtftimepoints_dict = {'date': date_lowtftimepoints,
                            'median VAF': [],
                            '# mutated genes': [],
                            'median VAF within mutated genes': [],
                            '# mutated genes TRUSTED': [],
                            'median VAF within mutated genes TRUSTED': []
                            }
    for date in date_lowtftimepoints:
        nmut = mutation_lowtftimepoints[
            (mutation_lowtftimepoints['date'] == date) & (mutation_lowtftimepoints['VAF'] != 0)].shape[0]
        nmuttrust = mutation_lowtftimepoints[
            (mutation_lowtftimepoints['date'] == date) & (mutation_lowtftimepoints['VAF'] != 0) & (
                    mutation_lowtftimepoints['TIERS'] == 'Trusted')].shape[0]
        medianvaf = np.median(mutation_lowtftimepoints[(mutation_lowtftimepoints['date'] == date)]['VAF'].values)
        medianvafn = np.median(mutation_lowtftimepoints[(mutation_lowtftimepoints['date'] == date) & (
                mutation_lowtftimepoints['VAF'] != 0)]['VAF'].values)
        medianvafntrust = np.median(mutation_lowtftimepoints[(mutation_lowtftimepoints['date'] == date) & (
                mutation_lowtftimepoints['VAF'] != 0) & (mutation_lowtftimepoints['TIERS'] == 'Trusted')][
                                        'VAF'].values)
        lowtftimepoints_dict['# mutated genes'].append(nmut)
        lowtftimepoints_dict['# mutated genes TRUSTED'].append(nmuttrust)
        lowtftimepoints_dict['median VAF'].append(medianvaf)
        lowtftimepoints_dict['median VAF within mutated genes'].append(medianvafn)
        lowtftimepoints_dict['median VAF within mutated genes TRUSTED'].append(medianvafntrust)
    lowtftimepoints_pd = pd.DataFrame.from_dict(lowtftimepoints_dict)
    lowtftimepoints_pd.set_index('date', inplace=True)
    return lowtftimepoints_pd


if __name__ == "__main__":
    import argparse
    from utils.config import Config
    from utils.viz import set_display_params
    # set working directory
    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    # # set config
    # config = Config("config/", "config_viz.yaml")
    # # example patient 986
    # patient = 986
    # # set display params
    # set_display_params(config)
    # # test functions
    # tf_df = get_tf(config, batch='all')
    # print(tf_df)
    # plot_patient_timeline(config, patient)
    # plot_patient_timeline(config, patient, treatment=True, mutations=True)
    # lowtftimepoints_pd = get_mutations_stats(config, patient)
    # print(lowtftimepoints_pd.dropna())

    parser = argparse.ArgumentParser()
    parser.add_argument('--config_path', default='config_viz.yaml', required=False)
    parser.add_argument('--patient', default=986, required=False)
    # parser.add_argument('--mutation_folder', default='None')
    # parser.add_argument('--prefix', default='CCG_226_')
    # parser.add_argument('--suffix', default='_reGeno.VEP.readable_tiers')
    args = parser.parse_args()
    print(args)

    # set config
    config = Config("config/", args.config_path)
    patient = args.patient
    # set display params
    set_display_params(config)
    tf_df = get_tf(config, batch='all')
    print(tf_df)
    plot_patient_timeline(config, patient)
    plot_patient_timeline(config, patient, treatment=True, mutations=True)
    mutations_stats_df = get_mutations_stats(config, patient)
    print(mutations_stats_df.dropna())





