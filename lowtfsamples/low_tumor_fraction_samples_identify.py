# Imports
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

from utils.config import Config
from utils.viz import *


def get_tf(config, batch='all'):
    if batch == 'deepwgs':
        tf_df = pd.read_csv(config.tumorfractionfiles.deepwgs, names=['sample_id', 'tumor_burden'])
    elif batch == 'batch1':
        tf_df = pd.read_csv(config.tumorfractionfiles.batch1, names=['sample_id', 'tumor_burden'])
    elif batch == 'batch2':
        tf_df = pd.read_csv(config.tumorfractionfiles.batch2, names=['sample_id', 'tumor_burden'])
    elif batch == 'all':
        tf_df = pd.concat(load_files([config.tumorfractionfiles[batch] for batch in config.batches]))
    else:
        raise ValueError("unknown batch {}. batch should be either 'batch1', 'batch2', 'deepwgs or 'all'".format(batch))
    tf_df['patient'] = [int(tf_df['sample_id'].iloc[i].split('-')[1].split('_')[0])
                        if tf_df['sample_id'].iloc[i].startswith('sortpanel')
                        else int(tf_df['sample_id'].iloc[i].split('_')[0]) for i in range(tf_df.shape[0])]
    tf_df['date'] = pd.to_datetime(tf_df['sample_id'].str.split('_').str[1], format='%d%m%y')
    return tf_df


def get_patient_tf_treatment(config, patient):
    # Tumor burden
    tf_df = get_tf(config, batch='all')
    tf_df['patient'] = [int(tf_df['sample_id'].iloc[i].split('-')[1].split('_')[0])
                        if tf_df['sample_id'].iloc[i].startswith('sortpanel')
                        else int(tf_df['sample_id'].iloc[i].split('_')[0]) for i in range(tf_df.shape[0])]
    tf_df['date'] = pd.to_datetime(tf_df['sample_id'].str.split('_').str[1], format='%d%m%y')
    tf_patient = tf_df[tf_df['patient'] == patient].sort_values('date')
    tf_patient['date'] = tf_patient['date'].astype(str)

    tumorburden_dates = tf_df[tf_df['patient'] == patient]['date'].sort_values().astype(str).unique()
    print(tumorburden_dates)
    print(tf_df[tf_df['patient'] == patient].shape[0])

    # Treatment file
    treatment_df = pd.read_csv(config.treatmentfile, sep='\t')
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


def plot_patient_timeline(config, patient, mutations=False):
    # Display paramaters
    #context = config.context  # talk or paper
    #if context == 'talk':
    #    sns.set(style="darkgrid", context="talk", rc={"lines.linewidth": 2, "legend.fontsize": 25})
    #    plt.style.use("dark_background")
    #    plt.rcParams.update({"grid.linewidth": 0.5, "grid.alpha": 0.5, 'font.size': 30})
    #    sns.set_palette("deep")
    #elif context == 'paper':
    #    sns.set(context="paper", font_scale=4, rc={"lines.linewidth": 2, "legend.fontsize": 25})
    #    plt.rcParams.update({'font.size': 35})
    #    sns.set_palette("deep")
    #else:
    #    raise ValueError("unknown context {}. Should be either 'talk' or 'paper'".format(config.context))
    lc = 'w' if config.context == 'talk' else 'k'

    # Patient treatment + tumor fraction infos
    df_patient, tf_patient, _ = get_patient_tf_treatment(config, patient)
    alldates = sorted(list(set(list(df_patient['date'].values) + list(tf_patient['date'].values))))
    print('alldates')
    print(len(alldates))
    print(alldates)
    tumorburden_dates = tf_patient[tf_patient['patient'] == patient]['date'].sort_values().astype(str).unique()

    # Plot patient timeline
    plt.figure()
    fig, ax2 = plt.subplots(figsize=(40, 10))
    # make a plot with different y-axis using second axis object
    sns.stripplot(y='treatment', x='date', hue='treatment', data=df_patient, s=10, ax=ax2)
    ax2.grid(False)
    plt.legend((), ())
    # twin object for two different y-axis on the sample plot
    ax = ax2.twinx()
    ele0 = ax.plot(df_patient['date'], df_patient['tumor_burden'], lc + '.', marker='s', markersize=10,
                   label='ichorCNA tumor burden')
    ax.plot(df_patient['date'][~df_patient['tumor_burden'].isna()],
            df_patient['tumor_burden'][~df_patient['tumor_burden'].isna()], lc + '-', linewidth=4)
    ax.plot(df_patient['date'][~df_patient['tumor_burden'].isna()],
            df_patient['tumor_burden'][~df_patient['tumor_burden'].isna()], lc + '.', marker='s', markersize=10)
    ax.set_ylabel('fraction (tumor burden or VAF)', fontsize=30)
    ax.set_ylim(-0.01, 1)
    fig.legend([ele0], ['ichorCNA tumor burden'], loc='upper left')
    ax.legend(loc='upper left')
    labels = [ad if ad in tumorburden_dates else '' for ad in alldates]
    ax2.set_xticklabels(labels, rotation=90) #fontsize=30)
    ax.set_xticklabels(labels, rotation=90) # fontsize=30)
    plt.title('Patient {}: VAF evolution across timepoints'.format(patient))
    if not os.path.exists(os.path.join(*config.outputpath, 'timeline_patient')):
        os.mkdir(os.path.join(*config.outputpath, 'timeline_patient'))
    if not os.path.exists(
            os.path.join(*config.outputpath, 'timeline_patient', 'timeline_patient_' + str(patient) + '.png')):
        plt.savefig(os.path.join(*config.outputpath, 'timeline_patient', 'timeline_patient_' + str(patient) + '.png'),
                    bbox_inches='tight')

    if mutations:
        date_lowtftimepoints = list(tf_patient[tf_patient['tumor_burden'] == 0]['date'].unique())
        if date_lowtftimepoints == []:
            print('no zero ichorCNA estimate tumor burden')
            print('min tumor burden is {}'.format(min(tf_patient['tumor_burden'])))
            if min(tf_patient['tumor_burden']) < 0.1:
                date_lowtftimepoints = list(
                    tf_patient[tf_patient['tumor_burden'] == min(tf_patient['tumor_burden'])]['date'].unique())
        print(date_lowtftimepoints)
        # 226 - panel
        mutation_df = pd.read_excel('../data/variant_calls/226 PANEL VARIANTS CLASSIFICSATION EXCEL/CCG_226_'+str(patient)+'_reGeno.VEP.readable_tiers.xls')
        print(mutation_df.columns)

        col = ['#CHROM', 'POS', 'REF', 'ALT', 'GENE', 'TIERS']
        targeted_lowtf = []
        for it in list(date_lowtftimepoints):
            aux1 = 'CCG_226_'+str(patient)+'.'+ it.split('-')[-1]+it.split('-')[1]+it.split('-')[0][-2:]+'.P'
            aux2 = 'CCG_226_'+str(patient)+'.'+ it.split('-')[-1]+it.split('-')[1]+it.split('-')[0][-2:]
            if sum(mutation_df.columns.str.contains(aux1)) == 1 or sum(mutation_df.columns.str.contains(aux2)) == 1:
                if sum(mutation_df.columns.str.contains(aux1)) == 1:
                    idx = mutation_df.columns.str.contains(aux1).tolist().index(True)
                elif sum(mutation_df.columns.str.contains(aux2)) == 1:
                    idx = mutation_df.columns.str.contains(aux2).tolist().index(True)
                col.append(mutation_df.columns[idx])
                targeted_lowtf.append(str(pd.to_datetime(mutation_df.columns[idx].split('.')[1], format='%d%m%y')).split(' ')[0])
        print(col)
        mutation_df = mutation_df[col]
        mutation_df.insert(loc = 6,
                           column = 'helper',
                           value = 'hello')
        mutation_lowtftimepoints_226 = (mutation_df.set_index(col[:6]+['helper'])
                                        .stack()
                                        .unstack(-2)
                                        .ffill(axis=1)
                                        .bfill(axis=1, downcast='infer')
                                        .add_prefix('new_')
                                        .reset_index()
                                        .rename({'level_6': 'date'}, axis=1))

        mutation_lowtftimepoints_226['date'] = mutation_lowtftimepoints_226['date'].str.split('.').str[1]
        mutation_lowtftimepoints_226['date'] = pd.to_datetime(mutation_lowtftimepoints_226['date'], format='%d%m%y').astype(str)
        foo2 = lambda x: pd.Series(float(x.split(' / ')[0])/float(x.split(' / ')[1].split(' = ')[0]))
        mutation_lowtftimepoints_226['VAF'] = mutation_lowtftimepoints_226['new_hello'].apply(foo2)
        mutation_lowtftimepoints_226.drop('new_hello', axis=1, inplace=True)
        mutation_lowtftimepoints_226 = mutation_lowtftimepoints_226[mutation_lowtftimepoints_226['TIERS'] == 'Trusted']
        print(mutation_lowtftimepoints_226.shape)
        print(targeted_lowtf)

        mutation_lowtftimepoints = mutation_lowtftimepoints_226
        #mutation_lowtftimepoints = pd.concat([mutation_lowtftimepoints_226, mutation_lowtftimepoints_MCP])
        print(mutation_lowtftimepoints.shape)
        print(mutation_lowtftimepoints['date'].unique())

        lowtftimepoints_dict = {'date': date_lowtftimepoints,
                                'median VAF': [],
                                '# mutated genes' : [],
                                'median VAF within mutated genes': [],
                                '# mutated genes TRUSTED': [],
                                'median VAF within mutated genes TRUSTED': []
                                }

        for date in date_lowtftimepoints:
            nmut = mutation_lowtftimepoints[(mutation_lowtftimepoints['date'] == date) & (mutation_lowtftimepoints['VAF'] != 0)].shape[0]
            nmuttrust = mutation_lowtftimepoints[(mutation_lowtftimepoints['date'] == date)  & (mutation_lowtftimepoints['VAF'] != 0)  & (mutation_lowtftimepoints['TIERS'] == 'Trusted')].shape[0]
            medianvaf = np.median(mutation_lowtftimepoints[(mutation_lowtftimepoints['date'] == date)]['VAF'].values)
            medianvafn = np.median(mutation_lowtftimepoints[(mutation_lowtftimepoints['date'] == date) & (mutation_lowtftimepoints['VAF'] != 0)]['VAF'].values)
            medianvafntrust = np.median(mutation_lowtftimepoints[(mutation_lowtftimepoints['date'] == date) & (mutation_lowtftimepoints['VAF'] != 0) & (mutation_lowtftimepoints['TIERS'] == 'Trusted')]['VAF'].values)
            lowtftimepoints_dict['# mutated genes'].append(nmut)
            lowtftimepoints_dict['# mutated genes TRUSTED'].append(nmuttrust)
            lowtftimepoints_dict['median VAF'].append(medianvaf)
            lowtftimepoints_dict['median VAF within mutated genes'].append(medianvafn)
            lowtftimepoints_dict['median VAF within mutated genes TRUSTED'].append(medianvafntrust)

        lowtftimepoints_pd = pd.DataFrame.from_dict(lowtftimepoints_dict)
        lowtftimepoints_pd.set_index('date', inplace=True)

        mutation_df_226 = pd.read_excel('../data/variant_calls/226 PANEL VARIANTS CLASSIFICSATION EXCEL/CCG_226_'+str(patient)+'_reGeno.VEP.readable_tiers.xls')

        col = ['#CHROM', 'POS', 'REF', 'ALT', 'GENE', 'TIERS']
        if mutation_df_226[mutation_df_226["TIERS"] == 'Trusted'].shape[0] != 0:
            mutation_df_226 = mutation_df_226[mutation_df_226["TIERS"] == 'Trusted']
        else:
            mutation_df_226 = mutation_df_226[mutation_df_226["TIERS"] == 'LowEvidence']
        print(mutation_df_226.shape)
        if mutation_df_226.shape[0] > 0:

            for c in list(mutation_df_226.columns[6:]):
                if c.startswith('CCG_226_'+str(patient)):
                    col.append(c)
            mutation_df_226 = mutation_df_226[col]
            mutation_df_226.insert(loc = 6,
                                   column = 'helper',
                                   value = 'hello')
            mutations_acrosstime_226 = (mutation_df_226.set_index(col[:6]+['helper'])
                                        .stack()
                                        .unstack(-2)
                                        .ffill(axis=1)
                                        .bfill(axis=1, downcast='infer')
                                        .add_prefix('new_')
                                        .reset_index()
                                        .rename({'level_6': 'date'}, axis=1))
            mutations_acrosstime_226['date'] = mutations_acrosstime_226['date'].str.split('.').str[1]
            mutations_acrosstime_226['date'] = pd.to_datetime(mutations_acrosstime_226['date'], format='%d%m%y').astype(str)
            foo2 = lambda x: pd.Series(float(x.split(' / ')[0])/float(x.split(' / ')[1].split(' = ')[0]))
            mutations_acrosstime_226['VAF'] = mutations_acrosstime_226['new_hello'].apply(foo2)
            mutations_acrosstime_226.drop('new_hello', axis=1, inplace=True)
            mutations_acrosstime_226 = mutations_acrosstime_226.pivot_table(values='VAF', index='GENE', columns='date', aggfunc='first')
            mutations_acrosstime_226 = mutations_acrosstime_226.T

            mutations_acrosstime = mutations_acrosstime_226

            xacrosstime = [i for i in mutations_acrosstime.index if i in df_patient['date'].values]

            eles = [ele0]
            collist = ['r', 'b', 'g', 'c', 'm', 'y', 'tab:pink', 'tab:orange', 'tab:olive', 'tab:purple', 'grey']
            if mutation_df_226["TIERS"].iloc[0] == 'Trusted':
                lstype = '-'
            elif mutation_df_226["TIERS"].iloc[0] == 'LowEvidence':
                lstype = '--'

            for gi, gene in enumerate(mutations_acrosstime.columns):
                yacrosstime = [m for i,m in enumerate(mutations_acrosstime[gene].values) if mutations_acrosstime.index[i] in df_patient['date'].values]
                elei = ax.plot(xacrosstime, yacrosstime, color=collist[gi%len(collist)], ls=lstype, linewidth=2, label=gene)
                ax.plot(xacrosstime, yacrosstime, color=collist[gi%len(collist)], ls=lstype, marker='D', markersize=10)
                eles.append(elei)

            labels = [ad if ad in tumorburden_dates else '' for ad in alldates]
            ax2.set_xticklabels(labels, rotation=90) #, fontsize=30)
            ax.set_xticklabels(labels, rotation=90) #, fontsize=30)

            print(tumorburden_dates)
            print(alldates)
            print(labels)
            print(date_lowtftimepoints)
            for dltbt in date_lowtftimepoints:
                ax.get_xticklabels()[labels.index(dltbt)].set_color('blue')
                ax2.get_xticklabels()[labels.index(dltbt)].set_color('blue')

            tf_deepwgs_df = get_tf(config, batch='deepwgs')
            listdeepwgs = list(pd.to_datetime(tf_deepwgs_df[tf_deepwgs_df['patient'] == patient]['date'], format='%d%m%y').astype(str).values)
            print(listdeepwgs)
            for ldw in listdeepwgs:
                ax.get_xticklabels()[labels.index(ldw)].set_color('red')
                ax2.get_xticklabels()[labels.index(ldw)].set_color('red')
    plt.show()


if __name__ == "__main__":
    # Config
    config = Config("../config/", "config_viz.yaml")
    patient = 986
    # tf_df = get_tf(config, batch='all')
    plot_patient_timeline(config, patient)
    plot_patient_timeline(config, patient, mutations=True)
