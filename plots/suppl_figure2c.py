# Imports

import os
import warnings

warnings.filterwarnings('ignore')

# set working directory
if not os.getcwd().endswith('cfdna_snv_benchmark'):
    os.chdir('../')
print('Current working directory: {}'.format(os.getcwd()))

from benchmark.table import *
from benchmark.metrics import *

# Config and Display paramaters

config = Config("config/", "config_viz.yaml")
set_display_params(config)

color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
alpha_dict = dict(zip(config.tissuebenchmark.fractions, [1-i*0.3 for i in range(len(config.tissuebenchmark.fractions))]))

print(color_dict)
print(alpha_dict)

df_table = get_call_table_tissue(config)
print(df_table.shape)
df_table.head()

for muttype in config.muttype:
    #for i, sample in enumerate(config.tissuebenchmark.samples):
    fig, ax = plt.subplots(figsize=(10,10))
    baseline_dict = {}
    for method in config.tissuebenchmark.methods:
        #for f in config.tissuebenchmark.fractions:
        #print(muttype, sample, f, method)
        vcf_sample = df_table[(df_table['mutation type'] == muttype)]
        df_sample_method = vcf_sample[[method+'_score', 'TRUTH']]
        df_sample_method[method + '_score'].fillna(0, inplace=True)
        precision, recall, thresholds = precision_recall_curve(df_sample_method['TRUTH'], df_sample_method[method + '_score'])
        f1 = f1_score(vcf_sample['TRUTH'], vcf_sample[method])
        estimator_name = method #if f == 1 else ''
        plot_pr_curve(precision, recall, estimator_name=estimator_name, f1_score=None, figax=(fig, ax), kwargs={'color':color_dict[method], 'alpha':1, 'lw':2})
    #baseline_dict[f] = len(vcf_sample['TRUTH'][vcf_sample['TRUTH']])/len(vcf_sample['TRUTH'])
    #plt.axhline(y=baseline_dict[f], c='k', ls='--', alpha=alpha_dict[f])
    handles, labels = plt.gca().get_legend_handles_labels()
    #list_lines = [Line2D([0], [0], color='black', alpha=alpha_dict[f], label='tumor purity = {:.2f}%'.format(round(100*config.tissuebenchmark.purities[i]*f, 2))) for f in config.tissuebenchmark.fractions]
    #legend_list = handles + list_lines + [Line2D([0], [0], color='black', ls='--', alpha=alpha_dict[f], label="baseline tf {:.2f}% = {:.2f}".format(round(100*config.tissuebenchmark.purities[i]*f, 2), baseline_dict[f])) for f in config.tissuebenchmark.fractions]
    # Creating legend with color box
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left")
    #plt.legend(bbox_to_anchor=(1, 1), loc="upper left", handles=legend_list)
    plt.title("Precision Recall curve for {} calling".format(muttype), pad=50)
    #plt.logx()
    plt.xlim([0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    if not os.path.exists(os.path.join(*config.outputpath, 'supplfigure2c')):
        os.mkdir(os.path.join(*config.outputpath, 'supplfigure2c'))
    plt.savefig(os.path.join(*config.outputpath, 'supplfigure2c', 'tissue_benchmark_'+muttype+'.svg'),
                bbox_inches='tight')