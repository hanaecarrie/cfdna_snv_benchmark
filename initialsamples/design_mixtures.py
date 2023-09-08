# design_mixtures.py

import os
import numpy as np
import pandas as pd

tfestimates_dict = {
    986: {'lowtf': ['300316', 0.15, [0.1, 0.2]],
          'hightf1': ['100215', 32.97, [27.10, 38.84]],
          'hightf2': ['261016', 30.48, [24.59, 36.37]],
          },
    1014: {'lowtf': ['090616', 2.85, [1.9, 3.8]],
          'hightf1': ['100816', 42.83, [40.65, 45.01]],
          'hightf2': ['100116', 33.00, [20.95, 45.05]],
          },
    123: {'lowtf': ['121115', 4.1, [2.73, 4.1]],
          'hightf1': ['310715', 60.95, [56.5, 65.4]],
          'hightf2': ['070116', 29.4, [27.8, 31.0]],
          }
}

for patient, patientdict in tfestimates_dict.items():
    print('######### ' + str(patient) + ' #########')
    lowtfsample, approxtflow, citflow = patientdict['lowtf']
    for hightf in ['hightf1', 'hightf2']:
        hightfsample, approxtfhigh, citfhigh = patientdict[hightf]
        print(citfhigh)
        print(str(patient) + '_' + hightfsample + '-' + str(patient) + '_' + lowtfsample)
        mixtures_df = pd.DataFrame()
        mixtures_df['cov high TF'] = [70, 70, 50, 30, 20, 10, 5, 70]
        mixtures_df['cov ultra low TF'] = [0, 80, 100, 120, 130, 140, 145, 230]
        mixtures_df['cov mixture'] = mixtures_df['cov high TF'] + mixtures_df['cov ultra low TF']
        mixtures_df['approx TF'] = (approxtfhigh * mixtures_df['cov high TF'] + approxtflow * mixtures_df['cov ultra low TF'])/mixtures_df['cov mixture']
        mixtures_df['approx TF'] = mixtures_df['approx TF'].round(1).astype(str) + '%'
        mixtures_df['min CI TF'] = (citfhigh[0] * mixtures_df['cov high TF'] + citflow[0] * mixtures_df['cov ultra low TF'])/mixtures_df['cov mixture']
        mixtures_df['max CI TF'] = (citfhigh[1] * mixtures_df['cov high TF'] + citflow[1] * mixtures_df['cov ultra low TF'])/mixtures_df['cov mixture']
        mixtures_df['CI mixture TF'] = mixtures_df['min CI TF'].round(1).astype(str)+ "-" + mixtures_df['max CI TF'].round(1).astype(str) + "%"
        mixtures_df['cov high TF'] = mixtures_df['cov high TF'].astype(str) + "x"
        mixtures_df['cov ultra low TF'] = mixtures_df['cov ultra low TF'].astype(str) + "x"
        mixtures_df['cov mixture'] = mixtures_df['cov mixture'].astype(str) + "%"
        mixtures_df.drop(['min CI TF', 'max CI TF'], axis=1, inplace=True)
        print(mixtures_df.head(10))

### Older code #TODO refrac

if not os.getcwd().endswith('cfdna_snv_benchmark'):
    os.chdir('../')
print('Current working directory: {}'.format(os.getcwd()))

tf_estimates_filepath = os.path.join('data', 'tumor_burden', 'ctDNA_burden_estimation_deepwgs_plasma_CRC_samples'
                                                             '.xlsx')
df_tf = pd.read_excel(tf_estimates_filepath, keep_default_na=False).set_index('sample name')
print(df_tf)
#df_tf.loc['1014_110116'] = ['1014', '1014_110116', 'CRC', 63.17, 0.14, 0.75, 0.50, 0.84, 0.14, 0.75, 0.50, 0.84, 0.14, 0.73, 0.50, 0.89, 0.56, 0.62, 0.10, np.nan, np.nan, np.nan]
# edit average coverage value with source bcbio Aquila
# bcbio_final/2015-07-31_XXX/mutiqc/report/metrics/XXX-T.bcbio.txt Avg_coverage
#df_tf.loc['986_100215'][3] = 77.68
#df_tf.loc['986_261016'][3] = 82.82
#df_tf.loc['1014_180816'][3] = 87.80
# studied samples
samples = ['986_100215', '986_261016', '1014_180816']  # '1014_110116',
# proposed ratios
mixseries = [['cov', 0], [70, 0], [70, 80], [50, 100], [30, 120], [20, 130], [10, 140], [5, 145], [70, 130],
             [70, 180], [70, 230]]

for sample in samples:
    cov_hightb = df_tf.loc[sample]['coverage']
    tf_tissuemethods_hightb = list(df_tf.loc[sample][4:16].values)  # hard coded column IDs for DNA methods
    tf_ichorcna_hightb = df_tf.loc[sample][18]  # hard coded column IDs for ichorcna
    print(cov_hightb, tf_tissuemethods_hightb, tf_ichorcna_hightb)
    # give same weight ichorCNA estimate VS DNA methods
    tf_list_hightb = tf_tissuemethods_hightb + len(tf_tissuemethods_hightb) * [tf_ichorcna_hightb]
    tf_list_hightb = np.multiply(tf_list_hightb, 100)
    print(tf_list_hightb)

    # 95% confidence interval
    cihigh = [np.mean(tf_list_hightb) - 1.96*np.std(tf_list_hightb)/np.sqrt(len(tf_list_hightb)),
              np.mean(tf_list_hightb) + 1.96*np.std(tf_list_hightb)/np.sqrt(len(tf_list_hightb))]
    if sample.split('_')[0] == '986':
        cilow = [0.03, 8.05]  # vaf in targeted VS ichorCNA estimate from deep WGS
    elif sample.split('_')[0] == '1014':
        cilow = [0.03, 5.07]  # vaf in targeted VS ichorCNA estimate from deep WGS
    else:
        raise ValueError('patient {} is not studied with low matched tumor burden'.format(sample.split('_')[0]))

    tfhigh = np.mean(tf_list_hightb)
    tflow = np.mean(cilow)

    res_df = pd.DataFrame(columns=['mean tf', 'ci tf min', 'mci tf max'])
    for ms in mixseries:
        covhigh, covlow = ms
        if covhigh == 'cov':
            covhigh = float(cov_hightb)
        tmin = (cihigh[0]*covhigh + cilow[0]*covlow)/(covhigh+covlow)
        tmax = (cihigh[1]*covhigh + cilow[1]*covlow)/(covhigh+covlow)
        tmean = (tfhigh*covhigh + tflow*covlow)/(covhigh+covlow)
        res_df.loc[str([covhigh, covlow])] = [np.round(tmean, 2), np.round(tmin, 2), np.round(tmax, 2)]
    print('#####')
    print(sample)
    print('#####')
    print(cov_hightb)
    print(np.round(tflow, 2), np.round(tfhigh, 2))
    print(res_df)
