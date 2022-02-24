import os
import numpy as np
import pandas as pd

if __name__ == "__main__":

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    tf_estimates_filepath = os.path.join('data', 'tumor_burden', 'ctDNA_burden_estimation_deepwgs_plasma_CRC_samples'
                                                                 '.xlsx')
    df_tf = pd.read_excel(tf_estimates_filepath, keep_default_na=False).set_index('sample name')

    # studied samples
    samples = ['986_100215', '986_261016', '1014_180816']
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
                covhigh = cov_hightb
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
