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