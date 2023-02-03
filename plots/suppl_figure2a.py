
import os
import warnings
warnings.filterwarnings('ignore')

# set working directory
if not os.getcwd().endswith('cfdna_snv_benchmark'):
    os.chdir('../')
print('Current working directory: {}'.format(os.getcwd()))

from utils.config import Config
from utils.viz import *
from utils.metricsseries import *
from utils.calltableseries import *
from utils.groundtruth import *

# Config and Display paramaters
config = Config("config/", "config_viz.yaml")
set_display_params(config)
print(config.methods)

# Chomosome
mixtureids = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T']
reload = False
save = True
fixedvars = ['coverage', 'ctdna']
filterparam = 'all'
markers = ['o', '^', 'X']
linestyles = ['-', '-', '-']
color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}

muttypes = ['snv', 'indel']
metrics = ['auprc', 'precision', 'recall']

chrom = 'all'

fixedvar = 'coverage'


if fixedvar == 'coverage':
    seriesorder = [(70, 0), (70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]
    xaxis = 'tumor burden'
elif fixedvar == 'ctdna':
    seriesorder = [(70, 0), (70, 80), (70, 180)]
    xaxis = 'coverage'

# for mixtureid in mixtureids:
mixtureid = 'CRC-123_310715-CW-T_CRC-123_121115-CW-T'
mixtureid = 'CRC-986_100215-CW-T_CRC-986_300316-CW-T'

print('############# {} ############'.format(mixtureid))
if mixtureid == 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T':
    chroms = [str(c) for c in range(1, 23) if c != 2 and c != 6 and c != 17 and c != 19 and c != 20 and c != 21]
    chroms = [str(c) for c in range(1, 9) if c != 2 and c != 6]
elif mixtureid == 'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
    chroms = [str(c) for c in range(1, 23) if c != 1 and c != 2 and c != 8 and c != 20 and c != 21 and c != 22]
else:
    chroms = [str(c) for c in range(1, 23) if c != 6 and c != 19 and c != 20]  # c !=1 and c!= 2 and
calltables = {'sampleid': [], 'tf': [], 'cov': [], 'snv': [], 'indel': [], 'snp': []}
aux_all = []
calltable_snv, aux = get_calltableseries(config, mixtureid, chroms, muttype='snv', filterparam=filterparam, reload=reload, save=save)
calltable_indel, aux = get_calltableseries(config, mixtureid, chroms, muttype='indel', filterparam=filterparam, reload=reload, save=save)
calltable_snp, aux = get_calltableseries(config, mixtureid, chroms, muttype='snp', filterparam=filterparam, reload=reload, save=save)
print(calltable_snv.shape, calltable_indel.shape, calltable_snp.shape)
print(aux)
plasmasample = '_'.join(mixtureid.split('_')[:2])
print(plasmasample)
healthysample = '_'.join(mixtureid.split('_')[2:])
print(healthysample)
calltables['snv'].append(calltable_snv)
calltables['indel'].append(calltable_indel)
calltables['snp'].append(calltable_snp)
calltables['sampleid'] = mixtureid
print(list(calltable_snv.columns))
calltables['tf'] = list(np.unique([cn.split('_')[0] if '_' in cn else '0' for cn in list(calltable_snv.columns)]).astype(float))
calltables['tf'].remove(0)
print(calltables['tf'])
calltables['snv'] = pd.concat(calltables['snv'])
calltables['indel'] = pd.concat(calltables['indel'])
calltables['snp'] = pd.concat(calltables['snp'])
dilutionseries = aux.T[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]].T

# for muttype in muttypes:
muttype = 'snv'
refsample = 'undiluted'
if muttype == 'snv':
    gtm = 4
else:  # elif muttype == 'indel':
    gtm = 2
print(max(aux['tf']))
if mixtureid ==  'CRC-986_100215-CW-T_CRC-986_300316-CW-T' and muttype == 'snv':
    gtm = 3
    refsample = 'tissue'
    calltablesseries  = generate_groundtruth(config, calltables[muttype], aux['tf'], ground_truth_method='tissue', muttype=muttype,
                                             matchedtissuepath=os.path.join('data', 'matchedtissue_ultradeep', '986_100215_T1-E', 'calls', '986_100215_T1-E_snv_calls_all.csv'))
    ### Missing 986 indels ###
else:
    calltablesseries = generate_groundtruth(config, calltables[muttype], aux['tf'], ground_truth_method=gtm, muttype=muttype,
                                            matchedtissuepath=None, methods=['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan', 'abemus', 'sinvict'])
if not os.path.exists(os.path.join(*config.outputpath, 'supplfigure2a')):
    os.mkdir(os.path.join(*config.outputpath, 'supplfigure2a'))
figure_curve_allchr(config, calltablesseries, dilutionseries, mixtureid, xy='pr', ground_truth_method=gtm,
                    refsample=refsample, muttype=muttype.upper(), methods=None, fixedvar=fixedvar, save=config.outputpath + ['supplfigure2a'])
