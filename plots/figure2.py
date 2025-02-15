# Imports
import os
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# set working directory
if not os.getcwd().endswith('cfdna_snv_benchmark'):
    os.chdir('../')
print('Current working directory: {}'.format(os.getcwd()))

from utils.config import Config
from utils.viz import set_display_params
from initialsamples.patient_timeline_analysis import plot_patient_timeline, get_mutations_stats
from initialsamples.pairedplots import paireplot

# Config and Display paramaters
config = Config("config/", "config_viz.yaml")
set_display_params(config)

# order of samples = 1) high tb sample 1, 2) high tb sample 2, 3) low tb sample
patientsample_dict = {
    #'1014': ['NCC_CRC-1014_180816-CW-T', 'NCC_CRC-1014_110116-CW-T', 'NCC_CRC-1014_090516-CW-T'],
    #'986': ['NCC_CRC-986_100215-CW-T', 'NCC_CRC-986_261016-CW-T', 'NCC_CRC-986_300316-CW-T'],
    #'123': ['NCC_CRC-123_310715-CW-T', 'NCC_CRC-123_070116-CW-T', 'NCC_CRC-123_121115-CW-T'],
    '412': ['NCC_BRA-412_240820', 'NCC_BRA-412_060220']
}
patients = list(patientsample_dict.keys())
print(patients)

targetbedhg19 = pd.read_csv(os.path.join(*config.extdatafolder, 'Cancer226-targets_hg19.bed'), sep='\t')
targetbedhg38 = pd.read_csv(os.path.join(*config.extdatafolder, 'Cancer226-targets_hg38.bed'), sep='\t')

targetbedhg19list = []
for i, r in targetbedhg19.iterrows():
    ra = r[2] + 1 - r[1]
    rb = targetbedhg38.iloc[i, 2]+1-targetbedhg38.iloc[i, 1]
    if r[2] <= r[1]:
        raise ValueError('issue with {}'.format(r))
    for p in range(r[1], r[2]+1):
        targetbedhg19list.append('{}_{}'.format(r[0][3:], p))
    if ra != rb:
        print(ra, rb, i)
        for a in range(rb - ra):
            targetbedhg19list.append('{}_{}'.format(r[0][3:], p))
print("length of target bed file is {}".format(len(targetbedhg19list)))
print("example bed file locus for nomenclature {}".format(targetbedhg19list[0]))

targetbedhg38list = []
for _, r in targetbedhg38.iterrows():
    for p in range(r[1], r[2]+1):
        targetbedhg38list.append('{}_{}'.format(r[0], p))
print("length of target bed file is {}".format(len(targetbedhg38list)))
print("example bed file locus for nomenclature {}".format(targetbedhg38list[0]))

if not os.path.exists(os.path.join(*config.outputpath, 'figure2')):
    os.mkdir(os.path.join(*config.outputpath, 'figure2'))

#####################################################################################################################
# Figure 2 top: Identify elligible patients
#####################################################################################################################

#for patient in patients:
#    res = plot_patient_timeline(config, int(patient), figsize=(30, 8), mutations=True, highlight='discovery', treatment=False,
#                                save=False, savepath=os.path.join(*config.outputpath, 'figure2'))
#    lowtftimepoints_pd = get_mutations_stats(config, patient)
#    lowtftimepoints_pd.dropna()

#####################################################################################################################
# Figure 2 bottom: Paired plots logscale and same scale
#####################################################################################################################

targetbedhg19 = pd.read_csv(os.path.join(*config.extdatafolder, 'Cancer226-targets_hg19.bed'), sep='\t')
targetbedhg38 = pd.read_csv(os.path.join(*config.extdatafolder, 'Cancer226-targets_hg38.bed'), sep='\t')

for patient in patients:
    for sc in ['samescale', 'logscale']:
        paireplot(config, patient, targetbedhg19, targetbedhg38, sc=sc, nbhightfsamples=1, save=True, savepath=None)
