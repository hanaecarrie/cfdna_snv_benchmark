
import io
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import warnings
from tqdm.notebook import tqdm
from sklearn.metrics import precision_recall_curve, f1_score, average_precision_score
warnings.filterwarnings('ignore')
from sklearn.metrics import confusion_matrix

# set working directory
if not os.getcwd().endswith('cfdna_snv_benchmark'):
    os.chdir('../')
print('Current working directory: {}'.format(os.getcwd()))

from utils.config import Config
from benchmark.groundtruth import *

# Config and Display paramaters

config = Config("config/", "config_viz.yaml")
print(config.methods)
print(config.methods_tissue)

calltablesseries = pd.read_csv(os.path.join('data', 'mixtures_wholegenome', 'mixtures_allchr', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T_snv_calls_all.csv'), index_col=0, memory_map=True)
aux = pd.read_csv(os.path.join('data', 'mixtures_wholegenome', 'mixtures_allchr', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T_tf_cov.csv'), index_col=0)
gtm = 5
muttype = 'snv'
calltablesseries = generate_groundtruth(config, calltablesseries, aux['tf'], ground_truth_method=gtm, muttype=muttype,
                                        matchedtissuepath=None, methods=['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan', 'varnet', 'abemus', 'sinvict'])
print(calltablesseries)

calltablesseries[calltablesseries['truth'] == True].to_csv(os.path.join('data', 'mixtures_wholegenome', 'mixtures_allchr', 'groundtruths', 'groundtruths_CRC-986_100215-CW-T_CRC-986_300316-CW-T_snv_calls_all.csv'))