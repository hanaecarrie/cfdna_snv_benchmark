import warnings
warnings.filterwarnings('ignore')

from utils.calltable import *


def get_calltableseries(config, mixtureid, chrom, muttype='snv', filterparam='PASS', reload=False, save=False):

    # Save table if do not exist and load tables
    calltables = {'sampleid': [], 'tf': [], 'cov': [], 'snv': [], 'indel': [], 'snp': []}
    mixturefolder = os.path.join(*config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+mixtureid)
    for mixturepath in [l for l in os.listdir(mixturefolder) if l.endswith('x') or l.endswith('T')]:
        print(mixturepath)
        if not os.path.exists(os.path.join(mixturefolder, mixturepath, 'calls', mixturepath+'_snv_calls_'+filterparam+'.csv')) or reload:
            calltable_snv, calltable_indel, calltable_snp = get_calltable(os.path.join(mixturefolder, mixturepath), config.methods, save=save, filter=filterparam)
        calltables['sampleid'].append(mixturepath)
        calltables['tf'].append(np.round(100*float(pd.read_csv(os.path.join(
            mixturefolder, mixturepath, 'estimated_tf_chr'+chrom+mixturepath[len(('mixture_chr'+chrom)):]+'.txt')).columns[0]), 4))
        calltables['cov'].append(np.round(float(pd.read_csv(os.path.join(
            mixturefolder, mixturepath, 'coverage_chr'+chrom+mixturepath[len(('mixture_chr'+chrom)):]+'.txt')).columns[0]), 4))
        calltable_snv = pd.read_csv(os.path.join(
            mixturefolder, mixturepath, 'calls', mixturepath+'_snv_calls_'+filterparam+'.csv'), index_col=0)
        calltable_indel = pd.read_csv(os.path.join(
            mixturefolder, mixturepath, 'calls', mixturepath+'_indel_calls_'+filterparam+'.csv'), index_col=0)
        calltable_snp = pd.read_csv(os.path.join(
            mixturefolder, mixturepath, 'calls', mixturepath+'_snp_calls_'+filterparam+'.csv'), index_col=0)
        calltables['snv'].append(calltable_snv)
        calltables['indel'].append(calltable_indel)
        calltables['snp'].append(calltable_snp)

    for mt in ['snv', 'indel', 'snp']:
        if not os.path.exists(os.path.join(mixturefolder, 'calls')):
            os.mkdir(os.path.join(mixturefolder, 'calls'))
        if not os.path.exists(os.path.join(mixturefolder, 'calls', mixtureid+'_'+mt+'_calls_'+filterparam+'.csv')) or reload:
            for ci, csnv in enumerate(calltables[mt]):
                cols = ['chrom', 'pos', 'ref', 'alt', 'type']
                for m in config.methods:
                    cols.append('{:.2f}_{}'.format(calltables['tf'][ci], m))
                    cols.append('{:.2f}_{}_score'.format(calltables['tf'][ci], m))
                for m in config.methods:
                    cols.append('{:.2f}_{}_altcov'.format(calltables['tf'][ci], m))
                    cols.append('{:.2f}_{}_totcov'.format(calltables['tf'][ci], m))
                    cols.append('{:.2f}_{}_vaf'.format(calltables['tf'][ci], m))
                csnv.columns = cols
            # ensure no duplicated index
            print(calltables[mt][0].loc[calltables[mt][0].index[calltables[mt][0].index.duplicated(keep=False)]].shape[0])
            # get call series
            calltablesseries = pd.concat([ct.set_index(['chrom', 'pos', 'ref', 'alt', 'type']) for ct in calltables[mt]], axis=1)
            calltablesseries.reset_index(inplace=True)
            calltablesseries['chrom_pos_ref_alt'] = calltablesseries['chrom'].astype('str').str.cat(calltablesseries['pos'].astype('str'), sep="_").str.cat(calltablesseries['ref'].astype('str'), sep='_').str.cat(calltablesseries['alt'].astype('str'), sep='_')
            calltablesseries.set_index('chrom_pos_ref_alt', inplace=True)
            print(calltablesseries.shape)
            calltablesseries.to_csv(os.path.join(mixturefolder, 'calls', mixtureid+'_'+mt+'_calls_'+filterparam+'.csv'))

    calltablesseries = pd.read_csv(os.path.join(mixturefolder, 'calls', mixtureid+'_'+muttype+'_calls_'+filterparam+'.csv'), index_col=0)

    return calltablesseries, calltables


if __name__ == "__main__":
    import os
    from utils.config import Config

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")

    chrom = '22'
    reload = True
    save = True
    filterparam = 'all'
    muttypes = ['snv', 'indel']
    mixtureids = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T']

    for muttype in muttypes:
        for mixtureid in mixtureids:
            print('############# {} {} ############'.format(mixtureid, muttype))
            calltablesseries, calltables = get_calltableseries(config, mixtureid, chrom, muttype, filterparam, reload, save)
            print(calltablesseries.head())
