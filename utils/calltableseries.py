import warnings
warnings.filterwarnings('ignore')

from utils.calltable import *


def get_calltableseries(config, mixtureid, chrom, muttype='snv', filterparam='PASS', reload=False, save=False):

    if chrom in [str(c) for c in range(1, 23)]:
        # Save table if do not exist and load tables
        calltables = {'sampleid': [], 'tf': [], 'cov': [], 'ichorcna': [], 'snv': [], 'indel': [], 'snp': []}
        mixturefolder = os.path.join(*config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+mixtureid)
        for mixturepath in [l for l in os.listdir(mixturefolder) if l.endswith('x') or l.endswith('T')]:
            print(mixturepath)
            if not os.path.exists(os.path.join(mixturefolder, mixturepath, 'calls', mixturepath+'_snv_calls_'+filterparam+'.csv')) or reload:
                calltable_snv, calltable_indel, calltable_snp = get_calltable(os.path.join(mixturefolder, mixturepath), config.methods, save=save, filter=filterparam)
            calltables['sampleid'].append(mixturepath)
            if mixtureid == 'CRC-123_310715-CW-T_CRC-123_121115-CW-T':
                tx, nx = int(mixturepath.split('_')[4][:-1]), int(mixturepath.split('_')[7][:-1])
                tf = 100 * 0.563 * tx / (tx + nx)
                calltables['tf'].append(tf)
                print(tf)
            else:
                calltables['tf'].append(np.round(100*float(pd.read_csv(os.path.join(
                    mixturefolder, mixturepath, 'estimated_tf_chr'+chrom+mixturepath[len(('mixture_chr'+chrom)):]+'.txt')).columns[0]), 4))
            print(os.path.exists(os.path.join(mixturefolder, mixturepath, 'ichorcna', mixturepath[len(('mixture_chr'+chrom)):]+'params.txt')))
            if os.path.exists(os.path.join(mixturefolder, mixturepath, 'ichorcna', mixturepath[len(('mixture_chr'+chrom)):]+'params.txt')):
                calltables['ichorcna'].append(np.round(100*float(os.system("cat " +os.path.join(
                    mixturefolder, mixturepath, 'ichorcna', mixturepath[len(('mixture_chr'+chrom)):]+'params.txt') + " | grep 'Tumor Fraction:'")), 4))
            else:
                calltables['ichorcna'].append(np.nan)
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
                calltables_aux = dict(calltables)
                calltables_aux.pop('snv')
                calltables_aux.pop('indel')
                calltables_aux.pop('snp')
                calltables_aux = pd.DataFrame.from_dict(calltables_aux)
                calltables_aux.set_index('sampleid', inplace=True)
                calltables_aux.to_csv(os.path.join(mixturefolder, 'calls', mixtureid+'_tf_cov.csv'))
        calltablesseries = pd.read_csv(os.path.join(mixturefolder, 'calls', mixtureid+'_'+muttype+'_calls_'+filterparam+'.csv'), index_col=0)
        calltablesaux = pd.read_csv(os.path.join(mixturefolder, 'calls', mixtureid+'_tf_cov.csv'), index_col=0)
        return calltablesseries, calltablesaux

    else:  # chrom == 'all'
        if chrom == 'all':
            chroms = [str(c) for c in range(1, 23)]
        else:
            chroms = chrom
        print(chroms)
        # get median tf and cov estimates
        caux_list = []
        for chrom in chroms:
            if os.path.exists(os.path.join(*config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+mixtureid, 'calls', mixtureid+'_tf_cov.csv')):
                caux = pd.read_csv(os.path.join(*config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+mixtureid, 'calls', mixtureid+'_tf_cov.csv'), index_col=0)
            else:
                _, caux = get_calltableseries(config, mixtureid, chrom, muttype=muttype, filterparam=filterparam, reload=reload, save=save)
            caux.index = [ci.split('_')[0] + '_'+'_'.join(ci.split('_')[2:]) for ci in list(caux.index)]
            caux_list.append(caux)
        caux_df = pd.DataFrame()
        caux_df['tf'] = pd.concat([c[['tf']] for c in caux_list], axis=1).median(axis=1)
        caux_df['cov'] = pd.concat([c[['cov']] for c in caux_list], axis=1).median(axis=1)
        if len(np.unique(caux_df['tf'].values)) < len(caux_df['tf'].values):
            raise ValueError('taking median of TFs gives duplicate values')
        calltablesserieschroms = []
        # get concatenated calltable
        for chrom in chroms:
            callfile = os.path.join(*config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+mixtureid, 'calls', mixtureid+'_'+muttype+'_calls_'+filterparam+'.csv')
            calltable = pd.read_csv(callfile, index_col=0)
            caux = pd.read_csv(os.path.join(*config.mixturefolder, 'mixtures_chr'+chrom, 'mixtures_chr'+chrom+'_'+mixtureid, 'calls', mixtureid+'_tf_cov.csv'), index_col=0)
            oldtf = caux['tf'].values
            newtf = [caux_df.loc[caux[caux['tf'] == otf].index[0].split('_')[0] + '_' + '_'.join(caux[caux['tf'] == otf].index[0].split('_')[2:]), 'tf'] for otf in oldtf]
            newcol = list(calltable.columns)
            for n, nc in enumerate(newcol):
                if nc.split('_')[0] in oldtf:
                    i = oldtf.index(nc.split('_')[0])
                    newcol[n] = nc.replace(oldtf[i], '{:.2f}'.format(newtf[i]))
            calltable.columns = newcol
            calltablesserieschroms.append(calltable)
        calltablesserieschroms = pd.concat(calltablesserieschroms, axis=0)
        if not os.path.exists(os.path.join(*config.mixturefolder, 'mixtures_allchr')):
            os.mkdir(os.path.join(*config.mixturefolder, 'mixtures_allchr'))
        calltablesserieschroms.to_csv(os.path.join(*config.mixturefolder, 'mixtures_allchr', mixtureid+'_'+muttype+'_calls_'+filterparam+'.csv'))
        calltablesserieschroms = pd.read_csv(os.path.join(*config.mixturefolder, 'mixtures_allchr', mixtureid+'_'+muttype+'_calls_'+filterparam+'.csv'), index_col=0)
        return calltablesserieschroms, caux_df


if __name__ == "__main__":
    import os
    from utils.config import Config

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")

    reload = True
    save = True
    filterparam = 'all'
    muttypes = ['snv', 'indel', 'snp']
    mixtureids = ['CRC-123_310715-CW-T_CRC-123_121115-CW-T', 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T']

    for muttype in muttypes:
        for mixtureid in mixtureids:
            for chrom in range(1, 23):
                #TODO fix call bugs
                if (mixtureid != 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T') or (chrom != 17 and chrom != 8):
                    chrom = str(chrom)
                    print('############# {} {} ############'.format(mixtureid, muttype))
                    calltablesseries, calltables = get_calltableseries(config, mixtureid, chrom, muttype, filterparam, reload, save)
                    print(calltablesseries.head())
