import warnings
warnings.filterwarnings('ignore')

from utils.calltable import *


def get_calltableseries(config, dilutionid, chrom, muttype='snv', filterparam='PASS', reload=False, save=False, diltype='mixture', concat='tf'):
    if diltype == 'mixture':
        confdilfolder = config.mixturefolder
    elif diltype == 'mixture_wes':
        confdilfolder = config.mixturefolderultradeep
    elif diltype == 'mixture_wgs':
        confdilfolder = config.mixturefolderwholegenome
    elif diltype == 'SEQC2':
        confdilfolder = config.mixturefolderSEQC2
    else: # spikein
        confdilfolder = config.spikeinfolder
    print(diltype, confdilfolder)
    if chrom in [str(c) for c in range(1, 23)] or diltype == 'mixture_wes' or diltype == 'SEQC2':
        if diltype.startswith('mixture_'):
            diltype = 'mixture'
        # Save table if do not exist and load tables
        calltables = {'sampleid': [], 'tf': [], 'vaf': [], 'cov': [], 'ichorcna': [], 'samplename' : [], 'snv': [], 'indel': [], 'snp': []}
        dilutionfolder = os.path.join(*confdilfolder, diltype+'s_chr' + chrom, diltype+'s_chr' + chrom +'_' + dilutionid)
        print(dilutionfolder)
        for dilutionpath in [l for l in os.listdir(dilutionfolder) if l.endswith('x') or l.endswith('T') or l.startswith('Sample')]:
            print(dilutionpath)
            print(reload, 'reload')
            if reload or not os.path.exists(os.path.join(dilutionfolder, dilutionpath, 'calls', dilutionpath+'_snv_calls_'+filterparam+'.csv')):
                calltable_snv, calltable_indel, calltable_snp = get_calltable(os.path.join(dilutionfolder, dilutionpath), config.methods, save=save, filter=filterparam)
            calltables['sampleid'].append(dilutionpath)
            if diltype == 'mixture':
                calltables['vaf'].append(np.nan)
                #TODO fix tb estimate 123 patient
                if dilutionid == 'CRC-123_310715-CW-T_CRC-123_121115-CW-T':
                    tx, nx = int(dilutionpath.split('_')[4][:-1]), int(dilutionpath.split('_')[7][:-1])
                    tf = 100 * 0.563 * tx / (tx + nx)
                    calltables['tf'].append(tf)
                    print(tf)
                else:
                    calltables['tf'].append(np.round(100*float(pd.read_csv(os.path.join(
                        dilutionfolder, dilutionpath, 'estimated_tf_chr'+chrom+dilutionpath[len((diltype+'_chr'+chrom)):]+'.txt')).columns[0]), 4))
                if os.path.exists(os.path.join(dilutionfolder, dilutionpath, 'ichorcna', dilutionpath[len((diltype+'_chr'+chrom)):]+'params.txt')):
                    calltables['ichorcna'].append(np.round(100*float(os.system("cat " +os.path.join(
                        dilutionfolder, dilutionpath, 'ichorcna', dilutionpath[len((diltype+'_chr'+chrom)):]+'params.txt') + " | grep 'Tumor Fraction:'")), 4)) #TODO fix
                else:
                    calltables['ichorcna'].append(np.nan)
                calltables['cov'].append(np.round(float(pd.read_csv(os.path.join(
                    dilutionfolder, dilutionpath, 'coverage_chr'+chrom+dilutionpath[len((diltype+'_chr'+chrom)):]+'.txt')).columns[0]), 4))
                calltables['samplename'].append(np.nan)
            elif diltype == 'SEQC2':
                calltables['samplename'].append(dilutionpath.split('_')[0])  # SampleDf or SampleEf
                calltables['tf'].append(np.nan)
                calltables['vaf'].append(np.nan)
                calltables['ichorcna'].append(np.nan)
                calltables['cov'].append(3000) #TODO fix
            else: # diltype spikein
                calltables['vaf'].append(float(dilutionpath.split('vaf')[1].split('_')[0]))
                calltables['tf'].append(np.nan)
                calltables['ichorcna'].append(np.nan)
                calltables['samplename'].append(np.nan)
                calltables['cov'].append(150) #TODO fix
            calltable_snv = pd.read_csv(os.path.join(
                dilutionfolder, dilutionpath, 'calls', dilutionpath+'_snv_calls_'+filterparam+'.csv'), index_col=0)
            calltable_indel = pd.read_csv(os.path.join(
                dilutionfolder, dilutionpath, 'calls', dilutionpath+'_indel_calls_'+filterparam+'.csv'), index_col=0)
            calltable_snp = pd.read_csv(os.path.join(
                dilutionfolder, dilutionpath, 'calls', dilutionpath+'_snp_calls_'+filterparam+'.csv'), index_col=0)
            calltables['snv'].append(calltable_snv)
            calltables['indel'].append(calltable_indel)
            calltables['snp'].append(calltable_snp)

            if concat == 'tf':
                for mt in ['snv', 'indel', 'snp']:
                    if not os.path.exists(os.path.join(dilutionfolder, 'calls')):
                        os.mkdir(os.path.join(dilutionfolder, 'calls'))
                    if not os.path.exists(os.path.join(dilutionfolder, 'calls', dilutionid + '_' + mt + '_calls_' + filterparam + '.csv')) or reload:
                        for ci, csnv in enumerate(calltables[mt]):
                            cols = ['chrom', 'pos', 'ref', 'alt', 'type']
                            if diltype == 'mixture':
                                varname = 'tf'
                            elif diltype == 'SEQC2':
                                varname = 'samplename'
                            else:
                                varname = 'vaf'
                            if diltype == 'SEQC2':
                                for m in config.methods:
                                    cols.append('{}_{}'.format(calltables[varname][ci], m))
                                    cols.append('{}_{}_score'.format(calltables[varname][ci], m))
                                for m in config.methods:
                                    cols.append('{}_{}_altcov'.format(calltables[varname][ci], m))
                                    cols.append('{}_{}_totcov'.format(calltables[varname][ci], m))
                                    cols.append('{}_{}_vaf'.format(calltables[varname][ci], m))
                            else:
                                for m in config.methods:
                                    cols.append('{:.2f}_{}'.format(calltables[varname][ci], m))
                                    cols.append('{:.2f}_{}_score'.format(calltables[varname][ci], m))
                                for m in config.methods:
                                    cols.append('{:.2f}_{}_altcov'.format(calltables[varname][ci], m))
                                    cols.append('{:.2f}_{}_totcov'.format(calltables[varname][ci], m))
                                    cols.append('{:.2f}_{}_vaf'.format(calltables[varname][ci], m))
                            #print(csnv.columns)
                            #print(cols)
                            csnv.columns = cols
                        # ensure no duplicated index
                        print(calltables[mt][0].loc[calltables[mt][0].index[calltables[mt][0].index.duplicated(keep=False)]].shape[0])
                        # get call series
                        calltablesseries = pd.concat([ct.set_index(['chrom', 'pos', 'ref', 'alt', 'type']) for ct in calltables[mt]], axis=1)
                        calltablesseries.reset_index(inplace=True)
                        calltablesseries['chrom_pos_ref_alt'] = calltablesseries['chrom'].astype('str').str.cat(calltablesseries['pos'].astype('str'), sep="_").str.cat(calltablesseries['ref'].astype('str'), sep='_').str.cat(calltablesseries['alt'].astype('str'), sep='_')
                        calltablesseries.set_index('chrom_pos_ref_alt', inplace=True)
                        print(calltablesseries.shape)
                        calltablesseries.to_csv(os.path.join(dilutionfolder, 'calls', dilutionid + '_' + mt + '_calls_' + filterparam + '.csv'))
                        calltables_aux = dict(calltables)
                        calltables_aux.pop('snv')
                        calltables_aux.pop('indel')
                        calltables_aux.pop('snp')
                        calltables_aux = pd.DataFrame.from_dict(calltables_aux)
                        calltables_aux.set_index('sampleid', inplace=True)
                        calltables_aux.sort_values(by='tf', ascending=False, inplace=True)
                        calltables_aux.to_csv(os.path.join(dilutionfolder, 'calls', dilutionid + '_tf_cov.csv'))
                calltablesseries = pd.read_csv(os.path.join(dilutionfolder, 'calls', dilutionid + '_' + muttype + '_calls_' + filterparam + '.csv'), index_col=0)
            elif concat == 'vaf':
                for mt in ['snv', 'indel', 'snp']:
                    if not os.path.exists(os.path.join(dilutionfolder, 'calls')):
                        os.mkdir(os.path.join(dilutionfolder, 'calls'))
                    if not os.path.exists(os.path.join(dilutionfolder, 'calls', dilutionid + '_' + mt + '_calls_' + filterparam + '_' + concat + '.csv')) or reload:
                        for ci, csnv in enumerate(calltables[mt]):
                            varname = 'tf' if diltype == 'mixture' else 'vaf'
                            print(varname, ci, mt)
                            print('{:.2f}'.format(calltables[varname][ci]))
                            calltables[mt][ci]['sampletf'] = '{:.2f}'.format(calltables[varname][ci])
                        # ensure no duplicated index
                        print(calltables[mt][0].loc[calltables[mt][0].index[calltables[mt][0].index.duplicated(keep=False)]].shape[0])
                        # get call series
                        calltablesseries = pd.concat(calltables[mt], axis=0)
                        calltablesseries.reset_index(inplace=True)
                        calltablesseries['chrom_pos_ref_alt'] = calltablesseries['chrom'].astype('str').str.cat(calltablesseries['pos'].astype('str'), sep="_").str.cat(calltablesseries['ref'].astype('str'), sep='_').str.cat(calltablesseries['alt'].astype('str'), sep='_')
                        calltablesseries.set_index('chrom_pos_ref_alt', inplace=True)
                        print(calltablesseries.shape)
                        calltablesseries.to_csv(os.path.join(dilutionfolder, 'calls', dilutionid + '_' + mt + '_calls_' + filterparam + '_' + concat + '.csv'))
                        calltables_aux = dict(calltables)
                        calltables_aux.pop('snv')
                        calltables_aux.pop('indel')
                        calltables_aux.pop('snp')
                        #print(calltables_aux)
                        calltables_aux = pd.DataFrame.from_dict(calltables_aux)
                        calltables_aux.set_index('sampleid', inplace=True)
                        calltables_aux.sort_values(by='tf', ascending=False, inplace=True)
        calltablesaux = pd.read_csv(os.path.join(dilutionfolder, 'calls', dilutionid + '_tf_cov.csv'), index_col=0)
        return calltablesseries, calltablesaux

    else:  # chrom == 'all'
        if diltype.startswith('mixture_'):
            diltype = 'mixture'
        if chrom == 'all':
            chroms = [str(c) for c in range(1, 23)]
        else:
            chroms = chrom
        print(chroms)
        # get median tf and cov estimates
        caux_list = []
        for chrom in chroms:
            if os.path.exists(os.path.join(*confdilfolder, diltype+'s_chr'+chrom, diltype+'s_chr' + chrom +'_' + dilutionid, 'calls', dilutionid + '_tf_cov.csv')) and not reload:
                caux = pd.read_csv(os.path.join(*confdilfolder, diltype+'s_chr' + chrom, diltype+'s_chr' + chrom +'_' + dilutionid, 'calls', dilutionid + '_tf_cov.csv'), index_col=0)
            else:
                _, caux = get_calltableseries(config, dilutionid, chrom, muttype=muttype, filterparam=filterparam, reload=reload, save=save, diltype=diltype, concat=concat)
            caux.index = [ci.split('_')[0] + '_'+ '_'.join(ci.split('_')[2:]) for ci in list(caux.index)]
            caux_list.append(caux)
        caux_df = pd.DataFrame()
        caux_df['tf'] = pd.concat([c[['tf']] for c in caux_list], axis=1).median(axis=1)
        caux_df['cov'] = pd.concat([c[['cov']] for c in caux_list], axis=1).median(axis=1)
        #print('new tb for columns names')
        if diltype == 'mixture':
            newtf_dict = dict(zip(list(caux_df['tf'].index), ['{:.2f}'.format(c) for c in caux_df['tf'].values]))
            if len(np.unique(caux_df['tf'].values)) < len(caux_df['tf'].values):
                raise ValueError('taking median of TFs gives duplicate values')
        calltablesserieschroms = []
        # get concatenated calltable
        for chrom in chroms:
            if concat == 'vaf':
                callfile = os.path.join(*confdilfolder, diltype+'s_chr' + chrom, diltype+'s_chr' + chrom +'_' + dilutionid, 'calls', dilutionid + '_' + muttype + '_calls_' + filterparam + '_' + concat + '.csv')
            else:
                callfile = os.path.join(*confdilfolder, diltype+'s_chr' + chrom, diltype+'s_chr' + chrom +'_' + dilutionid, 'calls', dilutionid + '_' + muttype + '_calls_' + filterparam + '.csv')
            calltable = pd.read_csv(callfile, index_col=0)
            if diltype == 'mixture':
                caux = pd.read_csv(os.path.join(*confdilfolder, diltype+'s_chr' + chrom, diltype+'s_chr' + chrom +'_' + dilutionid, 'calls', dilutionid + '_tf_cov.csv'), index_col=0)
                caux.index = [c.split('_')[0] + '_' +  '_'.join(c.split('_')[2:]) for c in list(caux.index)]
                oldtf_dict_inverted = dict(zip(['{:.2f}'.format(c) for c in caux['tf'].values], list(caux['tf'].index)))
                newcol = list(calltable.columns)
                for n, nc in enumerate(newcol):
                    if nc.split('_')[0] in list(oldtf_dict_inverted.keys()):
                        dilid = oldtf_dict_inverted[nc.split('_')[0]]
                        newcol[n] = nc.replace(nc.split('_')[0], newtf_dict[dilid])
                calltable.columns = newcol
            calltablesserieschroms.append(calltable)
        calltablesserieschroms = pd.concat(calltablesserieschroms, axis=0)
        if not os.path.exists(os.path.join(*confdilfolder, diltype+'s_allchr')):
            os.mkdir(os.path.join(*confdilfolder, diltype+'s_allchr'))
        if concat == 'vaf':
            calltablesserieschroms.to_csv(os.path.join(*confdilfolder, diltype+'s_allchr', dilutionid + '_' + muttype + '_calls_' + filterparam + '_' + concat + '.csv'))
            calltablesserieschroms = pd.read_csv(os.path.join(*confdilfolder, diltype+'s_allchr', dilutionid + '_' + muttype + '_calls_' + filterparam + '_' + concat + '.csv'), index_col=0)
        else:
            calltablesserieschroms.to_csv(os.path.join(*confdilfolder, diltype+'s_allchr', dilutionid + '_' + muttype + '_calls_' + filterparam + '.csv'))
            calltablesserieschroms = pd.read_csv(os.path.join(*confdilfolder, diltype+'s_allchr', dilutionid + '_' + muttype + '_calls_' + filterparam + '.csv'), index_col=0)
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
