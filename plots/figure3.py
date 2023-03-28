
if __name__ == "__main__":

    from utils.config import Config
    from featureanalysis.featureanalysis import parse_caller_feature
    from utils.config import Config
    from utils.viz import *
    from utils.table import *
    from utils.metrics import *
    from utils.calltable import *
    from utils.calltableseries import *
    from utils.groundtruth import *
    from utils.metricsseries import *

    import sklearn
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import RandomizedSearchCV
    from sklearn.model_selection import StratifiedKFold
    from sklearn.model_selection import cross_validate

    import seaborn as sns
    import matplotlib.pyplot as plt

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")

    # Chomosome
    mixtureids =  ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T']
    mixtureid = 'CRC-986_100215-CW-T_CRC-986_300316-CW-T'
    #mixtureid = 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T'
    #mixtureid = 'CRC-123_310715-CW-T_CRC-123_121115-CW-T'
    reload = False
    save = False
    fixedvars=['coverage', 'ctdna']
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
    print('############# {} ############'.format(mixtureid))
    if mixtureid ==  'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T':
        chroms = [str(c) for c in range(1,23) if c !=17 and c !=8 and c!=5 and c!=19 and c!=20 and c!=21 and c!=22]
    elif mixtureid ==  'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
        chroms = [str(c) for c in range(1,23) if c !=1 and c!= 2 and c !=8 and c!=20 and c!=21 and c!=22]
    else:
        chroms = [str(c) for c in range(1,23) if c !=1 and c!= 2 and c !=6 and c!=20 and c!=21]
    calltables = {'sampleid':[], 'tf':[], 'cov':[], 'snv':[], 'indel':[], 'snp':[]}
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
    calltables['tf'] = np.unique([cn.split('_')[0] for cn in list(calltable_snv.columns)])[:-5].astype(float)
    calltables['snv'] = pd.concat(calltables['snv'])
    calltables['indel'] = pd.concat(calltables['indel'])
    calltables['snp'] = pd.concat(calltables['snp'])
    dilutionseries = aux.T[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]].T
    #for muttype in muttypes:
    muttype = 'snv'
    refsample = 'undiluted'
    if muttype == 'snv':
        gtm = 4
    else:  # elif muttype == 'indel':
        gtm = 2
    print(max(aux['tf']))
    if mixtureid ==  'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
        gtm = 3
        refsample = 'tissue'
        calltablesseries  = generate_groundtruth(config, calltables[muttype], aux['tf'], ground_truth_method='tissue', muttype=muttype,
                                                 matchedtissuepath=os.path.join('data', 'matchedtissue', 'NCC_CRC-986_100215-T1W', 'calls', 'NCC_CRC-986_100215-T1W_snv_calls_PASS_exome.csv'))
    else:
        calltablesseries = generate_groundtruth(config, calltables[muttype], aux['tf'], ground_truth_method=gtm, muttype=muttype)


    for method in config.methods:
        color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
        callmethod_snv_list = []
        for serie in [(70, 0), (70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]:
            st, sn = serie
            #print(serie)
            for chrom in chroms:
                #print(chrom)
                callmethod_snv, _, _ = parse_caller_feature(
                    'data/mixtures/mixtures_chr'+chrom+'/mixtures_chr'+chrom+'_CRC-986_100215-CW-T_CRC-986_300316-CW-T/mixture_chr'+chrom+'_CRC-986_100215-CW-T_'+str(st)+'x_CRC-986_300316-CW-T_'+str(sn)+'x',
                    method, save=False)
                if method in ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']:
                    callmethod_snv.drop(['ANN'], axis=1, inplace=True)
                #print(callmethod_snv.shape)
                #print(list(callmethod_snv.columns))
                #print(callmethod_snv.head())
                callmethod_snv_list.append(callmethod_snv)
        callmethod_snv_allchrom = pd.concat(callmethod_snv_list)

        # callmethod_snv_allchrom = pd.concat([callmethod_snv_allchrom, calltablesseries['truth']], axis=1)
        callmethod_snv_allchrom = callmethod_snv_allchrom.join(calltablesseries[['truth']], how='outer')
        callmethod_snv_allchrom['truth'].fillna(False, inplace=True)
        callmethod_snv_allchrom['truth'].value_counts()
        #callmethod_snv_allchrom.head()

        features = list(callmethod_snv_allchrom.columns)
        print(features)
        for x in [method, method+'_score',  method+'_vaf', 'type', 'alt', 'chrom', 'pos', 'ref', method+'_totcov', method+'_altcov']:
            print(x)
            features.remove(x)
        print(features)

        #callmethod_snv_allchrom = pd.concat([callmethod_snv_allchrom, calltablesseries['truth']], axis=1)
        callmethod_snv_allchrom['truth'].fillna(False, inplace=True)
        callmethod_snv_allchrom['truth'].value_counts()

        callmethod_snv_allchrom[method].value_counts()
        callmethod_snv_allchrom[method][callmethod_snv_allchrom[method] == 'PASS'] = 1
        callmethod_snv_allchrom[method][callmethod_snv_allchrom[method] != 1] = 0

        forbiden_features = ['AF', 'OLD_VARIANT', 'TYPE', "cosmid_id", 'AN', "NS", 'CIGAR']
        score_features = ['ODDS', 'TLOD', 'SomaticEVS', 'SSF', 'SPV']
        print(features)

        Xlist = []
        for feature in features:
            print(feature)
            test = 0
            try:
                aux = callmethod_snv_allchrom[feature].fillna(0).astype(float)
                test = 1
            except:
                try:
                    aux = callmethod_snv_allchrom[feature].str.split(',').str[0].fillna(0).astype(float)
                    test = 1
                except:
                    test = 0
            if test == 1 and feature != 'truth' and feature not in forbiden_features:
                Xlist.append(aux)
        X = pd.concat(Xlist, axis=1)
        y = callmethod_snv_allchrom['truth'].fillna(0)
        print(X.shape, y.shape)

        feature_name_dict = {}
        for i, x in enumerate(list(X.columns)):
            feature_name_dict[i] = x
        feature_name_dict

        model = RandomForestClassifier(random_state=0)
        grid = {
            'bootstrap': [True, False],
            'max_depth':  [10, 25, 50, 75, 100, 200, None], #[10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
            'max_features': ['auto', 'sqrt'],
            'min_samples_leaf': [1, 2, 4],
            'min_samples_split': [2, 5, 10],
            'n_estimators': [200, 400, 600, 800, 1000] #, 1200, 1400, 1600, 1800, 2000]
        }
        # if model CV has already been done, read back best params
        if os.path.exists(os.path.join(os.path.join('data', 'featureanalysis', 'bestparams_'+method+'.csv'))):
            bestparams = pd.read_csv(os.path.join('data', 'featureanalysis', 'bestparams_'+method+'.csv'), index_col=0)
            bestparams = bestparams.to_dict()['0']
            model = model.set_params(**bestparams)
        else:
            gridsearch = RandomizedSearchCV(estimator=model, param_distributions=grid, n_jobs=1, scoring="roc_auc", n_iter=2, cv=5, verbose=6, random_state=0)
            gridsearch.fit(X, y)
            print(gridsearch.best_params_)
            gridsearch.best_params_.to_csv(os.path.join(os.getcwd(), 'featureanalysis', 'bestparams_'+method+'.csv'))
            model = model.set_params(**gridsearch.best_params_) # retrive best hyper params
        cv = StratifiedKFold(n_splits=10, shuffle = True, random_state=0) # prepare 10-fold CV
        cv_results = cross_validate(model, X, y, cv=cv, scoring=['roc_auc', 'average_precision'], return_estimator=True, verbose=6)
        plt.rc('font', size=20)
        fi_pd_list = []
        for idx, estimator in enumerate(cv_results['estimator']):
            fi_pd = pd.DataFrame(estimator.feature_importances_,
                                 columns=['importance']).sort_values('importance', ascending=False)
            fi_pd_list.append(fi_pd)
        feature_importances = pd.concat(fi_pd_list, axis=1)
        feature_importances.columns = ['estimator_'+str(i) for i in range(10)]
        feature_importances.index = feature_name_dict.values()
        feature_importances.sort_values(by='estimator_0', ascending=False, inplace=True)
        feature_importances = feature_importances[(feature_importances.T != 0).any()]
        g = sns.catplot(data=feature_importances.T, kind="box", color=color_dict[method],
                        height=6, aspect=3.5);
        g.map_dataframe(sns.stripplot,
                        palette=["k"], dodge=True)
        g.set_xticklabels(rotation=45, ha='right')
        plt.xlabel('Features')
        plt.ylabel('Feature importance in Random Forest models in CV')
        plt.title(method + ' feature importance analysis', pad=50)
        if not os.path.exists(os.path.join(*config.outputpath, 'figure3')):
            os.mkdir(os.path.join(*config.outputpath, 'figure3'))
        plt.savefig(os.path.join(*config.outputpath, 'figure3', 'featureimportance_'+method+'.svg'), bbox_inches='tight')
        plt.show()