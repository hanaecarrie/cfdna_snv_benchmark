
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
    from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV
    from sklearn.model_selection import cross_validate
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import average_precision_score
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    import pickle

    import seaborn as sns
    import matplotlib.pyplot as plt

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")

    for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict' if m!= 'freebayes']:

        print(method)

        # Prepare patient 1014 to study feature importance
        mixtureids = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T']
        mixtureids_train = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T']
        mixtureid_test = 'CRC-986_100215-CW-T_CRC-986_300316-CW-T'
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
            seriesorder = [(70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]
            xaxis = 'tumor burden'
        elif fixedvar == 'ctdna':
            seriesorder = [(70, 0), (70, 80), (70, 180)]
            xaxis = 'coverage'
        # print('############# {} ############'.format(mixtureid))
        X, y = pd.DataFrame(), pd.Series()
        print(X.shape)
        print(y.shape)
        for mixtureid in mixtureids_train:
            patientid = str(mixtureid.split('_')[0].split('-')[1])
            if mixtureid == 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T':
                chroms = [str(c) for c in range(1,23) if c != 17 and c != 8 and c != 5 and c != 19 and c != 20 and c != 21 and c != 22]
            elif mixtureid == 'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
                chroms = [str(c) for c in range(1,23) if c != 1 and c != 2 and c != 8 and c != 20 and c != 21 and c != 22]
            else:
                chroms = [str(c) for c in range(1,23) if c != 1 and c != 2 and c != 6 and c != 20 and c != 21]
            calltables = {'sampleid': [], 'tf': [], 'cov': [], 'snv': [], 'indel': [], 'snp': []}
            aux_all = []
            calltable_snv, aux = get_calltableseries(config, mixtureid, chroms, muttype='snv', filterparam=filterparam, reload=reload, save=save)
            #calltable_indel, aux = get_calltableseries(config, mixtureid, chroms, muttype='indel', filterparam=filterparam, reload=reload, save=save)
            #calltable_snp, aux = get_calltableseries(config, mixtureid, chroms, muttype='snp', filterparam=filterparam, reload=reload, save=save)
            #print(calltable_snv.shape, calltable_indel.shape, calltable_snp.shape)
            print(aux)
            calltables['snv'].append(calltable_snv)
            #calltables['indel'].append(calltable_indel)
            #calltables['snp'].append(calltable_snp)
            calltables['sampleid'] = mixtureid
            calltables['tf'] = np.unique([cn.split('_')[0] for cn in list(calltable_snv.columns)])[:-5].astype(float)
            calltables['snv'] = pd.concat(calltables['snv'])
            #calltables['indel'] = pd.concat(calltables['indel'])
            #calltables['snp'] = pd.concat(calltables['snp'])
            dilutionseries = aux.T[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]].T
            dilutionseries = pd.DataFrame(dilutionseries.iloc[0]).T
            print(dilutionseries)
            muttype = 'snv'
            refsample = 'undiluted'
            if muttype == 'snv':
                gtm = 5
            else:  # elif muttype == 'indel':
                gtm = 4
            print(max(aux['tf']))
            calltablesseries = generate_groundtruth(config, calltables[muttype], aux['tf'], ground_truth_method=gtm, muttype=muttype,
                                                    methods=['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan', 'varnet', 'abemus', 'sinvict'])
            print(calltablesseries['truth'].value_counts())
            print(calltablesseries.columns)
            print(calltablesseries.head())
            print(aux['tf'])
            aux = aux.loc[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]]
            aux = pd.DataFrame(aux.iloc[0]).T
            print(aux)
            tfl = list(aux['tf'].values)
            print(tfl)
            print(list(calltable_snv.columns))
            calltable_snv = calltable_snv[[c for c in list(calltable_snv.columns) if c in tfl or c in ['chrom', 'pos', 'ref', 'alt', 'type', 'sampletf']]]
            calltablesseries = calltablesseries[[c for c in list(calltablesseries.columns) if c in tfl or c in ['chrom', 'pos', 'ref', 'alt', 'type', 'sampletf', 'truth']]]

            # Extract features
            np.random.seed(1234)
            color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
            callmethod_snv_all_list = []
            for serie in [(70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]:
                st, sn = serie
                callmethod_snv_list = []
                for chrom in chroms:
                    callmethod_snv, _, _ = parse_caller_feature(
                        'data/mixtures/mixtures_chr'+chrom+'/mixtures_chr'+chrom+'_'+mixtureid + '/mixture_chr'+chrom+'_'+'_'.join(mixtureid.split('_')[:2])+'_'+str(st)+'x_'+'_'.join(mixtureid.split('_')[2:])+'_'+str(sn)+'x',
                        method, save=False)
                    if method in ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']:
                        callmethod_snv.drop(['ANN'], axis=1, inplace=True)
                    callmethod_snv_list.append(callmethod_snv)
                callmethod_snv_all = pd.concat(callmethod_snv_list)
                callmethod_snv_all = pd.concat([callmethod_snv_all, calltablesseries[['truth']]], axis=1)
                print(callmethod_snv_all['truth'].value_counts())
                callmethod_snv_all.index = [idx + '_' +str(st)+'x'+str(sn)+'x' for idx in list(callmethod_snv_all.index)]
                callmethod_snv_all_list.append(callmethod_snv_all)
            callmethod_snv_allchrom = pd.concat(callmethod_snv_all_list)
            print(callmethod_snv_allchrom['truth'].value_counts())
            callmethod_snv_allchrom['truth'].fillna(False, inplace=True)
            print(callmethod_snv_allchrom.columns)
            print(callmethod_snv_allchrom)
            features = list(callmethod_snv_allchrom.columns)
            print(features)
            for x in [method, method+'_score',  method+'_vaf', 'type', 'alt', 'chrom', 'pos', 'ref', method+'_totcov', method+'_altcov']:
                print(x)
                features.remove(x)
            print(features)
            callmethod_snv_allchrom['truth'].fillna(False, inplace=True)
            print(callmethod_snv_allchrom['truth'].value_counts())
            #print(callmethod_snv_allchrom[method].value_counts())
            #callmethod_snv_allchrom[method][callmethod_snv_allchrom[method] == 'PASS'] = 1
            #callmethod_snv_allchrom[method][callmethod_snv_allchrom[method] != 1] = 0
            forbiden_features = ['AF', 'OLD_VARIANT', 'TYPE', "cosmid_id", 'AN', "NS", 'CIGAR', 'END_POS_REF', 'DP', 'AF', 'AO', 'VD', 'N_refDepth', 'N_altDepth', 'relcov']
            score_features = ['ODDS', 'TLOD', 'SomaticEVS', 'SSF', 'SPV']
            print(features)
            Xlist = []
            for feature in features:
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
                    print(feature)
                else:
                    print(feature, 'not included')
            Xtmp = pd.concat(Xlist, axis=1)
            ytmp = callmethod_snv_allchrom['truth'].fillna(0)
            print(Xtmp.shape, ytmp.shape)
            if X.shape[0] == 0:
                X, y = Xtmp, ytmp
            else:
                X, y = pd.concat([X, Xtmp]), pd.concat([y, ytmp])
        print('mixtures_train size: X shape = {}, y shape = {}'.format(X.shape, y.shape))
        Xsave, ysave = X.copy(), y.copy()

        # train random forest to get feature importance on 2 train patients
        feature_name_dict = {}
        for i, x in enumerate(list(X.columns)):
            feature_name_dict[i] = x
        model = RandomForestClassifier(random_state=0)
        grid = {
            'bootstrap': [True, False],
            'max_depth':  [5, 10, 15, 20, 50], #[25, 50, 75, 100, 200, None 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
            'max_features': ['auto', 'sqrt'],
            'min_samples_leaf': [1, 2, 4],
            'min_samples_split': [2, 5, 10],
            'n_estimators': [200, 500, 1000] #400, 600, 800, 1200, 1400, 1600, 1800, 2000]
        }
        #TODO gridsearch
        if os.path.exists("data/featureanalysis/model_"+method+"_featureimportance_150x.pkl"):
            with open("data/featureanalysis/model_"+method+"_featureimportance_150x.pkl", "rb") as f:
                modelRFC = pickle.load(f)
        else:
            clf = RandomizedSearchCV(estimator=model, param_distributions=grid, n_iter=1, scoring=['average_precision', 'recall'], refit='average_precision', n_jobs=8, verbose=2)
            search = clf.fit(X, y)
            bestparams = search.best_params_
            print(bestparams)
            pd.Series(bestparams).to_csv("data/featureanalysis/bestparams_"+method+"_featureimportance_150x.csv")
            modelRFC = RandomForestClassifier(**bestparams, random_state=0)
            with open("data/featureanalysis/model_"+method+"_featureimportance_150x.pkl", "wb") as f:
                pickle.dump(modelRFC, f)
        if not os.path.exists('data/featureanalysis/cvresults_featureimportance_'+method+'_trainpatients_150x_RF.pkl'):
            cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=0)  # prepare 10-fold CV
            print(list(X.columns))
            print('is there nan? {}'.format(np.sum(np.isnan(X), axis=0)))
            print('is there only finite values? {}'.format(np.sum(~np.isfinite(X), axis=0)))
            if method == 'vardict':
                X['SOR'] = X['SOR'].astype('float32')
                X['SOR'].apply(np.nan_to_num)
            X = X.apply(np.nan_to_num)
            y = y.apply(np.nan_to_num)
            X, y = sklearn.utils.check_X_y(X, y)
            print('is there nan? {}'.format(np.sum(np.isnan(X), axis=0)))
            print('is there only finite values? {}'.format(np.sum(~np.isfinite(X), axis=0)))
            cv_results = cross_validate(modelRFC, X, y, cv=cv, scoring=['roc_auc', 'average_precision', 'recall'], return_estimator=True, verbose=6)
            # save cv results
            with open('data/featureanalysis/cvresults_featureimportance_'+method+'_trainpatients_150x_RF.pkl', 'wb') as fid:
                pickle.dump(cv_results, fid)
            #pd.DataFrame(cv_results).to_csv('data/featureanalysis/cvresults_featureimportance_'+method+'_trainpatients_150x_RF.csv')
        else:
            with open('data/featureanalysis/cvresults_featureimportance_'+method+'_trainpatients_150x_RF.pkl', 'rb') as fid:
                cv_results = pickle.load(fid)
            #cv_results = pd.read_csv('data/featureanalysis/cvresults_featureimportance_'+method+'_trainpatients_150x_RF.csv', index_col=0)
        plt.rc('font', size=20)
        fi_pd_list = []
        for idx, estimator in enumerate(cv_results['estimator']):
            fi_pd = pd.DataFrame(estimator.feature_importances_,
                                 columns=['importance']).sort_values('importance', ascending=False)
            fi_pd_list.append(fi_pd)
        feature_importances = pd.concat(fi_pd_list, axis=1)
        feature_importances.columns = ['estimator_'+str(i) for i in range(10)]
        feature_importances.index = feature_name_dict.values()
        print(feature_importances.columns)
        print(feature_importances.index)
        print(list(set((feature_importances[(feature_importances.T > 0.001).any()].index)) | set([l for l in list(feature_importances.index) if l in score_features])))
        feature_importances = feature_importances.loc[list(set((feature_importances[(feature_importances.T > 0.001).any()].index)) | set([l for l in list(feature_importances.index) if l in score_features]))]
        feature_importances['median'] = feature_importances.median(axis=1)
        feature_importances.sort_values(by='median', ascending=False, inplace=True)
        feature_importances.drop('median', axis=1, inplace=True)
        plt.figure(figsize=(20, 10))
        g = sns.barplot(data=feature_importances.T, color=color_dict[method])
        g.set_xticklabels(g.get_xticklabels(), rotation=90)
        plt.xlabel('Features')
        plt.ylabel('Feature importance in Random Forest models in CV')
        plt.title(method + ' feature importance analysis. CV ROC AUC (acc {:.2f} +/- {:.2f})'.format(np.mean(cv_results['test_roc_auc']), np.std(cv_results['test_roc_auc'])))
        feature_importances.T.to_csv('data/featureanalysis/barplot_featureimportance_'+method+'_'+patientid+'_150x_RF.csv')
        plt.savefig('data/featureanalysis/featureimportance_'+method+'_trainpatients_150x.svg')
        plt.savefig('data/featureanalysis/featureimportance_'+method+'_trainpatients_150x.png')
        plt.show()
        
        # Add top informative features
        c = 0
        features_selected = []
        feature_score = {'freebayes': 'ODDS', 'mutect2': 'TLOD', 'strelka2': 'SomaticEVS', 'vardict': 'SSF', 'varscan': 'SSC', 'smurf': 'FALSE.'}
        for fi in list(feature_importances.index):
            if fi not in ['DP', 'QSS', 'QSS_NT', 'AF', 'AO', 'VD', 'N_refDepth', 'N_altDepth', 'relcov'] and c < 5:
                features_selected.append(fi)
                c += 1
            if c == 4 and feature_score[method] not in features_selected:
                features_selected.append(feature_score[method])
                c += 1
        print(features_selected)

        # X = callmethod_snv_allchrom[list(features_selected)].fillna(0)
        X = Xsave[list(features_selected)]
        print(X.shape)
        if method == 'freebayes':
            if 'AC' in X.columns:
                X['AC'].loc[X['AC'].astype(str).str.contains(',') == True] = X['AC'].loc[X['AC'].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
                X['AC'] = X['AC'].astype(float)
        if method == 'mutect2':
            X['TLOD'].loc[X['TLOD'].astype(str).str.contains(',') == True] = X['TLOD'].loc[X['TLOD'].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
            X['TLOD'] = X['TLOD'].astype(float)
            for feat in ['NALOD', 'AS_SB_TABLE', 'FS', 'NLOD', 'MQRankSum']:
                if feat in X.columns:
                    X[feat].loc[X[feat].astype(str).str.contains(',') == True] = X[feat].loc[X[feat].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
                    X[feat] = X[feat].astype(float)
        # y = callmethod_snv_allchrom['truth']
        y = ysave
        print(y.value_counts())

        # check comma in X
        for col in list(X.columns):
            print(col)
            print(X[col])
            print(X[col].loc[X[col].astype(str).str.contains(',') == True])
            print(X[col].loc[X[col].astype(str).str.contains(',') == True].shape)
        # Build cfDNA specific model
        cfdnaspecific_model = RandomForestClassifier(random_state=0)
        grid = {
            'bootstrap': [True, False],
            'max_depth':  [5, 10, 15, 20, 50], #[25, 50, 75, 100, 200, None 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
            'max_features': ['auto', 'sqrt'],
            'min_samples_leaf': [1, 2, 4],
            'min_samples_split': [2, 5, 10],
            'n_estimators': [200, 500, 1000] #400, 600, 800, 1200, 1400, 1600, 1800, 2000]
        }
        #TODO gridsearch
        if os.path.exists("data/featureanalysis/cfdnaspecificmodel_"+method+"_986_150x.pkl"):
            with open("data/featureanalysis/cfdnaspecificmodel_"+method+"_986_150x.pkl", "rb") as f:
                cfdnaspecific_modeRFC = pickle.load(f)
        else:
            cfdnaspecific_clf = RandomizedSearchCV(estimator=cfdnaspecific_model, param_distributions=grid, n_iter=1, scoring=['average_precision', 'recall'], refit='average_precision', n_jobs=8, verbose=2)
            print(y.value_counts())
            search = cfdnaspecific_clf.fit(X, y)
            bestparams = search.best_params_
            print(bestparams)
            pd.Series(bestparams).to_csv("data/featureanalysis/bestparams_cfdnaspecificmodel_"+method+"_986_150x.csv")
            cfdnaspecific_modeRFC = RandomForestClassifier(**bestparams, random_state=0)
            with open("data/featureanalysis/cfdnaspecificmodel_"+method+"_986_150x.pkl", "wb") as f:
                pickle.dump(cfdnaspecific_modeRFC, f)
        cv = StratifiedKFold(n_splits=10, shuffle = True, random_state=0) # prepare 10-fold CV
        print('is there nan? {}'.format(np.sum(np.isnan(X), axis=0)))
        print('is there only finite values? {}'.format(np.sum(~np.isfinite(X), axis=0)))
        if method == 'vardict':
            if 'SOR' in list(X.columns):
                X['SOR'] = X['SOR'].astype('float32')
                X['SOR'] = X['SOR'].apply(np.nan_to_num)
        X = X.apply(np.nan_to_num)
        y = y.apply(np.nan_to_num)
        print(list(X.columns))
        X, y = sklearn.utils.check_X_y(X, y)
        print('is there nan? {}'.format(np.sum(np.isnan(X), axis=0)))
        print('is there only finite values? {}'.format(np.sum(~np.isfinite(X), axis=0)))
        cv_results = cross_validate(cfdnaspecific_modeRFC, X, y, cv=cv, scoring=['roc_auc', 'average_precision', 'recall'], return_estimator=True, verbose=6)
        cfdnaspecific_modeRFC.fit(X, y)
        ypred = cfdnaspecific_modeRFC.predict_proba(X)
        plt.hist(ypred)
        #print(pd.Series(ypred).describe())

        # Prepare patient 986 for test
        #mixtureid = 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T'
        mixtureid = 'CRC-986_100215-CW-T_CRC-986_300316-CW-T'
        #mixtureid = 'CRC-123_310715-CW-T_CRC-123_121115-CW-T'
        fixedvar = 'coverage'
        if fixedvar == 'coverage':
            seriesorder = [(70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]
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
        #calltable_indel, aux = get_calltableseries(config, mixtureid, chroms, muttype='indel', filterparam=filterparam, reload=reload, save=save)
        #calltable_snp, aux = get_calltableseries(config, mixtureid, chroms, muttype='snp', filterparam=filterparam, reload=reload, save=save)
        # print(calltable_snv.shape, calltable_indel.shape, calltable_snp.shape)
        print(aux)
        calltables['snv'].append(calltable_snv)
        #calltables['indel'].append(calltable_indel)
        #calltables['snp'].append(calltable_snp)
        calltables['sampleid'] = mixtureid
        calltables['tf'] = np.unique([cn.split('_')[0] for cn in list(calltable_snv.columns)])[:-5].astype(float)
        calltables['snv'] = pd.concat(calltables['snv'])
        #calltables['indel'] = pd.concat(calltables['indel'])
        #calltables['snp'] = pd.concat(calltables['snp'])
        dilutionseries = aux.T[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]].T
        dilutionseries = pd.DataFrame(dilutionseries.iloc[0]).T
        print(dilutionseries)
        muttype = 'snv'
        refsample = 'undiluted'
        if muttype == 'snv':
            gtm = 5
        else:  # elif muttype == 'indel':
            gtm = 4
        print(max(aux['tf']))
        calltablesseries = generate_groundtruth(config, calltables[muttype], aux['tf'], ground_truth_method=gtm, muttype=muttype,
                                                matchedtissuepath=None, methods=['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan', 'varnet', 'abemus', 'sinvict'])

        aux = aux.loc[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]]
        aux = pd.DataFrame(aux.iloc[0]).T
        print(aux)
        tfl = list(aux['tf'].values)
        print(tfl)
        print(str(tfl[0])[:3])
        calltable_snv = calltable_snv[[c for c in list(calltable_snv.columns) if c.startswith(str(tfl[0])[:3]) or c in ['chrom', 'pos', 'ref', 'alt', 'type', 'sampletf', 'truth']]]
        calltablesseries = calltablesseries[[c for c in list(calltablesseries.columns) if c.startswith(str(tfl[0])[:3]) or c in ['chrom', 'pos', 'ref', 'alt', 'type', 'sampletf', 'truth']]]
        print(list(calltable_snv.columns))
        np.random.seed(1234)
        color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
        callmethod_snv_list = []
        c = 0
        for serie in [(70, 80), (50, 100)]: #, (30, 120), (20, 130), (10, 140), (5, 145)]:
            if c < 1:
                st, sn = serie
                for chrom in chroms:
                    callmethod_snv, _, _ = parse_caller_feature(
                        'data/mixtures/mixtures_chr'+chrom+'/mixtures_chr'+chrom+'_CRC-986_100215-CW-T_CRC-986_300316-CW-T/mixture_chr'+chrom+'_CRC-986_100215-CW-T_'+str(st)+'x_CRC-986_300316-CW-T_'+str(sn)+'x',
                        method, save=False)
                    if method in ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']:
                        callmethod_snv.drop(['ANN'], axis=1, inplace=True)
                    callmethod_snv_list.append(callmethod_snv)
                c += 1
        callmethod_snv_allchrom = pd.concat(callmethod_snv_list)
        callmethod_snv_allchrom = callmethod_snv_allchrom.join(calltablesseries[['truth']], how='outer')
        callmethod_snv_allchrom['truth'].fillna(False, inplace=True)
        #features = list(callmethod_snv_allchrom.columns)
        #print(features)
        #for x in [method, method+'_score',  method+'_vaf', 'type', 'alt', 'chrom', 'pos', 'ref', method+'_totcov', method+'_altcov']:
        #    print(x)
        #    features.remove(x)
        #rint(features)
        #callmethod_snv_allchrom['truth'].fillna(False, inplace=True)
        #callmethod_snv_allchrom[method][callmethod_snv_allchrom[method] == 'PASS'] = 1
        #callmethod_snv_allchrom[method][callmethod_snv_allchrom[method] != 1] = 0
        print(callmethod_snv_allchrom['truth'].value_counts())
        # Plot PR curve naive vs cfDNA-specific on test patient 1014
        figax = plt.subplots(figsize=(12,12))
        ytrue = callmethod_snv_allchrom['truth']
        yscores = callmethod_snv_allchrom[method+'_score'].fillna(0)
        precision, recall, thresholds = precision_recall_curve(ytrue, yscores)
        print(average_precision_score(ytrue, yscores))
        plot_pr_curve(precision, recall, estimator_name=None, f1_score=None, figax=figax, kwargs={'color': config.colors[config.methods.index(method)], 'lw': 2}, plot='all')
        plt.ylim([0, 1])
        metric = 'maxrecallatleast0_03precision'
        threshold = float(metric.split('_')[1].split('precision')[0])/100
        print(threshold)
        p = np.where(precision >= threshold, 1, 0)
        rlist = recall * p
        dfresnaive = pd.DataFrame.from_dict({'precision': precision, 'recall': recall, 'AUPRC': average_precision_score(ytrue, yscores), 'maxrecallatleast0_03precision': np.max(rlist)})
        dfresnaive.to_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_986_150x_naive.csv')
        print('end '+method)
        ytrue = callmethod_snv_allchrom['truth']
        X = callmethod_snv_allchrom[list(features_selected)].fillna(0)
        if method == 'freebayes':
            if 'AC' in X.columns:
                X['AC'].loc[X['AC'].astype(str).str.contains(',') == True] = X['AC'].loc[X['AC'].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
                X['AC'] = X['AC'].astype(float)
        if method == 'mutect2':
            X['TLOD'].loc[X['TLOD'].astype(str).str.contains(',') == True] = X['TLOD'].loc[X['TLOD'].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
            X['TLOD'] = X['TLOD'].astype(float)
            for feat in ['NALOD', 'AS_SB_TABLE', 'FS', 'NLOD', 'MQRankSum']:
                if feat in X.columns:
                    X[feat].loc[X[feat].astype(str).str.contains(',') == True] = X[feat].loc[X[feat].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
                    X[feat] = X[feat].astype(float)
        yscores = cfdnaspecific_modeRFC.predict(X)
        precision, recall, thresholds = precision_recall_curve(ytrue, yscores)
        print(average_precision_score(ytrue, yscores))
        plot_pr_curve(precision, recall, estimator_name=None, f1_score=None, figax=figax, kwargs={'color': config.colors[config.methods.index(method)], 'lw': 5}, plot='all')
        plt.ylim([0, 1])
        plt.semilogx()
        plt.xlim([0.01, 1])
        plt.title(method)
        plt.show()
        metric = 'maxrecallatleast0_03precision'
        threshold = float(metric.split('_')[1].split('precision')[0])/100
        print(threshold)
        p = np.where(precision >= threshold, 1, 0)
        rlist = recall * p
        dfrescfdna = pd.DataFrame.from_dict({'precision': precision, 'recall': recall, 'AUPRC': average_precision_score(ytrue, yscores), 'maxrecallatleast0_03precision': np.max(rlist)})
        dfrescfdna.to_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_986_150x_cfdnaspecific.csv')


    # plot barplot feature importance
    mixtureid = 'CRC-123_310715-CW-T_CRC-123_121115-CW-T' # last patient used for training
    patientid = str(mixtureid.split('_')[0].split('-')[1])
    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(12, 14))
    for mi, method in enumerate([m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']):
        print(method)
        plt.subplot(3, 2, mi+1)
        feature_importances = pd.read_csv('data/featureanalysis/barplot_featureimportance_'+method+'_'+patientid+'_150x_RF.csv', index_col = 0)
        cvresults = pd.read_csv('data/featureanalysis/cvresults_featureimportance_'+method+'_'+patientid+'_150x_RF.csv', index_col = 0)
        g = sns.barplot(data=feature_importances, estimator=np.median, color=color_dict[method])
        g.set_xticklabels(g.get_xticklabels(), rotation=90)
        plt.xlabel('Features')
        plt.ylabel('Feature importance')
        plt.title('{} \n CV ROC AUC = {:.2f} +/- {:.2f}'.format(method, np.mean(cvresults['test_roc_auc']), np.std(cvresults['test_roc_auc'])))
        #plt.title(method + ' feature importance analysis. CV ROC AUC (acc {:.2f} +/- {:.2f})'.format(np.mean(cv_results['test_roc_auc']), np.std(cv_results['test_roc_auc'])))
    plt.savefig('data/featureanalysis/barplot_featureimportance_all_'+patientid+'_150x.svg')
    plt.savefig('data/featureanalysis/barplot_featureimportance_all_'+patientid+'_150x.png')
    plt.show()

    mixtureid = 'CRC-986_100215-CW-T_CRC-986_300316-CW-T'
    patientid = str(mixtureid.split('_')[0].split('-')[1])

    # plot PR curve overlay
    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
    figax = plt.subplots(figsize=(8, 8))
    for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:
        print(method)
        dfresnaive = pd.read_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_naive.csv', index_col=0)
        plot_pr_curve(dfresnaive['precision'], dfresnaive['recall'], estimator_name=None, f1_score=None,
                      figax=figax, kwargs={'color': config.colors[config.methods.index(method)], 'lw': 1}, plot='all')
        plt.ylim([0, 1])
        dfresnaive = pd.read_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_cfdnaspecific.csv', index_col=0)
        plot_pr_curve(dfresnaive['precision'], dfresnaive['recall'], estimator_name=method, f1_score=None,
                      figax=figax, kwargs={'color': config.colors[config.methods.index(method)], 'lw': 3}, plot='all')
        plt.ylim([0, 1])
        plt.semilogx()
        #plt.semilogy()
        plt.xlim([0.01, 1])
        plt.ylim([0.01, 1])
    plt.savefig('data/featureanalysis/PRcuve_all_'+patientid+'_150x.svg')
    plt.savefig('data/featureanalysis/PRcuve_all_'+patientid+'_150x.png')
    plt.show()


    # plot PR curve overlay
    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
    dfresnaivecfdna = pd.DataFrame(columns = ['method', 'model', 'AUPRC', 'maxrecallatleast0_03precision'])
    for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:
        print(method)
        dfresnaive = pd.read_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_naive.csv', index_col=0)
        dfresnaivecfdna = dfresnaivecfdna.append({'method': method, 'model': 'naive model', 'AUPRC': dfresnaive['AUPRC'].iloc[0], 'maxrecallatleast0_03precision': dfresnaive['maxrecallatleast0_03precision'].iloc[0]}, ignore_index=True)
        dfrescfdna = pd.read_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_cfdnaspecific.csv', index_col=0)
        dfresnaivecfdna = dfresnaivecfdna.append({'method': method, 'model': 'cfDNA-tuned model', 'AUPRC': dfrescfdna['AUPRC'].iloc[0], 'maxrecallatleast0_03precision': dfrescfdna['maxrecallatleast0_03precision'].iloc[0]}, ignore_index=True)
    print(dfresnaivecfdna)
    print(dfresnaivecfdna.shape)
    sns.catplot(data=dfresnaivecfdna, x="model", y="AUPRC", hue="method", capsize=.2, palette=color_dict, kind="point", height=6, aspect=.75)
    plt.savefig('data/featureanalysis/AUPRCgain_all_'+patientid+'_150x.svg')
    plt.savefig('data/featureanalysis/AUPRCgain_all_'+patientid+'_150x.png')
    sns.catplot(data=dfresnaivecfdna, x="model", y="maxrecallatleast0_03precision", hue="method", capsize=.2, palette=color_dict, kind="point", height=6, aspect=.75)
    plt.savefig('data/featureanalysis/maxrecallatleast0_03precisiongain_all_'+patientid+'_150x.svg')
    plt.savefig('data/featureanalysis/maxrecallatleast0_03precisiongain_all_'+patientid+'_150x.png')
    plt.show()

