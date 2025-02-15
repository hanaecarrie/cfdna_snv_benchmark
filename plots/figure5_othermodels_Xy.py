
if __name__ == "__main__":

    from utils.config import Config
    from featureanalysis.featureanalysis import parse_caller_feature
    from utils.config import Config
    from utils.viz import *
    from benchmark.table import *
    from benchmark.metrics import *
    from benchmark.calltable import *
    from benchmark.calltableseries import *
    from benchmark.groundtruth import *
    from benchmark.metricsseries import *

    import sklearn
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.model_selection import StratifiedKFold
    from sklearn.model_selection import cross_validate
    from sklearn.metrics import average_precision_score
    from sklearn.linear_model import LogisticRegression
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.neural_network import MLPClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    import os
    import pickle
    import random
    import numpy as np

    import seaborn as sns
    import matplotlib.pyplot as plt

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")

    mixtureids = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T'] # , 'BRA-412_240820-CW-T_BRA-412_060220-CW-T']
    mixtureids_train = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T',  'CRC-123_310715-CW-T_CRC-123_121115-CW-T'] # , 'BRA-412_240820-CW-T_BRA-412_060220-CW-T']
    mixtureid_test = 'CRC-986_100215-CW-T_CRC-986_300316-CW-T'

    mixture_sets = {'mixtureids_train': [['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T'],
                                         ['CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T'],
                                         ['CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T']],
                    'mixtureid_test': ['CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T']
                    }

    classifiers = {
                   #'LogisticRegressionLasso': Pipeline([('scaler', StandardScaler()), ('lasso', LogisticRegression(penalty='l1', solver='liblinear', random_state=0, verbose=True, n_jobs=12))], verbose=True),  #'LogisticRegressionLasso': LogisticRegression(penalty='l1', solver='liblinear', random_state=0),
                   #'LogisticRegressionRidge': Pipeline([('scaler', StandardScaler()), ('ridge', LogisticRegression(penalty='l2', solver='liblinear', random_state=0, verbose=True, n_jobs=12))], verbose=True),
                   #'RandomForest': RandomForestClassifier(n_estimators=100, random_state=0),
                   'LogisticRegressionElasticnet': Pipeline([('scaler', StandardScaler()), ('elasticnet', LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, random_state=0, n_jobs=12))], verbose=True), #'LogisticRegressionElasticnet': LogisticRegression(penalty='elasticnet', solver='saga', random_state=0),
                   #'DecisionTree': DecisionTreeClassifier(random_state=0),
                   'GradientBoosting': GradientBoostingClassifier(random_state=0),
                   #'MultiLayerPerceptron': Pipeline([('scaler', StandardScaler()), ('mlp', MLPClassifier(random_state=0))], verbose=True),  # 'MultiLayerPerceptron': MLPClassifier(random_state=0),
                   }

    # parameters
    reload = False
    save = False
    fixedvars=['coverage', 'ctdna']
    filterparam = 'all'
    markers = ['o', '^', 'X']
    linestyles = ['-', '-', '-']
    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
    metric = 'average_precision'
    chrom = 'all'
    muttype = 'snv' # indel not supported
    fixedvar = 'coverage'
    if fixedvar == 'coverage':
        seriesorder = [(70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]
        xaxis = 'tumor burden'
    elif fixedvar == 'ctdna':
        seriesorder = [(70, 0), (70, 80), (70, 180)]
        xaxis = 'coverage'
    refsample = 'undiluted'

    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}

    for msi in range(3):
        mixtureids_train = mixture_sets['mixtureids_train'][msi]
        mixtureid_test = mixture_sets['mixtureid_test'][msi]
        print(mixtureids_train)
        print(mixtureid_test)
        patientid = str(mixtureid_test.split('_')[0].split('-')[1])

        for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:
            print(method)

            for classifiername, classifiermodel in classifiers.items():
                print('############ '+ classifiername + ' ##########')

                if not os.path.exists('data/featureanalysis/othermodels/PRAUPRC_featureimportance_'+classifiername+'_'+method+'_' + patientid + '_150x_cfdnaspecific.csv'):

                    if not os.path.exists('data/featureanalysis/othermodels/Xtrain_'+method+'_without'+patientid+'.csv'):  # no Xtrain ytrain
                        X, y = pd.DataFrame(), pd.Series()
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
                            calltablesseries, aux = get_calltableseries(config, mixtureid, chroms, muttype='snv', filterparam=filterparam, reload=reload, save=save)
                            calltables['snv'].append(calltablesseries)
                            calltables['sampleid'] = mixtureid
                            calltables['tf'] = np.unique([cn.split('_')[0] for cn in list(calltablesseries.columns)])[:-5].astype(float)
                            calltables['snv'] = pd.concat(calltables['snv'])
                            dilutionseries = aux.T[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]].T
                            dilutionseries = pd.DataFrame(dilutionseries.iloc[0]).T
                            gtm = 5  # muttype = 'snv'
                            calltablesseries = generate_groundtruth(config, calltables[muttype], aux['tf'], ground_truth_method=gtm, muttype=muttype,
                                                                    methods=['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan', 'varnet', 'abemus', 'sinvict'])
                            print(calltablesseries['truth'].value_counts())
                            aux = aux.loc[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]]
                            print(aux)
                            calltablesseries = calltablesseries[[c for c in list(calltablesseries.columns) if c in ['chrom', 'pos', 'ref', 'alt', 'type', 'sampletf', 'truth']]]
                            print(calltablesseries)

                            # Extract features
                            np.random.seed(1234)
                            callmethod_snv_all_list = []
                            for serie in [(70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]:
                                st, sn = serie
                                callmethod_snv_list = []
                                for chrom in chroms:
                                    callmethod_snv, _, _ = parse_caller_feature(
                                        'data/mixtures/mixtures_chr'+chrom+'/mixtures_chr'+chrom+'_'+mixtureid + '/mixture_chr'+chrom+'_'+'_'.join(mixtureid.split('_')[:2])+'_'+str(st)+'x_'+'_'.join(mixtureid.split('_')[2:])+'_'+str(sn)+'x',
                                        method, save=False)
                                    if method in ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']:
                                        callmethod_snv.drop(['ANN'], axis=1, inplace=True) # annotation added by bcbio
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
                            features = list(callmethod_snv_allchrom.columns)

                            # Remove not applicable features
                            print(features)
                            for x in [method, method+'_score', method+'_vaf', 'type', 'alt', 'chrom', 'pos', 'ref', method+'_totcov', method+'_altcov']: #  let method_score in features
                                print(x)
                                features.remove(x)
                            print(features)
                            callmethod_snv_allchrom['truth'].fillna(False, inplace=True)
                            print(callmethod_snv_allchrom['truth'].value_counts())
                            forbiden_features = ['AF', 'OLD_VARIANT', 'TYPE', "cosmid_id", 'AN', "NS", 'CIGAR', 'DP', 'AF', 'AO', 'VD', 'N_refDepth', 'N_altDepth', 'relcov', 'END_POS_REF', 'REF_MFVdVs', 'ALT_MFVdVs', 'Sample_Name', 'cov', 'tf']
                            score_features = ['ODDS', 'TLOD', 'SomaticEVS', 'SSF', 'SPV']
                            Xlist = []
                            for feature in features:
                                if feature not in forbiden_features:
                                    test = 0
                                    callmethod_snv_allchrom[feature][callmethod_snv_allchrom[feature] == '.'] = np.nan
                                    try:
                                        aux = callmethod_snv_allchrom[feature].fillna(0).astype(float)
                                        test = 1
                                    except:
                                        try:
                                            aux = callmethod_snv_allchrom[feature].astype(str).str.split(',').str[0].fillna(0).astype(float)
                                            test = 1
                                        except:
                                            test = 0
                                    if test == 1 and feature != 'truth' and feature not in forbiden_features:
                                        Xlist.append(aux)
                                        print(feature)
                                    else:
                                        print(feature, 'not included')
                                else:
                                    print(feature, 'not included because forbidden')
                            Xtmp = pd.concat(Xlist, axis=1)
                            ytmp = callmethod_snv_allchrom['truth'].fillna(0)
                            print(Xtmp.shape, ytmp.shape)
                            if X.shape[0] == 0:
                                X, y = Xtmp, ytmp
                            else:
                                X, y = pd.concat([X, Xtmp]), pd.concat([y, ytmp])
                            if method == 'mutect2':
                                for f in ['N_ART_LOD', 'POP_AF', 'P_GERMLINE']:
                                    if f in list(X.columns):
                                        print(f, 'weird feature in mutect2')
                                        X.drop(f, axis=1, inplace=True)
                        print('mixtures_train size: X shape = {}, y shape = {}'.format(X.shape, y.shape))
                        print(list(X.columns))
                        listfeatures = X.columns
                        X.fillna(0, inplace=True)
                        if method == 'vardict':
                            X['SOR'] = X['SOR'].astype('float32')
                            mask = X['SOR'] != np.inf
                            X.loc[~mask, 'SOR'] = X.loc[mask, 'SOR'].max()
                        """
                        features_selected = X.columns
                        
                        # X = callmethod_snv_allchrom[list(features_selected)].fillna(0)
                        X = Xsave[list(features_selected)].fillna(0)
                        print(X.shape)
                        if method == 'freebayes':
                            for col in X.columns:
                                if not X[col].loc[X[col].astype(str).str.contains(',') == True].empty:
                                    print("################    " + col)
                                    X[col].loc[X[col].astype(str).str.contains(',') == True] = X[col].loc[X[col].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
                                    X[col] = X[col].astype(float)
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
                        if method == 'vardict':
                            X['SOR'] = X['SOR'].astype('float32')
                            mask = X['SOR'] != np.inf
                            X.loc[~mask, 'SOR'] = X.loc[mask, 'SOR'].max()
                        print(X.isna().sum())
                        print(X.max())
                        print(X.min())
                        # y = callmethod_snv_allchrom['truth']
                        y = ysave
                        print(y.value_counts())
                        print(X.shape, y.shape)
                        # check comma in X
                        for col in list(X.columns):
                            print(col)
                            print(X[col])
                            print(X[col].loc[X[col].astype(str).str.contains(',') == True])
                            print(X[col].loc[X[col].astype(str).str.contains(',') == True].shape)
                        """
                        # save X_train, y_train
                        mixtureid = mixtureid_test
                        patientid = str(mixtureid_test.split('_')[0].split('-')[1])
                        X.to_csv('data/featureanalysis/othermodels/Xtrain_'+method+'_without'+patientid+'.csv')
                        y.to_csv('data/featureanalysis/othermodels/ytrain_'+method+'_without'+patientid+'.csv')

                    # train model
                    mixtureid = mixtureid_test
                    patientid = str(mixtureid_test.split('_')[0].split('-')[1])
                    print('############ start model training ##########')
                    X = pd.read_csv('data/featureanalysis/othermodels/Xtrain_'+method+'_without'+patientid+'.csv', index_col=0, memory_map=True)
                    y = pd.read_csv('data/featureanalysis/othermodels/ytrain_'+method+'_without'+patientid+'.csv', index_col=0, memory_map=True)
                    print(X.shape)

                    print('############ donwsample training ##########')
                    tenpercentXshape = np.int(np.floor(X.shape[0]/10))
                    print(tenpercentXshape)
                    # sample only 10%
                    indexsample = random.sample(range(tenpercentXshape), tenpercentXshape)
                    X = X.iloc[indexsample]
                    y = y.iloc[indexsample]
                    print(X.shape, y.shape)
                    print('############ start model training ##########')
                    cfdnaspecific_modeRFC = classifiermodel
                    cfdnaspecific_modeRFC.fit(X, y)
                    features_train = X.columns
                    print(len(features_train))
                    print(features_train)
                    print('############ end model training ##########')

                    # TEST
                    if not os.path.exists('data/featureanalysis/othermodels/Xtest_'+method+'_'+patientid+'.csv'):  # no Xtrain ytrain

                        X, y = pd.DataFrame(), pd.Series()
                        patientid = str(mixtureid.split('_')[0].split('-')[1])
                        if mixtureid == 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T':
                            chroms = [str(c) for c in range(1,23) if c != 17 and c != 8 and c != 5 and c != 19 and c != 20 and c != 21 and c != 22]
                        elif mixtureid == 'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
                            chroms = [str(c) for c in range(1,23) if c != 1 and c != 2 and c != 8 and c != 20 and c != 21 and c != 22]
                        else:
                            chroms = [str(c) for c in range(1,23) if c != 1 and c != 2 and c != 6 and c != 20 and c != 21]
                        calltables = {'sampleid': [], 'tf': [], 'cov': [], 'snv': [], 'indel': [], 'snp': []}
                        aux_all = []
                        calltablesseries, aux = get_calltableseries(config, mixtureid, chroms, muttype='snv', filterparam=filterparam, reload=reload, save=save)
                        calltables['snv'].append(calltablesseries)
                        calltables['sampleid'] = mixtureid
                        calltables['tf'] = np.unique([cn.split('_')[0] for cn in list(calltablesseries.columns)])[:-5].astype(float)
                        calltables['snv'] = pd.concat(calltables['snv'])
                        dilutionseries = aux.T[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]].T
                        dilutionseries = pd.DataFrame(dilutionseries.iloc[0]).T
                        gtm = 5  # muttype = 'snv'
                        calltablesseries = generate_groundtruth(config, calltables[muttype], aux['tf'], ground_truth_method=gtm, muttype=muttype,
                                                                methods=['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan', 'varnet', 'abemus', 'sinvict'])
                        print(calltablesseries['truth'].value_counts())
                        aux = aux.loc[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]]
                        print(aux)
                        calltablesseries = calltablesseries[[c for c in list(calltablesseries.columns) if c in ['chrom', 'pos', 'ref', 'alt', 'type', 'sampletf', 'truth']]]
                        print(calltablesseries)

                        # Extract features
                        np.random.seed(1234)
                        callmethod_snv_all_list = []
                        for serie in [(70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]:
                            st, sn = serie
                            callmethod_snv_list = []
                            for chrom in chroms:
                                callmethod_snv, _, _ = parse_caller_feature(
                                    'data/mixtures/mixtures_chr'+chrom+'/mixtures_chr'+chrom+'_'+mixtureid + '/mixture_chr'+chrom+'_'+'_'.join(mixtureid.split('_')[:2])+'_'+str(st)+'x_'+'_'.join(mixtureid.split('_')[2:])+'_'+str(sn)+'x',
                                    method, save=False)
                                if method in ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']:
                                    callmethod_snv.drop(['ANN'], axis=1, inplace=True) # annotation added by bcbio
                                callmethod_snv_list.append(callmethod_snv)
                            callmethod_snv_all = pd.concat(callmethod_snv_list)
                            callmethod_snv_all = pd.concat([callmethod_snv_all, calltablesseries[['truth']]], axis=1)
                            print(callmethod_snv_all['truth'].value_counts())
                            callmethod_snv_all.index = [idx + '_' +str(st)+'x'+str(sn)+'x' for idx in list(callmethod_snv_all.index)]
                            callmethod_snv_all_list.append(callmethod_snv_all)
                        callmethod_snv_allchrom = pd.concat(callmethod_snv_all_list)
                        print(callmethod_snv_allchrom['truth'].value_counts())
                        callmethod_snv_allchrom['truth'].fillna(False, inplace=True)

                        # Keep only input features
                        X = callmethod_snv_allchrom[features_train]
                        for feature in features_train:
                            X[feature][X[feature] == '.'] = np.nan
                            try:
                                X[feature] = X[feature].fillna(0).astype(float)
                            except:
                                try:
                                    X[feature] = X[feature].astype(str).str.split(',').str[0].fillna(0).astype(float)
                                except:
                                    print('not included in different cases')
                        X.fillna(0, inplace=True)
                        if method == 'vardict':
                            X['SOR'] = X['SOR'].astype('float32')
                            mask = X['SOR'] != np.inf
                            X.loc[~mask, 'SOR'] = X.loc[mask, 'SOR'].max()
                        y = callmethod_snv_allchrom['truth'].fillna(0)
                        ydefault = callmethod_snv_allchrom[method+'_score'].fillna(0)
                        print(X.shape, y.shape, ydefault.shape)
                        print('mixtures_test size: X shape = {}, y shape = {}'.format(X.shape, y.shape))

                        # save X test and y true test
                        X.to_csv('data/featureanalysis/othermodels/Xtest_'+method+'_'+patientid+'.csv')
                        y.to_csv('data/featureanalysis/othermodels/ytrue_'+method+'_'+patientid+'.csv')
                        ydefault.to_csv('data/featureanalysis/othermodels/ydefault_'+method+'_'+patientid+'.csv')

                    # once Xtest and ytrue defined -> inference
                    print('############ load test set ##########')
                    X = pd.read_csv('data/featureanalysis/othermodels/Xtest_'+method+'_'+patientid+'.csv', index_col=0, memory_map=True)
                    ytrue = pd.read_csv('data/featureanalysis/othermodels/ytrue_'+method+'_'+patientid+'.csv', index_col=0, memory_map=True)
                    ydefault = pd.read_csv('data/featureanalysis/othermodels/ydefault_'+method+'_'+patientid+'.csv', index_col=0, memory_map=True)
                    print(X.shape)

                    figax = plt.subplots(figsize=(12,12))

                    # naive model
                    precision, recall, thresholds = precision_recall_curve(ytrue, ydefault)
                    print(average_precision_score(ytrue, ydefault))
                    metricscorenaive = average_precision_score(ytrue, ydefault)
                    plot_pr_curve(precision, recall, estimator_name='naive ' + method, f1_score=None, figax=figax, kwargs={'color': config.colors[config.methods.index(method)], 'lw': 2}, plot='all')
                    plt.ylim([0, 1])
                    metric = 'maxrecallatleast0_03precision'
                    threshold = float(metric.split('_')[1].split('precision')[0])/100
                    print(threshold)
                    p = np.where(precision >= threshold, 1, 0)
                    rlist = recall * p
                    dfresnaive = pd.DataFrame.from_dict({'precision': precision, 'recall': recall, 'AUPRC': average_precision_score(ytrue, ydefault), 'maxrecallatleast0_03precision': np.max(rlist)})
                    dfresnaive.to_csv('data/featureanalysis/othermodels/PRAUPRC_featureimportance_'+classifiername+'_'+method+'_'+patientid+'_150x_naive.csv')

                    # cfDNA-specific model
                    print('############ infer predictions on test set ##########')
                    yscores = cfdnaspecific_modeRFC.predict_proba(X)[:, 1]  # take predictions for class 1
                    precision, recall, thresholds = precision_recall_curve(ytrue, yscores)
                    print(average_precision_score(ytrue, yscores))
                    metricscorespecific = average_precision_score(ytrue, yscores)
                    plot_pr_curve(precision, recall, estimator_name='cfDNA-specific '+method, f1_score=None, figax=figax, kwargs={'color': config.colors[config.methods.index(method)], 'lw': 5}, plot='all')
                    metric = 'maxrecallatleast0_03precision'
                    threshold = float(metric.split('_')[1].split('precision')[0])/100
                    print(threshold)
                    p = np.where(precision >= threshold, 1, 0)
                    rlist = recall * p
                    dfrescfdna = pd.DataFrame.from_dict({'precision': precision, 'recall': recall, 'AUPRC': average_precision_score(ytrue, yscores), 'maxrecallatleast0_03precision': np.max(rlist)})
                    dfrescfdna.to_csv('data/featureanalysis/othermodels/PRAUPRC_featureimportance_'+classifiername+'_'+method+'_'+patientid+'_150x_cfdnaspecific.csv')
                    print(metricscorespecific - metricscorenaive)

                    plt.ylim([0, 1])
                    plt.semilogx()
                    plt.xlim([0.01, 1])
                    plt.title(method)
                    plt.savefig('data/featureanalysis/othermodels/PRcuvres_naive_finetuned_'+classifiername+'_'+method+'_'+patientid+'.png', bbox_inches='tight')
                    plt.savefig('data/featureanalysis/othermodels/PRcuvres_naive_finetuned_'+classifiername+'_'+method+'_'+patientid+'.svg', bbox_inches='tight')


                    """
                    # Build cfDNA specific model
                    cfdnaspecific_modeRFC = classifiermodel
            
                    mixtureid = mixtureid_test
                    patientid = str(mixtureid_test.split('_')[0].split('-')[1])
                    
                    if not os.path.exists('data/featureanalysis/ytrue_'+method+'_'+patientid+'.csv'):
                        # Prepare patient for test
                        mixtureid = mixtureid_test
                        patientid = str(mixtureid_test.split('_')[0].split('-')[1])
                        fixedvar = 'coverage'
                        if fixedvar == 'coverage':
                            seriesorder = [(70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]
                            xaxis = 'tumor burden'
                        elif fixedvar == 'ctdna':
                            seriesorder = [(70, 0), (70, 80), (70, 180)]
                            xaxis = 'coverage'
                        print('############# {} ############'.format(mixtureid))
                        if mixtureid == 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T':
                            chroms = [str(c) for c in range(1,23) if c !=17 and c !=8 and c!=5 and c!=19 and c!=20 and c!=21 and c!=22]
                        elif mixtureid == 'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
                            chroms = [str(c) for c in range(1,23) if c !=1 and c!= 2 and c !=8 and c!=20 and c!=21 and c!=22]
                        else:
                            chroms = [str(c) for c in range(1,23) if c !=1 and c!= 2 and c !=6 and c!=20 and c!=21]
                        calltables = {'sampleid':[], 'tf':[], 'cov':[], 'snv':[], 'indel':[], 'snp':[]}
                        aux_all = []
                        calltable_snv, aux = get_calltableseries(config, mixtureid, chroms, muttype='snv', filterparam=filterparam, reload=reload, save=save)
                        calltables['snv'].append(calltable_snv)
                        calltables['sampleid'] = mixtureid
                        calltables['tf'] = np.unique([cn.split('_')[0] for cn in list(calltable_snv.columns)])[:-5].astype(float)
                        calltables['snv'] = pd.concat(calltables['snv'])
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
                        tfl = list(aux['tf'].values)
                        print(str(tfl[0])[:3])
                        calltable_snv = calltable_snv[[c for c in list(calltable_snv.columns) if c.startswith(str(tfl[0])[:3]) or c in ['chrom', 'pos', 'ref', 'alt', 'type', 'sampletf', 'truth']]]
                        calltablesseries = calltablesseries[[c for c in list(calltablesseries.columns) if c.startswith(str(tfl[0])[:3]) or c in ['chrom', 'pos', 'ref', 'alt', 'type', 'sampletf', 'truth']]]
                        print(list(calltable_snv.columns))
                        np.random.seed(1234)
                        color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
                        callmethod_snv_list = []
                        c = 0
                        for serie in [(70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]:
                            if c < 1:
                                st, sn = serie
                                for chrom in chroms:
                                    callmethod_snv, _, _ = parse_caller_feature(
                                        'data/mixtures/mixtures_chr'+chrom+'/mixtures_chr'+chrom+'_'+mixtureid_test+'/mixture_chr'+chrom+'_'+'_'.join(mixtureid_test.split('_')[:2])+'_'+str(st)+'x_'+'_'.join(mixtureid_test.split('_')[2:])+'_'+str(sn)+'x',
                                        method, save=False)
                                    if method in ['freebayes', 'mutect2', 'strelka2', 'vardict', 'varscan']:
                                        callmethod_snv.drop(['ANN'], axis=1, inplace=True)
                                    callmethod_snv_list.append(callmethod_snv)
                                c += 1
                        callmethod_snv_allchrom = pd.concat(callmethod_snv_list)
                        callmethod_snv_allchrom = callmethod_snv_allchrom.join(calltablesseries[['truth']], how='outer')
                        callmethod_snv_allchrom['truth'].fillna(False, inplace=True)
                        print(callmethod_snv_allchrom['truth'].value_counts())
                        # Plot PR curve naive vs cfDNA-specific on test patient 1014
                        figax = plt.subplots(figsize=(12,12))
                        ytrue = callmethod_snv_allchrom['truth']
                        yscores = callmethod_snv_allchrom[method+'_score'].fillna(0)
                        # save prediction y_test naive model
                        yscores.to_csv('data/featureanalysis/yscores_naive'+method+'_'+patientid+'.csv')
                        ytrue.to_csv('data/featureanalysis/ytrue_'+method+'_'+patientid+'.csv')
                        precision, recall, thresholds = precision_recall_curve(ytrue, yscores)
                        print(average_precision_score(ytrue, yscores))
                        metricscorenaive = average_precision_score(ytrue, yscores)
                        plot_pr_curve(precision, recall, estimator_name='naive '+ method, f1_score=None, figax=figax, kwargs={'color': config.colors[config.methods.index(method)], 'lw': 2}, plot='all')
                        plt.ylim([0, 1])
                        metric = 'maxrecallatleast0_03precision'
                        threshold = float(metric.split('_')[1].split('precision')[0])/100
                        print(threshold)
                        p = np.where(precision >= threshold, 1, 0)
                        rlist = recall * p
                        dfresnaive = pd.DataFrame.from_dict({'precision': precision, 'recall': recall, 'AUPRC': average_precision_score(ytrue, yscores), 'maxrecallatleast0_03precision': np.max(rlist)})
                        dfresnaive.to_csv('data/featureanalysis/PRAUPRC_featureimportance_'+classifiername+'_'+method+'_'+patientid+'_150x_naive.csv')
                        print('end '+method)
                        ytrue = callmethod_snv_allchrom['truth']
                        print(list(callmethod_snv_allchrom.columns))
                        X = callmethod_snv_allchrom[list(features_selected)].fillna(0)
                        for feature in features_selected:
                            X[feature][X[feature] == '.'] = np.nan
                            try:
                                X[feature] = X[feature].fillna(0).astype(float)
                            except:
                                X[feature] = X[feature].astype(str).str.split(',').str[0].fillna(0).astype(float)
                        if method == 'freebayes':
                            for col in X.columns:
                                if not X[col].loc[X[col].astype(str).str.contains(',') == True].empty:
                                    print("################    " + col)
                                    X[col].loc[X[col].astype(str).str.contains(',') == True] = X[col].loc[X[col].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
                                    X[col] = X[col].astype(float)
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
                        if method == 'vardict':
                            X['SOR'] = X['SOR'].astype('float32')
                            mask = X['SOR'] != np.inf
                            X.loc[~mask, 'SOR'] = X.loc[mask, 'SOR'].max()
            
                        # save X test and y test
                        X.to_csv('data/featureanalysis/othermodels/Xtest_'+method+'_'+patientid+'.csv')
                        ytrue.to_csv('data/featureanalysis/othermodels/ytrue_'+method+'_'+patientid+'.csv')
            
                    # once Xtest and ytrue defined -> inference
                    X = pd.read_csv('data/featureanalysis/othermodels/Xtest_'+method+'_'+patientid+'.csv', index_col=0)
                    ytrue = pd.read_csv('data/featureanalysis/othermodels/ytrue_'+method+'_'+patientid+'.csv', index_col=0)
                    yscores = cfdnaspecific_modeRFC.predict_proba(X)[:,1]  # take predictions for class 1
                    precision, recall, thresholds = precision_recall_curve(ytrue, yscores)
                    print(average_precision_score(ytrue, yscores))
                    metricscorespecific = average_precision_score(ytrue, yscores)
                    plot_pr_curve(precision, recall, estimator_name='cfDNA-specific '+method, f1_score=None, figax=figax, kwargs={'color': config.colors[config.methods.index(method)], 'lw': 5}, plot='all')
                    plt.ylim([0, 1])
                    plt.semilogx()
                    plt.xlim([0.01, 1])
                    plt.title(method)
                    plt.savefig('data/featureanalysis/othermodels/PRcuvres_naive_finetuned_'+classifiername+'_'+method+'_'+patientid+'.png', bbox_inches='tight')
                    plt.savefig('data/featureanalysis/othermodels/PRcuvres_naive_finetuned_'+classifiername+'_'+method+'_'+patientid+'.svg', bbox_inches='tight')
                    metric = 'maxrecallatleast0_03precision'
                    threshold = float(metric.split('_')[1].split('precision')[0])/100
                    print(threshold)
                    p = np.where(precision >= threshold, 1, 0)
                    rlist = recall * p
                    dfrescfdna = pd.DataFrame.from_dict({'precision': precision, 'recall': recall, 'AUPRC': average_precision_score(ytrue, yscores), 'maxrecallatleast0_03precision': np.max(rlist)})
                    dfrescfdna.to_csv('data/featureanalysis/othermodels/PRAUPRC_featureimportance_'+classifiername+'_'+method+'_'+patientid+'_150x_cfdnaspecific.csv')
                    print(metricscorespecific - metricscorenaive)
                """
