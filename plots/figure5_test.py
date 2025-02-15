
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
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import StratifiedKFold
    from sklearn.model_selection import cross_validate
    from sklearn.metrics import average_precision_score
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    import os
    import pickle

    import seaborn as sns
    import matplotlib.pyplot as plt

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")

    mixtureids = ['CRC-1014_180816-CW-T_CRC-1014_090516-CW-T', 'CRC-986_100215-CW-T_CRC-986_300316-CW-T', 'CRC-123_310715-CW-T_CRC-123_121115-CW-T'] # , 'BRA-412_240820-CW-T_BRA-412_060220-CW-T']
    mixtureids_train = ['CRC-123_310715-CW-T_CRC-123_121115-CW-T', 'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T' ] # , 'BRA-412_240820-CW-T_BRA-412_060220-CW-T']
    mixtureid_test = 'CRC-986_100215-CW-T_CRC-986_300316-CW-T'
    """"
    method = 'smurf'
    #for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:

    print(method)

    # Prepare patient 1014 and 123 to study feature importance
    reload = False
    save = False
    fixedvars=['coverage', 'ctdna']
    filterparam = 'all'
    markers = ['o', '^', 'X']
    linestyles = ['-', '-', '-']
    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
    muttypes = ['snv', 'indel']
    metric = 'average_precision'
    chrom = 'all'
    fixedvar = 'coverage'
    if fixedvar == 'coverage':
        seriesorder = [(70, 80), (50, 100), (30, 120), (20, 130), (10, 140), (5, 145)]
        xaxis = 'tumor burden'
    elif fixedvar == 'ctdna':
        seriesorder = [(70, 0), (70, 80), (70, 180)]
        xaxis = 'coverage'

    X, y = pd.DataFrame(), pd.Series()
    for mixtureid in mixtureids_train:
        #mixtureid = mixtureid_test
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
        calltables['snv'].append(calltable_snv)
        calltables['sampleid'] = mixtureid
        calltables['tf'] = np.unique([cn.split('_')[0] for cn in list(calltable_snv.columns)])[:-5].astype(float)
        calltables['snv'] = pd.concat(calltables['snv'])
        dilutionseries = aux.T[['mixture_' + '_'.join(mixtureid.split('_')[:2]) + '_' + str(s[0]) + 'x_' + '_'.join(mixtureid.split('_')[2:4]) + '_' + str(s[1]) + 'x' for s in seriesorder]].T
        dilutionseries = pd.DataFrame(dilutionseries.iloc[0]).T
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
        #serie = (70, 80)
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
        # end for serie
        callmethod_snv_allchrom = pd.concat(callmethod_snv_all_list)
        print(callmethod_snv_allchrom['truth'].value_counts())
        callmethod_snv_allchrom['truth'].fillna(False, inplace=True)
        print(callmethod_snv_allchrom.columns)
        features = list(callmethod_snv_allchrom.columns)
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
                #try:
                #    aux = callmethod_snv_allchrom[feature].fillna(0).astype(float)
                #    test = 1
                #except:
                #    try:
                #        aux = callmethod_snv_allchrom[feature].str.split(',').str[0].fillna(0).astype(float)
                #        test = 1
                #    except:
                #        test = 0
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
    # end for loop mixtureids_train
    print('mixtures_train size: X shape = {}, y shape = {}'.format(X.shape, y.shape))
    Xsave, ysave = X.copy(), y.copy()
    print('################')
    print(list(X.columns))
    print('################')

    features_selected = X.columns
    if method == 'mutect2':
        features_selected = [f for f in features_selected if f not in ['N_ART_LOD', 'POP_AF', 'P_GERMLINE']]
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

    # Build cfDNA specific model
    cfdnaspecific_modeRFC = RandomForestClassifier(n_estimators=100, random_state=0)
    #cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=0)  # prepare 10-fold CV
    print('is there nan? {}'.format(np.sum(np.isnan(X), axis=0)))
    print('is there only finite values? {}'.format(np.sum(~np.isfinite(X), axis=0)))
    print(type(X))
    X = X.apply(np.nan_to_num)
    y = y.apply(np.nan_to_num)
    print(list(X.columns))
    X, y = sklearn.utils.check_X_y(X, y)
    print('is there nan? {}'.format(np.sum(np.isnan(X), axis=0)))
    print('is there only finite values? {}'.format(np.sum(~np.isfinite(X), axis=0)))
    #cv_results = cross_validate(cfdnaspecific_modeRFC, X, y, cv=cv, scoring=metric, return_estimator=True, verbose=6, n_jobs=8)
    cfdnaspecific_modeRFC.fit(X, y)

    mixtureid = mixtureid_test
    patientid = str(mixtureid_test.split('_')[0].split('-')[1])

    # feature importance
    feature_importances = pd.Series(cfdnaspecific_modeRFC.feature_importances_, index=features_selected)
    print(feature_importances)
    #print(list(set((feature_importances[(feature_importances.T > 0.001).any()].index)) | set([l for l in list(feature_importances.index) if l in score_features])))
    #feature_importances = feature_importances.loc[list(set((feature_importances[(feature_importances.T > 0.001).any()].index)) | set([l for l in list(feature_importances.index) if l in score_features]))]
    #feature_importances['median'] = feature_importances.median(axis=1)
    #feature_importances.sort_values(by='median', ascending=False, inplace=True)
    feature_importances.sort_values(ascending=False, inplace=True)
    #feature_importances.drop('median', axis=1, inplace=True)
    plt.figure(figsize=(20, 10))
    feature_importances = feature_importances.T
    feature_importances[feature_importances != 0].dropna().plot(kind='bar', color=color_dict[method])
    #g = sns.barplot(data=feature_importances, color=color_dict[method])
    #g.set_xticklabels(g.get_xticklabels(), rotation=90)
    plt.xlabel('Features')
    plt.ylabel('Feature importance in cfDNA tuned model')
    plt.title(method + ' feature importance analysis')
    # plt.title(method + ' feature importance analysis. CV ROC AUC (acc {:.2f} +/- {:.2f})'.format(np.mean(cv_results['test_roc_auc']), np.std(cv_results['test_roc_auc'])))
    feature_importances.to_csv('data/featureanalysis/barplot_featureimportance_'+method+'_cfDNAtuned_leftoutpatient'+patientid+'_150x_RF.csv')
    plt.savefig('data/featureanalysis/featureimportance_'+method+'_trainpatients_150x.svg')
    plt.savefig('data/featureanalysis/featureimportance_'+method+'_trainpatients_150x.png')
    #plt.show()


    # Prepare patient 986 for test
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
    if mixtureid ==  'CRC-1014_180816-CW-T_CRC-1014_090516-CW-T':
        chroms = [str(c) for c in range(1,23) if c !=17 and c !=8 and c!=5 and c!=19 and c!=20 and c!=21 and c!=22]
    elif mixtureid ==  'CRC-986_100215-CW-T_CRC-986_300316-CW-T':
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
    #serie = (70, 80)
    #serie = (70, 80)
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
    # serie for loop en
    callmethod_snv_allchrom = pd.concat(callmethod_snv_list)
    callmethod_snv_allchrom = callmethod_snv_allchrom.join(calltablesseries[['truth']], how='outer')
    callmethod_snv_allchrom['truth'].fillna(False, inplace=True)
    print(callmethod_snv_allchrom['truth'].value_counts())
    # Plot PR curve naive vs cfDNA-specific on test patient 1014
    figax = plt.subplots(figsize=(12,12))
    ytrue = callmethod_snv_allchrom['truth']
    yscores = callmethod_snv_allchrom[method+'_score'].fillna(0)
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
    dfresnaive.to_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_naive.csv')
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
    # if method == 'freebayes':
    #     if 'AC' in X.columns:
    #         X['AC'].loc[X['AC'].astype(str).str.contains(',') == True] = X['AC'].loc[X['AC'].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
    #         X['AC'] = X['AC'].astype(float)
    # if method == 'mutect2':
    #     X['TLOD'].loc[X['TLOD'].astype(str).str.contains(',') == True] = X['TLOD'].loc[X['TLOD'].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
    #     X['TLOD'] = X['TLOD'].astype(float)
    #     for feat in ['NALOD', 'AS_SB_TABLE', 'FS', 'NLOD', 'MQRankSum']:
    #         if feat in X.columns:
    #             X[feat].loc[X[feat].astype(str).str.contains(',') == True] = X[feat].loc[X[feat].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
    #             X[feat] = X[feat].astype(float)
    # if method == 'freebayes' or method == 'mutect2':
    #     for col in X.columns:
    #         if not X[col].loc[X[col].astype(str).str.contains(',') == True].empty:
    #             print("################    " + col)
    #             X[col].loc[X[col].astype(str).str.contains(',') == True] = X[col].loc[X[col].astype(str).str.contains(',') == True].astype(str).str.split(',').str[0]
    #             X[col] = X[col].astype(float)
    # if method == 'vardict':
    #     X['SOR'] = X['SOR'].astype('float32')
    #     mask = X['SOR'] != np.inf
    #     X.loc[~mask, 'SOR'] = X.loc[mask, 'SOR'].max()
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
    yscores = cfdnaspecific_modeRFC.predict_proba(X)[:,1]  # take predictions for class 1
    precision, recall, thresholds = precision_recall_curve(ytrue, yscores)
    print(average_precision_score(ytrue, yscores))
    metricscorespecific = average_precision_score(ytrue, yscores)
    plot_pr_curve(precision, recall, estimator_name='cfDNA-specific '+method, f1_score=None, figax=figax, kwargs={'color': config.colors[config.methods.index(method)], 'lw': 5}, plot='all')
    plt.ylim([0, 1])
    plt.semilogx()
    plt.xlim([0.01, 1])
    plt.title(method)
    plt.savefig('data/featureanalysis/PRcuvres_naive_finetuned_'+method+'_'+patientid+'.png', bbox_inches='tight')
    plt.savefig('data/featureanalysis/PRcuvres_naive_finetuned_'+method+'_'+patientid+'.svg', bbox_inches='tight')
    metric = 'maxrecallatleast0_03precision'
    threshold = float(metric.split('_')[1].split('precision')[0])/100
    print(threshold)
    p = np.where(precision >= threshold, 1, 0)
    rlist = recall * p
    dfrescfdna = pd.DataFrame.from_dict({'precision': precision, 'recall': recall, 'AUPRC': average_precision_score(ytrue, yscores), 'maxrecallatleast0_03precision': np.max(rlist)})
    dfrescfdna.to_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_cfdnaspecific.csv')
    print(metricscorespecific - metricscorenaive)
    """


    # # summary for one test patient
    # patientid = str(mixtureid_test.split('_')[0].split('-')[1])
    #
    # # plot PR curve overlay
    # color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
    # figax = plt.subplots(figsize=(8, 8))
    # for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:
    #     print(method)
    #     dfresnaive = pd.read_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_naive.csv', index_col=0)
    #     plot_pr_curve(dfresnaive['precision'], dfresnaive['recall'], estimator_name=None, f1_score=None,
    #                   figax=figax, kwargs={'color': config.colors[config.methods.index(method)], 'lw': 1}, plot='all')
    #     plt.ylim([0, 1])
    #     dfresnaive = pd.read_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_cfdnaspecific.csv', index_col=0)
    #     plot_pr_curve(dfresnaive['precision'], dfresnaive['recall'], estimator_name=method, f1_score=None,
    #                   figax=figax, kwargs={'color': config.colors[config.methods.index(method)], 'lw': 3}, plot='all')
    #     plt.ylim([0, 1])
    #     plt.semilogx()
    #     #plt.semilogy()
    #     plt.xlim([0.01, 1])
    #     plt.ylim([0.01, 1])
    # plt.savefig('data/featureanalysis/PRcuve_all_'+patientid+'_150x.svg')
    # plt.savefig('data/featureanalysis/PRcuve_all_'+patientid+'_150x.png')
    #
    # patientid = str(mixtureid_test.split('_')[0].split('-')[1])
    # # plot PR curve overlay
    # color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
    # dfresnaivecfdna = pd.DataFrame(columns = ['method', 'model', 'AUPRC', 'maxrecallatleast0_03precision'])
    # for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:
    #     print(method)
    #     dfresnaive = pd.read_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_naive.csv', index_col=0)
    #     dfresnaivecfdna = dfresnaivecfdna.append({'method': method, 'model': 'naive model', 'AUPRC': dfresnaive['AUPRC'].iloc[0], 'maxrecallatleast0_03precision': dfresnaive['maxrecallatleast0_03precision'].iloc[0]}, ignore_index=True)
    #     dfrescfdna = pd.read_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_cfdnaspecific.csv', index_col=0)
    #     dfresnaivecfdna = dfresnaivecfdna.append({'method': method, 'model': 'cfDNA-tuned model', 'AUPRC': dfrescfdna['AUPRC'].iloc[0], 'maxrecallatleast0_03precision': dfrescfdna['maxrecallatleast0_03precision'].iloc[0]}, ignore_index=True)
    # print(dfresnaivecfdna)
    # print(dfresnaivecfdna.shape)
    # sns.catplot(data=dfresnaivecfdna, x="model", y="AUPRC", hue="method", capsize=.2, palette=color_dict, kind="point", height=6, aspect=.75)
    # plt.savefig('data/featureanalysis/AUPRCgain_all_'+patientid+'_150x.svg')
    # plt.savefig('data/featureanalysis/AUPRCgain_all_'+patientid+'_150x.png')
    # sns.catplot(data=dfresnaivecfdna, x="model", y="maxrecallatleast0_03precision", hue="method", capsize=.2, palette=color_dict, kind="point", height=6, aspect=.75)
    # plt.savefig('data/featureanalysis/maxrecallatleast0_03precisiongain_all_'+patientid+'_150x.svg')
    # plt.savefig('data/featureanalysis/maxrecallatleast0_03precisiongain_all_'+patientid+'_150x.png')
    # plt.show()

    color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}

    
    # # plot PR curve overlay
    # for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:
    #     method = 'smurf'
    #     figax = plt.subplots(figsize=(8, 8))
    #     print(method)
    #     dfresnaiveall = pd.DataFrame()
    #     dfresspecificall = pd.DataFrame()
    #     for patientid in ['1014', '123', '986']:
    #         dfresnaive = pd.read_csv('data/featureanalysis/testpatient'+patientid+'/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_naive.csv', index_col=0)
    #         dfresnaive['patientid'] = patientid
    #         dfresspecific = pd.read_csv('data/featureanalysis/testpatient'+patientid+'PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_cfdnaspecific.csv', index_col=0)
    #         dfresspecific['patientid'] = patientid
    #         dfresnaiveall = pd.concat([dfresnaiveall, dfresnaive])
    #         dfresspecificall = pd.concat([dfresspecificall, dfresspecific])
    #
    #     npoints = 150
    #
    #     a1 = dfresnaiveall[dfresnaiveall['patientid'] == '1014']
    #     a2 = dfresnaiveall.loc[dfresnaiveall['patientid'] == '986']
    #     a3 = dfresnaiveall.loc[dfresnaiveall['patientid'] == '123']
    #
    #     min_a1_x, max_a1_x = min(a1['recall']), max(a1['recall'])
    #     min_a2_x, max_a2_x = min(a2['recall']), max(a2['recall'])
    #     min_a3_x, max_a3_x = min(a3['recall']), max(a3['recall'])
    #     new_a1_x = list(np.linspace(min_a1_x, max_a1_x, npoints))
    #     new_a2_x = list(np.linspace(min_a2_x, max_a2_x, npoints))
    #     new_a3_x = list(np.linspace(min_a3_x, max_a3_x, npoints))
    #     new_a1_y, new_a2_y, new_a3_y = [], [], []
    #     for i, nax in enumerate(new_a1_x[:-1]):
    #         new_a1_y.append(a1.loc[(a1['recall'] >= new_a1_x[i]) & (a1['recall'] < new_a1_x[i+1]) & (a1['recall'] !=0), 'precision'].mean())
    #     for i, nax in enumerate(new_a2_x[:-1]):
    #         new_a2_y.append(a2.loc[(a2['recall'] >= new_a2_x[i]) & (a2['recall'] < new_a2_x[i+1]) & (a2['recall'] !=0), 'precision'].mean())
    #     for i, nax in enumerate(new_a3_x[:-1]):
    #         new_a3_y.append(a3.loc[(a3['recall'] >= new_a3_x[i]) & (a3['recall'] < new_a3_x[i+1]) & (a3['recall'] !=0), 'precision'].mean())
    #     new_a1_x = new_a1_x[:-1]
    #     new_a2_x = new_a2_x[:-1]
    #     new_a3_x = new_a3_x[:-1]
    #     new_a1_x.insert(0, 0.01)
    #     new_a2_x.insert(0, 0.01)
    #     new_a3_x.insert(0, 0.01)
    #     new_a1_y.insert(0, [m for m in new_a1_y if str(m) != 'nan'][0])
    #     new_a2_y.insert(0, [m for m in new_a2_y if str(m) != 'nan'][0])
    #     new_a3_y.insert(0, [m for m in new_a3_y if str(m) != 'nan'][0])
    #     midx = [np.mean([new_a1_x[i], new_a2_x[i], new_a3_x[i]]) for i in range(npoints-1)]
    #     midy = [np.nanmean([new_a1_y[i], new_a2_y[i], new_a3_y[i]]) for i in range(npoints-1)]
    #     stdy = [np.std([new_a1_y[i], new_a2_y[i], new_a3_y[i]])/np.sqrt(3) for i in range(npoints-1)]
    #     #midx.insert(0, 0.01)
    #     #midy.insert(0, [m for m in midy if str(m) != 'nan'][0])
    #     midx.insert(len(midx), midx[midy.index([m for m in midy if str(m) != 'nan'][-1])])
    #     midy.insert(len(midy), 0)
    #     #stdy.insert(0, [m for m in stdy if str(m) != 'nan'][0])
    #     stdy.insert(len(stdy), [m for m in stdy if str(m) != 'nan'][-1])
    #     print(len(midx))
    #     print(len(stdy), len(midy))
    #     mid = pd.DataFrame.from_dict({'midx': midx, 'midy': midy, 'stdy': stdy})
    #     mid.dropna(inplace=True)
    #     fig, ax = plt.subplots(figsize=(8, 8))
    #     plt.plot(a1['recall'], a1['precision'], c='k')
    #     plt.plot(a2['recall'], a2['precision'], c='k')
    #     plt.plot(a3['recall'], a3['precision'], c='k')
    #     plt.plot(mid['midx'], mid['midy'], ls='--', c=config.colors[config.methods.index(method)])
    #     plt.plot(mid['midx'], mid['midy'] + mid['stdy'], ls='-', c=config.colors[config.methods.index(method)])
    #     plt.plot(mid['midx'], mid['midy'] - mid['stdy'], ls='-', c=config.colors[config.methods.index(method)])
    #     ax.fill_between(mid['midx'], mid['midy'] - mid['stdy'], mid['midy'] + mid['stdy'], color=config.colors[config.methods.index(method)], alpha=0.1)
    #     plt.ylim([0, 1])
    #     plt.semilogx()
    #     plt.xlim([0.01, 1])
    #     plt.ylim([0.01, 1])
    #     plt.show()
    #     mid = mid[mid['midx'] >= 0.01]
    #
    #     #mid.to_csv('data/featureanalysis/PRcurve_mean_std_naive_'+method+'.csv')
    #     figfull, axfull = plt.subplots(figsize=(8, 8))
    #     auprc_list = [a1['AUPRC'].unique(), a2['AUPRC'].unique(), a3['AUPRC'].unique()]
    #     plt.plot(mid['midx'], mid['midy'], ls='--', lw=1, c='grey', label='naive '+method + ' - AUPRC = {:.2f} ± {:.2f}'.format(np.mean(auprc_list), np.std(auprc_list)))
    #     axfull.fill_between(mid['midx'], mid['midy'] - mid['stdy'], mid['midy'] + mid['stdy'], alpha=0.2, color='grey') #  color=config.colors[config.methods.index(method)])
    #     plt.ylim([0, 1])
    #     plt.semilogx()
    #     plt.xlim([0.01, 1])
    #     plt.ylim([0.01, 1])
    #
    #
    #     a1 = dfresspecificall[dfresspecificall['patientid'] == '1014']
    #     a2 = dfresspecificall.loc[dfresspecificall['patientid'] == '986']
    #     a3 = dfresspecificall.loc[dfresspecificall['patientid'] == '123']
    #
    #     min_a1_x, max_a1_x = min(a1['recall']), max(a1['recall'])
    #     min_a2_x, max_a2_x = min(a2['recall']), max(a2['recall'])
    #     min_a3_x, max_a3_x = min(a3['recall']), max(a3['recall'])
    #     new_a1_x = list(np.linspace(min_a1_x, max_a1_x, npoints))
    #     new_a2_x = list(np.linspace(min_a2_x, max_a2_x, npoints))
    #     new_a3_x = list(np.linspace(min_a3_x, max_a3_x, npoints))
    #
    #     new_a1_y, new_a2_y, new_a3_y = [], [], []
    #     for i, nax in enumerate(new_a1_x[:-1]):
    #         new_a1_y.append(a1.loc[(a1['recall'] >= new_a1_x[i]) & (a1['recall'] < new_a1_x[i+1]) & (a1['recall'] !=0), 'precision'].mean())
    #     for i, nax in enumerate(new_a2_x[:-1]):
    #         new_a2_y.append(a2.loc[(a2['recall'] >= new_a2_x[i]) & (a2['recall'] < new_a2_x[i+1]) & (a2['recall'] !=0), 'precision'].mean())
    #     for i, nax in enumerate(new_a3_x[:-1]):
    #         new_a3_y.append(a3.loc[(a3['recall'] >= new_a3_x[i]) & (a3['recall'] < new_a3_x[i+1]) & (a3['recall'] !=0), 'precision'].mean())
    #     new_a1_x = new_a1_x[:-1]
    #     new_a2_x = new_a2_x[:-1]
    #     new_a3_x = new_a3_x[:-1]
    #     new_a1_x.insert(0, 0.01)
    #     new_a2_x.insert(0, 0.01)
    #     new_a3_x.insert(0, 0.01)
    #     new_a1_y.insert(0, [m for m in new_a1_y if str(m) != 'nan'][0])
    #     new_a2_y.insert(0, [m for m in new_a2_y if str(m) != 'nan'][0])
    #     new_a3_y.insert(0, [m for m in new_a3_y if str(m) != 'nan'][0])
    #     midx = [np.mean([new_a1_x[i], new_a2_x[i], new_a3_x[i]]) for i in range(npoints-1)]
    #     midy = [np.nanmean([new_a1_y[i], new_a2_y[i], new_a3_y[i]]) for i in range(npoints-1)]
    #     stdy = [np.std([new_a1_y[i], new_a2_y[i], new_a3_y[i]])/np.sqrt(3) for i in range(npoints-1)]
    #     midx.insert(len(midx), midx[midy.index([m for m in midy if str(m) != 'nan'][-1])])
    #     midy.insert(len(midy), 0)
    #     stdy.insert(len(stdy), [m for m in stdy if str(m) != 'nan'][-1])
    #     print(len(midx))
    #     print(len(stdy), len(midy))
    #     mid = pd.DataFrame.from_dict({'midx': midx, 'midy': midy, 'stdy': stdy})
    #     mid.dropna(inplace=True)
    #     # fig, ax = plt.subplots(figsize=(8, 8))
    #     # plt.plot(a1['recall'], a1['precision'], c='k')
    #     # plt.plot(a2['recall'], a2['precision'], c='k')
    #     # plt.plot(a3['recall'], a3['precision'], c='k')
    #     # plt.plot(mid['midx'], mid['midy'], ls='--', c=config.colors[config.methods.index(method)])
    #     # plt.plot(mid['midx'], mid['midy'] + mid['stdy'], ls='-', c=config.colors[config.methods.index(method)])
    #     # plt.plot(mid['midx'], mid['midy'] - mid['stdy'], ls='-', c=config.colors[config.methods.index(method)])
    #     # ax.fill_between(mid['midx'], mid['midy'] - mid['stdy'], mid['midy'] + mid['stdy'], color=config.colors[config.methods.index(method)], alpha=0.1)
    #     # plt.ylim([0, 1])
    #     # plt.semilogx()
    #     # plt.xlim([0.01, 1])
    #     # plt.ylim([0.01, 1])
    #     # plt.show()
    #     mid = mid[mid['midx'] >= 0.01]
    #
    #     auprc_list = [a1['AUPRC'].unique(), a2['AUPRC'].unique(), a3['AUPRC'].unique()]
    #     plt.plot(mid['midx'], mid['midy'], ls='-', lw=2, c=config.colors[config.methods.index(method)], label='cfDNA-tuned '+method + ' - AUPRC = {:.2f} ± {:.2f}'.format(np.mean(auprc_list), np.std(auprc_list)))
    #     axfull.fill_between(mid['midx'], mid['midy'] - mid['stdy'], mid['midy'] + mid['stdy'], color=config.colors[config.methods.index(method)], alpha=0.5)
    #     print(mid['midx'])
    #     #f_scores = np.linspace(0.1, 0.9, num=9)
    #     #for f_score in f_scores:
    #     #    x = np.linspace(0.01, 1)
    #     #    y = f_score * x / (2 * x - f_score)
    #     #    plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.5, lw=1)
    #     #    plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.6, y[45] + 0.02), color='grey')
    #     plt.legend()
    #     #plt.savefig('data/featureanalysis/PRcuve_'+method+'_allpatients_150x.png')
    #     #plt.savefig('data/featureanalysis/PRcuve_'+method+'_allpatients_150x.svg')
    #     plt.show()

    """   
    for method in [m for m in config.methods if m!='varnet' and m!= 'sinvict' and m!='abemus']:
        method = 'smurf'
        featureimp = []
        for mixtureid in mixtureids:
            patientid = str(mixtureid.split('_')[0].split('-')[1])
            fi = pd.read_csv('data/featureanalysis/barplot_featureimportance_'+method+'_cfDNAtuned_leftoutpatient'+str(patientid)+'_150x_RF.csv', index_col=0)
            fi.columns = [patientid]
            featureimp.append(fi)
        featureimp = pd.concat(featureimp, axis=1)
        featureimp['median'] = featureimp.median(axis=1)
        featureimp.sort_values(by='median', ascending=False, inplace=True)
        print(featureimp.T)
        featureimp = featureimp[featureimp['median'] > 0.001]
        plt.figure(figsize=(12, 10))
        g = sns.barplot(data=featureimp[['median']].T, color=color_dict[method])
        g.set_xticklabels(g.get_xticklabels(), rotation=90)
        featureimp.drop('median', axis=1, inplace=True)
        print(featureimp)
        for mi, mixtureid in enumerate(mixtureids):
            patientid = str(mixtureid.split('_')[0].split('-')[1])
            sns.stripplot(data=featureimp[[patientid]].T, marker=config.markers[mi], s=15, color='k', label=patientid, ax=g)
        plt.xlabel('Features')
        plt.ylabel('Feature importance in Random Forest models in CV')
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())
        plt.savefig('data/featureanalysis/featureimportance_'+method+'_all_150x.svg')
        plt.savefig('data/featureanalysis/featureimportance_'+method+'_all_150x.png')
        plt.show()
    """
    #
    #     fig, ax = plt.subplots(figsize=(8, 8))
    #     plt.plot(a1['recall'], a1['precision'], c='k')
    #     plt.plot(a2['recall'], a2['precision'], c='k')
    #     plt.plot(a3['recall'], a3['precision'], c='k')
    #     plt.plot(mid['midx'], mid['midy'], ls='--', c=config.colors[config.methods.index(method)])
    #     plt.plot(mid['midx'], mid['midy'] + mid['stdy'], ls='-', c=config.colors[config.methods.index(method)])
    #     plt.plot(mid['midx'], mid['midy'] - mid['stdy'], ls='-', c=config.colors[config.methods.index(method)])
    #     ax.fill_between(mid['midx'], mid['midy'] - mid['stdy'], mid['midy'] + mid['stdy'], color=config.colors[config.methods.index(method)], alpha=0.5)
    #     plt.ylim([0, 1])
    #     plt.semilogx()
    #     plt.xlim([0.01, 1])
    #     plt.ylim([0.01, 1])
    #     plt.show()
    #     mid.to_csv('data/featureanalysis/PRcurve_mean_std_cfDNAtuned_'+method+'.csv')

    # plot PR curve overlay
    # fig, ax = plt.subplots(figsize=(8, 8))
    # for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:
    #     print(method)
    #     mid = pd.read_csv('data/featureanalysis/PRcurve_mean_std_naive_'+method+'.csv', index_col=0)
    #     plt.plot(mid['midx'], mid['midy'], ls='-', lw=1.5, c=config.colors[config.methods.index(method)], label='naive '+method)
    #     ax.fill_between(mid['midx'], mid['midy'] - mid['stdy'], mid['midy'] + mid['stdy'], color=config.colors[config.methods.index(method)], alpha=0.1)
    #     mid = pd.read_csv('data/featureanalysis/PRcurve_mean_std_cfDNAtuned_'+method+'.csv', index_col=0)
    #     plt.plot(mid['midx'], mid['midy'], ls='-', lw=3, c=config.colors[config.methods.index(method)], label='cfDNA-tuned '+method)
    #     ax.fill_between(mid['midx'], mid['midy'] - mid['stdy'], mid['midy'] + mid['stdy'], color=config.colors[config.methods.index(method)], alpha=0.2)
    #     plt.ylim([0, 1])
    #     plt.semilogx()
    #     plt.xlim([0.01, 1])
    #     plt.ylim([0.01, 1])
    #     f_scores = np.linspace(0.1, 0.9, num=9)
    #     for f_score in f_scores:
    #         x = np.linspace(0.01, 1)
    #         y = f_score * x / (2 * x - f_score)
    #         plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.1, lw=1)
    #         plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.6, y[45] + 0.02), color='grey')
    # plt.show()
    """
    auprc_all = []
    maxprecision_all = []
    method_all = []
    patient_all = []
    model_all = []
    for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:
        for patientid in ['1014', '123', '986']:
            dfresnaive = pd.read_csv('data/featureanalysis/testpatient'+patientid+'/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_naive.csv', index_col=0)
            dfresspecific = pd.read_csv('data/featureanalysis/testpatient'+patientid+'/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_cfdnaspecific.csv', index_col=0)
            auprc_all.append(dfresnaive['AUPRC'].unique())
            maxprecision_all.append(dfresnaive['maxrecallatleast0_03precision'].unique())
            auprc_all.append(dfresspecific['AUPRC'].unique())
            maxprecision_all.append(dfresspecific['maxrecallatleast0_03precision'].unique())
            model_all = model_all + ['naive model', 'cfDNA-tuned model']
            patient_all = patient_all + [patientid, patientid]
            method_all = method_all + [method, method]

    metricspd = pd.DataFrame.from_dict({'patient': patient_all, 'model': model_all, 'AUPRC': auprc_all, 'sensitivity at 3% precision': maxprecision_all, 'caller': method_all})
    metricspd['AUPRC'] = metricspd['AUPRC'].astype(float)
    metricspd['sensitivity at 3% precision'] = metricspd['sensitivity at 3% precision'].astype(float)
    #sns.catplot(data=metricspd, x="model", y="AUPRC", hue="caller", capsize=.2, palette=color_dict, error_bars=('se', 1), kind="point", height=6, aspect=.75, dodge=True)
    #sns.catplot(data=metricspd, x="model", y="AUPRC", hue="caller", kind='point', error_bars=('se', 1), capsize=.2, palette=color_dict, height=8, aspect=0.75)
    fig, ax = plt.subplots(figsize=(5, 8))
    for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:
        tmpmetric = metricspd[metricspd['caller'] == method].groupby('model')[['AUPRC']].mean()
        tmpmetric = tmpmetric.loc[['naive model', 'cfDNA-tuned model']]
        print(tmpmetric)
        plt.plot(list(tmpmetric.index), tmpmetric['AUPRC'], ls='-', marker='o', markersize=10, lw=3, c=config.colors[config.methods.index(method)])
        tmpstd = metricspd[metricspd['caller'] == method].groupby('model')[['AUPRC']].std().apply(lambda x: x/np.sqrt(3))
        tmpstd = tmpstd.loc[['naive model', 'cfDNA-tuned model']]
        print(tmpstd)
        ax.fill_between(list(tmpmetric.index), tmpmetric['AUPRC'] - tmpstd['AUPRC'], tmpmetric['AUPRC'] + tmpstd['AUPRC'], color=config.colors[config.methods.index(method)], alpha=0.25)
    plt.ylabel('AUPRC')
    plt.savefig('data/featureanalysis/AUPRCgain_all_allpatients_150x.png')
    plt.savefig('data/featureanalysis/AUPRCgain_all_allpatients_150x.svg')
    plt.show()
    fig, ax = plt.subplots(figsize=(5, 8))
    for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:
        tmpmetric = metricspd[metricspd['caller'] == method].groupby('model')[['sensitivity at 3% precision']].mean()
        tmpmetric = tmpmetric.loc[['naive model', 'cfDNA-tuned model']]
        print(tmpmetric)
        plt.plot(list(tmpmetric.index), tmpmetric['sensitivity at 3% precision'], ls='-', marker='o', markersize=10, lw=3, c=config.colors[config.methods.index(method)])
        tmpstd = metricspd[metricspd['caller'] == method].groupby('model')[['sensitivity at 3% precision']].std().apply(lambda x: x/np.sqrt(3))
        tmpstd = tmpstd.loc[['naive model', 'cfDNA-tuned model']]
        print(tmpstd)
        ax.fill_between(list(tmpmetric.index), tmpmetric['sensitivity at 3% precision'] - tmpstd['sensitivity at 3% precision'], tmpmetric['sensitivity at 3% precision'] + tmpstd['sensitivity at 3% precision'], color=config.colors[config.methods.index(method)], alpha=0.25)
    plt.ylabel('sensitivity at 3% precision')
    plt.savefig('data/featureanalysis/maxrecallatleast0_03precisiongain_all_allpatients_150x.png')
    plt.savefig('data/featureanalysis/maxrecallatleast0_03precisiongain_all_allpatients_150x.svg')
    plt.show()
    """

    #     plt.ylim([0, 1])
    #     plt.semilogx()
    #     #plt.semilogy()
    #     plt.xlim([0.01, 1])
    #     plt.ylim([0.01, 1])
    #     #plt.savefig('data/featureanalysis/PRcuve_all_allpatients_150x.svg')
    #     #plt.savefig('data/featureanalysis/PRcuve_all_allpatients_150x.png')
    # plt.show()

        # patientid = str(mixtureid_test.split('_')[0].split('-')[1])
        # # plot PR curve overlay
        # color_dict = {config.methods[i]: config.colors[i] for i in range(len(config.methods))}
        # dfresnaivecfdna = pd.DataFrame(columns = ['method', 'model', 'AUPRC', 'maxrecallatleast0_03precision'])
        # for method in [m for m in config.methods if m != 'varnet' and m != 'abemus' and m != 'sinvict']:
        #     print(method)
        #     dfresnaive = pd.read_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_naive.csv', index_col=0)
        #     dfresnaivecfdna = dfresnaivecfdna.append({'method': method, 'model': 'naive model', 'AUPRC': dfresnaive['AUPRC'].iloc[0], 'maxrecallatleast0_03precision': dfresnaive['maxrecallatleast0_03precision'].iloc[0]}, ignore_index=True)
        #     dfrescfdna = pd.read_csv('data/featureanalysis/PRAUPRC_featureimportance_'+method+'_'+patientid+'_150x_cfdnaspecific.csv', index_col=0)
        #     dfresnaivecfdna = dfresnaivecfdna.append({'method': method, 'model': 'cfDNA-tuned model', 'AUPRC': dfrescfdna['AUPRC'].iloc[0], 'maxrecallatleast0_03precision': dfrescfdna['maxrecallatleast0_03precision'].iloc[0]}, ignore_index=True)
        # print(dfresnaivecfdna)
        # print(dfresnaivecfdna.shape)
        # sns.catplot(data=dfresnaivecfdna, x="model", y="AUPRC", hue="method", capsize=.2, palette=color_dict, kind="point", height=6, aspect=.75)
        # plt.savefig('data/featureanalysis/AUPRCgain_all_'+patientid+'_150x.svg')
        # plt.savefig('data/featureanalysis/AUPRCgain_all_'+patientid+'_150x.png')
        # sns.catplot(data=dfresnaivecfdna, x="model", y="maxrecallatleast0_03precision", hue="method", capsize=.2, palette=color_dict, kind="point", height=6, aspect=.75)
        # plt.savefig('data/featureanalysis/maxrecallatleast0_03precisiongain_all_'+patientid+'_150x.svg')
        # plt.savefig('data/featureanalysis/maxrecallatleast0_03precisiongain_all_'+patientid+'_150x.png')
        # plt.show()
