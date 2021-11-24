import gzip
import io
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import warnings
from sklearn.metrics import precision_recall_curve, roc_curve, f1_score, roc_auc_score, precision_score, recall_score, average_precision_score

warnings.filterwarnings("ignore")


def load_files(filenames):
    for filename in filenames:
        yield pd.read_csv(filename, names=['sample_id', 'tumor_burden'])


def read_vcf(path):
    if not os.path.exists(path) and os.path.exists(path + '.gz'):
        fp = open(path, "wb")
        with gzip.open(path + '.gz', "rb") as f:
            bindata = f.read()
        fp.write(bindata)
        fp.close()
    if os.path.exists(path):
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        res = pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})
        return res


def load_calls_from_vcf(vcf_path, methods, chrom='all'):
    sample = None
    for m in methods:
        vcf_path_m = vcf_path.split('ensemble')[0] + m + vcf_path.split('ensemble')[1]
        if os.path.exists(vcf_path_m) or os.path.exists(vcf_path_m + '.gz'):
            # print(vcf_path_m)
            res = read_vcf(vcf_path_m)
            res = res[res['FILTER'] == "PASS"]
            # res['callers'] = res['INFO'].apply(lambda x: pd.Series(x.split('CALLERS=')[1].split(';')[0]))
            res['type'] = np.nan
            res['type'][res['ALT'].str.len() - res['REF'].str.len() == 0] = 'SNV'
            res['type'][res['ALT'].str.len() - res['REF'].str.len() > 0] = 'INS'
            res['type'][res['ALT'].str.len() - res['REF'].str.len() < 0] = 'DEL'
            res['type'][res['ID'].str.contains('rs')] = 'SNP'
            # for m in methods:
            #    res[m] = res['INFO'].str.contains(m)
            info = res['INFO']
            res = res[['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'type']]
            res.rename(columns={'FILTER': m}, inplace=True)
            res[m] = True
            # sample = res[['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'type', *methods]]
            # for m in methods:
            if m == 'vardict':  # P-value /!\ opposite direction of score variation /!\
                res[m + '_score'] = [1-float(i.split('SSF=')[1].split(';')[0]) if 'SSF' in i else np.nan for i in info.tolist()]
            elif m == 'varscan':  # P-value
                res[m + '_score'] = [float(i.split('SSC=')[1].split(';')[0]) if 'SSC' in i else np.nan for i in info.tolist()]
            elif m == 'mutect2':  # logodds to probability score prob = exp(logTLOD)/(1+exp(logTLOD))
                res[m + '_score'] = [np.exp(np.log(float(i.split('TLOD=')[1].split(';')[0]))) / (1 + np.exp(np.log(float(i.split('TLOD=')[1].split(';')[0])))) if 'TLOD' in i else np.nan for i in info.to_list()]
            elif m == 'freebayes':  # logodds to probability score prob = exp(logODDS)/(1+exp(logODDS))
                res[m + '_score'] = [np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0]))) / (1 + np.exp(np.log(float(i.split('ODDS=')[1].split(';')[0])))) if 'ODDS' in i else np.nan for i in info.to_list()]
            elif m == 'strelka2':  # phred score to probability, prob = 1 - 10^(-SomaticEVS/10)
                res[m + '_score'] = [1 - (10 ** (-float(i.split('SomaticEVS=')[1].split(';')[0]) / 10)) if 'SomaticEVS' in i else np.nan for i in info.to_list()]
            res['CHROM_POS'] = res['CHROM'].astype('str').str.cat(res['POS'].astype('str'), sep="_")
            res.set_index('CHROM_POS', inplace=True)
            if sample is None:
                sample = res.copy()
            else:
                res.drop(['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'type'], axis=1, inplace=True)
                sample = sample.join(res, how='outer')  # pd.concat([sample, res], axis=1)
    if chrom != 'all' and sample is not None:
        print('select a single chrom = {} for analysis'.format(chrom))
        sample = sample[sample['CHROM'] == chrom]
        return sample
    else:
        print("sample is not present with path {}".format(vcf_path))
        return None


def load_calls_from_vcf_dilutionseries(dirpath, plasmasample, reference, dilutionseries, methods,
                                       prefix='dilution_chr22', chrom='all'):
    vcf_samples_dict = {}
    ci = 0
    dilutionseries_new = []
    for i, d in enumerate(dilutionseries):
        print('vcf_pd_' + str(ci), d)
        d0 = str(d[0]).replace('.', '_')
        d1 = str(d[1]).replace('.', '_')
        if d == (1, 0):
            ref = 'pooledhealthy'
        vcf_path = os.path.join(*dirpath, prefix + plasmasample + "_" + str(d[0]) + "_" + ref + "_" + str(d[1]),
                                prefix + plasmasample + "_" + d0 + "_" + ref + "_" + d1 + "-ensemble-annotated.vcf")
        if d == (1, 0):
            ref = reference

        vcf_sample = load_calls_from_vcf(vcf_path, methods, chrom=chrom)
        if vcf_sample is None:
            print("dilution {} is not present with path {}".format(d, vcf_path))
        else:
            vcf_samples_dict['sample_' + str(ci)] = vcf_sample
            ci += 1
            dilutionseries_new.append(d)
    return vcf_samples_dict, dilutionseries_new


def get_call_table(config, prefix, plasmasample, healthysample, dilutionseries, ground_truth_method=3, refsample='undiluted', chrom='22', muttype='SNV', tumorsample=None, vcf_ref_path=None):
    df_table = None
    methods = config.methods
    # tumor burden
    tb_dict = {}
    for i, d in enumerate(dilutionseries):
        tb_path = os.path.join(*config.dilutionfolder, "estimated_tf_chr22_"+plasmasample+"_"+str(dilutionseries[i][0])+"_"+healthysample+"_"+str(dilutionseries[i][1])+".txt")
        if d == (1, 0) and not os.path.exists(tb_path):
            tb_path = [os.path.join(*config.dilutionfolder, f) for f in os.listdir(os.path.join(*config.dilutionfolder)) if ("estimated_tf_chr22_"+plasmasample) and (f.endswith('_0.txt'))][0]
        tb_dict[str(dilutionseries[i])] = float(pd.read_csv(tb_path).columns[0])
    if vcf_ref_path is None:
        if refsample == 'tumor':
            if tumorsample is None:
                raise ValueError("no tumor sample passed while gt_sample='tumor'")
            vcf_ref_path = os.path.join(*config.bcbiofolder, tumorsample, tumorsample+"-ensemble-annotated.vcf")
        elif refsample == 'undiluted':
            vcf_ref_path = os.path.join(*config.bcbiofolder, prefix + plasmasample + "_1_pooledhealthy_0",
                                    prefix + plasmasample + "_1_pooledhealthy_0-ensemble-annotated.vcf")
    print(vcf_ref_path)
    if refsample == 'tumor':
        methods_ref = config.methods_tissue
    else:
        methods_ref = config.methods
    vcf_ref = load_calls_from_vcf(vcf_ref_path, methods_ref, chrom=chrom)
    print(vcf_ref.shape)
    vcf_ref = vcf_ref[vcf_ref['type'] == muttype]
    if type(ground_truth_method) == int:  # GROUND TRUTH = consensus across 3 callers in undiluted sample
        if int(ground_truth_method) <= len(methods_ref):
            y_true = pd.DataFrame(vcf_ref[methods_ref].T.sum())
            y_true.columns = ['truth']
            y_true['truth'][y_true['truth'] < ground_truth_method] = 0
            y_true = y_true.astype(bool)
            print(y_true.shape)
            df_table = y_true.copy()
            print(df_table.columns)
        else:
            print('number of common callers {} asked for consensus is too high. It should be <= {}'.format(ground_truth_method, len(methods_ref)))
    elif ground_truth_method == 'caller':  # GROUND TRUTH = calls detected by same method in ref sample
        df_table = vcf_ref[methods_ref]
        df_table.rename(columns={m: m + '_truth' for m in methods_ref}, inplace=True)
        df_table = df_table.astype(bool)
        print(df_table.columns)
    else:
        raise ValueError('unknown ground truth {}'.format(ground_truth_method))
    # Sample of interest
    for i, d in enumerate(dilutionseries):
        d0 = str(d[0]).replace('.', '_')
        d1 = str(d[1]).replace('.', '_')
        print('vcf_pd_'+str(i), d)
        vcf_path = os.path.join(*config.bcbiofolder, prefix+plasmasample+"_"+str(d[0])+"_"+healthysample+"_"+str(d[1]),
                                prefix+plasmasample+"_"+d0+"_"+healthysample+"_"+d1+"-ensemble-annotated.vcf")
        if d == (1, 0) and not os.path.exists(vcf_path):
            vcf_path_dir = [f for f in os.listdir(os.path.join(*config.bcbiofolder))
                        if (f.startswith(prefix+plasmasample+"_"+str(d[0])+"_")) and (f.endswith("_"+str(d[1])))]
            if vcf_path_dir:
                vcf_path = os.path.join(*config.bcbiofolder, vcf_path_dir[0], vcf_path_dir[0]+"-ensemble-annotated.vcf")
        vcf_sample = load_calls_from_vcf(vcf_path, methods, chrom=chrom)
        if vcf_sample is not None:
            vcf_sample = vcf_sample[vcf_sample['type'] == muttype]
            vcf_sample.rename(columns={m: str(round(100*tb_dict[str(d)], 3)) + '_' + m for m in methods}, inplace=True)
            vcf_sample.rename(columns={m+'_score': str(round(100*tb_dict[str(d)], 3)) + '_' + m + '_score' for m in methods}, inplace=True)
            colmerge = [str(round(100*tb_dict[str(d)], 3)) + '_' + m for m in methods] + [str(round(100*tb_dict[str(d)], 3)) + '_' + m + '_score' for m in methods]
            df_table = pd.concat([df_table, vcf_sample[colmerge]], axis=1)
    return df_table


def get_call_table_spikein(config, prefix, sample, bed_ref_path, dilutionseries, msstatus='MSS', chrom='22', muttype='SNV'):
    methods = config.methods
    y_true = pd.read_csv(bed_ref_path, sep='\t', header=None)
    y_true.columns = ['CHROM', 'POS', 'POSend', 'VAF', 'ALT']
    y_true['CHROM_POS'] = [str(y_true['CHROM'][i]) + '_' + str(y_true['POS'][i]) for i in range(y_true.shape[0])]
    y_true.set_index('CHROM_POS', inplace=True)
    y_true.drop(['CHROM', "POS", 'POSend'], axis=1, inplace=True)
    df_table = y_true.copy()

    for i, d in enumerate(dilutionseries):
        print(d)
        ds = str(d).replace('.', '_')
        vcf_path = os.path.join(*config.bcbiofolder, sample+prefix+msstatus+"_"+str(d)+"_"+muttype.lower(),
                                sample+prefix+msstatus+"_"+ds+"_"+muttype.lower()+"-ensemble-annotated.vcf")
        vcf_sample = load_calls_from_vcf(vcf_path, methods, chrom=chrom)
        if vcf_sample is not None:
            vcf_sample = vcf_sample[(vcf_sample['type'] == muttype) | (vcf_sample['type'] == 'SNP')]
            vcf_sample.rename(columns={m: 'vaf_' + str(d) + '_' + m for m in methods}, inplace=True)
            vcf_sample.rename(columns={m+'_score': 'vaf_' + str(d) + '_' + m + '_score' for m in methods}, inplace=True)
            colmerge = ['vaf_' + str(d) + '_' + m for m in methods] + ['vaf_' + str(d) + '_' + m + '_score' for m in methods]
            vcf_sample = vcf_sample.groupby(vcf_sample.index).first()  # keep first variant, different SNP variants may exist
            df_table = pd.concat([df_table, vcf_sample[colmerge]], axis=1)
    return df_table


def get_call_table_tissue(config):
    df_table = None
    for muttype in config.muttype:
        if muttype == 'snv':
            vcf_path = os.path.join(*config.tissuebenchmark.snv)
        elif muttype == 'indel':
            vcf_path = os.path.join(*config.tissuebenchmark.indel)
        else:
            raise ValueError("unknown muttype {}. Muttype should be 'snv' or 'indel'".format(muttype))
        print(vcf_path)
        vcf_samples = pd.read_csv(vcf_path, sep='\t')
        for i, sample in enumerate(config.tissuebenchmark.samples):
            print(sample)
            for f in config.tissuebenchmark.fractions:
                if f != 1:
                    print('_'.join(sample.split('_')[:-1]) + '_T' + str(int(round(100*f, 2))) + '_' + sample.split('_')[-1])
                    vcf_sample = vcf_samples[vcf_samples['Sample_Name'] == '_'.join(sample.split('_')[:-1]) + '_T' + str(int(round(100*f, 2))) + '_' + sample.split('_')[-1]]
                else:
                    vcf_sample = vcf_samples[vcf_samples['Sample_Name'] == sample]
                print(vcf_sample.shape)
                vcf_sample.reset_index(inplace=True)
                vcf_sample['CHROM_POS'] = vcf_sample['X.CHROM'].astype('str').str.cat(vcf_sample['POS'].astype('str'), sep="_")
                vcf_sample['mutect2_score'] = np.exp(np.log(vcf_sample['m2_TLOD']))/(1+np.exp(np.log(vcf_sample['m2_TLOD'])))
                vcf_sample['freebayes_score'] = np.exp(np.log(vcf_sample['f_ODDS']))/(1+np.exp(np.log(vcf_sample['f_ODDS'])))
                vcf_sample['strelka2_score'] = 1 - 10**(-vcf_sample['s2_SomaticEVS']/10)
                vcf_sample['varscan_score'] = vcf_sample['vs_SSC']
                vcf_sample['vardict_score'] = 1 - vcf_sample['vd_SSF']
                vcf_sample['sample'] = sample
                vcf_sample['purity'] = round(config.tissuebenchmark.purities[i]*f, 2)
                vcf_sample['mutation type'] = muttype
                vcf_sample.rename(columns={'FILTER_'+m[0].upper()+m[1:]: m for m in config.tissuebenchmark.methods}, inplace=True)
                vcf_columns = ['CHROM_POS', 'sample', 'purity', 'mutation type', 'TRUTH'] + config.tissuebenchmark.methods + [m+'_score' for m in config.tissuebenchmark.methods]
                if df_table is None:
                    df_table = vcf_sample[vcf_columns]
                else:
                    df_table = pd.concat([df_table, vcf_sample[vcf_columns]], ignore_index=True)
    return df_table


def count_mutations(samples, methods, samples_tf, mutationtypes=['all', 'INDEL', 'SNV', 'SNP'], threshold=0.5):
    for mutationtype in mutationtypes:
        if mutationtype not in ['all', 'INDEL', 'SNV', 'SNP']:
            raise ValueError(
                'mutation type {} is unknown. It should be in {}'.format(mutationtype, ['all', 'INDEL', 'SNV', 'SNP']))
        numbersnvs_pd = pd.DataFrame()
        for si, s in enumerate(samples):
            nb_snv = []
            for method in methods:
                if mutationtype == 'all':
                    nb_snv.append(s[s[method] > 0].shape[0])
                elif mutationtype == 'INDEL':
                    nb_snv.append(s[(s[method] > 0) & ((s['type'] == 'INS') | (s['type'] == 'DEL'))].shape[0])
                elif mutationtype == 'SNV':
                    nb_snv.append(s[(s[method] > 0) & (s['type'] == 'SNV')].shape[0])
                elif mutationtype == 'SNP':
                    nb_snv.append(s[(s[method] > 0) & (s['type'] == 'SNP')].shape[0])
                else:
                    raise ValueError('mutation type {} unknown'.format(mutationtype))
            if si == 0:
                numbersnvs_pd = pd.DataFrame.from_dict({'sample_' + str(si): nb_snv}).T
                numbersnvs_pd.columns = methods
            else:
                numbersnvs_pd.loc['sample_' + str(si)] = nb_snv
            numbersnvs_pd = numbersnvs_pd.rename(index=samples_tf)
        return numbersnvs_pd


if __name__ == "__main__":
    import os
    from utils.config import Config

    if not os.getcwd().endswith('cfdna_snv_benchmark'):
        os.chdir('../')
    print('Current working directory: {}'.format(os.getcwd()))

    config = Config("config/", "config_viz.yaml")
    df_table = get_call_table_tissue(config)
    print(df_table)