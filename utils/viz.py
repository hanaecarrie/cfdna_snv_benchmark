import io
import os
import numpy as np
import pandas as pd



def set_display_params(config):
    import seaborn as sns
    import matplotlib.pyplot as plt
    if config.context == 'talk':
        sns.set(style=config.context_talk.style, context=config.context_talk.context, rc={"lines.linewidth": config.context_talk.llw, "legend.fontsize": config.context_talk.llw})
        plt.style.use(config.context_talk.styleuse)
        plt.rcParams.update({"grid.linewidth": config.context_talk.glw, "grid.alpha": config.context_talk.galpha, 'font.size': config.context_talk.fs})
        sns.set_palette(config.context_talk.palette)
    elif config.context == 'paper':
        print(config.context_paper)
        sns.set(style=config.context_paper.style, context=config.context_paper.context, rc={"lines.linewidth": config.context_paper.llw, "legend.fontsize": config.context_paper.llw})
        plt.style.use(config.context_paper.styleuse)
        plt.rcParams.update({"grid.linewidth": config.context_paper.glw, "grid.alpha": config.context_paper.galpha, 'font.size': config.context_paper.fs})
        sns.set_palette(config.context_paper.palette)
    else:
        raise ValueError("unknown context {}. Should be either 'talk' or 'paper'".format(config.context))


def read_vcf(path):
    if not os.path.exists(path) and os.path.exists(path+'.gz'):
        fp = open(path, "wb")
        with gzip.open(path+'.gz', "rb") as f:
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
    if os.path.exists(vcf_path) or os.path.exists(vcf_path+'.gz'):
        res = read_vcf(vcf_path)
        res['callers'] = res['INFO'].apply(lambda x: pd.Series(x.split('CALLERS=')[1].split(';')[0]))
        res['type'] = np.nan
        res['type'][res['ALT'].str.len() - res['REF'].str.len() == 0] = 'SNV'
        res['type'][res['ALT'].str.len() - res['REF'].str.len() > 0] = 'INS'
        res['type'][res['ALT'].str.len() - res['REF'].str.len() < 0] = 'DEL'
        res['type'][res['ID'].str.contains('rs')] = 'SNP'
        for m in methods:
            res[m] = res['INFO'].str.contains(m)
        sample = res[['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'type', *methods]]
        for m in methods:
            if m == 'vardict': # P-value
                sample[m+'_score'] = [float(i.split('SSF=')[1].split(';')[0]) if 'SSF' in i else 0 for i in res['INFO']]
                # res['INFO'].apply(lambda x: pd.Series(x.split('SSF=')[1].split(';')[0]).astype(float) if 'SSF' in x else 0)
            elif m == 'varscan': # P-value
                sample[m+'_score'] = [float(i.split('SSC=')[1].split(';')[0])/100 if 'SSC' in i else 0 for i in res['INFO']]
            elif m == 'mutect2': # logodds to probability score prob = exp(logTLOD)/(1+exp(logTLOD))
                sample[m+'_score'] = [np.exp(float(i.split('TLOD=')[1].split(';')[0]))/(1+np.exp(float(i.split('TLOD=')[1].split(';')[0]))) if 'TLOD' in i else 0 for i in res['INFO']]
                # [float(i.split('TLOD=')[1].split(';')[0])/(1+float(i.split('TLOD=')[1].split(';')[0])) if 'TLOD' in i else 0 for i in res['INFO']]
            elif m == 'freebayes': # logodds to probability score prob = exp(logODDS)/(1+exp(logODDS))
                sample[m+'_score'] = [np.exp(float(i.split('ODDS=')[1].split(';')[0]))/(1+np.exp(float(i.split('ODDS=')[1].split(';')[0]))) if 'ODDS' in i else 0 for i in res['INFO']]
                # [float(i.split('ODDS=')[1].split(';')[0])/(1+float(i.split('ODDS=')[1].split(';')[0])) if 'ODDS' in i else 0 for i in res['INFO']] 
            elif m == 'strelka2': # phred score to probability, prob = 1 - 10^(-SomaticEVS/10)
                sample[m+'_score'] = [10**(-float(i.split('SomaticEVS=')[1].split(';')[0])/10) if 'SomaticEVS' in i else 0 for i in res['INFO']]
        sample['CHROM_POS'] = sample['CHROM'].astype('str').str.cat(sample['POS'].astype('str'),sep="_")
        sample.set_index('CHROM_POS', inplace = True)
        if chrom != 'all':
            print('select a single chrom = {} for analysis'.format(chrom))
            sample = sample[sample['CHROM'] == chrom]
        return sample
    else:
        print("sample is not present with path {}".format(vcf_path))
        return None
    

def load_calls_from_vcf_dilutionseries(dirpath, plasmasample, reference, dilutionseries, methods, prefix='dilution_chr22', chrom='all'):
    vcf_samples_dict = {}
    samples_dict = {}
    ci = 0
    dilutionseries_new = []
    for i, d in enumerate(dilutionseries):
        print('vcf_pd_'+str(ci), d)
        d0 = str(d[0]).replace('.', '_')
        d1 = str(d[1]).replace('.', '_')
        if d == (1,0):
            ref = 'pooledhealthy'
        vcf_path = os.path.join(*dirpath, prefix+plasmasample+"_"+str(d[0])+"_"+ref+"_"+str(d[1]),
                                                            prefix+plasmasample+"_"+d0+"_"+ref+"_"+d1+"-ensemble-annotated.vcf") 
        if d == (1,0):
            ref = reference
        
        vcf_sample = load_calls_from_vcf(vcf_path, methods, chrom=chrom)
        if vcf_sample is None:
            print("dilution {} is not present with path {}".format(d, vcf_path))
        else:
            vcf_samples_dict['sample_'+str(ci)] = vcf_sample
            ci += 1
            dilutionseries_new.append(d)      
    return vcf_samples_dict, dilutionseries_new


def count_mutations(samples, methods, samples_tf, mutationtypes=['all', 'INDEL', 'SNV', 'SNP'], threshold=0.5):
    for mutationtype in mutationtypes:
        if mutationtype not in ['all', 'INDEL', 'SNV', 'SNP']:
            raise ValueError('mutation type {} is unknown. It should be in {}'.format(mutationtype, ['all', 'INDEL', 'SNV', 'SNP']))
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
                numbersnvs_pd = pd.DataFrame.from_dict({'sample_'+ str(si): nb_snv}).T
                numbersnvs_pd.columns = methods
            else:
                numbersnvs_pd.loc['sample_'+ str(si)] = nb_snv
            numbersnvs_pd = numbersnvs_pd.rename(index=samples_tf)
        return numbersnvs_pd
    
    
def get_pr_table():
    #get PR table of values
    """
    PR.table = function (data, value, truth, verbose=FALSE, decreasing=F) {
  
      df = data[,c(value, truth)]
      colnames(df) = c('value','truth')

      thresholds = sort(unique(df$value))
      thresholds = thresholds[1:length(thresholds)]

      if(length(thresholds)>10000) {
        seq.idx = seq(from=thresholds[1], to=length(thresholds), length.out = 10000)
        thresholds = thresholds[seq.idx]
      }

      table = data.frame(matrix(nrow = length(thresholds),ncol = 7))
      colnames(table) = c('Threshold', 'TP','FP','FN','Precision','Recall','F1')
      table$Threshold = thresholds

      print(paste('Total thresholds:', length(thresholds)))

      for (i in 1:length(thresholds)) {

        if (verbose==T) {
          print(thresholds[i])
        }

        df$predict = NA

        if (decreasing==F){
        df[df$value>=thresholds[i],'predict'] = T 
        } else {
          df[df$value<=thresholds[i],'predict'] = T 
        }
        df[is.na(df$predict),'predict'] = F 

        table$TP[i] = nrow(df[df$predict==T & df$truth==T,])
        table$FP[i] = nrow(df[df$predict==T & df$truth==F,])
        table$FN[i] = nrow(df[df$predict==F & df$truth==T,])

      }

      table$Precision = table$TP/(table$TP + table$FP)
      table$Recall = table$TP/(table$TP + table$FN)
      table$F1 = (2*table$Precision*table$Recall)/(table$Precision + table$Recall)

      return(table)
    }
    """
    pass


def load_files(filenames):
    for filename in filenames:
        yield pd.read_csv(filename, names=['sample_id', 'tumor_burden'])

