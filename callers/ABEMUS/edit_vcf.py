# design_mixtures.py

import os
import numpy as np
import pandas as pd
import subprocess
from itertools import chain


def chainer(s):
    """Returns list from series of comma-separated strings.

    Parameters
    ----------
    s : pandas Series
        series of string items containing commas

    Returns
    -------
    list
        a list containing items splitted where the comma used to be
    """

    return list(chain.from_iterable(s.str.split(',')))


def edit_vcf(vcfpath, save=True):
    """Edits a vcf file by replacing comma-separated mutliple alternate variants to (ex: G->A,T)
    to simplified entry (ex: G->A only).

    Parameters
    ----------
    vcfpath : string
        path to the vcf file to edit.
    save : bool
        whether to save the edited vcf (default is True).

    Returns
    -------
    res : pandas DataFrame
        a dataframe of the edited vcf without any comma in 'ALT'.
    """

    print(vcfpath)
    nrowstoskip = int(subprocess.check_output("grep -n -m 1 '#CHROM' " + vcfpath + " | cut -f1 -d:", shell=True))
    print(nrowstoskip)
    vcf_df = pd.read_csv(vcfpath, skiprows=nrowstoskip-1, sep='\t', memory_map=True)
    print(vcf_df.head())
    print(vcf_df.shape)
    # calculate lengths of splits
    lens = vcf_df['ALT'].str.split(',').map(len)
    print(lens)
    # create new dataframe, repeating or chaining as appropriate
    res = pd.DataFrame({'#CHROM': np.repeat(vcf_df['#CHROM'], lens),
                        'POS': np.repeat(vcf_df['POS'], lens),
                        'ID': np.repeat(vcf_df['ID'], lens),
                        'REF': np.repeat(vcf_df['REF'], lens),
                        'ALT': chainer(vcf_df['ALT']),
                        'QUAL': np.repeat(vcf_df['QUAL'], lens),
                        'FILTER': np.repeat(vcf_df['FILTER'], lens),
                        'INFO': np.repeat(vcf_df['INFO'], lens)})
    res.drop_duplicates(subset=['#CHROM', 'POS'], inplace=True)  # keep first occurrence
    res = res[(res['REF'] == 'A') | (res['REF'] == 'C') | (res['REF'] == 'G') | (res['REF'] == 'T')]  # ensure ref is single base
    res = res[(res['ALT'] == 'A') | (res['ALT'] == 'C') | (res['ALT'] == 'G') | (res['ALT'] == 'T')]  # ensure alt is single base
    print(res)
    if save:
        if os.path.exists(vcfpath[:-4]+'_edited.vcf'):
            print("{} already exists. No overwritting, saving is skipped.".format(vcfpath[:-4]+'_edited.vcf'))
        else:
            res.to_csv(vcfpath[:-4]+'_edited.vcf', sep='\t', header=False, index=False)
    return res


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcfpath')
    args = parser.parse_args()
    print(args)

    vcfpath = args.vcfpath
    res = edit_vcf(vcfpath, save=True)
