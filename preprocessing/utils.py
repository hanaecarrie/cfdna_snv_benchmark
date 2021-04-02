import io
import pandas as pd


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    if path.endswith('.gz'):
        res = pd.read_csv(io.StringIO(''.join(lines[:])), compression='gzip',
                      dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                             'QUAL': str, 'FILTER': str, 'INFO': str}, sep='\t')
    else:
        res = pd.read_csv(io.StringIO(''.join(lines[:])),
                          dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                                 'QUAL': str, 'FILTER': str, 'INFO': str}, sep='\t')
    return res