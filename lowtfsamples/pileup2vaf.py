import os
import re
import pandas as pd
import numpy as np
from collections import Counter

def removehatascii(listex): # "\^[\x00-\x7F|\s]
	postoremove = [x.start() for x in re.finditer("\^.", listex)] + [x.end()-1  for x in re.finditer("\^.", listex)]
	#print(postoremove)
	#postoremove = list(itertools.chain(*postoremove))
	#print(postoremove)
	listex = [lex for li, lex in enumerate(listex) if li not in postoremove]
	#print(''.join([lex for lex in listex]))
	return ''.join([lex for lex in listex])

def getindelpatternid(listex):
	patterns = set()
	for li, lx in enumerate(listex):
		if lx == '+' or lx == '-':
			c = 1
			while re.match('[0-9]', listex[li+1+c]) and li+1+c < len(listex) :
				c += 1
			nl = int(listex[li+1:li+1+c])
			pattern = listex[li:(li+len(str(nl))+nl+1)]
			if pattern not in patterns:
				patterns.add(pattern)
	return patterns


def getindelpattern(x):
	res = 0
	pileup, patterns, indelpatterns = x[0], x[1], x[2]
	for pattern in patterns:
		res += pileup.count(pattern)
	x[2] = res
	return x

def getotherpattern(listex):
	patterns = []
	listres = listex[:]
	for li, lx in enumerate(listex):
		if lx == '+' or lx == '-': 
			c = 1
			while re.match('[0-9]', listex[li+1+c]) and li+1+c < len(listex) :
				c += 1
			nl = int(listex[li+1:li+1+c])
			pattern = listex[li:(li+len(str(nl))+nl+1)]
			listres = listres.replace(pattern, '')
	return str(listres)

def getoccurencemostcommonindelpattern(x):
	res = 0
	pileup, patterns = x[0], x[1]
	for pattern in patterns:
		occ = pileup.count(pattern)
		res = occ if occ > res else res
	x[1] = res
	return x
	

def pileup2vaf(pileuppath):
	outputpath = os.path.join(os.path.dirname(pileuppath), os.path.basename(pileuppath).replace('pileup.txt', '')+'vaf.txt')
	print(outputpath)
	pileup_df = pd.read_csv(pileuppath, sep='\t', memory_map=True)
	print(pileup_df.shape)
	pileup_df = pileup_df.iloc[:, :-1] # remove last column about base quality
	pileup_df.columns = ['chrom', 'pos', 'ref', 'totcov', 'pileup']
	pileup_df = pileup_df[pileup_df['pileup'] != '0'] # drop comment rows	
	pileup_df['totcov'] =  pileup_df['totcov'].astype(int)
	pileup_df['pileuporiginal'] = pileup_df['pileup']
	pileup_df['pileup'] = pileup_df['pileup'].str.replace('.', '', regex=False) # base match in F strand
	pileup_df['pileup'] = pileup_df['pileup'].str.replace(',', '', regex=False) # base match in R strand
	pileup_df['pileup'] = pileup_df['pileup'].str.replace('>', '', regex=False) # ref base skipped in F strand (due to CIGAR 'N')
	pileup_df['pileup'] = pileup_df['pileup'].str.replace('<', '', regex=False) # ref base skipped in R strand (due to CIGAR 'N')
	pileup_df['pileup'] = pileup_df['pileup'].str.replace('$', '', regex=False) # last position covered by read
	pileup_df['pileup'] = pileup_df['pileup'].apply(removehatascii) # first position covered by read followed by mapping quality in ASCII character
	pileup_df['pileup'] = pileup_df['pileup'].str.upper()
	#TODO get a clean way to remove ASCII characters
	pileup_df['pileup'] = pileup_df['pileup'].str.replace(']', '', regex=False)
	pileup_df['pileup'] = pileup_df['pileup'].str.replace('Q', '', regex=False)
	pileup_df['pileup'] = pileup_df['pileup'].str.replace('!', '', regex=False)
	pileup_df['pileup'] = pileup_df['pileup'].str.replace('^', '', regex=False)
	pileup_df['altcov'] = np.nan
	pileup_df['othercov'] =  np.nan
	# case no mismatch
	pileup_df.loc[pileup_df['pileup'] == '', 'altcov'] = 0 
	pileup_df.loc[pileup_df['pileup'] == '', 'othercov'] = 0
	# case single alternate base
	pileup_df['aux'] = pileup_df['pileup'].apply(lambda x: len(Counter(x).keys()))
	pileup_df.loc[pileup_df['aux'] == 1, 'altcov'] =  pileup_df['pileup'].str.len() 
	pileup_df.loc[((pileup_df['aux'] > 1) & (~pileup_df['pileup'].str.contains('+', regex=False)) & (~pileup_df['pileup'].str.contains('-', regex=False))), 'altcov'] = pileup_df.loc[((pileup_df['aux'] > 1) & (~pileup_df['pileup'].str.contains('+', regex=False)) & (~pileup_df['pileup'].str.contains('-', regex=False))), 'pileup'].apply(lambda x: Counter(x).most_common(1)[0][1])
	pileup_df['othercov'] =  pileup_df['pileup'].str.len() - pileup_df['altcov']
	# case indels other than ref deletion * or #
	pileup_df['patterns'] = np.nan
	pileup_df['indelpatterns'] = np.nan
	pileup_df['otherpatterns'] = ''
	pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), 'patterns'] =  pileup_df['pileup'].apply(getindelpatternid)
	pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), ['pileup', 'patterns', 'indelpatterns']] =  pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), ['pileup', 'patterns', 'indelpatterns']].apply(getindelpattern, axis=1)
	pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), 'otherpatterns'] =  pileup_df['pileup'].apply(getotherpattern)
	pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), 'othercov'] =  pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), 'otherpatterns'].apply(lambda x: sum(list(Counter(x).values())) if x != '' else 0) +  pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), 'indelpatterns']
	pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), ['pileup', 'patterns']] =   pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), ['pileup', 'patterns']].apply(getoccurencemostcommonindelpattern, axis=1)
	pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), 'otherpatterns'] = pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), 'otherpatterns'].apply(lambda x: Counter(x).most_common(1)[0][1] if x != '' else 0)
	pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), 'altcov'] = pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), ['patterns', 'otherpatterns']].max(axis=1)
	pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), 'othercov'] =  pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), 'othercov']  -  pileup_df.loc[((pileup_df['pileup'].str.contains('+', regex=False)) | (pileup_df['pileup'].str.contains('-', regex=False))), 'altcov']
	pileup_df['vaf'] =  pileup_df['altcov']/pileup_df['totcov']
	pileup_df['vaf other variants'] = pileup_df['othercov']/ pileup_df['totcov']
	pileup_df['vaf all variants'] = (pileup_df['altcov']+pileup_df['othercov'])/ pileup_df['totcov']
	output_df = pileup_df[['chrom', 'pos', 'ref', 'totcov', 'pileup', 'altcov', 'othercov', 'vaf' , 'vaf other variants', 'vaf all variants']]
	output_df.to_csv(outputpath)
	return output_df

if __name__== "__main__":
	pileuppath = '/mnt/projects/huangwt/wgs/hanae/targeted/NCC_CRC-1014_090516-CW-T_deepWGS_pileup.txt'
	output_df = pileup2vaf(pileuppath)
	output_df = output_df[output_df['pileup'] != '']
	#print(output_df['pileup'].value_counts().head(50))
	print(output_df.loc[(((output_df['pileup'].str.contains('+', regex=False)) | (output_df['pileup'].str.contains('-', regex=False))) & (output_df['pileup'].str.len() > 10))].head(50))
	#aux = output_df[output_df['pileup'].str.contains('-')][['pileup', 'pileuporiginal']]
	#for i, ri in aux.iterrows():
#		print(ri)
#		print(ri['pileup'], ri['pileuporiginal'])
