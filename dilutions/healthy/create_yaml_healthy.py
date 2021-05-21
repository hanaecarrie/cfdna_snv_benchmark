
import os
import pandas as pd

cohorts_names_paths = {'GIS1':'/mnt/projects/zwpoh/cfDNA/bulk/ccg/ccg_batch1/lpwgs/hg19_bam/healthy/bam/*.bam', \
'GIS2': '/mnt/projects/zwpoh/cfDNA/bulk/ccg/ccg_batch2/lpwgs/hg19_bam/healthy/*.bam', \
'CUHK': '/mnt/projects/zhug/cfDNA/data-cuhk/rmduplic-hkC*.bam', \
'CUHKCD': '/mnt/projects/zhug/cfDNA/data-cuhk-cd/healthy-*.bam', \
'GRAZ': '/mnt/projects/zhug/cfDNA/data-graz-control/healthy-*.bam'}

tsv_dict = {}

for cohort_name, cohort_path in cohorts_names_paths.items():
	dirname, filepattern = os.path.split(cohort_path)
	filepattern = filepattern.split('*')[0]
	if cohort_name == 'CUHK':
		yamlconfig = 'illumina'
	else:
		 yamlconfig = 'standard'
	yamlconfigfile = '/mnt/projects/skanderupamj/wgs/data/cfdna.crc/mix.samples/without.gatk.hap/tumor-only-'+yamlconfig+'.yaml'
	for file in [f for f in os.listdir(dirname) if f.startswith(filepattern) and f.endswith('.bam')]:
		print(file)
		tsv_dict['healthy_'+cohort_name+'_'+os.path.splitext(os.path.basename(file))[0]] = ['symlink::'+dirname, file, 'GIS_SRA', 'NA', 'MIX', 'NA', 'NA', yamlconfigfile]
	
tsv_pd = pd.DataFrame.from_dict(tsv_dict).T
tsv_pd.reset_index(inplace=True)
print(tsv_pd.head())

with open("/mnt/projects/carriehc/cfDNA/cfSNV/benchmark/patients_healthy.tsv", 'w') as write_tsv:
	write_tsv.write(tsv_pd.to_csv(sep='\t', index=False, header=False))


