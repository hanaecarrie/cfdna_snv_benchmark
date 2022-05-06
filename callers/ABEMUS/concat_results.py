import os
import pandas as pd
import sys

if __name__ == "__main__":
    dilutionseriesfolder = sys.argv[1]
    outdir = sys.argv[2]
    print(dilutionseriesfolder)
    print(outdir)
    if not os.path.exists(os.path.join(outdir, 'results')):
        os.mkdir(os.path.join(outdir, 'results'))
    for dil in [d for d in os.listdir(dilutionseriesfolder) if d.endswith('T') or d.endswith('x')]:
        print(dil)
        if not os.path.exists(os.path.join(outdir, 'results', dil)):
            os.mkdir(os.path.join(outdir, 'results', dil))
        nchunks = len([f for f in os.listdir(outdir) if f.startswith('chunk_')])
        for j in range(1,4):
            print(j)
            df_list = [pd.read_csv(os.path.join(outdir, 'chunk_'+str(i).rjust(2, '0'), 'Results',  dil, 'pmtab_F'+str(j)+'_'+dil+'.tsv'), sep="\t") for i in range(nchunks) if os.path.exists(os.path.join(outdir, 'chunk_'+str(i).rjust(2, '0'), 'Results',  dil, 'pmtab_F'+str(j)+'_'+dil+'.tsv'))]
            if df_list:
                res_df = pd.concat(df_list)
                print(res_df.head())
                res_df.to_csv(os.path.join(outdir, 'results',  dil, 'pmtab_F'+str(j)+'_'+dil+'.csv'))
            else:
                print('no results to concat for pmtab_F'+str(j)+'_'+dil)
        print('optimal R')
        df_list = [pd.read_csv(os.path.join(outdir, 'chunk_'+str(i).rjust(2, '0'), 'Results',  dil, 'pmtab_F3_optimalR_'+dil+'.tsv'), sep="\t") for i in range(nchunks) if os.path.exists(os.path.join(outdir, 'chunk_'+str(i).rjust(2, '0'), 'Results',  dil, 'pmtab_F3_optimalR_'+dil+'.tsv'))]
        if df_list:
            res_df = pd.concat(df_list)
            print(res_df.head())
            res_df.to_csv(os.path.join(outdir, 'results',  dil, 'pmtab_F3_optimalR_'+dil+'.csv'))

