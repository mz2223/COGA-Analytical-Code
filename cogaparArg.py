import pandas as pd
import numpy as np
import os
import time
import multiprocessing as mp
import sys



numprocesses = 15

startchr = int(sys.argv[1])
endchr = int(sys.argv[2])
# startchr = input('Start CHR:')
# endchr = input('End CHR:')
dfcoga = pd.read_csv('/N/project/mmge_audprs/coga/snps/chr1-22.csv',dtype=str)
# print(dfcoga)
# dfall = pd.DataFrame()
dfall = None
d = './'
start_time = time.time()


def process_file(outpath):
    print(outpath)    
    df = pd.read_csv(outpath)

    split_columns = df['snp'].str.split(r'[._]', expand=True)
    df[['CHR', 'BP', 'C', 'D']] = split_columns.iloc[:, :4]
    df['CHR'] = df['CHR'].str[1:]

    df['EA'] = df['snp'].str.split('_', expand=True)[1]
    df['NEA'] = np.where(df['EA'] == df['D'], df['C'], np.where(df['EA'] == df['C'], df['D'], np.nan))
    df['EAF'] = (df['n1'] + 2*df['n2']) / (2*df['n0'] + 2*df['n1'] + 2*df['n2'])
    # print(df)
    # print(df.shape)
    dfm = pd.merge(df, dfcoga, how='left',on=['CHR','BP'])
    # print(dfm)
    match = ((dfm['EA'] == dfm['Effect.Allele']) & (dfm['NEA'] == dfm['non_Effect.Allele'])) | ((dfm['EA'] == dfm['non_Effect.Allele']) & (dfm['NEA'] == dfm['Effect.Allele']))
    dfm = dfm[match]
    dfm['N'] = dfm['n0'] + dfm['n1'] + dfm['n2']
    dfm = dfm.drop(['nd0','nd1','nd2','n0','n1','n2','miss.0','miss.1','miss.diff.p','chisq','C','D','df','model','remark'], axis=1)       

    dfm = dfm.drop(['Effect.Allele','non_Effect.Allele','snp'], axis=1)
    return dfm


outpaths = []
for chr in range(int(startchr),int(endchr)+1):
    print('chr =',chr)
    directory = f'{d}chr{chr}/'
    ls = os.listdir(directory)
    for outfile in ls:
        outpath = directory + outfile
        outpaths.append(outpath)
print(outpaths)
# cl = mp.Pool(processes=mp.cpu_count()) 
cl = mp.Pool(processes=numprocesses) 
results = cl.map(process_file, outpaths)

cl.close()
cl.join()
dfall = pd.concat(results, ignore_index=True)

dfall.to_csv(d + f'allout{startchr}-{endchr}.csv', index=False)
end_time = time.time()
time_taken_seconds = end_time - start_time
time_taken_minutes = time_taken_seconds / 60
print(f"You took {time_taken_minutes:.2f} minutes.")
