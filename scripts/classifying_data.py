#!/usr/bin/env python3

import pandas as pd
import numpy as np

df = pd.read_csv('class_test.tsv', sep='\t', index_col=0)
df.fillna(0, inplace=True)
df.map_av.astype(int)

#print(df.head(100))

classes = [
    (df['reads'] <= 10) & (df['map_av'] <= 10) & (df['phred_av'] <= 7),
    (df['map_av'] <= 10) & (df['phred_av'] <= 7),
    (df['reads'] <= 10) & (df['phred_av'] <= 7),
    (df['reads'] <= 10) & (df['map_av'] <= 10),
    df['reads'] <= 10,
    df['map_av'] <= 10,
    df['phred_av'] <= 7      
]

outputs = [
    'depth, mapq, phred',
    'mapq, phred',
    'depth, phred',
    'depth, mapq',
    'depth',
    'mapq',
    'phred' 
]

df['dark_reason'] = np.select(classes, outputs)
    
#print(df.head(100))

#Now trying to add in the merging functionality 

df.rename(columns = {'pos':'start'}, inplace=True) #renaming position 'start' and adding an 'end' column
df['end'] = df['start']
    
df = df[['chrom', 'start', 'end', 'dark_reason']] #Dropping other columns, the BED file just needs chromosome, start, and end position.

#print(df.dtypes)

df = df.astype({"start":'int64', "end":'int64'}) #Needed to include this because for some reason start and end are just "objects" before this - really weird considering that the code below works perfectly fine in isolation without the need to convert anything

print(df.head(30))

with open('dark_regions.bed', 'w') as fout:
    for chrom, tmp_df in df.groupby('chrom'):
        df_end = tmp_df[~((tmp_df['end'].shift(0) == tmp_df['end'].shift(-1)-1))]
        df_start = tmp_df[~((tmp_df['start'].shift(0) == tmp_df['start'].shift(+1)+1))]
        for start, end in zip(df_start['start'], df_end['end']):
            tmp_tmp = tmp_df[(tmp_df.start >= start) & (tmp_df.start <= end)]
            if 'depth, mapq, phred' in tmp_tmp.dark_reason:
                labels = 'depth, mapq, phred'
            else:
                labels = set(' '.join(list(tmp_df[(tmp_df.start >= start) & (tmp_df.start <= end)].dark_reason)).replace(',','').split(' '))
            #print (tmp_tmp.dark_reason)
            #print (set(tmp_tmp.dark_reason))
            #labels = set(' '.join(list(tmp_df[(tmp_df.start >= start) & (tmp_df.start <= end)].dark_reason)).replace(',','').split(' '))          
            print(chrom, start, end, labels)