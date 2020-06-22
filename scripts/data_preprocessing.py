#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv('multi_chrom_raw_data.tsv', sep='\t', index_col=0)

df.rename(columns = {'pos':'start'}, inplace=True) #renaming position 'start' and adding an 'end' column
df['end'] = df['start']

df['annot'] = 'NA'
    
df = df[['chrom', 'start', 'end']] #Dropping other columns, the BED file just needs chromosome, start, and end position.

    #print(df.dtypes)

df = df.astype({"start":'int64', "end":'int64'}) #Needed to include this because for some reason start and end are just "objects" before this - really weird considering that the code below works perfectly fine in isolation without the need to convert anything

df.to_csv('mc_processed.tsv', sep='\t')





def class_dark(df):
    for row in 