#!/usr/bin/env python3

import pandas as pd
import sys
import time

df = pd.read_csv('mc_processed.tsv', sep='\t', index_col=0)

time1 = time.time()

with open('TEST_file.bed', 'w') as fout:
    for chrom, tmp_df in df.groupby('chrom'):
        df_end = tmp_df[~((tmp_df['end'].shift(0) == tmp_df['end'].shift(-1)-1))]
        df_start = tmp_df[~((tmp_df['start'].shift(0) == tmp_df['start'].shift(+1)+1))]
        for start, end in zip(df_start['start'], df_end['end']):
            print(chrom, start, end)
            fout.write(chrom + '\t' + str(start) + '\t' + str(end) + '\n')




time2 = time.time()

print('liam time', time2-time1)


