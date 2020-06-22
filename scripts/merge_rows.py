#!/usr/bin/env python3

import pandas as pd
import sys
import time

#This was solved by posing the following question: https://stackoverflow.com/questions/61490531/merging-pandas-dataframe-rows-to-create-ranges/61493472#61493472

def findEnd(df, index):
    while index < len(df)-1:
        if(df.iloc[index]['end']+1) == df.iloc[index+1]['start']:
            index+=1
        else: return(df.iloc[index]['end'], index)
    return (df.iloc[index]['end'], index)


df = pd.read_csv('chr22_test_processed.tsv', sep='\t', index_col=0)

time1 = time.time()

lst = []
i = 0
genLen = len(df)
#Traverse entire dataframe
while i < genLen:  
    #Check if we have at least one more row
    if i < genLen-1: 
        #Check the next row is the same chrom
        if(df.iloc[i]['chrom'] == df.iloc[i+1]['chrom']):
            start = df.iloc[i]['start']
            end,i = findEnd(df,i)
            lst.append([df.iloc[i]['chrom'],start,end])
        else:
            #if the next row is a different 
            lst.append(list(df.iloc[i]))
    elif i == genLen -1:
        lst.append(list(df.iloc[i]))
    i+=1

chrom = pd.DataFrame(lst,columns=['chrom','start','end'])

#chrom.to_csv('SO_dark_regions.tsv', sep='\t')
print(chrom)

time2 = time.time()
print('Time to join rows, SO method:', time2 - time1)
