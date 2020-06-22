#!/usr/bin/env python3

import pandas as pd
import pysam
from Bio import SeqIO
import numpy as np


#def chunk_ref(ref, size, chroms):

#    '''
#    Split ref into chunks for threading
#    '''

#    chunks = []
#    #size = size #just for testing, put back to 1mb
##    total_len = 0
#    for record in SeqIO.parse(ref, 'fasta'):
#        if record.id in chroms:
#            for i in range(0, len(record.seq), size):
#                if len(record.seq) > i + size:
#                    end = i + size
#                else:#end of chrom
##                    end = len(record.seq)
#                chunks.append((record.id, i, i+size))
#                
#    return chunks

chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']

#chunks = chunk_ref('/reference/genomes/GRCh38_no_alt_analysis_set/indexes/BWAKIT_0.7.12/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa', 1000000, chroms)

def chrom_lens(ref, chroms):

    chroml = []

    for record in SeqIO.parse(ref, 'fasta'):
        if record.id in chroms:
            chroml.append(len(record.seq))

    return chroml

chroml = chrom_lens('/reference/genomes/GRCh38_no_alt_analysis_set/indexes/BWAKIT_0.7.12/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa', chroms)

#print(chroml)

#Ok the above works, I can incorporate that into the main script pretty easily

df = pd.read_csv('multi_chrom_bigtest_raw_data.tsv', sep='\t', index_col=0)

#print(df.head())


drs = []

for chrom, tmp_df in df.groupby('chrom'):
    #print(len(tmp_df))
    drs.append(len(tmp_df))



#This above now works too, I can incorporate both of these into the existing script pretty easily

#This for loop below works but only every 6th value is the one I want

#for chrom, chromlen in chroml:
#    for chrom, drl in drs:
#        print(drl / chromlen * 100)


coverage = [x/y*100 for x, y in zip(drs, chroml)]
#for chrom, value in zip(chroms, coverage):
#    print(chrom, value)

#print(drs / chroml * 100)

#Writing stats

with open('TEST' + '_summary_stats.txt', 'w') as fout: #Will need to change to args.output
    fout.write('#summary statistics for' + str(chroms) + '\n')
    fout.write('chrom' + '\t' + 'dark_region_coverage' + '\t' + 'largest_DR' + '\t' + 'mean_DR_size' + '\t' + 'DR_size_sd' + '\n')
    for chrom, value in zip(chroms, coverage):
        fout.write(chrom + '\t' + str(value) + '%' + '\n')




