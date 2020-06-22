#!/usr/bin/env python3
from multiprocessing import Pool, TimeoutError
from glob import glob
import gzip
import pysam
from Bio import SeqIO
from collections import Counter, defaultdict
import scipy.stats as stats
import operator
import pandas as pd
import numpy as np
import argparse
import sys
import time
from itertools import chain

def main ():

    parser = argparse.ArgumentParser(description='dark region finder (returns raw data (TSV), bed file with dark regions and annotations, and summary stats (txt file))')

    required = parser.add_argument_group(
            'Required',
            'ref, pileup, and output location')

    required.add_argument(
            '-r',
            '--ref',
            type=str,
            help='ref 37 or 38 [38]',
            default = '/reference/genomes/GRCh38_no_alt_analysis_set/indexes/BWAKIT_0.7.12/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')

    required.add_argument(
            '-f',
            '--input',
            type=str,
            help='pileup')

    #required.add_argument(
    #        '-d',
    #        '--depth',
    #        type=int,
    #        default=10,
    #        help='read depth to filter on [default 10]')
    
    #required.add_argument(
    #        '-p',
    #        '--phred',
    #        type=int,
    #        default=20,
    #        help='phred quality to filter on [default 20]')

    #required.add_argument(
    #        '-m',
    #        '--map',
    #        type=int,
    #        default=10,
    #        help='mapping quality to filter on [default 10]')

    #required.add_argument(
    #        '-o',
    #        '--output',
    #        type=str,
    #        help= 'location to write output')

    optional = parser.add_argument_group(
            'Optional',
            'threads, chroms, window size')

    optional.add_argument(
            '-x',
            '--chrom',
            type=str,
            help='Which chromosomes to query. comma,separated,list or [all]',
            default='all')

    optional.add_argument(
            '-t',
            '--threads',
            type=int,
            help='Threads [5]',
            default=5)

    optional.add_argument(
            '-w',
            '--window_size',
            type=int,
            help='To use threading the genome is chunked into chunks of window size. the bigger the better but more RAM needed [1000000]',
            default=1000000)

    args = parser.parse_args()

    start1 = time.time()

    '''
    1. Chunking reference genome to split input apart for mulitprocessing
    '''

    if args.chrom == 'all':
        chroms = ['chr' + str(i+1) for i in range(22)] + ['chrX', 'chrY']
    else:
        chroms = args.chrom.split(',')
    chunks = chunk_ref(args.ref, args.window_size, chroms)[20:25] #REMEMBER TO RESET ONCE TESTING IS DONE!!!!!!
    print(f'len chunks {len(chunks)}')

    #print(chunks)
    #chroms = str(chroms) #Need to make chroms into strings so that it's accepted by find_DRF
    #print(type(chroms))
    #print(chroms)

    #chroml = chrom_lens(args.ref, chroms) #Just temporarily adding this function in here, would like to incorporate into the chunk_ref function, see below

    #end1 = time.time()
    #print('Chunking time:', end1 - start1)

    d={}
    with Pool(processes=args.threads) as pool:
        tmp = [(args, chunk) for chunk in chunks]
        res = pool.map(find_DRF, tmp)
        #print(res)
        for sub in res:
            print(sub)
            d = {**d, **sub}#needs work
        
    #    print(d)

        #print(len(d))
    
    #end2 = time.time()
    #print('Time to aggregate data:', end2 - end1)

'''

#############################

    result = defaultdict(lambda: defaultdict(int))
    keys = ['seq', 'pos', 'ref', 'reads', 'res', 'qual'] 

    tabixfile = pysam.TabixFile(args.input)
    pos = (line for line in tabixfile.fetch('chr22', 16000000, 17000000))
    list_line = (s.split() for s in pos)
    pile = (dict(zip(keys, data)) for data in list_line)

    key = (pi['seq'] + ':' + pi['pos'] for pi in pile) #Should I create keys here or wait until I am in the main for loop?
    
    res_str = ((pi['seq'], pi['pos'], pi['reads'], pi['res'], pi['qual']) for pi in pile)

    for seq, pos, reads, res, qual in res_str:
        #print(seq, pos, res, qual)

        key = str(seq + ':' + pos) #Should I generate the keys here or before I enter the for loop?
        
        result[key]['reads'] = reads

        for score in qual:
            #print(score)
            tmp_l = []
            tmp_l.append(ord(score)-33)
            result[key]['phred_av'] = sum(tmp_l)/len(tmp_l) 
            
        for i, nuc in enumerate(res):
            if nuc =='^':
                qual = ord(res[i+1])-33
                #print(qual)
                for j in range(int(pos), int(pos)+150):
                    #print(j)
                    key =  str(seq) + ':' + str(j)
                    result[key]['map_av'] += qual

    result = dict(result)
    print(result)
    #print(result, '\n')

'''
########################

    



def find_DRF(tup):

    args, chunk = tup
    chrom, start, end = chunk
    
    '''
    This function extracts the data from the pileup that is necessary for later filtering, it also assigns mapq values to each position by finding the last 
    mapq score and then giving that to all positions within 150 bp upstream of that position. This function returns a dictionary with each position and the relevant
    information at that position.

    *important to note here that 'find_DRF' is a bit of a misnomer, as that step is actually performed later, the dictionary that is returned from this function
    contains ALL positions, not just those that are 'dark'
    '''

    result = defaultdict(lambda: defaultdict(int))
    keys = ['seq', 'pos', 'ref', 'reads', 'res', 'qual'] 

    tabixfile = pysam.TabixFile(args.input)
    pos = (line for line in tabixfile.fetch(chrom, start, end))
    list_line = (s.split() for s in pos)
    pile = (dict(zip(keys, data)) for data in list_line)

    key = (pi['seq'] + ':' + pi['pos'] for pi in pile) #Should I create keys here or wait until I am in the main for loop?
    
    res_str = ((pi['seq'], pi['pos'], pi['reads'], pi['res'], pi['qual']) for pi in pile)

    for seq, pos, reads, res, qual in res_str:
        #print(seq, pos, res, qual)

        key = str(seq + ':' + pos) #Should I generate the keys here or before I enter the for loop?
        
        result[key]['reads'] = reads

        for score in qual:
            #print(score)
            tmp_l = []
            tmp_l.append(ord(score)-33)
            result[key]['phred_av'] = sum(tmp_l)/len(tmp_l) 
            
        for i, nuc in enumerate(res):
            if nuc =='^':
                qual = ord(res[i+1])-33
                #print(qual)
                for j in range(int(pos), int(pos)+150):
                    #print(j)
                    key =  str(seq) + ':' + str(j)
                    result[key]['map_av'] += qual


    return dict(result)

#    result = defaultdict(lambda: defaultdict(int))
#    keys = ['seq', 'pos', 'ref', 'reads', 'res', 'qual'] 
#    tabixfile = pysam.TabixFile(args.input)

#    pos = (line for line in tabixfile.fetch(chrom, start, end))
#    #pile = dict(zip(keys, pos.split()))

    


#    for pos in tabixfile.fetch(chrom, start, end):
#        pile = dict(zip(keys, pos.split()))
#        #print (chunk, 'pile',pile)
#        tmp_l = []
#        for elem in pile.get('qual'):
#                tmp_l.append(ord(elem)-33)
#        key = pile.get('seq') +':' + pile.get('pos')
#        result[key]['phred_av'] = sum(tmp_l)/len(tmp_l)
#        result[key]['reads'] = int(pile.get('reads'))
#        res_str = pile.get('res')
#        for i, nuc in enumerate(res_str):
#                if nuc =='^':
#                        qual = ord(res_str[i+1])-33
#                        read_start = int(pile.get('pos'))
#                        for j in range(read_start, read_start+150):
#                                key = pile.get('seq') +':' + str(j)
#                                result[key]['map_av'] += qual

    #print (chunk, result)
#    return dict(result)


def chunk_ref(ref, size, chroms):

    '''
    Split ref into chunks for threading
    '''

    chunks = []
    #chroml = [] #Can't include getting chrom lengths here because it messes up multiprocessing
    #size = size #just for testing, put back to 1mb
    total_len = 0
    for record in SeqIO.parse(ref, 'fasta'):
        if record.id in chroms:
            #chroml.append(len(record.seq))
            for i in range(0, len(record.seq), size):
                if len(record.seq) > i + size:
                    end = i + size
                else:#end of chrom
                    end = len(record.seq)
                chunks.append((record.id, i, i+size))
                
    return chunks #, chroml

if __name__ == '__main__':
    main()
    