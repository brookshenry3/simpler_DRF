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
import functools
import logging
import struct

def main ():

    parser = argparse.ArgumentParser(description='multi-threaded dark region finder for short read sequencing data (returns raw data (TSV), bed file with dark regions and annotations, and summary stats (txt file))')

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

    required.add_argument(
            '-d',
            '--depth',
            type=int,
            default=10,
            help='read depth to filter on [default 10]')
    
    required.add_argument(
            '-p',
            '--phred',
            type=int,
            default=20,
            help='phred quality to filter on [default 20]')

    required.add_argument(
            '-m',
            '--map',
            type=int,
            default=10,
            help='mapping quality to filter on [default 10]')

    required.add_argument(
            '-o',
            '--output',
            type=str,
            help= 'location to write output')

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

    patch_mp_connection_bpo_17560()

    start1 = time.time()

    '''
    1. Chunking reference genome to split input apart for mulitprocessing
    '''

    if args.chrom == 'all':
        chroms = ['chr' + str(i+1) for i in range(22)] + ['chrX', 'chrY']
    else:
        chroms = args.chrom.split(',')
    chunks = chunk_ref(args.ref, args.window_size, chroms) #[20:40] #Comment out last bit for full test
    print(f'len chunks {len(chunks)}')

    #print(chunks)
    chroms = str(chroms) #Need to make chroms into strings so that it's accepted by find_DRF
    #print(type(chroms))
    #print(chroms)

    #chroml = chrom_lens(args.ref, chroms) #Just temporarily adding this function in here, would like to incorporate into the chunk_ref function, see below

    end1 = time.time()
    print('Chunking time:', end1 - start1)

    '''
    2. Parsing the input pileup file and just extracting the data that will be used to identify dark regions
    '''

    d={}
    with Pool(processes=args.threads) as pool:
        tmp = [(args, chunk) for chunk in chunks]
        res = pool.map(find_DRF, tmp)
        #print(res)
        for sub in res:
                d = {**d, **sub}#needs work
        
        #print(d)

        #print(len(d))
    
    end2 = time.time()
    print('Time to aggregate data:', end2 - end1)

    '''
    3. Converting the dictionary generated in the last step into a pandas dataframe for easier manipulation 
    & reordering/cleaning things up for filtering and merging
    '''

    df = pd.DataFrame.from_dict(d, orient = 'index')

    #print (df)

    df['chrom_pos'] = df.index #creating a new column called 'chrom_pos' from the index, this will be split in the next step and then the index will eventually be returned to just integer values

    #print(df)

    df[['chrom', 'pos']] = df.chrom_pos.str.split(':', expand=True) #splitting chrom_pos into two new data columns 
    df = df[['chrom', 'pos', 'reads', 'map_av', 'phred_av']] #reordering the columns and dropping the ones that I don't need

    end3 = time.time()
    print('Time to convert to DF:', end3 - end2)

    #print(df.get('chrom'))

    '''
    4. Finding/filtering for dark regions, and saving outputs to TSV and BED files
    '''

    with Pool(processes=args.threads) as pool:
        tup = [(args, df, chrom) for chrom in chroms]
        #print(tmp)
        res = pool.map(big_one, tup)
        #print(res)
        #for sub in res:
        #    d = {**d, **sub}#needs work

    end4 = time.time()
    print('Time for big_one:', end4 - end3)

    

#Functions 

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

    for pos in tabixfile.fetch(chrom, start, end):
        pile = dict(zip(keys, pos.split()))
        #print (chunk, 'pile',pile)
        tmp_l = []
        for elem in pile.get('qual'):
                tmp_l.append(ord(elem)-33)
        key = pile.get('seq') +':' + pile.get('pos')
        result[key]['phred_av'] = sum(tmp_l)/len(tmp_l)
        result[key]['reads'] = int(pile.get('reads'))
        res_str = pile.get('res')
        for i, nuc in enumerate(res_str):
                if nuc =='^':
                        qual = ord(res_str[i+1])-33
                        read_start = int(pile.get('pos'))
                        for j in range(read_start, read_start+150):
                                key = pile.get('seq') +':' + str(j)
                                result[key]['map_av'] += qual

    #print (chunk, result)
    return dict(result)


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


def chrom_lens(ref, chroms):

    '''    
    Getting the chromosome lengths for use in summary statistics. 
    I'd like to incorporate this into chunk_ref but I get the following error when I do include it in chunk ref:
    
    multiprocessing.pool.RemoteTraceback: 
    
    """
    Traceback (most recent call last):
    File "/software/python/Python-3.6.1/lib/python3.6/multiprocessing/pool.py", line 119, in worker
    result = (True, func(*args, **kwds))
    File "/software/python/Python-3.6.1/lib/python3.6/multiprocessing/pool.py", line 44, in mapstar
    return list(map(*args))
    File "simpler_DRF.py", line 275, in find_DRF
    chrom, start, end = chunk
    ValueError: not enough values to unpack (expected 3, got 1)
    """

    '''

    chroml = []

    for record in SeqIO.parse(ref, 'fasta'):
        if record.id in chroms:
            chroml.append(len(record.seq))

    return chroml


def big_one(tup):

    args, df, chromsome = tup

    
    #print(df)
    #df.reset_index(drop=True, inplace=True)

    '''
    4.1: Filtering to find positions that are considered 'dark' based on depth, mapq, and phred score (add in further filters here?).
    & saving 'raw' data into a tsv that shows all three of the values for the above filters at each position
    '''

    #Filtering
    df = df[(df['reads'] <= args.depth) | (df['map_av'] <= args.map) | (df['phred_av'] <= args.phred)] #The actual filtering step to find 'dark' regions

    #end3 = time.time()
    #print('Time to filter:', end3 - end2)

    df.reset_index(drop=True, inplace=True) #resetting the index so that it's just integers, this is important for the iteration that occurs later on in step 5

    #Saving the raw data file
    #print('Saving raw data...')
    df.to_csv(args.output + '_raw_data.tsv', sep='\t') #for full run add in: , compression='gzip'

    #end4 = time.time()
    #print('Time to save raw data:', end4 - end3)

    #Getting length of dark regions below for use in summary stats (see section 7)

    #drs = []
    #for chrom, tmp_df in df.groupby('chrom'):
    #    #print(len(tmp_df))
    #    drs.append(len(tmp_df))


    '''
    4.2: Classifying dark regions

    Here each postion in the above dataframe is being classified according to why it is dark.
    There is probably a better/more dynamic way to do this...
    Might want this step to go BEFORE the raw data is saved 

    '''

    #print('Classifying positions...')

    df.fillna(0, inplace=True)

    classes = [
        (df['reads'] <= args.depth) & (df['map_av'] <= args.map) & (df['phred_av'] <= args.phred),
        (df['map_av'] <= args.map) & (df['phred_av'] <= args.phred),
        (df['reads'] <= args.depth) & (df['phred_av'] <= args.phred),
        (df['reads'] <= args.depth) & (df['map_av'] <= args.map),
        df['reads'] <= args.depth,
        df['map_av'] <= args.map,
        df['phred_av'] <= args.phred      
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

    df['dark_class'] = np.select(classes, outputs)

    #end5 = time.time()
    #print('Time to classify positions:', end5 - end4)

    '''
    4.3: Joining rows that contain adjacent genomic positons & saving file
    
    Thanks goes to Liam for finding a solution to the joining rows problem
    '''

    #print('Joining adjacent rows and creating BED file...')

    df.rename(columns = {'pos':'start'}, inplace=True) #renaming position 'start' and adding an 'end' column
    df['end'] = df['start']
    
    df = df[['chrom', 'start', 'end', 'dark_class']] #Dropping other columns, the BED file just needs chromosome, start, and end position.

    #print(df.dtypes)

    df = df.astype({"start":'int64', "end":'int64'}) #Needed to include this because for some reason start and end are just "objects" before this - really weird considering that the code below works perfectly fine in isolation without the need to convert anything

    with open(args.output + '_dark_regions.bed', 'w') as fout:
        for chrom, tmp_df in df.groupby('chrom'):
            df_end = tmp_df[~((tmp_df['end'].shift(0) == tmp_df['end'].shift(-1)-1))]
            df_start = tmp_df[~((tmp_df['start'].shift(0) == tmp_df['start'].shift(+1)+1))]
            for start, end in zip(df_start['start'], df_end['end']):
                tmp_tmp = tmp_df[(tmp_df.start >= start) & (tmp_df.start <= end)]
                if 'depth, mapq, phred' in tmp_tmp.dark_class:
                    labels = 'depth, mapq, phred'
                else:
                    labels = set(' '.join(list(tmp_df[(tmp_df.start >= start) & (tmp_df.start <= end)].dark_class)).replace(',','').split(' '))
                #!!!! FINISH FIXING THIS PART USING THE EXAMPLE classifying_data.py
                #labels = set(' '.join(list(tmp_df[(tmp_df.start >= start) & (tmp_df.start <= end)].dark_class)).replace(',','').split(' '))      #Ok this works to get the annotations but seems to be extremely slow    
                fout.write(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + str(labels) + '\n')
                #print(chrom, start, end, labels)

    #end6 = time.time()
    #print('Time to join rows & save file:', end6 - end5)

logger = logging.getLogger()

def patch_mp_connection_bpo_17560():
    """Apply PR-10305 / bpo-17560 connection send/receive max size update

    See the original issue at https://bugs.python.org/issue17560 and 
    https://github.com/python/cpython/pull/10305 for the pull request.

    This only supports Python versions 3.3 - 3.7, this function
    does nothing for Python versions outside of that range.

    """
    patchname = "Multiprocessing connection patch for bpo-17560"
    if not (3, 3) < sys.version_info < (3, 8):
        logger.info(
            patchname + " not applied, not an applicable Python version: %s",
            sys.version
        )
        return

    from multiprocessing.connection import Connection

    orig_send_bytes = Connection._send_bytes
    orig_recv_bytes = Connection._recv_bytes
    if (
        orig_send_bytes.__code__.co_filename == __file__
        and orig_recv_bytes.__code__.co_filename == __file__
    ):
        logger.info(patchname + " already applied, skipping")
        return

    @functools.wraps(orig_send_bytes)
    def send_bytes(self, buf):
        n = len(buf)
        if n > 0x7fffffff:
            pre_header = struct.pack("!i", -1)
            header = struct.pack("!Q", n)
            self._send(pre_header)
            self._send(header)
            self._send(buf)
        else:
            orig_send_bytes(self, buf)

    @functools.wraps(orig_recv_bytes)
    def recv_bytes(self, maxsize=None):
        buf = self._recv(4)
        size, = struct.unpack("!i", buf.getvalue())
        if size == -1:
            buf = self._recv(8)
            size, = struct.unpack("!Q", buf.getvalue())
        if maxsize is not None and size > maxsize:
            return None
        return self._recv(size)

    Connection._send_bytes = send_bytes
    Connection._recv_bytes = recv_bytes

    logger.info(patchname + " applied")





if __name__ == '__main__':
    main()
    