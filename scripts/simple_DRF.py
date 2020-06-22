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
import argparse
#import cProfile
import sys

#what portion of pile < 20x
#what portion of coding regions is < 20x
#how much is lost due to soft clipping?

def main ():

    parser = argparse.ArgumentParser(description='dark region finder')

    required = parser.add_argument_group(
            'Required',
            'ref, pileup, and output location')

    required.add_argument(
            '-r',
            '--ref',
            type=str,
            help='ref 37 or 38 [38]',
            default = '/reference/genomes/GRCh38_no_alt_analysis_set/indexes/BWAKIT_0.7.12/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')

    #required.add_argument(
    #        '-g',
    #        '--genes',
    #        type=str,
    #        help='Gene coords [38]',
    #        default = '/working/lab_nicw/kevinH/projects/dark_regions/GRCh38.bedish')

    required.add_argument(
            '-n',
            '--pile',
            type=str,
            help='pileup')

    required.add_argument(
            '-o',
            '--output',
            type=str,
            help= 'location to write output')

    optional = parser.add_argument_group(
            'Optional',
            'threads and chroms')

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
            help='To use threading the genome is chunked into chunks of window size. the bigger the better but more RAM needed [32000]',
            default=32000)

    args = parser.parse_args()


    if args.chrom == 'all':
        chroms = ['chr' + str(i+1) for i in range(22)] + ['chrX', 'chrY']
    else:
        chroms = args.chrom.split(',')
    chunks = chunk_ref(args, chroms)
    print (f'len chunks {len(chunks)}')

    p_count = 0
    total = 0
    cp_count = 0
    ctotal = 0

    with Pool(processes=args.threads) as pool:
        tmp = [(args, chunk) for chunk in chunks]
        res = pool.map(doit, tmp)
        for p_c, t, cp_c, ct in res:
            p_count += p_c
            total += t
            cp_count += cp_c
            ctotal += ct

    print(d)

    #print (p_count, total, cp_count, ctotal)
    with open(args.output + '_counts.txt','w') as fout:
        fout.write('\t'.join(['p_count', 'total', 'cp_count', 'ctotal']) + '\n')
        fout.write('\t'.join(list(map(str, [p_count, total, cp_count, ctotal]))) + '\n')


def doit(tup):
    
    args, genomic_region = tup
    return genomic_region.run(args)

class GenomicRegion:

    def __init__(self, chrom, start, end, seq, pile):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.seq = seq
        self.pile = pile
        #self.hugo_positions = set([])

    def name(self):

        print (f'{str(self.chrom)}:{str(self.start)}:{str(self.end)}')


    def run(self, args):

        #self.find_hugo(args)
        pile_d, coding_pile_d = self.parse_tab(pysam.TabixFile(self.pile))

        p_count, total = self.less_than_20x(pile_d)
        cp_count, ctotal = self.less_than_20x(coding_pile_d)

        print('gg',self.chrom, self.start, self.end, p_count, total, cp_count, ctotal)
        return p_count, total, cp_count, ctotal

    def less_than_20x(self, pile_d, lower_bound = 5):

        #assert len(tumour_d) == len(normal_d)
        p_count = 0
        
        for pos in pile_d:
            if pile_d.get(pos) < lower_bound:
                p_count +=1

        return p_count, len(pile_d)

    def parse_tab(self, tabixfile):


        keys = ['seq', 'pos', 'ref', 'reads', 'res', 'qual']
        d = {}
        hd = {}
        for pos in tabixfile.fetch(self.chrom, self.start, self.end):
            tmp_d = dict(zip(keys, pos.split()))
            d[tmp_d.get('pos')] = int(tmp_d.get('reads'))
            #if int(tmp_d.get('pos')) in self.hugo_positions:
            #    hd[tmp_d.get('pos')] = int(tmp_d.get('reads'))

        return d, hd

    #def find_hugo(self, args):

    #    chunk_positions = set(range(self.start, self.end))
    #    with open(args.genes, 'r') as fin:
    #        fin.readline()
    #        fin.readline()
    #        for line in fin:
    #            bits = line.strip().split('\t')
    #            #print (bits, self.chrom)
    #            if bits[1] == self.chrom:
    #                #print (bits)
    #                for exon_start, exon_end in zip(map(int, bits[5].split(',')[:-1]), map(int, bits[6].split(',')[:-1])):
    #                    #print ('gg',exon_start, exon_end, exon_end in chunk_positions, len(chunk_positions))
    #                    if exon_start in chunk_positions or exon_end in chunk_positions:
    #                        #print ('jjj',exon_start, exon_end, exon_end in chunk_positions, len(chunk_positions))
    #                        self.hugo_positions = self.hugo_positions.union(set(range(exon_start, exon_end)))


def chunk_ref(args, chroms):

    '''
    Split ref into 1Mb chunks for threading
    '''
    #make server mode where just takes one chrom and scatters with WDL

    chunks = []
    size = args.window_size #just for testing, put back to 1mb
    total_len = 0
    for record in SeqIO.parse(args.ref, 'fasta'):
        if record.id in chroms:
            for i in range(0, len(record.seq), size):
                if len(record.seq) > i + size:
                    end = i + size
                else:#end of chrom
                    end = len(record.seq)
                chunks.append(GenomicRegion(
                    str(record.id),
                    i,
                    end,
                    str(record.seq)[i:end],
                    args.pile))
                total_len+=len(str(record.seq)[i:end])

    assert total_len == sum([genomic_region.end - genomic_region.start for genomic_region in chunks])

    return chunks

if __name__ == '__main__':
    main()
    #cProfile.run('main()')
