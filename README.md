# simpler_DRF
A python script to find genomic dark regions in short- and long-read sequencing data






Options:

usage: simpler_DRF.py [-h] [-r REF] [-s METHOD] [-f INPUT] [-d DEPTH]
                      [-p PHRED] [-m MAP] [-o OUTPUT] [-x CHROM] [-t THREADS]
                      [-w WINDOW_SIZE]

dark region finder (returns raw data (TSV), bed file with dark regions and
annotations, and summary stats (txt file))

optional arguments:
  -h, --help            show this help message and exit

Required:
  ref, pileup, and output location

  -r REF, --ref REF     ref 37 or 38 [38]
  -s METHOD, --method METHOD
                        sequence type [short] or [long], need to set depth,
                        mapq, and phred filters for the given method
  -f INPUT, --input INPUT
                        pileup
  -d DEPTH, --depth DEPTH
                        read depth to filter on [default 10]
  -p PHRED, --phred PHRED
                        phred quality to filter on [default 20]
  -m MAP, --map MAP     mapping quality to filter on [default 10]
  -o OUTPUT, --output OUTPUT
                        location to write output

Optional:
  threads, chroms, window size

  -x CHROM, --chrom CHROM
                        Which chromosomes to query. comma,separated,list or
                        [all]
  -t THREADS, --threads THREADS
                        Threads [5]
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        To use threading the genome is chunked into chunks of
                        window size. the bigger the better but more RAM needed
                        [1000000]
