# simpler_DRF
A python script to find genomic dark regions in short- and long-read sequencing data

simpler_DRF.py can be run from the command line, and is inteneded for use in a cluster/HPC environment. **Due to memory issues it is advisable to run simpler_DRF on a chromosome-by-chromosome basis.** simpler_DRF.py supports multithreading, however performance gains are marginal as only the first step of the program is multithreaded. Generally 3 threads/cores is enough. 

As an input simpler_DRF.py takes the following:

* a .pileup file that has been bgzipped and indexed
* a reference genome fasta file
* the sequencing method used:
 * Use "short" for short-read methods such as Illumina HiSeq, etc.
 * Use "long" for long-read methods such as ONT MinION/PromethION, etc.
* filtering criteria:
 * read depth (defaults to 10)
 * PHRED score (defaults to 20)
 * MAPQ (defaults to 10, do not include this option if sequencing method used is "long")
* Location to write output file


simpler_DRF will return 3 files: 

* A TSV file containing the raw data (i.e. all positions that simpler_DRF identifies as 'dark' according to filtering criteria used)
* A file in BED format listing the genomic coordinates that have been identified as dark, as well as an additional column containing the reason they were identified as dark
* A txt file containing the total % dark regions for the given chromosome(s)

A workflow with simpler_DRF could look like:

short- or long-read BAM files > create .pileup file using samtools mpileup > bgzip and index .pileup file > run simpler_DRF.py > run exploratory analysis in DRF_summary.Rmd file. 




Options:
```
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
```
