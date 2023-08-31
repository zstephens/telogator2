#!/usr/bin/env python
import argparse
import gzip

from source.tg_reader import TG_Reader
from source.tg_util   import get_file_type, RC


def main(raw_args=None):
    parser = argparse.ArgumentParser(description='get_kmer_hits.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    parser.add_argument('-i', type=str, required=True,  metavar='input.fa',     nargs='*',              help="* Input reads (fa / fa.gz / fq / fq.gz / bam)")
    parser.add_argument('-o', type=str, required=True,  metavar='output.fa',                            help="* Output reads (fa / fa.gz / fq / fq.gz)")
    parser.add_argument('-k', type=str, required=False, metavar='TTAGGGTTAGGG', default='TTAGGGTTAGGG', help="Kmer to search")
    parser.add_argument('-m', type=int, required=False, metavar='5',            default=5,              help="Minimum hit count to output read")
    args = parser.parse_args()
    #
    IN_READS  = args.i
    OUT_READS = args.o
    KMER      = args.k
    KMER_RC   = RC(KMER)
    MIN_HITS  = args.m
    #
    total_reads = 0
    total_bases = 0
    output_reads = 0
    output_bases = 0
    #
    output_type = get_file_type(OUT_READS)
    if output_type[1]:
        f_out = gzip.open(OUT_READS, 'wt')
    else:
        f_out = open(OUT_READS, 'w')
    #
    for ifn in IN_READS:
        input_type = get_file_type(ifn)
        if input_type[0] not in ['fasta', 'fastq', 'bam']:
            print('Error: input must be fasta, fastq, or bam')
            exit(1)
        if output_type[0] not in ['fasta', 'fastq']:
            print('Error: output must be fasta or fastq')
            exit(1)
        if input_type[0] == 'fasta' and output_type[0] == 'fastq':
            print('Error: input is fasta and output is fastq')
            exit(1)
        if input_type[0] == 'bam' and output_type[0] == 'fastq':
            print('Warning: input is bam and output is fastq')
        #
        my_reader = TG_Reader(ifn)
        hit_count_dict = {}
        while True:
            (my_name, my_rdat, my_qdat) = my_reader.get_next_read()
            if not my_name:
                break
            total_reads += 1
            total_bases += len(my_rdat)
            count_fwd = my_rdat.count(KMER)
            count_rev = my_rdat.count(KMER_RC)
            if count_fwd not in hit_count_dict:
                hit_count_dict[count_fwd] = 0
            hit_count_dict[count_fwd] += 1
            if count_rev not in hit_count_dict:
                hit_count_dict[count_rev] = 0
            hit_count_dict[count_rev] += 1
            if count_fwd >= MIN_HITS or count_rev >= MIN_HITS:
                if output_type[0] == 'fastq':
                    f_out.write(f'@{my_name}\n{my_rdat}\n+\n{my_qdat}\n')
                elif output_type[0] == 'fasta':
                    f_out.write(f'>{my_name}\n{my_rdat}\n')
                output_reads += 1
                output_bases += len(my_rdat)
        my_reader.close()
    #
    f_out.close()
    #
    ####for k in sorted(hit_count_dict.keys()):
    ####    print(k, hit_count_dict[k])
    print('get_reads_with_kmer_hits.py summary:')
    print(f' - input reads:  {total_reads} ({total_bases} bp)')
    print(f' - output reads: {output_reads} ({output_bases} bp)')


if __name__ == '__main__':
    main()
