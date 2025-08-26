#!/usr/bin/env python
import argparse
import numpy as np
import pathlib
import pysam

from source.tg_kmer import get_telomere_kmer_density, get_telomere_regions, read_kmer_tsv
from source.tg_plot import plot_tel_signal
from source.tg_ref  import ReferenceFasta
from source.tg_util import exists_and_is_nonzero, posmax, RC

SUBTEL_BUFF = 500000
TEL_SEARCH  = 50000
READ_TYPE   = 'hifi'
DEFAULT_CHR = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

TEL_WINDOW_SIZE   = 100
P_VS_Q_AMP_THRESH = 0.500
MIN_TEL_SCORE     = 100

FASTA_WIDTH = 60


def main(raw_args=None):
    parser = argparse.ArgumentParser(description='make_telogator_ref.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    parser.add_argument('-i', type=str, required=True,                 metavar='* input.fa',   help="* Input T2T reference")
    parser.add_argument('-o', type=str, required=True,                 metavar='* output.fa',  help="* Output Telogator reference")
    parser.add_argument('-s', type=str, required=True,                 metavar='* sampname',   help="Sample name (prepends contig names)")
    parser.add_argument('-c', type=str, required=False, default='',    metavar='chr1,chr2,..', help="Sorted list of contigs (comma-delimited)")
    parser.add_argument('-r', type=str, required=False, default='',    metavar='chr1,chr2,..', help="Remap sorted contig names to this (comma-delimited)")
    parser.add_argument('-k', type=str, required=False, default='',    metavar='kmers.tsv',    help="Telomere kmers")
    parser.add_argument('-m', type=int, required=False, default=0,     metavar='0',            help="Minimum required TL at contig ends")
    parser.add_argument('--add-tel',    required=False, default=False, action='store_true',    help="Include masked tels as separate contigs")
    parser.add_argument('--plot',       required=False, default=False, action='store_true',    help="Output tel signal plots")
    args = parser.parse_args()

    IN_REF     = args.i
    OUT_REF    = args.o
    CHR_LIST   = args.c
    CHR_REMAP  = args.r
    KMER_FILE  = args.k
    MIN_END_TL = args.m
    SAMP_NAME  = args.s
    APPEND_TEL = args.add_tel
    MAKE_PLOTS = args.plot

    if CHR_LIST:
        SORTED_CHR = CHR_LIST.split(',')
        if CHR_REMAP:
            remap_names = CHR_REMAP.split(',')
            if len(SORTED_CHR) != len(remap_names):
                print('Error: number of contigs specified in -c and -r must be the same')
                exit(1)
            CHR_REMAP_DICT = {SORTED_CHR[n]:remap_names[n] for n in range(len(SORTED_CHR))}
    else:
        SORTED_CHR = DEFAULT_CHR
        CHR_REMAP_DICT = {SORTED_CHR[n]:SORTED_CHR[n] for n in range(len(SORTED_CHR))}

    if SAMP_NAME != '':
        if '_' in SAMP_NAME:
            print('Error: sample name cannot include an underscore')
            exit(1)
        SAMP_NAME += '_'

    OUT_DIR = '/'.join(OUT_REF.split('/')[:-1]) + '/'

    if KMER_FILE == '':
        print('using default telomere kmers.')
        sim_path  = pathlib.Path(__file__).resolve().parent
        KMER_FILE = str(sim_path) + '/resources/kmers.tsv'
    elif exists_and_is_nonzero(KMER_FILE):
        fn_suffix = KMER_FILE.split('/')[-1]
        print('using user-specified kmer list:', fn_suffix)
    else:
        print('Error: kmer list not found')
        exit(1)
    (KMER_METADATA, KMER_ISSUBSTRING, CANONICAL_STRINGS) = read_kmer_tsv(KMER_FILE, READ_TYPE)
    [KMER_LIST, KMER_COLORS, KMER_LETTER, KMER_FLAGS]    = KMER_METADATA
    SIGNAL_KMERS     = [n for n in KMER_LIST]
    SIGNAL_KMERS_REV = [RC(n) for n in KMER_LIST]

    REF_FASTA = ReferenceFasta(IN_REF)

    out_sequences = []
    tel_sequences = []

    for chr_name in SORTED_CHR:
        refseq = REF_FASTA.get_refseq(chr_name)
        if refseq is not None:
            print(f'{chr_name} {len(refseq)}bp')
            p_arm = refseq[:SUBTEL_BUFF]
            q_arm = refseq[len(refseq)-SUBTEL_BUFF:]
            p_arm_telsearch = p_arm[:TEL_SEARCH]
            q_arm_telsearch = q_arm[len(q_arm)-TEL_SEARCH:]
            both_dat = [(p_arm_telsearch, 'p'), (q_arm_telsearch, 'q')]
            #
            pq_tlens = []
            for (subtel_seq, pq) in both_dat:
                (td_p_e0, td_p_e1) = get_telomere_kmer_density(subtel_seq, SIGNAL_KMERS,     TEL_WINDOW_SIZE)
                (td_q_e0, td_q_e1) = get_telomere_kmer_density(subtel_seq, SIGNAL_KMERS_REV, TEL_WINDOW_SIZE)
                (p_vs_q_power, tel_regions) = get_telomere_regions(td_p_e0, td_p_e1, td_q_e0, td_q_e1, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH)
                #
                my_score = []
                my_tel_len = None
                tl_vals = [0,0]
                if pq == 'p':
                    for i in range(len(tel_regions)):
                        my_s = (tel_regions[i][1] - tel_regions[i][0])
                        if tel_regions[i][2] == 'p':
                            my_score.append(my_s)
                        elif tel_regions[i][2] == 'q':
                            my_score.append(my_s)
                        elif tel_regions[i][2] is None:
                            my_score.append(-my_s)
                    cum_score = np.cumsum(my_score)
                    max_i     = posmax(cum_score)
                    if cum_score[max_i] >= MIN_TEL_SCORE:
                        my_tel_len = int(tel_regions[max_i][1] + TEL_WINDOW_SIZE/2)
                        tl_vals = [my_tel_len, 0]
                elif pq == 'q':
                    for i in range(len(tel_regions)):
                        my_s = (tel_regions[i][1] - tel_regions[i][0])
                        if tel_regions[i][2] == 'q':
                            my_score.append(my_s)
                        elif tel_regions[i][2] == 'p':
                            my_score.append(my_s)
                        elif tel_regions[i][2] is None:
                            my_score.append(-my_s)
                    cum_score = np.cumsum(my_score[::-1])[::-1]
                    max_i     = posmax(cum_score)
                    if cum_score[max_i] >= MIN_TEL_SCORE:
                        my_tel_len = int(tel_regions[-1][1] - tel_regions[max_i][0] + TEL_WINDOW_SIZE/2)
                        tl_vals = [0, my_tel_len]
                if my_tel_len is None:  # if we got this far then I guess the contig had no tel repeats
                    my_tel_len = 0
                pq_tlens.append((pq, my_tel_len))
                #
                if MAKE_PLOTS:
                    density_data = [td_p_e0, td_p_e1, td_q_e0, td_q_e1, p_vs_q_power, len(subtel_seq), TEL_WINDOW_SIZE]
                    plot_title = f'{chr_name}{pq} {my_tel_len} bp'
                    plot_fn    = f'{OUT_DIR}{SAMP_NAME}telsignal_{chr_name}{pq}.png'
                    plot_tel_signal(density_data, plot_title, plot_fn, tl_vals=tl_vals)
            #
            for pqt in pq_tlens:
                if pqt[1] < MIN_END_TL:
                    print(f' -- {SAMP_NAME}{chr_name}{pqt[0]} masked tel = {pqt[1]}bp [skipped]')
                else:
                    print(f' -- {SAMP_NAME}{chr_name}{pqt[0]} masked tel = {pqt[1]}bp')
                    if pqt[0] == 'p':
                        out_sequences.append((f'{SAMP_NAME}{CHR_REMAP_DICT[chr_name]}p', 'N'*pqt[1] + p_arm[pqt[1]:]))
                        tel_sequences.append((f'{SAMP_NAME}tel-{CHR_REMAP_DICT[chr_name]}p', p_arm[:pqt[1]]))
                    elif pqt[0] == 'q':
                        out_sequences.append((f'{SAMP_NAME}{CHR_REMAP_DICT[chr_name]}q', q_arm[:len(q_arm)-pqt[1]] + 'N'*pqt[1]))
                        tel_sequences.append((f'{SAMP_NAME}tel-{CHR_REMAP_DICT[chr_name]}q', q_arm[len(q_arm)-pqt[1]:]))
    tel_sequences.append(('tel-TAACCC', 'TAACCC'*10000))

    with open(OUT_REF, 'w') as f:
        for n in out_sequences:
            f.write(f'>{n[0]}\n')
            for j in range(0, len(n[1]), FASTA_WIDTH):
                f.write(f'{n[1][j:j+FASTA_WIDTH]}\n')
        if APPEND_TEL:
            for n in tel_sequences:
                f.write(f'>{n[0]}\n')
                for j in range(0, len(n[1]), FASTA_WIDTH):
                    f.write(f'{n[1][j:j+FASTA_WIDTH]}\n')
    pysam.faidx(OUT_REF)


if __name__ == '__main__':
    main()
