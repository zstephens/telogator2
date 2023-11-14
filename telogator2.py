#!/usr/bin/env python
import argparse
import copy
import gzip
import numpy as np
import pathlib
import pickle
import pysam
import subprocess
import sys
import time

from source.tg_kmer   import get_nonoverlapping_kmer_hits, get_telomere_base_count, read_kmer_tsv
from source.tg_muscle import check_muscle_version
from source.tg_plot   import plot_kmer_hits, tel_len_violin_plot
from source.tg_reader import quick_grab_all_reads_nodup, TG_Reader
from source.tg_tel    import get_allele_tsv_dat, get_terminating_tl
from source.tg_tvr    import cluster_consensus_tvrs, cluster_tvrs, convert_colorvec_to_kmerhits, make_tvr_plots
from source.tg_util   import exists_and_is_nonzero, get_downsample_inds, get_file_type, LEXICO_2_IND, makedir, parse_read, RC, rm, strip_paths_from_string

TEL_WINDOW_SIZE = 100
P_VS_Q_AMP_THRESH = 0.5
#
DUMMY_TEL_MAPQ = 60
# how much subtel should we use for de-clustering alleles? [min_size, max_size]
SUBTEL_CLUSTER_SIZE = [1500, 3000]
#
MIN_TEL_UNANCHORED = 500
NONTEL_EDGE_UNANCHORED = 200
# if >30% this fraction of reads have no terminating tel, skip the cluster
TERM_TEL_ZERO_FRAC = 0.30
# if the nontel content of reads exceeds 190 for >49% of reads, skip this cluster
NONTEL_END_FILT_PARAMS = (190, 0.49)
#
MAX_QUAL_SCORE = 60
ANCHORING_ASSIGNMENT_FRAC = 0.20


def main(raw_args=None):
    parser = argparse.ArgumentParser(description='Telogator2', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    parser.add_argument('-i', type=str, required=True,  metavar='input.fa',     nargs='*',      help="* Input reads (fa / fa.gz / fq / fq.gz / bam)")
    parser.add_argument('-o', type=str, required=True,  metavar='output/',                      help="* Path to output directory")
    parser.add_argument('-k', type=str, required=False, metavar='kmers.tsv',    default='',     help="Telomere kmers file")
    parser.add_argument('-t', type=str, required=False, metavar='telogator.fa', default='',     help="Telogator reference fasta")
    parser.add_argument('-r', type=str, required=False, metavar='hifi',         default='hifi', help="Read type: hifi / clr / ont")
    parser.add_argument('-l', type=int, required=False, metavar='4000',         default=4000,   help="Minimum read length")
    parser.add_argument('-s', type=int, required=False, metavar='1000',         default=1000,   help="Minimum subtelomere anchor size")
    parser.add_argument('-c', type=int, required=False, metavar='5',            default=5,      help="Minimum hits to canonical kmer")
    parser.add_argument('-n', type=int, required=False, metavar='3',            default=3,      help="Minimum number of reads per cluster")
    parser.add_argument('-m', type=str, required=False, metavar='p75',          default='p75',  help="Method for computing chr TL: mean / median p75 / max")
    parser.add_argument('-d', type=int, required=False, metavar='nreads',       default=-1,     help="Downsample to this many telomere reads")
    parser.add_argument('-p', type=int, required=False, metavar='4',            default=4,      help="Number of processes to use")
    #
    parser.add_argument('-ti', type=float, required=False, metavar='0.180', default=0.180, help="Treecut value: initial TVR clustering")
    parser.add_argument('-tp', type=float, required=False, metavar='0.120', default=0.120, help="Treecut value: merging TVRs with similar prefixes")
    parser.add_argument('-tt', type=float, required=False, metavar='0.250', default=0.250, help="Treecut value: cluster refinement [TVR]")
    parser.add_argument('-ts', type=float, required=False, metavar='0.250', default=0.250, help="Treecut value: cluster refinement [SUBTEL]")
    #
    parser.add_argument('--plot-filt-tvr',  required=False, action='store_true', default=False, help="Plot denoised TVR instead of raw signal")
    parser.add_argument('--debug-print',    required=False, action='store_true', default=False, help="[DEBUG] Print extra info for each read as its processed")
    parser.add_argument('--debug-npy',      required=False, action='store_true', default=False, help="[DEBUG] Save .npy files and use existing .npy files")
    parser.add_argument('--debug-noplot',   required=False, action='store_true', default=False, help="[DEBUG] Do not regenerate plots that already exist")
    parser.add_argument('--debug-nosubtel', required=False, action='store_true', default=False, help="[DEBUG] Skip cluster refinement step that uses subtels")
    parser.add_argument('--debug-noanchor', required=False, action='store_true', default=False, help="[DEBUG] Do not align reads or do any anchoring")
    parser.add_argument('--fast-aln',       required=False, action='store_true', default=False, help="Use faster but less accurate pairwise alignment")
    parser.add_argument('--fast-filt',      required=False, action='store_true', default=False, help="Remove interstitial telomere reads earlier")
    #
    parser.add_argument('-afa-x', type=int, required=False, metavar='15000', default=15000, help="[all_final_alleles.png] X axis max")
    parser.add_argument('-afa-t', type=int, required=False, metavar='1000',  default=1000,  help="[all_final_alleles.png] X axis tick steps")
    parser.add_argument('-afa-a', type=int, required=False, metavar='1000',  default=1000,  help="[all_final_alleles.png] Min ATL to include allele")
    #
    parser.add_argument('--muscle',    type=str, required=False, metavar='muscle',    default='', help="/path/to/muscle")
    parser.add_argument('--minimap2',  type=str, required=False, metavar='minimap2',  default='', help="/path/to/minimap2")
    parser.add_argument('--winnowmap', type=str, required=False, metavar='winnowmap', default='', help="/path/to/winnowmap")
    parser.add_argument('--pbmm2',     type=str, required=False, metavar='pbmm2',     default='', help="/path/to/pbmm2")
    #
    args = parser.parse_args()
    #

    INPUT_ALN     = args.i
    OUT_DIR       = args.o
    KMER_FILE     = args.k
    TELOGATOR_REF = args.t
    #
    READ_TYPE           = args.r
    MINIMUM_READ_LEN    = args.l
    MINIMUM_CANON_HITS  = args.c
    MIN_SUBTEL_BUFF     = args.s
    MIN_READS_PER_PHASE = args.n
    ALLELE_TL_METHOD    = args.m
    DOWNSAMPLE_READS    = args.d
    NUM_PROCESSES       = args.p
    #
    TREECUT_INITIAL       = args.ti
    TREECUT_PREFIXMERGE   = args.tp
    TREECUT_REFINE_SUBTEL = args.ts
    TREECUT_REFINE_TVR    = args.tt
    #
    FINAL_PLOTTING_XMAX  = args.afa_x
    FINAL_PLOTTING_XTICK = args.afa_t
    MIN_ATL_FOR_FINAL_PLOTTING = args.afa_a

    # executables
    #
    MUSCLE_EXE    = args.muscle
    MINIMAP2_EXE  = args.minimap2
    WINNOWMAP_EXE = args.winnowmap
    PBMM2_EXE     = args.pbmm2

    # debug params
    #
    DONT_OVERWRITE_PLOTS = args.debug_noplot     # (True = don't replot figures if they already exist)
    ALWAYS_REPROCESS     = not(args.debug_npy)   # (True = don't write out .npy matrices, always recompute)
    SKIP_SUBTEL_REFINE   = args.debug_nosubtel
    SKIP_ANCHORING       = args.debug_noanchor
    PLOT_FILT_CVECS      = args.plot_filt_tvr
    PRINT_DEBUG          = args.debug_print
    FAST_ALIGNMENT       = args.fast_aln
    FAST_FILTERING       = args.fast_filt

    # check input
    #
    any_pickle = False
    for ifn in INPUT_ALN:
        if exists_and_is_nonzero(ifn) is False:
            print('Error: input not found:')
            print(ifn)
            exit(1)
        input_type = get_file_type(ifn)[0]
        if input_type not in ['fasta', 'fastq', 'bam', 'pickle']:
            print('Error: input must be fasta, fastq, bam, or pickle')
            exit(1)
        if input_type == 'pickle':
            any_pickle = True
    if len(INPUT_ALN) >= 2 and any_pickle:
        print('Error: if multiple inputs are specified, none can be pickle')
        exit(1)

    # prep output directory
    #
    if OUT_DIR[-1] != '/':
        OUT_DIR += '/'
    OUT_CLUST_DIR = OUT_DIR + 'clust_dat/'
    OUT_CDIR_INIT = OUT_CLUST_DIR + '00_initial/'
    OUT_CDIR_TVR  = OUT_CLUST_DIR + '01_tvr/'
    OUT_CDIR_SUB  = OUT_CLUST_DIR + '02_subtel/'
    OUT_CDIR_FIN  = OUT_CLUST_DIR + '03_final/'
    makedir(OUT_DIR)
    makedir(OUT_CLUST_DIR)
    for d in [OUT_CDIR_INIT, OUT_CDIR_TVR, OUT_CDIR_SUB, OUT_CDIR_FIN]:
        makedir(d)
        makedir(d + 'dendro/')
        makedir(d + 'fa/')
        makedir(d + 'npy/')
        makedir(d + 'results/')

    TELOMERE_READS = OUT_DIR + 'tel_reads.fa.gz'
    OUT_ALLELE_TL = OUT_DIR + 'tlens_by_allele.tsv'
    OUT_PICKLE_UNANCHORED = OUT_DIR + 'unanchored-dat.p'
    OUT_UNANCHORED_SUBTELS = OUT_DIR + 'unanchored_subtels.fa.gz'
    ALIGNED_SUBTELS = OUT_DIR + 'subtel_aln'
    VIOLIN_ATL = OUT_DIR + 'violin_atl.png'
    FINAL_TVRS = OUT_DIR + 'all_final_alleles.png'

    RAND_SHUFFLE = 3
    if FAST_ALIGNMENT:
        RAND_SHUFFLE = 1

    #
    # check exes
    #
    if len(MUSCLE_EXE) == 0:
        print('Error: --muscle not specified')
        exit(1)
    check_muscle_version(MUSCLE_EXE)
    #
    if len([n for n in [MINIMAP2_EXE, WINNOWMAP_EXE, PBMM2_EXE] if len(n)]) == 0:
        print('Error: no aligner was specified')
        exit(1)
    elif len([n for n in [MINIMAP2_EXE, WINNOWMAP_EXE, PBMM2_EXE] if len(n)]) > 1:
        print('Error: multiple aligners specified, use only one of the following:')
        print(' --minimap2')
        print(' --winnowmap')
        print(' --pbmm2')
        exit(1)
    if len(MINIMAP2_EXE):
        ALIGNER_EXE = MINIMAP2_EXE
        WHICH_ALIGNER = 'minimap2'
    if len(WINNOWMAP_EXE):
        ALIGNER_EXE = WINNOWMAP_EXE
        WHICH_ALIGNER = 'winnowmap'
    if len(PBMM2_EXE):
        ALIGNER_EXE = PBMM2_EXE
        WHICH_ALIGNER = 'pbmm2'

    #
    # reference seq
    #
    WINNOWMAP_K15 = None
    if TELOGATOR_REF != '':
        fn_suffix = TELOGATOR_REF.split('/')[-1]
        print('using user-specified subtel reference:', fn_suffix)
    else:
        print('using default telogator subtel reference')
        sim_path  = pathlib.Path(__file__).resolve().parent
        TELOGATOR_REF = str(sim_path) + '/resources/telogator-ref.fa.gz'
        WINNOWMAP_K15 = str(sim_path) + '/resources/telogator-ref-k15.txt'

    #
    # parse telomere kmer data
    #
    if KMER_FILE != '':
        fn_suffix = KMER_FILE.split('/')[-1]
        print('using user-specified kmer list:', fn_suffix)
    else:
        print('using default telomere kmers')
        sim_path  = pathlib.Path(__file__).resolve().parent
        KMER_FILE = str(sim_path) + '/resources/kmers.tsv'
    (KMER_METADATA, KMER_ISSUBSTRING, CANONICAL_STRINGS) = read_kmer_tsv(KMER_FILE, READ_TYPE)
    [KMER_LIST, KMER_COLORS, KMER_LETTER, KMER_FLAGS] = KMER_METADATA
    KMER_LIST_REV         = [RC(n) for n in KMER_LIST]
    CANONICAL_STRINGS_REV = [RC(n) for n in CANONICAL_STRINGS]
    KMER_INITIAL = CANONICAL_STRINGS[0] + CANONICAL_STRINGS[0]
    KMER_INITIAL_RC = CANONICAL_STRINGS_REV[0] + CANONICAL_STRINGS_REV[0]

    #
    # param lists that are needed in both anchored and unanchored mode
    #
    gatd_params = [ALLELE_TL_METHOD, MIN_READS_PER_PHASE]
    mtp_params  = [KMER_METADATA, KMER_COLORS, MIN_READS_PER_PHASE, PLOT_FILT_CVECS, DUMMY_TEL_MAPQ, DONT_OVERWRITE_PLOTS]
    blank_chr   = 'chrBq'
    fake_chr    = 'chrUq'
    fake_pos    = 0

    #
    # write cmd out, for debugging purposes
    #
    with open(OUT_DIR + 'cmd.txt', 'w') as f:
        stripped_strings = [strip_paths_from_string(n) for n in sys.argv]
        f.write(' '.join(stripped_strings) + '\n')

    #
    # begin unanchored telogator2
    #
    ALLELE_TEL_DAT = []
    #
    print('beginning telogator2...')
    #
    if any_pickle:
        print('getting unanchored read dat from pickle file...')
        with open(INPUT_ALN[0], 'rb') as f:
            my_pickle = pickle.load(f)
        kmer_hit_dat = my_pickle['kmer-hit-dat']
        all_tvrtel_seq = my_pickle['all-tvrtel-seq']
        all_subtel_seq = my_pickle['all-subtel-seq']
        all_terminating_tl = my_pickle['all-terminating-tl']
        all_nontel_end = my_pickle['all-nontel-end']
        print(' -', len(kmer_hit_dat), 'reads')
    #
    else:
        all_readcount = 0
        tel_readcount = 0
        sys.stdout.write(f'getting reads with at least {MINIMUM_CANON_HITS} matches to {KMER_INITIAL}...')
        sys.stdout.flush()
        tt = time.perf_counter()
        total_bp_all = 0
        total_bp_tel = 0
        with gzip.open(TELOMERE_READS, 'wt') as f:
            for ifn in INPUT_ALN:
                my_reader = TG_Reader(ifn, verbose=False)
                while True:
                    (my_name, my_rdat, my_qdat) = my_reader.get_next_read()
                    if not my_name:
                        break
                    count_fwd = my_rdat.count(KMER_INITIAL)
                    count_rev = my_rdat.count(KMER_INITIAL_RC)
                    all_readcount += 1
                    total_bp_all += len(my_rdat)
                    if count_fwd >= MINIMUM_CANON_HITS or count_rev >= MINIMUM_CANON_HITS:
                        f.write(f'>{my_name}\n{my_rdat}\n')
                        tel_readcount += 1
                        total_bp_tel += len(my_rdat)
                my_reader.close()
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()
        print(f' - {all_readcount} --> {tel_readcount} reads')
        print(f' - ({total_bp_all} bp) --> ({total_bp_tel} bp)')
        if tel_readcount <= 0:
            print('Error: No telomere reads found, stopping here...')
            exit(1)
        #
        sys.stdout.write(f'filtering by read length (>{MINIMUM_READ_LEN}bp)...')
        sys.stdout.flush()
        tt = time.perf_counter()
        (all_read_dat, readcount_len_filtered) = quick_grab_all_reads_nodup(TELOMERE_READS, min_len=MINIMUM_READ_LEN)
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()
        print(f' - {len(all_read_dat)+readcount_len_filtered} --> {len(all_read_dat)} reads')
        if len(all_read_dat) <= 0:
            print('Error: No telomere reads remaining, stopping here...')
            exit(1)
        #
        sys.stdout.write('getting telomere repeat composition...')
        sys.stdout.flush()
        tt = time.perf_counter()
        kmer_hit_dat = []
        all_tvrtel_seq = []
        all_subtel_seq = []
        all_terminating_tl = []
        all_nontel_end = []
        gtt_params = [KMER_LIST, KMER_LIST_REV, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH]
        num_starting_reads = len(all_read_dat)
        for (my_rnm, my_rdat, my_qdat) in all_read_dat:
            tel_bc_fwd = get_telomere_base_count(my_rdat, CANONICAL_STRINGS, mode=READ_TYPE)
            tel_bc_rev = get_telomere_base_count(my_rdat, CANONICAL_STRINGS_REV, mode=READ_TYPE)
            # put everything into q orientation
            if tel_bc_fwd > tel_bc_rev:
                my_rdat = RC(my_rdat)
                if my_qdat is not None:
                    my_qdat = my_qdat[::-1]
            # make sure read actually ends in telomere (remove interstitial telomere regions now, if desired)
            # - removing interstitial tel reads now is less accurate than keeping them and removing them after clustering
            (my_terminating_tel, my_nontel_end) = get_terminating_tl(my_rdat, 'q', gtt_params)
            # some paranoid bounds checking
            my_terminating_tel = min(my_terminating_tel, len(my_rdat))
            my_nontel_end = min(my_nontel_end, len(my_rdat))
            #
            if FAST_FILTERING and my_terminating_tel < MIN_TEL_UNANCHORED:
                continue
            if FAST_FILTERING and my_nontel_end > NONTEL_EDGE_UNANCHORED:
                continue
            # too little subtel sequence?
            if MIN_SUBTEL_BUFF > 0 and len(my_rdat) < my_terminating_tel + MIN_SUBTEL_BUFF:
                continue
            my_subtel_end = max(len(my_rdat)-my_terminating_tel-MIN_SUBTEL_BUFF, 0)
            my_teltvr_seq = my_rdat[my_subtel_end:]
            # if there's no terminating tel at all (and no fast-filt), then lets pretend entire read is tvr+tel then
            # - so that we can use this read for removing interstitial tel regions later
            # - though we're probably just as well off removing them entirely at this point, I don't know for sure.
            if len(my_teltvr_seq) == 0:
                my_teltvr_seq = my_rdat
                out_tvrtel_seq = my_rdat
                out_subtel_seq = ''
            else:
                out_tvrtel_seq = my_rdat[len(my_rdat)-my_terminating_tel:]
                out_subtel_seq = my_rdat[:len(my_rdat)-my_terminating_tel]
            #
            kmer_hit_dat.append([get_nonoverlapping_kmer_hits(my_teltvr_seq, KMER_LIST_REV, KMER_ISSUBSTRING),
                                 len(my_teltvr_seq),   # atb, lets pretend entire read is tel
                                 0,                    # my_dbta
                                 'q',                  # my_type
                                 my_rnm.split(' ')[0], # my_rnm
                                 DUMMY_TEL_MAPQ,       # my_mapq
                                 None])                # out_fasta_dat
            all_tvrtel_seq.append(out_tvrtel_seq)
            all_subtel_seq.append(out_subtel_seq)
            all_terminating_tl.append(my_terminating_tel)
            all_nontel_end.append(my_nontel_end)
        fast_filt_str = ''
        if FAST_FILTERING:
            fast_filt_str = ' [fast-filt applied]'
        num_ending_reads = len(kmer_hit_dat)
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads' + fast_filt_str + '\n')
        sys.stdout.flush()
        if num_ending_reads <= 0:
            print('Error: No telomere reads remaining, stopping here...')
            exit(1)
        #
        if DOWNSAMPLE_READS > 0 and len(kmer_hit_dat) > DOWNSAMPLE_READS:
            sys.stdout.write('downsampling reads...')
            sys.stdout.flush()
            tt = time.perf_counter()
            num_starting_reads = len(kmer_hit_dat)
            del_keys = get_downsample_inds(len(kmer_hit_dat), len(kmer_hit_dat) - DOWNSAMPLE_READS)
            for di in del_keys:
                del kmer_hit_dat[di]
                del all_tvrtel_seq[di]
                del all_subtel_seq[di]
                del all_terminating_tl[di]
                del all_nontel_end[di]
            num_ending_reads = len(kmer_hit_dat)
            sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
            sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads\n')
            sys.stdout.flush()
        #
        f = open(OUT_PICKLE_UNANCHORED, 'wb')
        pickle.dump({'kmer-hit-dat':kmer_hit_dat,
                     'all-tvrtel-seq':all_tvrtel_seq,
                     'all-subtel-seq':all_subtel_seq,
                     'all-terminating-tl':all_terminating_tl,
                     'all-nontel-end':all_nontel_end}, f)
        f.close()
        #
    my_rlens = [len(all_subtel_seq[n]) + len(all_tvrtel_seq[n]) for n in range(len(all_subtel_seq))]
    my_rnames = [n[4] for n in kmer_hit_dat]
    #
    # [1] INITIAL CLUSTERING
    #
    sys.stdout.write('initial clustering of all reads...')
    sys.stdout.flush()
    tt = time.perf_counter()
    init_dendrogram_fn  = OUT_CDIR_INIT + 'dendro/' + 'dendrogram.png'
    init_dend_prefix_fn = OUT_CDIR_INIT + 'dendro/' + 'dendrogram-prefixmerge.png'
    init_dist_matrix_fn = OUT_CDIR_INIT + 'npy/'    + 'dist-matrix.npy'
    init_dist_prefix_fn = OUT_CDIR_INIT + 'npy/'    + 'dist-matrix-prefixmerge.npy'
    init_consensus_fn   = OUT_CDIR_INIT + 'fa/'     + 'consensus.fa'
    #
    if ALWAYS_REPROCESS:
        init_dist_matrix_fn = None
        init_dist_prefix_fn = None
        init_consensus_fn = None
    #
    read_clust_dat = cluster_tvrs(kmer_hit_dat, KMER_METADATA, fake_chr, fake_pos, TREECUT_INITIAL, TREECUT_PREFIXMERGE,
                                  aln_mode='ds',
                                  alignment_processes=NUM_PROCESSES,
                                  rand_shuffle_count=RAND_SHUFFLE,
                                  dist_in=init_dist_matrix_fn,
                                  dist_in_prefix=init_dist_prefix_fn,
                                  fig_name=init_dendrogram_fn,
                                  fig_prefix_name=init_dend_prefix_fn,
                                  save_msa=init_consensus_fn,
                                  muscle_exe=MUSCLE_EXE,
                                  PRINT_DEBUG=PRINT_DEBUG)
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()
    n_clusters = len([n for n in read_clust_dat[0] if len(n) >= MIN_READS_PER_PHASE])
    print(f' - {n_clusters} clusters formed (>= {MIN_READS_PER_PHASE} supporting reads)')
    #
    if True:
        sys.stdout.write('plotting initial clusters...')
        sys.stdout.flush()
        tt = time.perf_counter()
        for clust_i in range(len(read_clust_dat[0])):
            current_clust_inds = read_clust_dat[0][clust_i]
            if len(current_clust_inds) < MIN_READS_PER_PHASE:
                continue
            khd_subset = [copy.deepcopy(kmer_hit_dat[n]) for n in current_clust_inds]
            clustdat_to_plot = [[list(range(len(current_clust_inds)))],
                                [read_clust_dat[1][clust_i]],
                                [read_clust_dat[2][clust_i]],
                                [read_clust_dat[3][n] for n in current_clust_inds],
                                [read_clust_dat[4][clust_i]],
                                [read_clust_dat[5][clust_i]],
                                [read_clust_dat[6][n] for n in current_clust_inds],
                                [read_clust_dat[7][clust_i]]]
            plot_fn_reads = OUT_CDIR_INIT + 'results/' + 'reads_' + str(clust_i).zfill(3) + '.png'
            if DONT_OVERWRITE_PLOTS and exists_and_is_nonzero(plot_fn_reads):
                pass
            else:
                make_tvr_plots(khd_subset, clustdat_to_plot, fake_chr, fake_pos, plot_fn_reads, None, mtp_params)
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()
    #
    # [2] TVR CLUSTERING
    #
    sys.stdout.write('refining clusters [TVR]...')
    sys.stdout.flush()
    tt = time.perf_counter()
    clust_num = 0
    blank_inds = []
    fail_clusters = []
    fail_blank = []
    clusters_with_tvrs = []
    for clust_i in range(len(read_clust_dat[0])):
        current_clust_inds = read_clust_dat[0][clust_i]
        # cluster filter: not enough reads
        if len(current_clust_inds) < MIN_READS_PER_PHASE:
            continue
        # cluster filter: blank tvr, add to blank_inds and move on
        # - also copying the terminating-tel filters from below to avoid interstitial tel repeats
        if read_clust_dat[7][clust_i] <= 0:
            clust_failed = False
            term_tel = [all_terminating_tl[n] for n in current_clust_inds]
            term_zero_frac = len([n for n in term_tel if n <= 0.0]) / len(term_tel)
            if term_zero_frac > TERM_TEL_ZERO_FRAC:
                clust_failed = True
            nt_end = [all_nontel_end[n] for n in current_clust_inds]
            nontel_long_frac = len([n for n in nt_end if n > NONTEL_END_FILT_PARAMS[0]]) / len(nt_end)
            if nontel_long_frac > NONTEL_END_FILT_PARAMS[1]:
                clust_failed = True
            if clust_failed:
                fail_blank.append((term_zero_frac, nontel_long_frac, [n for n in current_clust_inds]))
                continue
            blank_inds.extend(current_clust_inds)
            continue
        khd_subset = [copy.deepcopy(kmer_hit_dat[n]) for n in current_clust_inds]
        zfcn = str(clust_num).zfill(3)
        telcompplot_fn = OUT_CDIR_TVR + 'results/' + 'reads_'                   + zfcn + '.png'
        telcompcons_fn = OUT_CDIR_TVR + 'results/' + 'consensus_'               + zfcn + '.png'
        dendrogram_fn  = OUT_CDIR_TVR + 'dendro/'  + 'dendrogram_'              + zfcn + '.png'
        dend_prefix_fn = OUT_CDIR_TVR + 'dendro/'  + 'dendrogram-prefixmerge_'  + zfcn + '.png'
        dist_matrix_fn = OUT_CDIR_TVR + 'npy/'     + 'dist-matrix_'             + zfcn + '.npy'
        dist_prefix_fn = OUT_CDIR_TVR + 'npy/'     + 'dist-matrix-prefixmerge_' + zfcn + '.npy'
        consensus_fn   = OUT_CDIR_TVR + 'fa/'      + 'consensus_'               + zfcn + '.fa'
        #
        if ALWAYS_REPROCESS:
            dist_matrix_fn = None
            dist_prefix_fn = None
            consensus_fn = None
        #
        subset_clustdat = cluster_tvrs(khd_subset, KMER_METADATA, fake_chr, fake_pos, TREECUT_REFINE_TVR, TREECUT_PREFIXMERGE,
                                       aln_mode='ds',
                                       alignment_processes=NUM_PROCESSES,
                                       rand_shuffle_count=RAND_SHUFFLE,
                                       dist_in=dist_matrix_fn,
                                       dist_in_prefix=dist_prefix_fn,
                                       fig_name=dendrogram_fn,
                                       fig_prefix_name=dend_prefix_fn,
                                       save_msa=consensus_fn,
                                       muscle_exe=MUSCLE_EXE,
                                       PRINT_DEBUG=PRINT_DEBUG)
        #
        make_tvr_plots(khd_subset, subset_clustdat, fake_chr, fake_pos, telcompplot_fn, telcompcons_fn, mtp_params)
        #
        for sci,subclust_inds in enumerate(subset_clustdat[0]):
            subclust_read_inds = [current_clust_inds[n] for n in subclust_inds]
            if len(subclust_read_inds) < MIN_READS_PER_PHASE:
                continue
            #
            # filters to try and remove interstitial telomere repeats
            #
            clust_failed = False
            term_tel = [all_terminating_tl[n] for n in subclust_read_inds]
            term_zero_frac = len([n for n in term_tel if n <= 0.0]) / len(term_tel)
            if term_zero_frac > TERM_TEL_ZERO_FRAC:
                clust_failed = True
            nt_end = [all_nontel_end[n] for n in subclust_read_inds]
            nontel_long_frac = len([n for n in nt_end if n > NONTEL_END_FILT_PARAMS[0]]) / len(nt_end)
            if nontel_long_frac > NONTEL_END_FILT_PARAMS[1]:
                clust_failed = True
            if clust_failed:
                fail_clusters.append((term_zero_frac, nontel_long_frac, [n for n in subclust_read_inds]))
                continue
            #
            # are any of the subclusters blank? --> add to blank inds
            #
            subclust_tvr_len = len(subset_clustdat[4][sci])
            if subclust_tvr_len <= 0:
                blank_inds.extend(subclust_read_inds)
                continue
            #
            # pass all checks? --> add to list for subsequent subtel clustering
            #
            clusters_with_tvrs.append([n for n in subclust_read_inds])
        #
        clust_num += 1
    have_blank = False
    if len(blank_inds) >= MIN_READS_PER_PHASE:
        have_blank = True
        clusters_with_tvrs.append([n for n in blank_inds])
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()
    print(f' - {len(clusters_with_tvrs)} clusters with tvrs (does not yet include blanks)')
    print(f' - {len(fail_clusters) + len(fail_blank)} clusters removed for not termintating in tel)')
    #
    # [3] SUBTEL CLUSTERING
    #
    if SKIP_SUBTEL_REFINE:
        print('skipping subtel cluster refinement...')
        final_clustered_read_inds = []
        for sci,subclust_read_inds in enumerate(clusters_with_tvrs):
            if sci == len(clusters_with_tvrs) - 1 and have_blank:
                my_chr = blank_chr
            else:
                my_chr = fake_chr
            final_clustered_read_inds.append((my_chr, [n for n in subclust_read_inds]))
    else:
        sys.stdout.write('refining clusters [SUBTEL]...')
        sys.stdout.flush()
        tt = time.perf_counter()
        final_clustered_read_inds = []
        for sci,subclust_read_inds in enumerate(clusters_with_tvrs):
            if sci == len(clusters_with_tvrs) - 1 and have_blank:
                my_chr = blank_chr
            else:
                my_chr = fake_chr
            subtel_sizes = [len(all_subtel_seq[n]) for n in subclust_read_inds]
            subtel_size = max(min(min(subtel_sizes), SUBTEL_CLUSTER_SIZE[1]), SUBTEL_CLUSTER_SIZE[0])
            my_subtels = [all_subtel_seq[n][-subtel_size:] for n in subclust_read_inds]
            #
            zfcn = str(sci).zfill(3)
            subtel_dendro_fn = OUT_CDIR_SUB + 'dendro/' + 'dendrogram_'  + zfcn + '.png'
            subtel_dist_fn   = OUT_CDIR_SUB + 'npy/'    + 'dist-matrix_' + zfcn + '.npy'
            my_dendro_title  = zfcn
            subtel_labels    = None
            #
            if ALWAYS_REPROCESS:
                subtel_dist_fn = None
            #
            subtel_clustdat = cluster_consensus_tvrs(my_subtels, KMER_METADATA, TREECUT_REFINE_SUBTEL,
                                                     alignment_processes=NUM_PROCESSES,
                                                     aln_mode='ms',
                                                     gap_bool=(False,False),
                                                     rand_shuffle_count=RAND_SHUFFLE,
                                                     adjust_lens=False,
                                                     dist_in=subtel_dist_fn,
                                                     dendro_name=subtel_dendro_fn,
                                                     samp_labels=subtel_labels,
                                                     linkage_method='ward',
                                                     normalize_dist_matrix=False,
                                                     job=(1,1),
                                                     dendrogram_title=my_dendro_title,
                                                     dendrogram_height=8)
            #
            for sci,subclust_inds in enumerate(subtel_clustdat):
                subsubclust_read_inds = [subclust_read_inds[n] for n in subclust_inds]
                if len(subsubclust_read_inds) < MIN_READS_PER_PHASE:
                    continue
                final_clustered_read_inds.append((my_chr, [n for n in subsubclust_read_inds]))
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()
        print(f' - {len(final_clustered_read_inds)} clusters')
    #
    #
    #
    sys.stdout.write('final reclustering to get TVRs and TLs for each allele...')
    sys.stdout.flush()
    tt = time.perf_counter()
    clust_num = 0
    subtels_out = []
    allele_outdat = []
    allele_consensus = []
    for final_clustdat in final_clustered_read_inds:
        (my_chr, current_clust_inds) = final_clustdat
        khd_subset = [copy.deepcopy(kmer_hit_dat[n]) for n in current_clust_inds]
        rlens_subset = [my_rlens[n] for n in current_clust_inds]
        zfcn = str(clust_num).zfill(3)
        telcompplot_fn = OUT_CDIR_FIN + 'results/' + 'reads_'                   + zfcn + '.png'
        telcompcons_fn = OUT_CDIR_FIN + 'results/' + 'consensus_'               + zfcn + '.png'
        dendrogram_fn  = OUT_CDIR_FIN + 'dendro/'  + 'dendrogram_'              + zfcn + '.png'
        dend_prefix_fn = OUT_CDIR_FIN + 'dendro/'  + 'dendrogram-prefixmerge_'  + zfcn + '.png'
        dist_matrix_fn = OUT_CDIR_FIN + 'npy/'     + 'dist-matrix_'             + zfcn + '.npy'
        dist_prefix_fn = OUT_CDIR_FIN + 'npy/'     + 'dist-matrix-prefixmerge_' + zfcn + '.npy'
        consensus_fn   = OUT_CDIR_FIN + 'fa/'      + 'consensus_'               + zfcn + '.fa'
        #
        if ALWAYS_REPROCESS:
            dist_matrix_fn = None
            dist_prefix_fn = None
            consensus_fn = None
        #
        solo_clustdat = cluster_tvrs(khd_subset, KMER_METADATA, my_chr, fake_pos, 100., TREECUT_PREFIXMERGE,
                                     aln_mode='ds',
                                     alignment_processes=NUM_PROCESSES,
                                     rand_shuffle_count=1,
                                     dist_in=dist_matrix_fn,
                                     dist_in_prefix=dist_prefix_fn,
                                     fig_name=dendrogram_fn,
                                     fig_prefix_name=dend_prefix_fn,
                                     save_msa=consensus_fn,
                                     muscle_exe=MUSCLE_EXE,
                                     PRINT_DEBUG=PRINT_DEBUG)
        #
        my_subtels = []
        for i,sri in enumerate(solo_clustdat[0][0]):
            read_i = current_clust_inds[sri]
            subtel_used_in_tvr = MIN_SUBTEL_BUFF - solo_clustdat[2][0][i]
            if subtel_used_in_tvr >= 0:
                my_subtels.append(all_subtel_seq[read_i][:len(all_subtel_seq[read_i])-subtel_used_in_tvr])
            else:
                my_subtels.append(all_subtel_seq[read_i])
        for i,subtel_seq in enumerate(my_subtels):
            subtels_out.append((f'cluster-{clust_num}_read-{i}_{my_rnames[read_i]}', subtel_seq))
        #
        make_tvr_plots(khd_subset, solo_clustdat, my_chr, fake_pos, telcompplot_fn, telcompcons_fn, mtp_params)
        #
        allele_outdat.extend(get_allele_tsv_dat(khd_subset, solo_clustdat, my_chr, fake_pos, rlens_subset, gatd_params))
        allele_consensus.append(solo_clustdat[4][0])
        #
        clust_num += 1
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()

    #
    # fill out allele id fields (there are no multimapped alleles at this point)
    # -- print out some TVR stats to console
    #
    for i,atdat in enumerate(allele_outdat):
        ALLELE_TEL_DAT.append([n for n in atdat])
        ALLELE_TEL_DAT[-1][3] = str(i)
    num_unique_alleles = len(set([n[3] for n in ALLELE_TEL_DAT]))
    num_blank_alleles  = len([n[0] for n in ALLELE_TEL_DAT if n[0] == blank_chr])
    print(f' - {num_unique_alleles} unique alleles')
    print(f' - {num_blank_alleles} / {num_unique_alleles} have blank TVRs')

    #
    #
    #
    if SKIP_ANCHORING is False:
        #
        # write out subtels
        #
        with gzip.open(OUT_UNANCHORED_SUBTELS, 'wt') as f:
            for n in subtels_out:
                if len(n[1]) > 0:
                    f.write(f'>{n[0]}\n{n[1]}\n')

        #
        # subtel alignment
        #
        sys.stdout.write('aligning subtels to reference...')
        sys.stdout.flush()
        tt = time.perf_counter()
        cmd = ''
        aln_log = ALIGNED_SUBTELS + '.log'
        aln_sam = ALIGNED_SUBTELS + '.sam'
        aln_bam = ALIGNED_SUBTELS + '.bam'
        #
        if WHICH_ALIGNER == 'minimap2':
            aln_params = ''
            if READ_TYPE == 'ont':
                aln_params = '-ax map-ont -Y'
            elif READ_TYPE == 'hifi':
                aln_params = '-ax map-hifi -Y'
            cmd = ALIGNER_EXE + ' ' + aln_params + ' -o ' + aln_sam + ' ' + TELOGATOR_REF + ' ' + OUT_UNANCHORED_SUBTELS
        #
        elif WHICH_ALIGNER == 'winnowmap':
            aln_params = ''
            if READ_TYPE == 'ont':
                aln_params = '-ax map-ont -Y'
            elif READ_TYPE == 'hifi':
                aln_params = '-ax map-pb -Y'
            cmd = ALIGNER_EXE + ' -W ' + WINNOWMAP_K15 + ' ' + aln_params + ' -o ' + aln_sam + ' ' + TELOGATOR_REF + ' ' + OUT_UNANCHORED_SUBTELS
        #
        elif WHICH_ALIGNER == 'pbmm2':
            cmd = ALIGNER_EXE + ' align ' + TELOGATOR_REF + ' ' + OUT_UNANCHORED_SUBTELS + ' ' + aln_bam + ' --preset HiFi --sort'
        if len(cmd):
            with open(aln_log, 'w') as f:
                try:
                    subprocess.check_output(cmd.split(' '), stderr=f, text=True)
                except subprocess.CalledProcessError as exc:
                    print('Error: alignment command returned an error:', exc.returncode)
                    print(exc.output)
                    exit(1)
        if exists_and_is_nonzero(aln_sam) is False and exists_and_is_nonzero(aln_bam) is False:
            print()
            print('Error: subtel alignment failed. Check the log:')
            print(f'{aln_log}')
            exit(1)
        #
        if WHICH_ALIGNER != 'pbmm2':
            pysam.sort("-o", aln_bam, aln_sam)
            pysam.index(aln_bam)
            rm(aln_sam)
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()
        if exists_and_is_nonzero(aln_bam) is False:
            print('Error: failed to create subtel bam')
            exit(1)

        #
        # get anchors from subtel alignment
        #
        sys.stdout.write('getting anchors from alignment...')
        sys.stdout.flush()
        tt = time.perf_counter()
        ALIGNMENTS_BY_RNAME = {}
        samfile = pysam.AlignmentFile(aln_bam, "rb")
        refseqs = samfile.references
        for aln in samfile.fetch(until_eof=True):
            sam_line = str(aln).split('\t')
            # pysam weirdness:
            my_ref_ind = sam_line[2].replace('#','')
            if my_ref_ind.isdigit():
                sam_line[2] = refseqs[int(my_ref_ind)]
            elif my_ref_ind == '-1':
                sam_line[2] = refseqs[-1]
            else:
                sam_line[2] = my_ref_ind
            [rnm, ref_key, pos, read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat] = parse_read(sam_line)
            if rnm not in ALIGNMENTS_BY_RNAME:
                ALIGNMENTS_BY_RNAME[rnm] = []
            ALIGNMENTS_BY_RNAME[rnm].append([read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat])
        samfile.close()
        #
        top_alns_by_cluster = {}
        for readname in ALIGNMENTS_BY_RNAME:
            sorted_choices = []
            for i,aln in enumerate(ALIGNMENTS_BY_RNAME[readname]):
                dist_to_end = len(aln[7]) - max(aln[0], aln[1])
                read_span = abs(aln[1] - aln[0])
                my_mapq = aln[6]
                sorted_choices.append((read_span, my_mapq, -dist_to_end, i))
            sorted_choices = sorted(sorted_choices, reverse=True)
            top_i = sorted_choices[0][3]
            #print(readname, sorted_choices)
            my_cluster = int(readname.split('_')[0][8:])
            if my_cluster not in top_alns_by_cluster:
                top_alns_by_cluster[my_cluster] = []
            top_alns_by_cluster[my_cluster].append((readname, copy.deepcopy(ALIGNMENTS_BY_RNAME[readname][top_i])))
        #
        for clustnum in sorted(top_alns_by_cluster.keys()):
            #print(clustnum)
            chr_arm_scores = {}
            anchors_by_ref = {}
            for n in top_alns_by_cluster[clustnum]:
                my_ref = n[1][2]
                my_refspan = abs(n[1][3] - n[1][4])
                my_anchor = n[1][4] # second refpos is always anchor coord?
                my_orr = n[1][5]
                my_mapq = n[1][6]
                #print(my_ref, my_refspan, my_anchor, my_orr, my_mapq)
                if my_ref != '*':
                    if my_ref not in chr_arm_scores:
                        chr_arm_scores[my_ref] = 0.
                        anchors_by_ref[my_ref] = []
                    chr_arm_scores[my_ref] += my_refspan * ((my_mapq + 1) / (MAX_QUAL_SCORE + 1))
                    anchors_by_ref[my_ref].append(my_anchor)
            total_score = sum(chr_arm_scores.values())
            my_anchors = []
            if total_score > 0:
                sorted_chr_scores = sorted([(chr_arm_scores[k]/total_score, k) for k in chr_arm_scores], reverse=True)
                #print(sorted_chr_scores)
                my_anchors = [n[1] for n in sorted_chr_scores if n[0] >= ANCHORING_ASSIGNMENT_FRAC]
            else:
                pass # all unmapped
            if len(my_anchors):
                anchor_pos = []
                ref_builds = []
                my_chrs    = []
                chrs_encountered_so_far = {}
                for n in my_anchors:
                    my_chr  = n.split('_')[1]
                    my_samp = n.split('_')[0]
                    my_pos  = str(int(np.median(anchors_by_ref[n])))
                    if my_chr in chrs_encountered_so_far:
                        continue
                    chrs_encountered_so_far[my_chr] = True
                    anchor_pos.append(my_pos)
                    ref_builds.append(my_samp)
                    my_chrs.append(my_chr)
                ALLELE_TEL_DAT[clustnum][0] = ','.join(my_chrs)    # assign chr
                ALLELE_TEL_DAT[clustnum][1] = ','.join(anchor_pos) # assign pos
                ALLELE_TEL_DAT[clustnum][2] = ','.join(ref_builds) # assign ref builds
        #
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()

        #
        # resort out allele dat by chr & pos
        #
        sorted_ad_inds = sorted([(LEXICO_2_IND[n[0].split(',')[0][:-1]], n[0].split(',')[0][-1], int(n[1].split(',')[0]), i) for i,n in enumerate(ALLELE_TEL_DAT)])
        ALLELE_TEL_DAT = [ALLELE_TEL_DAT[n[3]] for n in sorted_ad_inds]

    #
    # write out allele TLs
    #
    print('writing allele TL results to "' + OUT_ALLELE_TL.split('/')[-1] + '"...')
    with open(OUT_ALLELE_TL, 'w') as f:
        f.write('#chr' + '\t' +
                'position' + '\t' +
                'ref_samp' + '\t' +
                'allele_id' + '\t' +
                'TL_' + ALLELE_TL_METHOD + '\t' +
                'read_TLs' + '\t' +
                'read_lengths' + '\t' +
                'read_mapq' + '\t' +
                'tvr_len' + '\t' +
                'tvr_consensus' + '\t' +
                'supporting_reads' + '\n')
        for n in ALLELE_TEL_DAT:
            f.write('\t'.join(n) + '\n')

    #
    # violin plots
    #
    if SKIP_ANCHORING is False:
        atl_by_arm = {}
        for atd in ALLELE_TEL_DAT:
            my_chr = atd[0].split(',')[0]
            if my_chr in [fake_chr, blank_chr]:
                my_chr = 'unanchored'
            if my_chr not in atl_by_arm:
                atl_by_arm[my_chr] = []
            atl_by_arm[my_chr].extend([int(n) for n in atd[5].split(',')])
        vparams = {'norm_by_readcount':True,
                   'include_unanchored':True,
                   'p_ymax':20000,
                   'q_ymax':20000,
                   'y_step':5000,
                   'y_label':'<-- q      ATL      p -->'}
        tel_len_violin_plot(atl_by_arm, VIOLIN_ATL, custom_plot_params=vparams)

    #
    # plot all final tvrs
    #
    tvrs_to_plot = []
    tvr_labels_to_plot = []
    clustdat_to_plot = [[[]], [[]], [[]], []]
    current_i = 0
    for atd in ALLELE_TEL_DAT:
        my_id = int(atd[3])
        my_rep_atl = int(atd[4])
        my_max_atl = max([int(n) for n in atd[5].split(',')])
        my_tvr_len = int(atd[8])
        if my_max_atl < MIN_ATL_FOR_FINAL_PLOTTING:
            continue
        if int(atd[8]) <= 0:
            my_annot = str(atd[3]) + 'b'
        else:
            my_annot = str(atd[3])
        tvrs_to_plot.append(allele_consensus[my_id][:my_tvr_len+my_rep_atl])
        tvr_labels_to_plot.append(f'({my_annot}) {atd[0]}')
        clustdat_to_plot[0][0].append(current_i)
        clustdat_to_plot[1][0].append(60)
        clustdat_to_plot[2][0].append(0)
        clustdat_to_plot[3].append(0)
        current_i += 1
    redrawn_tvrs = convert_colorvec_to_kmerhits(tvrs_to_plot, KMER_METADATA)
    clust_khd = []
    for i in range(len(redrawn_tvrs)):
        clust_khd.append([redrawn_tvrs[i], len(tvrs_to_plot[i]), 0, 'FWD', tvr_labels_to_plot[i], MAX_QUAL_SCORE, None])
    custom_plot_params = {'xlim':[0,FINAL_PLOTTING_XMAX],
                          'xstep':FINAL_PLOTTING_XTICK,
                          'custom_title':'',
                          'number_label_rows':False,
                          'font.size':16,
                          'font.weight':'bold'}
    plot_kmer_hits(clust_khd, KMER_COLORS, '', 0, FINAL_TVRS, clust_dat=clustdat_to_plot, plot_params=custom_plot_params)


if __name__ == '__main__':
    main()
