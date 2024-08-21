#!/usr/bin/env python
import argparse
import copy
import gzip
import numpy as np
import pathlib
import pysam
import subprocess
import sys
import time

from source.tg_align  import get_nucl_consensus, quick_compare_tvrs
from source.tg_kmer   import get_canonical_letter, get_nonoverlapping_kmer_hits, get_telomere_base_count, read_kmer_tsv
from source.tg_plot   import convert_colorvec_to_kmerhits, make_tvr_plots, plot_kmer_hits, readlen_plot, tel_len_violin_plot
from source.tg_reader import quick_grab_all_reads, quick_grab_all_reads_nodup, TG_Reader
from source.tg_tel    import get_allele_tsv_dat, get_terminating_tl, merge_allele_tsv_dat
from source.tg_tvr    import cluster_consensus_tvrs, cluster_tvrs
from source.tg_util   import annotate_interstitial_tel, check_aligner_exe, dir_exists, exists_and_is_nonzero, get_downsample_inds, get_file_type, LEXICO_2_IND, makedir, mv, parse_read, RC, rm, strip_paths_from_string

TEL_WINDOW_SIZE = 100
P_VS_Q_AMP_THRESH = 0.5
DUMMY_TEL_MAPQ = 60
# how much subtel should we use for de-clustering alleles? [min_size, max_size]
SUBTEL_CLUSTER_SIZE = [1500, 3000]
# if >30% this fraction of reads have no terminating tel, skip the cluster
TERM_TEL_ZERO_FRAC = 0.30
# interstitial telomere filtering parameters:
# - if the terminating content of reads is >190bp of non-telomere sequence for >49% of reads, skip this cluster
NONTEL_END_FILT_PARAMS = (190, 0.49)
# treecut values for collapsing multiple (presumably false positive) clusters together
COLLAPSE_TVR_THRESH = 0.050
COLLAPSE_SUB_THRESH = 0.050 # applies to non-normalized distance matrix
#
MAX_QUAL_SCORE = 60
ANCHORING_ASSIGNMENT_FRAC = 0.20


def main(raw_args=None):
    parser = argparse.ArgumentParser(description='Telogator2', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    parser.add_argument('-i', type=str, required=True,  metavar='',             nargs='*',      help="* Input reads (fa / fa.gz / fq / fq.gz / bam)")
    parser.add_argument('-o', type=str, required=True,  metavar='output/',                      help="* Path to output directory")
    parser.add_argument('-k', type=str, required=False, metavar='kmers.tsv',    default='',     help="Telomere kmers file")
    parser.add_argument('-t', type=str, required=False, metavar='telogator.fa', default='',     help="Telogator reference fasta")
    parser.add_argument('-r', type=str, required=False, metavar='hifi',         default='hifi', help="Read type: hifi / ont")
    parser.add_argument('-l', type=int, required=False, metavar='4000',         default=4000,   help="Minimum read length")
    parser.add_argument('-s', type=int, required=False, metavar='1000',         default=1000,   help="Minimum subtelomere anchor size")
    parser.add_argument('-c', type=int, required=False, metavar='8',            default=8,      help="Minimum hits to tandem canonical kmer")
    parser.add_argument('-n', type=int, required=False, metavar='3',            default=3,      help="Minimum number of reads per cluster")
    parser.add_argument('-m', type=str, required=False, metavar='p75',          default='p75',  help="Method for computing chr TL: mean / median / p75 / max")
    parser.add_argument('-d', type=int, required=False, metavar='-1',           default=-1,     help="Downsample to this many telomere reads")
    parser.add_argument('-p', type=int, required=False, metavar='4',            default=4,      help="Number of processes to use")
    #
    parser.add_argument('-ti', type=float, required=False, metavar='0.250', default=0.250, help="[TREECUT] initial TVR clustering")
    parser.add_argument('-tt', type=float, required=False, metavar='0.250', default=0.250, help="[TREECUT] cluster refinement [TVR]")
    parser.add_argument('-ts', type=float, required=False, metavar='0.250', default=0.250, help="[TREECUT] cluster refinement [SUBTEL]")
    #
    parser.add_argument('-afa-x', type=int, required=False, metavar='15000', default=15000, help="[ALL_FINAL_ALLELES.PNG] X axis max")
    parser.add_argument('-afa-t', type=int, required=False, metavar='1000',  default=1000,  help="[ALL_FINAL_ALLELES.PNG] X axis tick steps")
    parser.add_argument('-afa-a', type=int, required=False, metavar='100',   default=100,   help="[ALL_FINAL_ALLELES.PNG] Min ATL to include allele")
    #
    parser.add_argument('-va-y', type=int, required=False, metavar='20000', default=20000, help="[VIOLIN_ATL.PNG] Y axis max")
    parser.add_argument('-va-t', type=int, required=False, metavar='5000',  default=5000,  help="[VIOLIN_ATL.PNG] Y axis tick steps")
    parser.add_argument('-va-p', type=int, required=False, metavar='2',     default=2,     help="[VIOLIN_ATL.PNG] ploidy. i.e. number of alleles per arm")
    #
    parser.add_argument('--plot-filt-tvr',  required=False, action='store_true', default=False, help="[DEBUG] Plot denoised TVR instead of raw signal")
    parser.add_argument('--plot-all-reads', required=False, action='store_true', default=False, help="[DEBUG] Plot all tel reads during initial clustering")
    parser.add_argument('--plot-signals',   required=False, action='store_true', default=False, help="[DEBUG] Plot tel signals for all reads")
    parser.add_argument('--dont-reprocess', required=False, action='store_true', default=False, help="[DEBUG] Use existing intermediary files (for redoing failed runs)")
    parser.add_argument('--debug-noplot',   required=False, action='store_true', default=False, help="[DEBUG] Do not regenerate plots that already exist")
    parser.add_argument('--debug-realign',  required=False, action='store_true', default=False, help="[DEBUG] Do not redo subtel alignment if bam exists")
    parser.add_argument('--debug-nosubtel', required=False, action='store_true', default=False, help="[DEBUG] Skip subtel cluster refinement")
    parser.add_argument('--debug-noanchor', required=False, action='store_true', default=False, help="[DEBUG] Do not align reads or do any anchoring")
    parser.add_argument('--fast-aln',       required=False, action='store_true', default=False, help="[PERFORMANCE] Use faster but less accurate pairwise alignment")
    #
    parser.add_argument('--init-filt',    type=int, required=False, metavar='', nargs=2, default=(-1,-1), help="[PERFORMANCE] Apply term-tel filters earlier")
    parser.add_argument('--collapse-hom', type=int, required=False, metavar='',          default=1000,    help="[PERFORMANCE] Merge alleles mapped within this bp of each other")
    #
    parser.add_argument('--minimap2',  type=str, required=False, metavar='exe',    default='', help="/path/to/minimap2")
    parser.add_argument('--winnowmap', type=str, required=False, metavar='exe',    default='', help="/path/to/winnowmap")
    parser.add_argument('--pbmm2',     type=str, required=False, metavar='exe',    default='', help="/path/to/pbmm2")
    parser.add_argument('--ref',       type=str, required=False, metavar='ref.fa', default='', help="Reference filename (only needed if input is cram)")
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
    CRAM_REF_FILE       = args.ref
    #
    TREECUT_INITIAL       = args.ti
    TREECUT_REFINE_SUBTEL = args.ts
    TREECUT_REFINE_TVR    = args.tt
    #
    VIOLIN_YMAX   = args.va_y
    VIOLIN_TICK   = args.va_t
    VIOLIN_PLOIDY = args.va_p
    #
    FINAL_PLOTTING_XMAX  = args.afa_x
    FINAL_PLOTTING_XTICK = args.afa_t
    MIN_ATL_FOR_FINAL_PLOTTING = args.afa_a

    # executables
    #
    MINIMAP2_EXE  = args.minimap2
    WINNOWMAP_EXE = args.winnowmap
    PBMM2_EXE     = args.pbmm2

    # debug params
    #
    DONT_OVERWRITE_PLOTS = args.debug_noplot          # (True = don't replot figures if they already exist)
    ALWAYS_REPROCESS     = not(args.dont_reprocess)   # (True = don't write out .npz matrices, always recompute)
    NO_SUBTEL_REALIGN    = args.debug_realign
    SKIP_SUBTEL_REFINE   = args.debug_nosubtel
    SKIP_ANCHORING       = args.debug_noanchor
    PLOT_ALL_INITIAL     = args.plot_all_reads
    PLOT_FILT_CVECS      = args.plot_filt_tvr
    PLOT_TEL_SIGNALS     = args.plot_signals
    FAST_ALIGNMENT       = args.fast_aln
    INIT_FILTERING_TUPLE = args.init_filt
    COLLAPSE_HOM_BP      = args.collapse_hom
    #
    PLOT_SEPARATE_INITIAL = False # turn on for debugging, otherwise is mostly a waste of time
    #
    INIT_FILTERING = False
    MIN_TEL_UNANCHORED = 500
    NONTEL_EDGE_UNANCHORED = 200
    if INIT_FILTERING_TUPLE[0] >= 0 and INIT_FILTERING_TUPLE[1] >= 0:
        INIT_FILTERING = True
        MIN_TEL_UNANCHORED = INIT_FILTERING_TUPLE[0]
        NONTEL_EDGE_UNANCHORED = INIT_FILTERING_TUPLE[1]

    # check input
    #
    any_cram = False
    for ifn in INPUT_ALN:
        if exists_and_is_nonzero(ifn) is False:
            print('Error: input not found:')
            print(ifn)
            exit(1)
        input_type = get_file_type(ifn)[0]
        if input_type not in ['fasta', 'fastq', 'bam', 'cram']:
            print('Error: input must be fasta, fastq, bam, or cram')
            exit(1)
        if input_type == 'cram':
            any_cram = True
    if any_cram and CRAM_REF_FILE == '':
        print('Error: cram input requires reference fasta via --ref')
        exit(1)

    # check parameters (TODO: make this more comprehensive)
    #
    if READ_TYPE not in ['hifi', 'ont']:
        print('Error: -r must be hifi or ont')
        exit(1)
    if MIN_READS_PER_PHASE < 1:
        print('Error: -n must be >= 1')
        exit(1)

    # prep output directory
    #
    if OUT_DIR[-1] != '/':
        OUT_DIR += '/'
    OUT_CLUST_DIR = OUT_DIR + 'temp/'
    OUT_CDIR_INIT = OUT_CLUST_DIR + '01_initial/'
    OUT_CDIR_TVR  = OUT_CLUST_DIR + '02_tvr/'
    OUT_CDIR_SUB  = OUT_CLUST_DIR + '03_subtel/'
    OUT_CDIR_COL  = OUT_CLUST_DIR + '04_collapse/'
    OUT_CDIR_FIN  = OUT_CLUST_DIR + '05_final/'
    OUT_QC_DIR    = OUT_DIR + 'qc/'
    makedir(OUT_DIR)
    if dir_exists(OUT_DIR) is False:
        print('Error: could not create output directory')
        exit(1)
    makedir(OUT_CLUST_DIR)
    makedir(OUT_QC_DIR)
    for d in [OUT_CDIR_INIT, OUT_CDIR_TVR, OUT_CDIR_SUB, OUT_CDIR_FIN]:
        makedir(d)
        makedir(d + 'dendro/')
        makedir(d + 'fa/')
        makedir(d + 'npz/')
        makedir(d + 'results/')
    makedir(OUT_CDIR_COL)
    TEL_SIGNAL_DIR = OUT_DIR + 'tel_signal/'
    if PLOT_TEL_SIGNALS:
        makedir(TEL_SIGNAL_DIR)
    UNCOMPRESSED_TELOGATOR_REF = OUT_CLUST_DIR + 'telogator-ref.fa'

    TELOMERE_READS = OUT_CLUST_DIR + 'tel_reads.fa.gz'
    OUT_ALLELE_TL = OUT_DIR + 'tlens_by_allele.tsv'
    ####OUT_PICKLE_UNANCHORED = OUT_CLUST_DIR + 'unanchored-dat.p'
    OUT_UNANCHORED_SUBTELS = OUT_CLUST_DIR + 'unanchored_subtels.fa.gz'
    ALIGNED_SUBTELS = OUT_CLUST_DIR + 'subtel_aln'
    VIOLIN_ATL = OUT_DIR + 'violin_atl.png'
    FINAL_TVRS = OUT_DIR + 'all_final_alleles.png'
    READLEN_NPZ = OUT_CLUST_DIR + 'readlens.npz'
    QC_READLEN = OUT_QC_DIR + 'qc_readlens.png'
    QC_CMD     = OUT_QC_DIR + 'cmd.txt'
    QC_STATS   = OUT_QC_DIR + 'stats.txt'

    RAND_SHUFFLE = 3
    if FAST_ALIGNMENT:
        RAND_SHUFFLE = 1

    #
    # check exes
    #
    if len([n for n in [MINIMAP2_EXE, WINNOWMAP_EXE, PBMM2_EXE] if len(n)]) == 0:
        print('Warning: no aligner was specified, trying "minimap2"...')
        MINIMAP2_EXE = 'minimap2'
    elif len([n for n in [MINIMAP2_EXE, WINNOWMAP_EXE, PBMM2_EXE] if len(n)]) > 1:
        print('Error: multiple aligners specified, use only one of the following:')
        print(' --minimap2')
        print(' --winnowmap')
        print(' --pbmm2')
        exit(1)
    if len(MINIMAP2_EXE):
        ALIGNER_EXE = MINIMAP2_EXE
        WHICH_ALIGNER = 'minimap2'
    elif len(WINNOWMAP_EXE):
        ALIGNER_EXE = WINNOWMAP_EXE
        WHICH_ALIGNER = 'winnowmap'
    elif len(PBMM2_EXE):
        ALIGNER_EXE = PBMM2_EXE
        WHICH_ALIGNER = 'pbmm2'
    if check_aligner_exe(ALIGNER_EXE) is False:
        print(f'Error: {WHICH_ALIGNER} executable not found.')
        exit(1)

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
    canonical_letter = get_canonical_letter(KMER_FILE, READ_TYPE)
    if canonical_letter is None:
        print('Error: kmer list does not contain any kmer marked as canonical')
        exit(1)

    #
    # param lists that are needed in both anchored and unanchored mode
    #
    gatd_params = [ALLELE_TL_METHOD, MIN_READS_PER_PHASE]
    mtp_params  = [KMER_METADATA, KMER_COLORS, MIN_READS_PER_PHASE, PLOT_FILT_CVECS, DUMMY_TEL_MAPQ, DONT_OVERWRITE_PLOTS]
    blank_chr   = 'chrBq'
    unclust_chr = 'chrUq'
    fake_chr    = 'chrFq'
    fake_pos    = 0

    #
    # write cmd out, for debugging purposes
    #
    with open(QC_CMD, 'w') as f:
        stripped_strings = [strip_paths_from_string(n) for n in sys.argv]
        f.write(' '.join(stripped_strings) + '\n')

    #
    # begin!
    #
    ALLELE_TEL_DAT = []
    #
    if ALWAYS_REPROCESS or exists_and_is_nonzero(TELOMERE_READS) is False or exists_and_is_nonzero(READLEN_NPZ) is False:
        sys.stdout.write(f'getting reads with at least {MINIMUM_CANON_HITS} matches to {KMER_INITIAL}...')
        sys.stdout.flush()
        tt = time.perf_counter()
        all_readcount = 0
        tel_readcount = 0
        sup_readcount = 0
        total_bp_all = 0
        total_bp_tel = 0
        readlens_all = []
        readlens_tel = []
        with gzip.open(TELOMERE_READS+'.temp', 'wt') as f:
            for ifn in INPUT_ALN:
                my_reader = TG_Reader(ifn, verbose=False, ref_fasta=CRAM_REF_FILE)
                while True:
                    (my_name, my_rdat, my_qdat, my_issup) = my_reader.get_next_read()
                    if not my_name:
                        break
                    if my_issup:
                        sup_readcount += 1
                        continue
                    count_fwd = my_rdat.count(KMER_INITIAL)
                    count_rev = my_rdat.count(KMER_INITIAL_RC)
                    all_readcount += 1
                    my_rlen = len(my_rdat)
                    total_bp_all += my_rlen
                    readlens_all.append(my_rlen)
                    if count_fwd >= MINIMUM_CANON_HITS or count_rev >= MINIMUM_CANON_HITS:
                        f.write(f'>{my_name}\n{my_rdat}\n')
                        tel_readcount += 1
                        total_bp_tel += my_rlen
                        readlens_tel.append(my_rlen)
                my_reader.close()
        mv(TELOMERE_READS+'.temp', TELOMERE_READS) # temp file as to not immediately overwrite tel_reads.fa.gz if it's the input
        np.savez_compressed(READLEN_NPZ, readlen_all=np.array(readlens_all,dtype='<i8'), readlen_tel=np.array(readlens_tel,dtype='<i8'))
        del readlens_all
        del readlens_tel
        #
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()
        if sup_readcount:
            print(f' - skipped {sup_readcount} supplementary alignments')
        print(f' - {all_readcount} --> {tel_readcount} reads')
        print(f' - ({total_bp_all} bp) --> ({total_bp_tel} bp)')
        if tel_readcount <= 0:
            print('Error: No telomere reads found, stopping here...')
            exit(1)
    else:
        print('found telomere reads from a previous run, using them instead of reprocessing input file')
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
    tel_signal_plot_num = 0
    for (my_rnm, my_rdat, my_qdat) in all_read_dat:
        tel_bc_fwd = get_telomere_base_count(my_rdat, CANONICAL_STRINGS, mode=READ_TYPE)
        tel_bc_rev = get_telomere_base_count(my_rdat, CANONICAL_STRINGS_REV, mode=READ_TYPE)
        # put everything into q orientation
        if tel_bc_fwd > tel_bc_rev:
            my_rdat = RC(my_rdat)
            if my_qdat is not None:
                my_qdat = my_qdat[::-1]
        if PLOT_TEL_SIGNALS:
            tel_signal_plot_dat = (my_rnm, TEL_SIGNAL_DIR + 'signal_' + str(tel_signal_plot_num).zfill(5) + '.png')
            tel_signal_plot_num += 1
        else:
            tel_signal_plot_dat = None
        # make sure read actually ends in telomere (remove interstitial telomere regions now, if desired)
        # - removing interstitial tel reads now is less accurate than keeping them and removing them after clustering
        (my_terminating_tel, my_nontel_end) = get_terminating_tl(my_rdat, 'q', gtt_params, telplot_dat=tel_signal_plot_dat)
        # some paranoid bounds checking
        my_terminating_tel = min(my_terminating_tel, len(my_rdat))
        my_nontel_end = min(my_nontel_end, len(my_rdat))
        #
        if INIT_FILTERING and my_terminating_tel < MIN_TEL_UNANCHORED:
            continue
        if INIT_FILTERING and my_nontel_end > NONTEL_EDGE_UNANCHORED:
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
    if INIT_FILTERING:
        fast_filt_str = ' [init-filt applied]'
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
    my_rlens = [len(all_subtel_seq[n]) + len(all_tvrtel_seq[n]) for n in range(len(all_subtel_seq))]
    my_rnames = [n[4] for n in kmer_hit_dat]

    #
    # [1] INITIAL CLUSTERING
    #

    sys.stdout.write('initial clustering of all reads...')
    sys.stdout.flush()
    tt = time.perf_counter()
    init_dendrogram_fn  = OUT_CDIR_INIT + 'dendro/' + 'dendrogram.png'
    init_dist_matrix_fn = OUT_CDIR_INIT + 'npz/'    + 'dist-matrix.npz'
    if ALWAYS_REPROCESS:
        rm(init_dist_matrix_fn)
    clustering_only = True
    if PLOT_ALL_INITIAL or PLOT_SEPARATE_INITIAL:
        clustering_only = False
    #
    read_clust_dat = cluster_tvrs(kmer_hit_dat, KMER_METADATA, fake_chr, fake_pos, TREECUT_INITIAL,
                                  aln_mode='ds',
                                  alignment_processes=NUM_PROCESSES,
                                  rand_shuffle_count=RAND_SHUFFLE,
                                  dist_in=init_dist_matrix_fn,
                                  fig_name=init_dendrogram_fn,
                                  clustering_only=clustering_only)
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()
    n_clusters = len([n for n in read_clust_dat[0] if len(n) >= MIN_READS_PER_PHASE])
    n_reads = sum([len(n) for n in read_clust_dat[0] if len(n) >= MIN_READS_PER_PHASE])
    print(f' - {n_clusters} clusters ({n_reads} reads)')
    print(f' - (>= {MIN_READS_PER_PHASE} reads per cluster)')
    #
    if PLOT_ALL_INITIAL:
        sys.stdout.write('plotting all tel reads in a single big plot...')
        sys.stdout.flush()
        tt = time.perf_counter()
        init_telcompplot_fn = OUT_CDIR_INIT + 'results/' + 'all_reads.png'
        init_telcompcons_fn = OUT_CDIR_INIT + 'results/' + 'all_consensus.png'
        make_tvr_plots(kmer_hit_dat, read_clust_dat, fake_chr, fake_pos, init_telcompplot_fn, init_telcompcons_fn, mtp_params, dpi=100)
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()
    #
    if PLOT_SEPARATE_INITIAL:
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
    # we're going to use these distances for all the remaining steps instead of recomputing it each time
    init_dist_matrix = np.load(init_dist_matrix_fn)['dist']

    #
    # [1.5] INTERSTITIAL TEL FILTER
    #

    sys.stdout.write('removing interstitial telomeres...')
    sys.stdout.flush()
    tt = time.perf_counter()
    unclustered_inds = []
    pass_clust_inds_list = []
    fail_clust_inds_list = []
    for clust_i in range(len(read_clust_dat[0])):
        current_clust_inds = read_clust_dat[0][clust_i]
        # terminating-tel filters (to prevent interstitial tels from getting into fail/blank sets)
        term_tel = [all_terminating_tl[n] for n in current_clust_inds]
        term_zero_frac = len([n for n in term_tel if n <= 0.0]) / len(term_tel)
        nt_end = [all_nontel_end[n] for n in current_clust_inds]
        nontel_long_frac = len([n for n in nt_end if n > NONTEL_END_FILT_PARAMS[0]]) / len(nt_end)
        passed_termtel = True
        if term_zero_frac > TERM_TEL_ZERO_FRAC or nontel_long_frac > NONTEL_END_FILT_PARAMS[1]:
            passed_termtel = False
        # cluster filter: normal vs. interstitial
        if passed_termtel:
            if len(current_clust_inds) < MIN_READS_PER_PHASE:
                unclustered_inds.extend(current_clust_inds)
            else:
                pass_clust_inds_list.append(read_clust_dat[0][clust_i])
        else:
            fail_clust_inds_list.append(read_clust_dat[0][clust_i])
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()
    print(f' - {len(fail_clust_inds_list)} clusters removed for not ending in enough tel sequence')
    print(f' - ({len(unclustered_inds)} tel reads unclustered)')

    #
    # [2] TVR CLUSTERING
    #

    sys.stdout.write('refining clusters [TVR]...')
    sys.stdout.flush()
    tt = time.perf_counter()
    clust_num = 0
    fail_tvr_clust = []
    clusters_with_tvrs = []
    all_consensus_tvrs = []
    blank_inds = []
    n_reads = 0
    if len(unclustered_inds) >= MIN_READS_PER_PHASE:
        pass_clust_inds_list += [unclustered_inds]
    for current_clust_inds in pass_clust_inds_list:
        my_chr = fake_chr
        if current_clust_inds == unclustered_inds:
            my_chr = unclust_chr
        current_clust_inds = sorted(current_clust_inds)
        khd_subset = [copy.deepcopy(kmer_hit_dat[n]) for n in current_clust_inds]
        zfcn = str(clust_num).zfill(3)
        telcompplot_fn = OUT_CDIR_TVR + 'results/' + 'reads_'       + zfcn + '.png'
        telcompcons_fn = OUT_CDIR_TVR + 'results/' + 'consensus_'   + zfcn + '.png'
        dendrogram_fn  = OUT_CDIR_TVR + 'dendro/'  + 'dendrogram_'  + zfcn + '.png'
        dist_matrix_fn = OUT_CDIR_TVR + 'npz/'     + 'dist-matrix_' + zfcn + '.npz'
        consensus_fn   = OUT_CDIR_TVR + 'fa/'      + 'consensus_'   + zfcn + '.fa'
        if ALWAYS_REPROCESS:
            rm(consensus_fn)
        my_dist_matrix = init_dist_matrix[:,current_clust_inds]
        my_dist_matrix = my_dist_matrix[current_clust_inds,:]
        np.savez_compressed(dist_matrix_fn, dist=my_dist_matrix)
        #
        subset_clustdat = cluster_tvrs(khd_subset, KMER_METADATA, my_chr, fake_pos, TREECUT_REFINE_TVR,
                                       aln_mode='ds',
                                       alignment_processes=NUM_PROCESSES,
                                       rand_shuffle_count=RAND_SHUFFLE,
                                       dist_in=dist_matrix_fn,
                                       fig_name=dendrogram_fn,
                                       save_msa=consensus_fn)
        #
        make_tvr_plots(khd_subset, subset_clustdat, my_chr, fake_pos, telcompplot_fn, telcompcons_fn, mtp_params)
        #
        for sci,subclust_inds in enumerate(subset_clustdat[0]):
            subclust_read_inds = [current_clust_inds[n] for n in subclust_inds]
            my_tvr_len = subset_clustdat[7][sci]
            # blank tvr? --> add to blank inds
            if my_tvr_len <= 0:
                blank_inds.extend(subclust_read_inds)
                continue
            # readcount filter
            if len(subclust_read_inds) < MIN_READS_PER_PHASE:
                continue
            n_reads += len(subclust_read_inds)
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
                fail_tvr_clust.append((term_zero_frac, nontel_long_frac, [n for n in subclust_read_inds]))
                continue
            #
            # pass all checks? --> add to list for subsequent subtel clustering
            #
            clusters_with_tvrs.append((my_chr, [n for n in subclust_read_inds]))
            all_consensus_tvrs.append(subset_clustdat[4][sci][:my_tvr_len])
        #
        clust_num += 1
    #
    if len(blank_inds) >= MIN_READS_PER_PHASE:
        clusters_with_tvrs.append((blank_chr, [n for n in blank_inds]))
    #
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()
    print(f' - {len(clusters_with_tvrs)} clusters with tvrs ({n_reads} reads)')
    print(f' - {len(fail_tvr_clust)} clusters removed for not ending in enough tel sequence')
    print(f' - ({len(blank_inds)} tel reads with blank tvr)')

    #
    # [3] SUBTEL CLUSTERING
    #

    all_consensus_tvr_subtel_pairs = []
    if SKIP_SUBTEL_REFINE:
        print('skipping subtel cluster refinement...')
        final_clustered_read_inds = []
        for sci,(my_chr,subclust_read_inds) in enumerate(clusters_with_tvrs):
            final_clustered_read_inds.append((my_chr, [n for n in subclust_read_inds]))
    else:
        sys.stdout.write('refining clusters [SUBTEL]...')
        sys.stdout.flush()
        tt = time.perf_counter()
        final_clustered_read_inds = []
        n_reads = 0
        for sci,(my_chr,subclust_read_inds) in enumerate(clusters_with_tvrs):
            subtel_sizes = [len(all_subtel_seq[n]) for n in subclust_read_inds]
            subtel_size = max(min(min(subtel_sizes), SUBTEL_CLUSTER_SIZE[1]), SUBTEL_CLUSTER_SIZE[0])
            my_subtels = [all_subtel_seq[n][-subtel_size:] for n in subclust_read_inds]
            #
            zfcn = str(sci).zfill(3)
            subtel_dendro_fn = OUT_CDIR_SUB + 'dendro/' + 'dendrogram_'  + zfcn + '.png'
            subtel_dist_fn   = OUT_CDIR_SUB + 'npz/'    + 'dist-matrix_' + zfcn + '.npz'
            subtel_muscle_fn = OUT_CDIR_SUB + 'fa/'     + 'consensus_'   + zfcn
            my_dendro_title  = zfcn
            subtel_labels    = None
            if ALWAYS_REPROCESS:
                rm(subtel_dist_fn)
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
                                                     dendrogram_title=my_dendro_title,
                                                     dendrogram_height=8)
            #
            if ALWAYS_REPROCESS or exists_and_is_nonzero(subtel_muscle_fn + '.fa') is False:
                with open(subtel_muscle_fn + '.fa', 'w') as f:
                    for sci_s,subclust_inds in enumerate(subtel_clustdat):
                        subsubclust_read_inds = [subclust_read_inds[n] for n in subclust_inds]
                        if len(subsubclust_read_inds) < MIN_READS_PER_PHASE:
                            continue
                        zfcn_s = str(sci_s).zfill(3)
                        clustered_subtels = [my_subtels[n] for n in subclust_inds]
                        subtel_consensus_seq = get_nucl_consensus(clustered_subtels)
                        f.write(f'>subtel-consensus_{zfcn}_{zfcn_s}\n')
                        f.write(f'{subtel_consensus_seq}\n')
            my_cons_subtels = [n[1] for n in quick_grab_all_reads(subtel_muscle_fn + '.fa')]
            for cons_subtel in my_cons_subtels:
                if my_chr != blank_chr:
                    all_consensus_tvr_subtel_pairs.append((all_consensus_tvrs[sci], cons_subtel))
                else:
                    all_consensus_tvr_subtel_pairs.append((None, cons_subtel))
            for sci_s,subclust_inds in enumerate(subtel_clustdat):
                subsubclust_read_inds = [subclust_read_inds[n] for n in subclust_inds]
                if len(subsubclust_read_inds) < MIN_READS_PER_PHASE:
                    continue
                n_reads += len(subsubclust_read_inds)
                final_clustered_read_inds.append((my_chr, [n for n in subsubclust_read_inds]))
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()
        print(f' - {len(final_clustered_read_inds)} clusters ({n_reads} reads)')

    #
    # [4] COLLAPSE CLUSTERS WITH VERY SIMILAR CONSENSUS TVRS + VERY SIMILAR CONSENSUS SUBTELS
    #     -- this is to mitigate false positives due to alleles that were split in initial clustering
    #     -- since we have consensus sequences at this point, comparisons should be easier
    #

    sys.stdout.write(f'checking for duplicate alleles...')
    sys.stdout.flush()
    tt = time.perf_counter()
    n_alleles_before_collapsing = len(final_clustered_read_inds)
    # [4a]: TVRs
    clust1 = []
    tvrs_to_compare = [n[0] for n in all_consensus_tvr_subtel_pairs if n[0] is not None]
    if tvrs_to_compare:
        collapse_dendro_fn = OUT_CDIR_COL + 'tvrs_dendrogram.png'
        collapse_dist_fn   = OUT_CDIR_COL + 'tvrs_dist-matrix.npz'
        collapse_fig_fn    = OUT_CDIR_COL + 'tvrs_plot.png'
        collapse_fa_fn     = OUT_CDIR_COL + 'tvrs.fa'
        with open(collapse_fa_fn, 'w') as f:
            for seq_i, seq in enumerate(tvrs_to_compare):
                f.write(f'>tvr_{seq_i}\n{seq}\n')
        if ALWAYS_REPROCESS:
            rm(collapse_dist_fn)
        clust1 = cluster_consensus_tvrs(tvrs_to_compare, KMER_METADATA, COLLAPSE_TVR_THRESH,
                                        dist_in=collapse_dist_fn,
                                        fig_name=collapse_fig_fn,
                                        dendro_name=collapse_dendro_fn,
                                        aln_mode='ms',
                                        gap_bool=(True,False),      # allow grouping prefixes
                                        adjust_lens=True,
                                        rand_shuffle_count=3,
                                        linkage_method='single',
                                        normalize_dist_matrix=True,
                                        alignment_processes=NUM_PROCESSES,
                                        dendrogram_xlim=[1,0])
    # [4b]: subtels in clusters with TVRs
    clust2 = []
    subs_to_compare = [n[1] for n in all_consensus_tvr_subtel_pairs if n[0] is not None]
    if subs_to_compare:
        collapse_dendro_fn = OUT_CDIR_COL + 'subtels_dendrogram.png'
        collapse_dist_fn   = OUT_CDIR_COL + 'subtels_dist-matrix.npz'
        collapse_fig_fn    = OUT_CDIR_COL + 'subtels_plot.png'
        collapse_fa_fn     = OUT_CDIR_COL + 'subtels.fa'
        with open(collapse_fa_fn, 'w') as f:
            for seq_i, seq in enumerate(subs_to_compare):
                f.write(f'>subtel_{seq_i}\n{seq}\n')
        if ALWAYS_REPROCESS:
            rm(collapse_dist_fn)
        clust2 = cluster_consensus_tvrs(subs_to_compare, KMER_METADATA, COLLAPSE_SUB_THRESH,
                                        dist_in=collapse_dist_fn,
                                        fig_name=collapse_fig_fn,
                                        dendro_name=collapse_dendro_fn,
                                        aln_mode='ms',
                                        gap_bool=(False,False),
                                        adjust_lens=True,
                                        rand_shuffle_count=3,
                                        linkage_method='single',
                                        normalize_dist_matrix=False,
                                        alignment_processes=NUM_PROCESSES,
                                        dendrogram_xlim=[1,0])
    # [4c]: subtels in clusters without TVRs
    clust3 = []
    subs_to_compare = [n[1] for n in all_consensus_tvr_subtel_pairs if n[0] is None]
    if subs_to_compare:
        collapse_dendro_fn = OUT_CDIR_COL + 'subtels-blanktvr_dendrogram.png'
        collapse_dist_fn   = OUT_CDIR_COL + 'subtels-blanktvr_dist-matrix.npz'
        collapse_fig_fn    = OUT_CDIR_COL + 'subtels-blanktvr_plot.png'
        collapse_fa_fn     = OUT_CDIR_COL + 'subtels-blanktvr.fa'
        with open(collapse_fa_fn, 'w') as f:
            for seq_i, seq in enumerate(subs_to_compare):
                f.write(f'>subtel-blank_{seq_i}\n{seq}\n')
        if ALWAYS_REPROCESS:
            rm(collapse_dist_fn)
        clust3 = cluster_consensus_tvrs(subs_to_compare, KMER_METADATA, COLLAPSE_SUB_THRESH,
                                        dist_in=collapse_dist_fn,
                                        fig_name=collapse_fig_fn,
                                        dendro_name=collapse_dendro_fn,
                                        aln_mode='ms',
                                        gap_bool=(False,False),
                                        adjust_lens=True,
                                        rand_shuffle_count=3,
                                        linkage_method='single',
                                        normalize_dist_matrix=False,
                                        alignment_processes=NUM_PROCESSES,
                                        dendrogram_xlim=[1,0])
        clust3 = [[n + len(tvrs_to_compare) for n in l3] for l3 in clust3]
    #
    intersections = []
    for l1 in clust1:
        for l2 in clust2:
            lint = [n for n in l1 if n in l2]
            if len(lint) > 1:
                intersections.append(lint)
    for l3 in clust3:
        if len(l3) > 1:
            intersections.append([n + len(tvrs_to_compare) for n in l3])
    if intersections:
        new_final_inds = []
        final_inds_seen = {}
        for i in range(len(final_clustered_read_inds)):
            if i in final_inds_seen:
                continue
            intersection_ind = None
            for j in range(len(intersections)):
                if i in intersections[j]:
                    intersection_ind = j
                    break
            if intersection_ind is not None:
                new_final_inds.append((final_clustered_read_inds[intersections[intersection_ind][0]][0],
                                       sum([final_clustered_read_inds[n][1] for n in intersections[intersection_ind]], [])))
                for k in intersections[intersection_ind]:
                    final_inds_seen[k] = True
            else:
                new_final_inds.append(final_clustered_read_inds[i])
                final_inds_seen[i] = True
        final_clustered_read_inds = new_final_inds
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()
    if len(final_clustered_read_inds) != n_alleles_before_collapsing:
        print(f' - {n_alleles_before_collapsing} clusters --> {len(final_clustered_read_inds)} clusters')
    else:
        print(' - no duplicates found')

    #
    # [5] FINAL RE-CLUSTERING
    #

    sys.stdout.write('final reclustering to get TVRs and TLs for each allele...')
    sys.stdout.flush()
    tt = time.perf_counter()
    clust_num = 0
    plot_num = 0
    subtels_out = []
    allele_outdat = []
    allele_consensus = []
    n_final_clusters_added = 0
    n_final_clusters_removed = 0
    for (my_chr, current_clust_inds) in final_clustered_read_inds:
        current_clust_inds = sorted(current_clust_inds)
        khd_subset = [copy.deepcopy(kmer_hit_dat[n]) for n in current_clust_inds]
        rlens_subset = [my_rlens[n] for n in current_clust_inds]
        zfcn = str(plot_num).zfill(3)
        telcompplot_fn = OUT_CDIR_FIN + 'results/' + 'reads_'       + zfcn + '.png'
        telcompcons_fn = OUT_CDIR_FIN + 'results/' + 'consensus_'   + zfcn + '.png'
        dendrogram_fn  = OUT_CDIR_FIN + 'dendro/'  + 'dendrogram_'  + zfcn + '.png'
        dist_matrix_fn = OUT_CDIR_FIN + 'npz/'     + 'dist-matrix_' + zfcn + '.npz'
        consensus_fn   = OUT_CDIR_FIN + 'fa/'      + 'consensus_'   + zfcn + '.fa'
        if ALWAYS_REPROCESS:
            rm(consensus_fn)
        my_dist_matrix = init_dist_matrix[:,current_clust_inds]
        my_dist_matrix = my_dist_matrix[current_clust_inds,:]
        np.savez_compressed(dist_matrix_fn, dist=my_dist_matrix)
        #
        if my_chr == unclust_chr:
            # if this cluster is from unclust_chr, we're allowed to scrutinize the TVRs
            solo_clustdat = cluster_tvrs(khd_subset, KMER_METADATA, my_chr, fake_pos, TREECUT_REFINE_TVR,
                                         aln_mode='ds',
                                         alignment_processes=NUM_PROCESSES,
                                         rand_shuffle_count=RAND_SHUFFLE,
                                         dist_in=dist_matrix_fn,
                                         fig_name=dendrogram_fn,
                                         save_msa=consensus_fn)
        else:
            # otherwise we're just reclustering to get tl boundaries so treecut can be set arbitrarily high
            solo_clustdat = cluster_tvrs(khd_subset, KMER_METADATA, my_chr, fake_pos, 100.,
                                         aln_mode='ms',
                                         alignment_processes=NUM_PROCESSES,
                                         rand_shuffle_count=1,
                                         dist_in=dist_matrix_fn,
                                         fig_name=dendrogram_fn,
                                         save_msa=consensus_fn)
        #
        my_tsvdat = get_allele_tsv_dat(khd_subset, solo_clustdat, my_chr, fake_pos, rlens_subset, gatd_params)
        prior_clustnum = clust_num
        for sci,subclust_inds in enumerate(solo_clustdat[0]):
            subclust_read_inds = [current_clust_inds[n] for n in subclust_inds]
            if len(subclust_read_inds) < MIN_READS_PER_PHASE:
                continue
            my_subtels = []
            for i,read_i in enumerate(subclust_read_inds):
                subtel_used_in_tvr = MIN_SUBTEL_BUFF - solo_clustdat[2][sci][i]
                if subtel_used_in_tvr >= 0:
                    my_subtels.append((read_i, all_subtel_seq[read_i][:len(all_subtel_seq[read_i])-subtel_used_in_tvr]))
                else:
                    my_subtels.append((read_i, all_subtel_seq[read_i]))
            for i,(read_i,subtel_seq) in enumerate(my_subtels):
                subtels_out.append((f'cluster-{len(allele_outdat)}_read-{i}_{my_rnames[read_i]}', subtel_seq))
            allele_outdat.append(my_tsvdat[sci])
            allele_consensus.append(solo_clustdat[4][sci])
            clust_num += 1
        if clust_num == prior_clustnum:
            n_final_clusters_removed += 1
        elif clust_num - prior_clustnum > 1:
            n_final_clusters_added += clust_num - prior_clustnum - 1
        #
        make_tvr_plots(khd_subset, solo_clustdat, my_chr, fake_pos, telcompplot_fn, telcompcons_fn, mtp_params)
        plot_num += 1
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
    #num_blank_alleles  = len([n[0] for n in ALLELE_TEL_DAT if n[0] == blank_chr])
    print(f' - {num_unique_alleles} unique alleles')
    #print(f' - {num_blank_alleles} / {num_unique_alleles} alleles with blank tvrs')
    print(f' - ({n_final_clusters_removed} clusters removed during tvr reanalysis)')
    print(f' - ({n_final_clusters_added} clusters added during tvr reanalysis)')

    #
    # SUBTEL ALIGNMENT
    #
    if SKIP_ANCHORING is False:
        #
        with gzip.open(OUT_UNANCHORED_SUBTELS, 'wt') as f:
            for n in subtels_out:
                if len(n[1]) > 0:
                    f.write(f'>{n[0]}\n{n[1]}\n')
        #
        aln_log = ALIGNED_SUBTELS + '.log'
        aln_sam = ALIGNED_SUBTELS + '.sam'
        aln_bam = ALIGNED_SUBTELS + '.bam'
        if exists_and_is_nonzero(aln_bam) and NO_SUBTEL_REALIGN:
            pass
        else:
            sys.stdout.write('aligning subtels to reference...')
            sys.stdout.flush()
            tt = time.perf_counter()
            cmd = ''
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
                with gzip.open(TELOGATOR_REF,'rt') as f_tr: # pbmm2 wants an uncompressed reference
                    with open(UNCOMPRESSED_TELOGATOR_REF,'w') as f_utr:
                        for line in f_tr:
                            f_utr.write(line)
                pysam.faidx(UNCOMPRESSED_TELOGATOR_REF)
                cmd = ALIGNER_EXE + ' align ' + UNCOMPRESSED_TELOGATOR_REF + ' ' + OUT_UNANCHORED_SUBTELS + ' ' + aln_bam + ' --preset HiFi --sort'
            if len(cmd):
                with open(aln_log, 'w') as f:
                    try:
                        subprocess.check_output(cmd.split(' '), stderr=f, text=True)
                    except subprocess.CalledProcessError as exc:
                        print('Error: alignment command returned an error:', exc.returncode)
                        print()
                        print(exc.output)
                        print()
                        print('Here is the command we tried to run:')
                        print()
                        print(cmd)
                        print()
                        exit(1)
            if exists_and_is_nonzero(aln_sam) is False and exists_and_is_nonzero(aln_bam) is False:
                print()
                print('Error: subtel alignment failed. Check the log:')
                print(f'{aln_log}')
                exit(1)
            #
            if WHICH_ALIGNER == 'pbmm2':
                rm(UNCOMPRESSED_TELOGATOR_REF)
                rm(UNCOMPRESSED_TELOGATOR_REF+'.fai')
            else:
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
            [rnm, pos, read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat] = parse_read(sam_line)
            if rnm not in ALIGNMENTS_BY_RNAME:
                ALIGNMENTS_BY_RNAME[rnm] = []
            ALIGNMENTS_BY_RNAME[rnm].append([read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat])
        samfile.close()
        #
        top_alns_by_cluster = {}
        best_mapq_by_readname = {}
        for readname in ALIGNMENTS_BY_RNAME:
            original_readname = '_'.join(readname.split('_')[2:])
            sorted_choices = []
            for i,aln in enumerate(ALIGNMENTS_BY_RNAME[readname]):
                dist_to_end = len(aln[7]) - max(aln[0], aln[1])
                read_span = abs(aln[1] - aln[0])
                my_mapq = aln[6]
                sorted_choices.append((read_span, my_mapq, -dist_to_end, i))
                if original_readname not in best_mapq_by_readname:
                    best_mapq_by_readname[original_readname] = my_mapq
                elif my_mapq > best_mapq_by_readname[original_readname]:
                    best_mapq_by_readname[original_readname] = my_mapq
            sorted_choices = sorted(sorted_choices, reverse=True)
            top_i = sorted_choices[0][3]
            my_cluster = int(readname.split('_')[0][8:])
            if my_cluster not in top_alns_by_cluster:
                top_alns_by_cluster[my_cluster] = []
            top_alns_by_cluster[my_cluster].append((readname, copy.deepcopy(ALIGNMENTS_BY_RNAME[readname][top_i])))
        #
        for clustnum in sorted(top_alns_by_cluster.keys()):
            chr_arm_scores = {}
            anchors_by_ref = {}
            for n in top_alns_by_cluster[clustnum]:
                my_ref = n[1][2]
                my_refspan = abs(n[1][3] - n[1][4])
                my_anchor = n[1][4] # second refpos is always anchor coord?
                #my_orr = n[1][5]
                my_mapq = n[1][6]
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
                    my_pos  = int(np.median(anchors_by_ref[n]))
                    if my_chr in chrs_encountered_so_far:
                        continue
                    chrs_encountered_so_far[my_chr] = True
                    anchor_pos.append(str(my_pos))
                    ref_builds.append(my_samp)
                    my_chrs.append(my_chr)
                ALLELE_TEL_DAT[clustnum][0] = ','.join(my_chrs)    # assign chr
                ALLELE_TEL_DAT[clustnum][1] = ','.join(anchor_pos) # assign pos
                ALLELE_TEL_DAT[clustnum][2] = ','.join(ref_builds) # assign ref builds
                # annotate alignments to known interstitial telomere regions
                any_interstitial = any([annotate_interstitial_tel(ref_builds[n]+'_'+my_chrs[n], int(anchor_pos[n])) for n in range(len(my_chrs))])
                if any_interstitial:
                    ALLELE_TEL_DAT[clustnum][3] += 'i'
                # assign mapq values
                clust_readnames = ALLELE_TEL_DAT[clustnum][10].split(',')
                clust_mapq = []
                for og_readname in clust_readnames:
                    if og_readname in best_mapq_by_readname:
                        clust_mapq.append(best_mapq_by_readname[og_readname])
                    else:
                        clust_mapq.append(0)
                ALLELE_TEL_DAT[clustnum][7] = ','.join([str(n) for n in clust_mapq])
        #
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()

        #
        # resort out allele dat by chr & pos
        #
        sorted_ad_inds = sorted([(LEXICO_2_IND[n[0].split(',')[0][:-1]], n[0].split(',')[0][-1], int(n[1].split(',')[0]), i) for i,n in enumerate(ALLELE_TEL_DAT)])
        ALLELE_TEL_DAT = [ALLELE_TEL_DAT[n[3]] for n in sorted_ad_inds]

    #
    # collapse homozygous alleles
    #
    if COLLAPSE_HOM_BP > 0:
        sys.stdout.write(f'merging alleles mapped within {COLLAPSE_HOM_BP}bp of each other...')
        sys.stdout.flush()
        tt = time.perf_counter()
        n_alleles_before_collapsing = len(ALLELE_TEL_DAT)
        while True:
            del_ind = None
            for i in range(len(ALLELE_TEL_DAT)):
                for j in range(i+1, len(ALLELE_TEL_DAT)):
                    if all([ALLELE_TEL_DAT[i][0].split(',')[0] == ALLELE_TEL_DAT[j][0].split(',')[0],
                            ALLELE_TEL_DAT[i][2].split(',')[0] == ALLELE_TEL_DAT[j][2].split(',')[0],
                            abs(int(ALLELE_TEL_DAT[i][1].split(',')[0]) - int(ALLELE_TEL_DAT[j][1].split(',')[0])) <= COLLAPSE_HOM_BP]):
                        my_dist = 0.
                        tvr_i, tvr_j = ALLELE_TEL_DAT[i][9], ALLELE_TEL_DAT[j][9]
                        if len(tvr_i) and len(tvr_j):
                            my_dist = quick_compare_tvrs(tvr_i, tvr_j)
                        # handling blank tvrs
                        elif len(tvr_i) and len(tvr_j) == 0:
                            my_dist = quick_compare_tvrs(tvr_i, canonical_letter*len(tvr_i))
                        elif len(tvr_i) == 0 and len(tvr_j):
                            my_dist = quick_compare_tvrs(canonical_letter*len(tvr_j), tvr_j)
                        if my_dist <= COLLAPSE_TVR_THRESH:
                            which_ind, merged_dat = merge_allele_tsv_dat(ALLELE_TEL_DAT[i], ALLELE_TEL_DAT[j], ALLELE_TL_METHOD)
                            if which_ind == 0:
                                ALLELE_TEL_DAT[i] = merged_dat
                                del_ind = j
                            else:
                                ALLELE_TEL_DAT[j] = merged_dat
                                del_ind = i
                            break
                if del_ind is not None:
                    break
            if del_ind is not None:
                del ALLELE_TEL_DAT[del_ind]
            else:
                break
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()
        print(f' - {n_alleles_before_collapsing} --> {len(ALLELE_TEL_DAT)} alleles')

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
    # violin plot
    #
    if SKIP_ANCHORING is False:
        atl_by_arm = [{} for n in range(max(1,VIOLIN_PLOIDY))]
        for atd in ALLELE_TEL_DAT:
            my_chr = atd[0].split(',')[0]
            my_tvr_len = int(atd[8])
            if my_chr in [fake_chr, blank_chr, unclust_chr]:
                my_chr = 'unanchored'
            for i in range(len(atl_by_arm)):
                if my_chr not in atl_by_arm[i]:
                    atl_by_arm[i][my_chr] = [my_tvr_len + int(n) for n in atd[5].split(',')]
                    break
        vparams = {'norm_by_readcount':False,
                   'include_unanchored':False,
                   'p_ymax':VIOLIN_YMAX,
                   'q_ymax':VIOLIN_YMAX,
                   'y_step':VIOLIN_TICK,
                   'y_label':'<-- q      ATL      p -->'}
        tel_len_violin_plot(atl_by_arm, VIOLIN_ATL, custom_plot_params=vparams)

    #
    # plot all final tvrs
    #
    tvrs_to_plot = []
    tvr_labels_to_plot = []
    clustdat_to_plot = [[[]], [[]], [[]], []]
    current_i = 0
    num_alleles_too_short = 0
    num_alleles_interstitial = 0
    num_alleles_unmapped = 0
    num_pass_alleles = 0
    readlens_final = []
    for atd in ALLELE_TEL_DAT:
        my_id = atd[3]
        while my_id[-1].isdigit() is False:
            my_id = my_id[:-1]
        my_id = int(my_id)
        my_rep_atl = int(atd[4])
        my_max_atl = max([int(n) for n in atd[5].split(',')])
        my_tvr_len = int(atd[8])
        if atd[0] in [fake_chr, blank_chr, unclust_chr]:
            num_alleles_unmapped += 1
            continue
        if my_max_atl < MIN_ATL_FOR_FINAL_PLOTTING:
            num_alleles_too_short += 1
            continue
        if 'i' in atd[3]:
            num_alleles_interstitial += 1
            continue
        num_pass_alleles += 1
        if int(atd[8]) <= 0:
            my_annot = str(atd[3]) + 'b'
        else:
            my_annot = str(atd[3])
        readlens_final.extend([int(n) for n in atd[6].split(',')])
        tvrs_to_plot.append(allele_consensus[my_id][:my_tvr_len+my_rep_atl])
        tvr_labels_to_plot.append(f'({my_annot}) {atd[0]}')
        clustdat_to_plot[0][0].append(current_i)
        clustdat_to_plot[1][0].append(60)
        clustdat_to_plot[2][0].append(0)
        clustdat_to_plot[3].append(0)
        current_i += 1
    #
    in_npz = np.load(READLEN_NPZ)
    readlens_all = in_npz['readlen_all']
    readlens_tel = in_npz['readlen_tel']
    readlen_plot(readlens_all, readlens_tel, readlens_final, QC_READLEN)
    rm(READLEN_NPZ)
    #
    if len(tvrs_to_plot):
        redrawn_tvrs = convert_colorvec_to_kmerhits(tvrs_to_plot, KMER_METADATA)
        clust_khd = []
        for i in range(len(redrawn_tvrs)):
            clust_khd.append([redrawn_tvrs[i], len(tvrs_to_plot[i]), 0, 'FWD', tvr_labels_to_plot[i], MAX_QUAL_SCORE, None])
        custom_plot_params = {'xlim':[0,FINAL_PLOTTING_XMAX],
                              'xstep':FINAL_PLOTTING_XTICK,
                              'custom_title':'',
                              'number_label_rows':False,
                              'font.size':16,
                              'font.weight':'normal'}
        plot_kmer_hits(clust_khd, KMER_COLORS, '', 0, FINAL_TVRS, clust_dat=clustdat_to_plot, plot_params=custom_plot_params)
    else:
        print('no alleles to plot.')
    print(f' - {num_pass_alleles} final telomere alleles')
    print(f' - {num_alleles_unmapped} (unmapped)')
    print(f' - {num_alleles_too_short} (atl < {MIN_ATL_FOR_FINAL_PLOTTING} bp)')
    print(f' - {num_alleles_interstitial} (interstitial)')
    with open(QC_STATS, 'w') as f:
        f.write(f' - {num_pass_alleles} final telomere alleles\n')
        f.write(f' - {num_alleles_unmapped} (unmapped)\n')
        f.write(f' - {num_alleles_too_short} (atl < {MIN_ATL_FOR_FINAL_PLOTTING} bp)\n')
        f.write(f' - {num_alleles_interstitial} (interstitial)\n')


if __name__ == '__main__':
    main()
