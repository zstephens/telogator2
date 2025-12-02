#!/usr/bin/env python
import argparse
import copy
import gzip
import importlib.resources as ir
import numpy as np
import os
import pysam
import random
import subprocess
import sys
import time

from source.tg_align  import quick_compare_tvrs
from source.tg_kmer   import get_canonical_letter, read_kmer_tsv
from source.tg_plot   import convert_colorvec_to_kmerhits, make_tvr_plots, plot_fusion, plot_kmer_hits, plot_some_tvrs, readlen_plot, tel_len_violin_plot
from source.tg_reader import quick_grab_all_reads_nodup, TG_Reader
from source.tg_tel    import get_allele_tsv_dat, get_tel_repeat_comp_parallel, merge_allele_tsv_dat
from source.tg_tvr    import cluster_consensus_tvrs, cluster_tvrs, quick_get_tvrtel_lens
from source.tg_util   import annotate_interstitial_tel, check_aligner_exe, dir_exists, exists_and_is_nonzero, get_downsample_inds, get_file_type, LEXICO_2_IND, makedir, mv, parse_read, RC, rm, strip_paths_from_string, BLANK_CHR, UNCLUST_CHR, UNCLUST_POS

TEL_WINDOW_SIZE = 100
P_VS_Q_AMP_THRESH = 0.5
MIN_TEL_SCORE = 100
#
DUMMY_TEL_MAPQ = 60
MAX_QUAL_SCORE = 60
ANCHORING_ASSIGNMENT_FRAC = 0.150
#
HUMAN_GENOME_BP = 3100000000


def main(raw_args=None):
    parser = argparse.ArgumentParser(description='Telogator', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    #
    parser.add_argument('--version', action='version', version='%(prog)s 2.2.3')
    #
    parser.add_argument('-i', type=str, required=True,  metavar='',             nargs='*',      help="* Input reads (fa / fa.gz / fq / fq.gz / bam)")
    parser.add_argument('-o', type=str, required=True,  metavar='output/',                      help="* Path to output directory")
    parser.add_argument('-k', type=str, required=False, metavar='kmers.tsv',    default='',     help="Telomere kmers file")
    parser.add_argument('-t', type=str, required=False, metavar='telogator.fa', default='',     help="Telogator reference fasta")
    parser.add_argument('-r', type=str, required=False, metavar='ont',          default='ont',  help="Read type: hifi / ont")
    parser.add_argument('-l', type=int, required=False, metavar='4000',         default=4000,   help="Minimum read length")
    parser.add_argument('-c', type=int, required=False, metavar='8',            default=8,      help="Minimum hits to tandem canonical kmer")
    parser.add_argument('-n', type=int, required=False, metavar='3',            default=3,      help="Minimum number of reads per cluster")
    parser.add_argument('-m', type=str, required=False, metavar='p75',          default='p75',  help="Method for choosing ATL: mean / median / p75 / max")
    parser.add_argument('-d', type=int, required=False, metavar='-1',           default=-1,     help="Downsample to this many telomere reads")
    parser.add_argument('-p', type=int, required=False, metavar='4',            default=4,      help="Number of processes to use")
    #
    parser.add_argument('--filt-tel',    type=int, required=False, metavar='400',  default=400,  help="[FILTERING] Remove reads that end in < this much tel")
    parser.add_argument('--filt-nontel', type=int, required=False, metavar='100',  default=100,  help="[FILTERING] Remove reads that end in > this much non-tel")
    parser.add_argument('--filt-sub',    type=int, required=False, metavar='1000', default=1000, help="[FILTERING] Remove reads that end in < this much subtel")
    #
    parser.add_argument('-t0', type=float, required=False, metavar='0.200', default=0.200, help="[TREECUT] TVR clustering (iteration 0)")
    parser.add_argument('-t1', type=float, required=False, metavar='0.150', default=0.150, help="[TREECUT] TVR clustering (iteration 1)")
    parser.add_argument('-t2', type=float, required=False, metavar='0.100', default=0.100, help="[TREECUT] TVR clustering (iteration 2)")
    parser.add_argument('-tc', type=float, required=False, metavar='0.050', default=0.050, help="[TREECUT] TVR clustering (collapse)")
    parser.add_argument('-ts', type=float, required=False, metavar='0.200', default=0.200, help="[TREECUT] Subtel cluster refinement")
    parser.add_argument('-th', type=float, required=False, metavar='0.050', default=0.050, help="[TREECUT] Collapsing aligned alleles")
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
    parser.add_argument('--plot-t1clust',   required=False, action='store_true', default=False, help="[DEBUG] Plot tel reads during TVR clustering (iteration 1)")
    parser.add_argument('--plot-t2clust',   required=False, action='store_true', default=False, help="[DEBUG] Plot tel reads during TVR clustering (iteration 2)")
    parser.add_argument('--plot-finclust',  required=False, action='store_true', default=False, help="[DEBUG] Plot tel reads during final clustering")
    parser.add_argument('--plot-signals',   required=False, action='store_true', default=False, help="[DEBUG] Plot tel signals for all reads")
    parser.add_argument('--plot-fusions',   required=False, action='store_true', default=False, help="[DEBUG] Plot tel fusions for all reads")
    parser.add_argument('--dont-reprocess', required=False, action='store_true', default=False, help="[DEBUG] Use existing intermediary files (for redoing failed runs)")
    parser.add_argument('--debug-replot',   required=False, action='store_true', default=False, help="[DEBUG] Regenerate plots that already exist")
    parser.add_argument('--debug-nosubtel', required=False, action='store_true', default=False, help="[DEBUG] Skip subtel cluster refinement")
    parser.add_argument('--debug-noanchor', required=False, action='store_true', default=False, help="[DEBUG] Do not align reads or do any anchoring")
    parser.add_argument('--debug-telreads', required=False, action='store_true', default=False, help="[DEBUG] Stop immediately after extracting tel reads")
    parser.add_argument('--debug-progress', required=False, action='store_true', default=False, help="[DEBUG] Print progress to screen during init matrix computation")
    parser.add_argument('--debug-dendro',   required=False, action='store_true', default=False, help="[DEBUG] Plot dendrogram during init clustering")
    parser.add_argument('--debug-collapse', required=False, action='store_true', default=False, help="[DEBUG] Plot TVR distances during final collapse")
    parser.add_argument('--fast-aln',       required=False, action='store_true', default=False, help="[PERFORMANCE] Use faster but less accurate pairwise alignment")
    #
    parser.add_argument('--collapse-hom',  type=int, required=False, metavar='',        default=1000, help="Merge similar alleles mapped <= this bp from each other")
    parser.add_argument('--minimap2',      type=str, required=False, metavar='exe',     default='',   help="/path/to/minimap2")
    parser.add_argument('--pbmm2',         type=str, required=False, metavar='exe',     default='',   help="/path/to/pbmm2")
    parser.add_argument('--winnowmap',     type=str, required=False, metavar='exe',     default='',   help="/path/to/winnowmap")
    parser.add_argument('--winnowmap-k15', type=str, required=False, metavar='k15.txt', default='',   help="high freq kmers file (only needed for winnowmap)")
    parser.add_argument('--ref',           type=str, required=False, metavar='ref.fa',  default='',   help="Reference filename (only needed if input is cram)")
    parser.add_argument('--rng',           type=int, required=False, metavar='-1',      default=-1,   help="RNG seed value")
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
    MIN_READS_PER_PHASE = args.n
    ALLELE_TL_METHOD    = args.m
    DOWNSAMPLE_READS    = args.d
    NUM_PROCESSES       = args.p
    CRAM_REF_FILE       = args.ref
    RNG_SEED            = args.rng
    #
    TREECUT_TVR_ITER_0    = args.t0
    TREECUT_TVR_ITER_1    = args.t1
    TREECUT_TVR_ITER_2    = args.t2
    TREECUT_TVR_COLLAPSE  = args.tc
    TREECUT_REFINE_SUBTEL = args.ts
    TREECUT_HOM_COLLAPSE  = args.th
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
    WINNOWMAP_K15 = args.winnowmap_k15
    PBMM2_EXE     = args.pbmm2

    # debug params
    #
    DONT_OVERWRITE_PLOTS  = not(args.debug_replot)     # (True = don't replot figures if they already exist)
    ALWAYS_REPROCESS      = not(args.dont_reprocess)   # (True = don't write out .npz matrices, always recompute)
    SKIP_SUBTEL_REFINE    = args.debug_nosubtel
    SKIP_ANCHORING        = args.debug_noanchor
    STOP_AFTER_TELREADS   = args.debug_telreads
    PRINT_PROGRESS        = args.debug_progress
    PLOT_INIT_DENDRO      = args.debug_dendro
    DEBUG_COLLAPSE        = args.debug_collapse
    PLOT_TVR_CLUST_ITER_1 = args.plot_t1clust
    PLOT_TVR_CLUST_ITER_2 = args.plot_t2clust
    PLOT_FINAL_CLUST      = args.plot_finclust
    PLOT_FILT_CVECS       = args.plot_filt_tvr
    PLOT_TEL_SIGNALS      = args.plot_signals
    PLOT_TEL_FUSIONS      = args.plot_fusions
    FILT_TERM_TEL         = args.filt_tel         # reads must terminate in at least this much tel sequence
    FILT_TERM_NONTEL      = args.filt_nontel      # remove reads that have more than this much terminating nontel sequence
    FILT_TERM_SUBTEL      = args.filt_sub         # reads must terminate in at least this much subtel sequence
    FAST_ALIGNMENT        = args.fast_aln
    COLLAPSE_HOM_BP       = args.collapse_hom

    MIN_INTERSTITIAL_TL = 100
    MIN_FUSION_ANCHOR = 1000

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
    if len(set(INPUT_ALN)) < len(INPUT_ALN):
        print('Error: an input file was provided multiple times')
        exit(1)
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
    OUT_CLUST_DIR    = OUT_DIR + 'temp/'
    OUT_DIR_INIT     = OUT_CLUST_DIR + '01_initial/'
    OUT_DIR_UNCLUST  = OUT_CLUST_DIR + '02_unclust/'
    OUT_DIR_COLLAPSE = OUT_CLUST_DIR + '03_collapse/'
    OUT_DIR_SUBTEL   = OUT_CLUST_DIR + '04_subtel/'
    OUT_DIR_FINAL    = OUT_CLUST_DIR + '05_final/'
    OUT_DIR_ANCHOR   = OUT_CLUST_DIR + '06_anchoring/'
    OUT_DIR_COL2     = OUT_CLUST_DIR + '07_collapse/'
    OUT_QC_DIR       = OUT_DIR + 'qc/'
    makedir(OUT_DIR)
    if dir_exists(OUT_DIR) is False:
        print('Error: could not create output directory')
        exit(1)
    makedir(OUT_CLUST_DIR)
    makedir(OUT_QC_DIR)
    for d in [OUT_DIR_INIT, OUT_DIR_UNCLUST, OUT_DIR_SUBTEL, OUT_DIR_FINAL]:
        makedir(d)
        makedir(d + 'dendro/')
        makedir(d + 'fa/')
        makedir(d + 'npz/')
        makedir(d + 'results/')
    makedir(OUT_DIR_COLLAPSE)
    makedir(OUT_DIR_ANCHOR)
    makedir(OUT_DIR_COL2)

    # optional dirs
    TEL_SIGNAL_DIR = OUT_DIR + 'tel_signal/'
    if PLOT_TEL_SIGNALS:
        makedir(TEL_SIGNAL_DIR)
    TEL_FUSION_DIR = OUT_DIR + 'tel_fusion/'
    if PLOT_TEL_FUSIONS:
        makedir(TEL_FUSION_DIR)

    TELOMERE_READS = OUT_CLUST_DIR + 'tel_reads.fa.gz'
    VIOLIN_ATL = OUT_DIR + 'violin_atl.png'
    FINAL_TVRS = OUT_DIR + 'all_final_alleles.png'
    FINAL_TSV  = OUT_DIR + 'tlens_by_allele.tsv'
    READLEN_NPZ = OUT_QC_DIR + 'readlens.npz'
    QC_READLEN  = OUT_QC_DIR + 'qc_readlens.png'
    QC_CMD      = OUT_QC_DIR + 'cmd.txt'
    QC_STATS    = OUT_QC_DIR + 'stats.txt'
    QC_RNG      = OUT_QC_DIR + 'rng.txt'
    UNCOMPRESSED_TELOGATOR_REF = OUT_CLUST_DIR + 'telogator-ref.fa'
    UNANCHORED_SUBTELS = OUT_DIR_ANCHOR + 'unanchored_subtels.fa.gz'
    ALIGNED_SUBTELS    = OUT_DIR_ANCHOR + 'subtel_aln'

    # rng
    if RNG_SEED <= -1:
        RNG_SEED = random.randint(1, 99999999)
    random.seed(RNG_SEED)
    np.random.seed(RNG_SEED)
    np.random.default_rng(RNG_SEED)
    with open(QC_RNG, 'w') as f:
        f.write(str(RNG_SEED) + '\n')

    # only compare up to this much tvr / subtel sequence for clustering
    TVR_TRUNCATE_ITER_0 = 1000
    TVR_TRUNCATE_ITER_1 = 2000
    TVR_TRUNCATE_ITER_2 = 3000
    SUBTEL_TRUNCATE = 3000
    RAND_SHUFFLE_ITER_0 = 2
    RAND_SHUFFLE_ITER_1 = 3
    RAND_SHUFFLE_ITER_2 = 3
    TVR_REFINE_ITERATIONS = [1,2]
    TVR_ITERATIONS_TO_PLOT = [1]*PLOT_TVR_CLUST_ITER_1 + [2]*PLOT_TVR_CLUST_ITER_2
    if FAST_ALIGNMENT:
        RAND_SHUFFLE_ITER_0 = 1
        RAND_SHUFFLE_ITER_1 = 1
        RAND_SHUFFLE_ITER_2 = 1
        TVR_REFINE_ITERATIONS = [2]

    # TVR clustering parameters by iteration
    tvr_clust_params = [None, None, None]
    tvr_clust_params[0] = {'aln_mode':'ds',
                           'tvr_truncate':TVR_TRUNCATE_ITER_0,
                           'rand_shuffle_count':RAND_SHUFFLE_ITER_0,
                           'tree_cut':TREECUT_TVR_ITER_0,
                           'clustering_only':True,
                           'num_processes':NUM_PROCESSES,
                           'rng_seed':RNG_SEED,
                           'print_matrix_progress':PRINT_PROGRESS}
    #
    tvr_clust_params[1] = {'aln_mode':'ds',
                           'tvr_truncate':TVR_TRUNCATE_ITER_1,
                           'rand_shuffle_count':RAND_SHUFFLE_ITER_1,
                           'tree_cut':TREECUT_TVR_ITER_1,
                           'clustering_only':not(PLOT_TVR_CLUST_ITER_1),
                           'num_processes':NUM_PROCESSES,
                           'rng_seed':RNG_SEED,
                           'print_matrix_progress':False}
    #
    tvr_clust_params[2] = {'aln_mode':'ds',
                           'tvr_truncate':TVR_TRUNCATE_ITER_2,
                           'rand_shuffle_count':RAND_SHUFFLE_ITER_2,
                           'tree_cut':TREECUT_TVR_ITER_2,
                           'clustering_only':False,
                           'num_processes':NUM_PROCESSES,
                           'rng_seed':RNG_SEED,
                           'print_matrix_progress':False}
    #
    collapse_clust_params = {'aln_mode':'ms',
                             'gap_bool':(True,False),      # allow grouping prefixes
                             'adjust_lens':False,
                             'rand_shuffle_count':3,
                             'tree_cut':TREECUT_TVR_COLLAPSE,
                             'linkage_method':'single',
                             'normalize_dist_matrix':True,
                             'num_processes':NUM_PROCESSES,
                             'rng_seed':RNG_SEED,
                             'overwrite_figures':not(DONT_OVERWRITE_PLOTS),
                             'dendrogram_xlim':[1.1,0]}
    #
    subtel_clust_params = {'aln_mode':'ms',
                           'tvr_truncate':SUBTEL_TRUNCATE,
                           'gap_bool':(True,False),
                           'adjust_lens':False,
                           'rand_shuffle_count':3,
                           'tree_cut':TREECUT_REFINE_SUBTEL,
                           'linkage_method':'ward',
                           'normalize_dist_matrix':True,
                           'num_processes':NUM_PROCESSES,
                           'rng_seed':RNG_SEED,
                           'dendrogram_height':8,
                           'overwrite_figures':not(DONT_OVERWRITE_PLOTS)}
    #
    final_clust_params = {'aln_mode':'ds',
                          'tvr_truncate':1000,
                          'rand_shuffle_count':1,
                          'tree_cut':100.,                 # arbitrarily high
                          'clustering_only':False,
                          'num_processes':NUM_PROCESSES,
                          'rng_seed':RNG_SEED}

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
    # find subtel reference seq
    #
    if TELOGATOR_REF != '':
        fn_suffix = TELOGATOR_REF.split('/')[-1]
        if not os.path.exists(TELOGATOR_REF):
            # try to find it using importlib.resources
            resource_path = str(ir.files("source").joinpath(TELOGATOR_REF))
            if os.path.exists(resource_path):
                TELOGATOR_REF = resource_path
            else:
                raise FileNotFoundError(f'Provided subtel reference (-t) could not be found: {TELOGATOR_REF} or at {resource_path}')
        print(f'using user-specified subtel reference: {fn_suffix}')
        print(f' - found at {TELOGATOR_REF}')
    else:
        print('using default subtel reference')
        TELOGATOR_REF = str(ir.files("source").joinpath('resources/telogator-ref.fa.gz'))
        WINNOWMAP_K15 = str(ir.files("source").joinpath('resources/telogator-ref-k15.txt'))
    if WHICH_ALIGNER == 'winnowmap' and WINNOWMAP_K15 == '':
        print('Error: winnowmap was chosen as aligner but no --winnowmap-k15 was provided.')
        exit(1)

    #
    # parse telomere kmer data
    #
    if KMER_FILE != '':
        fn_suffix = KMER_FILE.split('/')[-1]
        if not os.path.exists(KMER_FILE):
            # try to find it using importlib.resources
            kmer_path = str(ir.files("source").joinpath(KMER_FILE))
            if os.path.exists(kmer_path):
                KMER_FILE = kmer_path
            else:
                raise FileNotFoundError(f'Provided telomere kmers (-k) could not be found: {KMER_FILE} or at {kmer_path}')
        print(f'using user-specified telomere kmers: {fn_suffix}')
        print(f' - found at {KMER_FILE}')
    else:
        print('using default telomere kmers')
        KMER_FILE = str(ir.files("source").joinpath('resources/kmers.tsv'))
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
    gtrc_params = [READ_TYPE, DUMMY_TEL_MAPQ,
                   KMER_LIST, KMER_LIST_REV, KMER_ISSUBSTRING, CANONICAL_STRINGS, CANONICAL_STRINGS_REV,
                   TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH,
                   MIN_TEL_SCORE, FILT_TERM_TEL, FILT_TERM_NONTEL, FILT_TERM_SUBTEL,
                   PLOT_TEL_SIGNALS, TEL_SIGNAL_DIR,
                   MIN_INTERSTITIAL_TL, MIN_FUSION_ANCHOR]

    #
    # write cmd out, for debugging purposes
    #
    with open(QC_CMD, 'w') as f:
        stripped_strings = [strip_paths_from_string(n) for n in sys.argv]
        f.write(' '.join(stripped_strings) + '\n')

    #=====================================================#
    #
    # [0a] GET TELOMERE READS
    #
    #=====================================================#

    ALLELE_TEL_DAT = []

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
                    if not my_rdat:
                        # this can happen in aligned bam if sequence is not present (because read is multimapped maybe?)
                        continue
                    if my_issup:
                        sup_readcount += 1
                        continue
                    count_fwd = my_rdat.count(KMER_INITIAL)
                    count_rev = 0
                    if count_fwd < MINIMUM_CANON_HITS:  # don't need to do this is we already have enough fwd hits
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
        print(f' - ({total_bp_all} bp) [{int(total_bp_all/HUMAN_GENOME_BP + 0.5)}x] --> ({total_bp_tel} bp)')
        if tel_readcount <= 0:
            print('Error: No telomere reads found, stopping here...')
            exit(1)
    else:
        print('found telomere reads from a previous run, using them instead of reprocessing input file')
    #
    if STOP_AFTER_TELREADS:
        print('--debug-telreads was specified, stopping now.')
        exit(0)
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

    #=====================================================#
    #
    # [0b] GET TELOMERE REPEAT COMPOSITION
    #
    #=====================================================#

    sys.stdout.write('getting telomere repeat composition...')
    sys.stdout.flush()
    tt = time.perf_counter()
    num_starting_reads = len(all_read_dat)
    #
    gtrc_results = get_tel_repeat_comp_parallel(all_read_dat, gtrc_params, max_workers=NUM_PROCESSES)
    [kmer_hit_dat,
     all_terminating_tl,
     all_nontel_end,
     reads_removed_term_tel,
     reads_removed_term_unk,
     reads_removed_term_sub,
     interstitial_subtels,
     interstitial_khd_by_readname,
     interstitial_khd_rev_by_readname] = gtrc_results
    #
    num_ending_reads = len(kmer_hit_dat)
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()
    print(f' - {num_starting_reads} --> {num_ending_reads} reads')
    if FILT_TERM_TEL > 0:
        print(f' - ({reads_removed_term_tel} reads removed for ending in < {FILT_TERM_TEL}bp tel sequence)')
    if FILT_TERM_NONTEL > 0:
        print(f' - ({reads_removed_term_unk} reads removed for ending in > {FILT_TERM_NONTEL}bp non-tel sequence)')
    if FILT_TERM_SUBTEL > 0:
        print(f' - ({reads_removed_term_sub} reads removed for ending in < {FILT_TERM_SUBTEL}bp subtel sequence)')
    #
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
            del all_terminating_tl[di]
            del all_nontel_end[di]
        num_ending_reads = len(kmer_hit_dat)
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads\n')
        sys.stdout.flush()
    #
    my_rlens = [len(n[6]) for n in kmer_hit_dat]
    my_rnames = [n[4] for n in kmer_hit_dat]

    #=====================================================#
    #
    # [1a] INITIAL CLUSTERING
    #
    #=====================================================#

    sys.stdout.write('initial clustering of all reads [iteration 0]...')
    sys.stdout.flush()
    if PRINT_PROGRESS:
        print()
    tt = time.perf_counter()
    tvr_clust_params[0]['fig_name'] = None
    tvr_clust_params[0]['dist_in'] = f'{OUT_DIR_INIT}npz/i0_dist-matrix.npz'
    if PLOT_INIT_DENDRO:
        tvr_clust_params[0]['fig_name'] = f'{OUT_DIR_INIT}dendro/i0_dendrogram.png'
    if ALWAYS_REPROCESS:
        rm(tvr_clust_params[0]['dist_in'])
    #
    init_clust_dat = cluster_tvrs(kmer_hit_dat, KMER_METADATA,  **tvr_clust_params[0])
    #
    init_clusters = []
    unclustered_inds = []
    for clust_inds in init_clust_dat[0]:
        if len(clust_inds) >= MIN_READS_PER_PHASE:
            init_clusters.append(clust_inds)
        else:
            unclustered_inds.extend(clust_inds)
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()
    print(f' - {len(init_clusters)} clusters ({sum([len(n) for n in init_clusters])} reads)')
    print(f' - {len(unclustered_inds)} unclustered reads')

    #=====================================================#
    #
    # [1b] TVR CLUSTER REFINEMENT
    #
    #=====================================================#

    def tvr_refinement(init_clusters, temp_dir, which_chr, message, print_stats=False):
        #
        out_unclustered_inds = []
        out_blank_inds = []
        out_clusters_with_tvrs = []
        out_all_consensus_tvrs = []
        #
        current_iteration_clusters = init_clusters
        for tvr_clust_iter in TVR_REFINE_ITERATIONS:
            sys.stdout.write(f'{message} [iteration {tvr_clust_iter}]...')
            sys.stdout.flush()
            if PRINT_PROGRESS:
                print()
            tt = time.perf_counter()
            next_iteration_clusters = []
            for clust_i, current_clust_inds in enumerate(current_iteration_clusters):
                khd_subset = [copy.deepcopy(kmer_hit_dat[n]) for n in current_clust_inds]
                zfcn = str(clust_i).zfill(3)
                tvr_clust_params[tvr_clust_iter]['fig_name'] = f'{temp_dir}dendro/i{tvr_clust_iter}_dendrogram_{zfcn}.png'
                tvr_clust_params[tvr_clust_iter]['dist_in']  = f'{temp_dir}npz/i{tvr_clust_iter}_dist-matrix_{zfcn}.npz'
                tvr_clust_params[tvr_clust_iter]['save_msa'] = f'{temp_dir}fa/i{tvr_clust_iter}_consensus_{zfcn}.fa'
                if tvr_clust_iter == TVR_REFINE_ITERATIONS[-1]:
                    tvr_clust_params[tvr_clust_iter]['clustering_only'] = False
                if ALWAYS_REPROCESS:
                    rm(tvr_clust_params[tvr_clust_iter]['dist_in'])
                    rm(tvr_clust_params[tvr_clust_iter]['save_msa'])
                #
                init_refine_clust_dat = cluster_tvrs(khd_subset, KMER_METADATA, **tvr_clust_params[tvr_clust_iter])
                #
                telcompplot_fn = f'{temp_dir}results/i{tvr_clust_iter}_reads_{zfcn}.png'
                telcompcons_fn = f'{temp_dir}results/i{tvr_clust_iter}_consensus_{zfcn}.png'
                if tvr_clust_iter in TVR_ITERATIONS_TO_PLOT:
                    if DONT_OVERWRITE_PLOTS is False or (exists_and_is_nonzero(telcompplot_fn) is False and exists_and_is_nonzero(telcompcons_fn) is False):
                        make_tvr_plots(khd_subset, init_refine_clust_dat, which_chr, UNCLUST_POS, telcompplot_fn, telcompcons_fn, mtp_params)
                #
                num_pass_subclust_reads = []
                for sub_clust in init_refine_clust_dat[0]:
                    next_clust = [current_clust_inds[n] for n in sub_clust]
                    if len(next_clust) >= MIN_READS_PER_PHASE:
                        next_iteration_clusters.append([n for n in next_clust])
                        num_pass_subclust_reads.append(len(next_clust))
                    else:
                        out_unclustered_inds.extend(next_clust)
                #
                if PRINT_PROGRESS:
                    print(f' - processed cluster {zfcn} ({len(current_clust_inds)} reads) --> {len(num_pass_subclust_reads)} clusters ({sum(num_pass_subclust_reads)} reads)')
                #
                # if we're on the last iteration output the clusters for subsequent analyses
                #
                if tvr_clust_iter == TVR_REFINE_ITERATIONS[-1]:
                    for sci,sub_clust in enumerate(init_refine_clust_dat[0]):
                        next_clust = [current_clust_inds[n] for n in sub_clust]
                        my_tvr_len = init_refine_clust_dat[7][sci]
                        if my_tvr_len <= 0:
                            out_blank_inds.extend(next_clust)
                            continue
                        if len(next_clust) < MIN_READS_PER_PHASE:
                            continue
                        tvrtel_lens = [init_refine_clust_dat[1][n] for n in sub_clust]
                        my_subtels = []
                        for tti,sri in enumerate(next_clust):
                            my_subtel_len = len(kmer_hit_dat[sri][6])-tvrtel_lens[tti]
                            my_subtel_seq_rev = kmer_hit_dat[sri][6][:my_subtel_len][::-1]
                            my_subtels.append(my_subtel_seq_rev)
                        out_clusters_with_tvrs.append((which_chr, [n for n in next_clust], my_subtels))
                        out_all_consensus_tvrs.append(init_refine_clust_dat[4][sci][:my_tvr_len])
            #
            current_iteration_clusters = next_iteration_clusters
            sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
            sys.stdout.flush()
            if print_stats:
                print(f' - {len(current_iteration_clusters)} clusters ({sum([len(n) for n in current_iteration_clusters])} reads)')
                print(f' - {len(out_blank_inds)} reads with blank TVR')
                print(f' - {len(out_unclustered_inds)} unclustered reads')
        #
        return [out_unclustered_inds,
                out_blank_inds,
                out_clusters_with_tvrs,
                out_all_consensus_tvrs]

    [new_unclustered_inds,
     blank_inds,
     clusters_with_tvrs,
     all_consensus_tvrs] = tvr_refinement(init_clusters, OUT_DIR_INIT, UNCLUST_CHR, 'initial clustering of all reads', print_stats=True)

    # combine unclustered reads from initial clustering with refinement iterations
    unclustered_inds += new_unclustered_inds

    #=====================================================#
    #
    # [2] SALVAGE UNCLUSTERED READS
    #
    #=====================================================#

    [salvaged_unclustered_inds,
     salvaged_blank_inds,
     salvaged_clusters_with_tvrs,
     salvaged_all_consensus_tvrs] = tvr_refinement([unclustered_inds], OUT_DIR_UNCLUST, UNCLUST_CHR, 'reanalyzing unclustered reads')
    print(f' - {len(salvaged_clusters_with_tvrs)} new candidate clusters found')

    # combine salvaged inds
    blank_inds += salvaged_blank_inds
    clusters_with_tvrs += salvaged_clusters_with_tvrs
    all_consensus_tvrs += salvaged_all_consensus_tvrs

    #=====================================================#
    #
    # [3] COLLAPSE CLUSTERS WITH SIMILAR TVRS
    #
    #=====================================================#

    sys.stdout.write(f'collapsing duplicates...')
    sys.stdout.flush()
    tt = time.perf_counter()
    collapse_clust_params['dendro_name'] = OUT_DIR_COLLAPSE + 'tvrs_dendrogram.png'
    collapse_clust_params['dist_in']     = OUT_DIR_COLLAPSE + 'tvrs_dist-matrix.npz'
    collapse_clust_params['fig_name']    = OUT_DIR_COLLAPSE + 'tvrs_plot.png'
    if ALWAYS_REPROCESS:
        rm(collapse_clust_params['dist_in'])
    #
    clustered_consensus_tvrs = cluster_consensus_tvrs(all_consensus_tvrs, KMER_METADATA, **collapse_clust_params)
    #
    temp_clusters_with_tvrs = copy.deepcopy(clusters_with_tvrs)
    clusters_with_tvrs = []
    for tvr_clust in clustered_consensus_tvrs:
        clusters_with_tvrs.append((UNCLUST_CHR, sum([temp_clusters_with_tvrs[n][1] for n in tvr_clust], []), sum([temp_clusters_with_tvrs[n][2] for n in tvr_clust], [])))
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()
    if len(temp_clusters_with_tvrs) == len(clusters_with_tvrs):
        print(' - no duplicates found')
    else:
        print(f' - {len(temp_clusters_with_tvrs)} clusters --> {len(clusters_with_tvrs)} clusters')
    del temp_clusters_with_tvrs

    # get subtels for reads with blank tvrs
    if len(blank_inds) >= MIN_READS_PER_PHASE:
        khd_subset = [copy.deepcopy(kmer_hit_dat[n]) for n in blank_inds]
        tvrtel_lens = quick_get_tvrtel_lens(khd_subset, KMER_METADATA)
        my_subtels = []
        for tti,sri in enumerate(blank_inds):
            my_subtel_len = len(kmer_hit_dat[sri][6])-tvrtel_lens[tti]
            my_subtel_seq_rev = kmer_hit_dat[sri][6][:my_subtel_len][::-1]
            my_subtels.append(my_subtel_seq_rev)
        clusters_with_tvrs.append((BLANK_CHR, [n for n in blank_inds], my_subtels))
        print('appending cluster of reads with blank TVR...')
        print(f' - {len(clusters_with_tvrs)} clusters')

    #=====================================================#
    #
    # [4] SUBTEL CLUSTER REFINEMENT
    #
    #=====================================================#

    if SKIP_SUBTEL_REFINE:
        print('skipping subtel cluster refinement...')
        final_clustered_read_inds = []
        for sci, (my_chr, subclust_read_inds, my_subtels) in enumerate(clusters_with_tvrs):
            final_clustered_read_inds.append((my_chr, subclust_read_inds, my_subtels))
    else:
        sys.stdout.write('refining clusters [SUBTEL]...')
        sys.stdout.flush()
        if PRINT_PROGRESS:
            print()
        tt = time.perf_counter()
        final_clustered_read_inds = []
        n_reads = 0
        for sci, (my_chr, subclust_read_inds, my_subtels) in enumerate(clusters_with_tvrs):
            zfcn = str(sci).zfill(3)
            subtel_clust_params['dendro_name'] = OUT_DIR_SUBTEL + 'dendro/' + 'dendrogram_'  + zfcn + '.png'
            subtel_clust_params['dist_in']     = OUT_DIR_SUBTEL + 'npz/'    + 'dist-matrix_' + zfcn + '.npz'
            subtel_clust_params['samp_labels'] = None
            subtel_clust_params['dendrogram_title'] = zfcn
            if ALWAYS_REPROCESS:
                rm(subtel_clust_params['dist_in'])
            #
            subtel_clustdat = cluster_consensus_tvrs(my_subtels, KMER_METADATA, **subtel_clust_params)
            #
            num_pass_subclust_reads = []
            for subclust_inds in subtel_clustdat:
                subsubclust_read_inds = [subclust_read_inds[n] for n in subclust_inds]
                if len(subsubclust_read_inds) >= MIN_READS_PER_PHASE:
                    final_clustered_read_inds.append((my_chr, [n for n in subsubclust_read_inds], [my_subtels[n] for n in subclust_inds]))
                    num_pass_subclust_reads.append(len(subsubclust_read_inds))
                    n_reads += len(subsubclust_read_inds)
            #
            if PRINT_PROGRESS:
                print(f' - processed cluster {zfcn} ({len(subclust_read_inds)} reads) --> {len(num_pass_subclust_reads)} clusters ({sum(num_pass_subclust_reads)} reads)')
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()
        print(f' - {len(final_clustered_read_inds)} clusters ({n_reads} reads)')

    #=====================================================#
    #
    # [5] FINAL RECLUSTERING
    #
    #=====================================================#

    sys.stdout.write('final reclustering to get TVRs and TLs for each allele...')
    sys.stdout.flush()
    if PRINT_PROGRESS:
        print()
    tt = time.perf_counter()

    subtels_out = []
    allele_outdat = []
    allele_consensus = []
    for plot_num, (my_chr, current_clust_inds, rev_subtels) in enumerate(final_clustered_read_inds):
        current_clust_inds = sorted(current_clust_inds)
        khd_subset = [copy.deepcopy(kmer_hit_dat[n]) for n in current_clust_inds]
        rlens_subset = [my_rlens[n] for n in current_clust_inds]
        #
        zfcn = str(plot_num).zfill(3)
        final_clust_params['fig_name'] = f'{OUT_DIR_FINAL}dendro/dendrogram_{zfcn}.png'
        final_clust_params['dist_in']  = f'{OUT_DIR_FINAL}npz/dist-matrix_{zfcn}.npz'
        final_clust_params['save_msa'] = f'{OUT_DIR_FINAL}fa/consensus_{zfcn}.fa'
        final_clust_params['my_chr']   = my_chr
        if ALWAYS_REPROCESS:
            rm(final_clust_params['dist_in'])
            rm(final_clust_params['save_msa'])
        #
        solo_clustdat = cluster_tvrs(khd_subset, KMER_METADATA, **final_clust_params)
        #
        my_tsvdat = get_allele_tsv_dat(khd_subset, solo_clustdat, UNCLUST_CHR, UNCLUST_POS, rlens_subset, gatd_params)
        for sci,subclust_inds in enumerate(solo_clustdat[0]):
            subclust_read_inds = [current_clust_inds[n] for n in subclust_inds]
            for i,sri in enumerate(subclust_read_inds):
                subtels_out.append((f'cluster-{len(allele_outdat)}_read-{i}_{my_rnames[sri]}', rev_subtels[subclust_inds[i]][::-1]))
            allele_outdat.append(my_tsvdat[sci])
            allele_consensus.append(solo_clustdat[4][sci])
        #
        telcompplot_fn = f'{OUT_DIR_FINAL}results/reads_{zfcn}.png'
        telcompcons_fn = f'{OUT_DIR_FINAL}results/consensus_{zfcn}.png'
        if PLOT_FINAL_CLUST:
            if DONT_OVERWRITE_PLOTS is False or (exists_and_is_nonzero(telcompplot_fn) is False and exists_and_is_nonzero(telcompcons_fn) is False):
                make_tvr_plots(khd_subset, solo_clustdat, my_chr, UNCLUST_POS, telcompplot_fn, telcompcons_fn, mtp_params)
        #
        if PRINT_PROGRESS:
            print(f' - processed cluster {zfcn} ({len(current_clust_inds)} reads)')
        #
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()

    # fill out allele id fields (there are no multimapped alleles at this point)
    for i,atdat in enumerate(allele_outdat):
        ALLELE_TEL_DAT.append([n for n in atdat])
        ALLELE_TEL_DAT[-1][3] = str(i)

    #=====================================================#
    #
    # [6] SUBTEL ANCHORING
    #
    #=====================================================#

    if SKIP_ANCHORING is False:
        aln_log = ALIGNED_SUBTELS + '.log'
        aln_sam = ALIGNED_SUBTELS + '.sam'
        aln_bam = ALIGNED_SUBTELS + '.bam'
        if ALWAYS_REPROCESS or exists_and_is_nonzero(aln_bam) is False:
            sys.stdout.write('aligning subtels to reference...')
            sys.stdout.flush()
            tt = time.perf_counter()
            #
            with gzip.open(UNANCHORED_SUBTELS, 'wt') as f:
                for n in subtels_out:
                    if len(n[1]) > 0:
                        f.write(f'>{n[0]}\n{n[1]}\n')
                for n in interstitial_subtels:
                    if len(n[1]) > 0:
                        f.write(f'>{n[0]}\n{n[1]}\n')
            #
            cmd = ''
            #
            if WHICH_ALIGNER == 'minimap2':
                aln_params = ''
                if READ_TYPE == 'ont':
                    aln_params = '-ax map-ont -Y'
                elif READ_TYPE == 'hifi':
                    aln_params = '-ax map-hifi -Y'
                cmd = ALIGNER_EXE + ' ' + aln_params + ' -o ' + aln_sam + ' ' + TELOGATOR_REF + ' ' + UNANCHORED_SUBTELS
            #
            elif WHICH_ALIGNER == 'winnowmap':
                aln_params = ''
                if READ_TYPE == 'ont':
                    aln_params = '-ax map-ont -Y'
                elif READ_TYPE == 'hifi':
                    aln_params = '-ax map-pb -Y'
                cmd = ALIGNER_EXE + ' -W ' + WINNOWMAP_K15 + ' ' + aln_params + ' -o ' + aln_sam + ' ' + TELOGATOR_REF + ' ' + UNANCHORED_SUBTELS
            #
            elif WHICH_ALIGNER == 'pbmm2':
                with gzip.open(TELOGATOR_REF,'rt') as f_tr: # pbmm2 wants an uncompressed reference
                    with open(UNCOMPRESSED_TELOGATOR_REF,'w') as f_utr:
                        for line in f_tr:
                            f_utr.write(line)
                pysam.faidx(UNCOMPRESSED_TELOGATOR_REF)
                cmd = ALIGNER_EXE + ' align ' + UNCOMPRESSED_TELOGATOR_REF + ' ' + UNANCHORED_SUBTELS + ' ' + aln_bam + ' --preset HiFi --sort'
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
        print('getting anchors from alignment...')
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
        interstial_anchors = {}
        for readname in ALIGNMENTS_BY_RNAME:
            #
            # parse interstitial telomere region anchors (possible fusions)
            #
            if readname[:12] == 'interstitial':
                original_readname = '_'.join(readname.split('_')[2:])
                original_readlen = int(readname.split('_')[1])
                if original_readname not in interstial_anchors:
                    interstial_anchors[original_readname] = [None, None]
                sorted_choices = []
                #print(readname)
                for i,aln in enumerate(ALIGNMENTS_BY_RNAME[readname]):
                    my_refbuild = aln[2]
                    if my_refbuild == '*':
                        continue
                    # check for known interstitial regions
                    my_refpos = None
                    l_vs_r = None
                    if readname[:17] == 'interstitial-left':
                        l_vs_r = 0
                        my_refpos = aln[4]
                    elif readname[:18] == 'interstitial-right':
                        l_vs_r = 1
                        my_refpos = aln[3]
                    #print(aln[:7], annotate_interstitial_tel(my_refbuild, my_refpos))
                    if annotate_interstitial_tel(my_refbuild, my_refpos) is False:
                        read_span = abs(aln[1] - aln[0])
                        my_mapq = aln[6]
                        sorted_choices.append((my_mapq, read_span, (len(aln[7]), original_readlen), l_vs_r,
                                               my_refbuild, my_refpos, copy.deepcopy(aln[:7]))) # sorting first by mapq, then by span
                if sorted_choices:
                    sorted_choices = sorted(sorted_choices, reverse=True)
                    top_choice = sorted_choices[0]
                    interstial_anchors[original_readname][top_choice[3]] = top_choice
            #
            # parse subtel anchors
            #
            elif readname[:7] == 'cluster':
                original_readname = '_'.join(readname.split('_')[2:])
                sorted_choices = []
                for i,aln in enumerate(ALIGNMENTS_BY_RNAME[readname]):
                    dist_to_end = len(aln[7]) - max(aln[0], aln[1])
                    read_span = abs(aln[1] - aln[0])
                    my_mapq = aln[6]
                    sorted_choices.append((read_span, my_mapq, -dist_to_end, i)) # sorting first by span, then by mapq
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
        # annotate possible tel fusions
        #
        tel_fusion_plot_num = 0
        for readname,anchor_dat in interstial_anchors.items():
            if anchor_dat[0] is not None and anchor_dat[1] is not None:
                # require both anchors have mapq > 0
                if anchor_dat[0][0] > 0 and anchor_dat[1][0] > 0:
                    #print(readname, readname in interstitial_khd_by_readname)
                    #print(anchor_dat[0])
                    #print(anchor_dat[1])
                    if PLOT_TEL_FUSIONS:
                        plot_fn = TEL_FUSION_DIR + 'fusion_' + str(tel_fusion_plot_num).zfill(5) + '.png'
                        khd_tuple = (interstitial_khd_by_readname[readname], interstitial_khd_rev_by_readname[readname])
                        plot_fusion(readname, anchor_dat[0], anchor_dat[1], khd_tuple, KMER_METADATA, plot_fn)
                    tel_fusion_plot_num += 1
        #
        # subtel assignment
        #
        for clustnum in sorted(top_alns_by_cluster.keys()):
            chr_arm_scores = {}
            anchors_by_ref = {}
            readnames_by_ref = {}
            for (readname,aln) in top_alns_by_cluster[clustnum]:
                my_ref = aln[2]
                my_refspan = abs(aln[3] - aln[4])
                my_anchor = aln[4] # second refpos is always anchor coord?
                #my_orr = aln[5]
                my_mapq = aln[6]
                if my_ref != '*':
                    if my_ref not in chr_arm_scores:
                        chr_arm_scores[my_ref] = 0.
                        anchors_by_ref[my_ref] = []
                        readnames_by_ref[my_ref] = []
                    chr_arm_scores[my_ref] += my_refspan * ((my_mapq + 1) / (MAX_QUAL_SCORE + 1))
                    anchors_by_ref[my_ref].append(my_anchor)
                    readnames_by_ref[my_ref].append(readname)
            total_score = sum(chr_arm_scores.values())
            #
            my_anchors = []
            if total_score > 0:
                sorted_chr_scores = sorted([(chr_arm_scores[k]/total_score, k) for k in chr_arm_scores], reverse=True)
                my_anchors = [n[1] for n in sorted_chr_scores if n[0] >= ANCHORING_ASSIGNMENT_FRAC]
            else:
                pass # all unmapped
            #
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

        # resort out allele dat by chr & atl & pos
        sorted_ad_inds = sorted([(LEXICO_2_IND[n[0].split(',')[0][:-1]], n[0].split(',')[0][-1], -1 * int(n[4]), int(n[1].split(',')[0]), i) for i,n in enumerate(ALLELE_TEL_DAT)])
        ALLELE_TEL_DAT = [ALLELE_TEL_DAT[n[4]] for n in sorted_ad_inds]

    #=====================================================#
    #
    # [7] COLLAPSE SIMILAR ALLELES MAPPED NEAR EACH OTHER
    #
    #=====================================================#

    if COLLAPSE_HOM_BP > 0:
        sys.stdout.write(f'merging alleles mapped within {COLLAPSE_HOM_BP}bp of each other...')
        sys.stdout.flush()
        tt = time.perf_counter()
        n_alleles_before_collapsing = len(ALLELE_TEL_DAT)
        plot_num = 0
        while True:
            del_ind = None
            for i in range(len(ALLELE_TEL_DAT)):
                for j in range(i+1, len(ALLELE_TEL_DAT)):
                    if all([ALLELE_TEL_DAT[i][0].split(',')[0] == ALLELE_TEL_DAT[j][0].split(',')[0],
                            ALLELE_TEL_DAT[i][2].split(',')[0] == ALLELE_TEL_DAT[j][2].split(',')[0],
                            abs(int(ALLELE_TEL_DAT[i][1].split(',')[0]) - int(ALLELE_TEL_DAT[j][1].split(',')[0])) <= COLLAPSE_HOM_BP]):
                        my_dist = TREECUT_HOM_COLLAPSE * 2
                        tvr_i, tvr_j = ALLELE_TEL_DAT[i][9], ALLELE_TEL_DAT[j][9]
                        if len(tvr_i) and len(tvr_j):
                            my_dist = quick_compare_tvrs(tvr_i, tvr_j)
                        if DEBUG_COLLAPSE:
                            print(tvr_i)
                            print(tvr_j)
                            print()
                            print(f'DEBUG: ({ALLELE_TEL_DAT[i][3]}) vs. ({ALLELE_TEL_DAT[j][3]}): {my_dist} <= {TREECUT_HOM_COLLAPSE} : {my_dist <= TREECUT_HOM_COLLAPSE}')
                            print()
                        if my_dist <= TREECUT_HOM_COLLAPSE:
                            which_ind, merged_dat = merge_allele_tsv_dat(ALLELE_TEL_DAT[i], ALLELE_TEL_DAT[j], ALLELE_TL_METHOD)
                            if which_ind == 0:
                                ALLELE_TEL_DAT[i] = merged_dat
                                del_ind = j
                            else:
                                ALLELE_TEL_DAT[j] = merged_dat
                                del_ind = i
                            plot_fn_col2 = OUT_DIR_COL2 + f'merge_{str(plot_num).zfill(3)}.png'
                            if exists_and_is_nonzero(plot_fn_col2) is False or DONT_OVERWRITE_PLOTS is False:
                                plot_some_tvrs([tvr_i, tvr_j],
                                               [' '.join(ALLELE_TEL_DAT[i][:3]), ' '.join(ALLELE_TEL_DAT[j][:3])],
                                               KMER_METADATA,
                                               plot_fn_col2,
                                               custom_plot_params={'custom_title':f'tvr dist = {my_dist:.3f}'})
                            plot_num += 1
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

    #=====================================================#
    #
    # [8] WRITE OUT THE FINAL RESULTS
    #
    #=====================================================#

    # change chrB to chrU if any blank clusters are unanchored, to simplify output
    for i in range(len(ALLELE_TEL_DAT)):
        ALLELE_TEL_DAT[i][0] = ALLELE_TEL_DAT[i][0].replace(BLANK_CHR, UNCLUST_CHR)

    print('writing allele TL results to "' + FINAL_TSV.split('/')[-1] + '"...')
    with open(FINAL_TSV, 'w') as f:
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
            if my_chr in [UNCLUST_CHR, BLANK_CHR]:
                my_chr = 'unanchored'
            if 'i' in atd[3]:
                continue # don't include interstitial alleles in violin
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
    readcount_fail_final_filters = [0, 0, 0]
    for atd in ALLELE_TEL_DAT:
        my_id = atd[3].split(',')[0]
        while my_id[-1].isdigit() is False:
            my_id = my_id[:-1]
        my_id = int(my_id)
        my_rep_atl = int(atd[4])
        #my_max_atl = max([int(n) for n in atd[5].split(',')])
        my_tvr_len = int(atd[8])
        n_reads = len(atd[6].split(','))
        if atd[0] in [UNCLUST_CHR, BLANK_CHR]:
            num_alleles_unmapped += 1
            readcount_fail_final_filters[0] += n_reads
            continue
        if my_tvr_len + my_rep_atl < MIN_ATL_FOR_FINAL_PLOTTING:
            num_alleles_too_short += 1
            readcount_fail_final_filters[1] += n_reads
            continue
        if 'i' in atd[3]:
            num_alleles_interstitial += 1
            readcount_fail_final_filters[2] += n_reads
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
    np.savez_compressed(READLEN_NPZ, readlen_all=np.array(readlens_all,dtype='<i8'), readlen_tel=np.array(readlens_tel,dtype='<i8'), readlen_tel_final=np.array(readlens_final,dtype='<i8'))
    readlen_plot(readlens_all, readlens_tel, readlens_final, QC_READLEN)
    del readlens_all
    del readlens_tel
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

    #
    # stats
    #
    print(f' - {num_pass_alleles} final telomere alleles ({len(readlens_final)} reads)')
    print(f' - {num_alleles_unmapped} alleles unmapped ({readcount_fail_final_filters[0]} reads)')
    print(f' - {num_alleles_too_short} alleles with atl < {MIN_ATL_FOR_FINAL_PLOTTING} bp ({readcount_fail_final_filters[1]} reads)')
    print(f' - {num_alleles_interstitial} interstitial alelles ({readcount_fail_final_filters[2]} reads)')
    with open(QC_STATS, 'w') as f:
        f.write(f' - {num_pass_alleles} final telomere alleles ({len(readlens_final)} reads)\n')
        f.write(f' - {num_alleles_unmapped} alleles unmapped ({readcount_fail_final_filters[0]} reads)\n')
        f.write(f' - {num_alleles_too_short} alleles with atl < {MIN_ATL_FOR_FINAL_PLOTTING} bp ({readcount_fail_final_filters[1]} reads)\n')
        f.write(f' - {num_alleles_interstitial} interstitial alelles ({readcount_fail_final_filters[2]} reads)\n')


if __name__ == '__main__':
    main()
