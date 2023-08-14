#!/usr/bin/env python
import argparse
import copy
import multiprocessing
import numpy as np
import pathlib
import pickle
import pysam
import sys
import time

from source.tg_kmer   import get_telomere_composition, read_kmer_tsv
from source.tg_muscle import check_muscle_version
from source.tg_plot   import tel_len_violin_plot
from source.tg_tel    import choose_tl_from_observations, get_allele_tsv_dat, get_double_anchored_tels, get_tels_below_canonical_thresh, parallel_anchored_tel_job, parallel_filtering_job
from source.tg_tvr    import cluster_consensus_tvrs, cluster_tvrs, make_tvr_plots
from source.tg_util   import cluster_list, exists_and_is_nonzero, LEXICO_2_IND, makedir, parse_read, RC

# hardcoded parameters
TEL_WINDOW_SIZE = 100
# toss reads if too much of the non-telomere sequence couldn't be aligned anywhere
MAXIMUM_UNEXPLAINED_FRAC = 0.7
# toss reads if the median |p_vs_q| of the non-telomere sequence is above this value
MAX_NONTEL_MEDIAN_KMER_DENSITY = 0.25
#
MIN_DOUBLE_ANCHOR_LEN   = 1000
MIN_DOUBLE_ANCHOR_READS = 3
#
MAX_NOISE_BASES = 5000
MAX_NOISE_FRAC  = 0.60
#
DUMMY_TEL_MAPQ = 60
#
ORR_SWAP = {'FWD':'REV', 'REV':'FWD', 'p':'q', 'q':'p'}
# how much subtel should we use for de-clustering alleles? [min_size, max_size]
SUBTEL_CLUSTER_SIZE = [1500, 3000]

#
# ANCHORING_STRATEGY = 'largest'     - anchor tels onto largest non-tel alignment
# ANCHORING_STRATEGY = 'closest'     - anchor tels onto nearest non-tel alignment
#
# MATCH_TRIM_STRATEGY = 'largest'    - prioritize largest alignments when trimming overlaps
# MATCH_TRIM_STRATEGY = 'mapq'       - prioritize alignments with highest MAPQ when trimming overlaps
# MATCH_TRIM_STRATEGY = 'none'       - do not trim overlaps (use at your own risk for mm2 BAMs)
#


def main(raw_args=None):
    parser = argparse.ArgumentParser(description='Telogator2', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    parser.add_argument('-i',  type=str,   required=True,  metavar='input.bam',                      help="* Long reads aligned to t2t-telogator-ref.fa")
    parser.add_argument('-o',  type=str,   required=True,  metavar='output/',                        help="* Path to output directory")
    parser.add_argument('-r',  type=str,   required=False, metavar='hifi',        default='hifi',    help="Read type: hifi / clr / ont")
    parser.add_argument('-k',  type=str,   required=False, metavar='kmers.tsv',   default='',        help="Telomere kmers to search for")
    parser.add_argument('-t',  type=float, required=False, metavar='0.5',         default=0.5,       help="Telomere signal threshold (0-1)")
    #
    parser.add_argument('-fl', type=int,   required=False, metavar='5000',        default=5000,      help="Initial read filtering: Min read length")
    parser.add_argument('-ft', type=int,   required=False, metavar='200',         default=200,       help="Initial read filtering: Min tel bases")
    parser.add_argument('-ff', type=float, required=False, metavar='0.9',         default=0.9,       help="Anchored-tel filtering: Max frac of read that can be tel")
    parser.add_argument('-fp', type=float, required=False, metavar='0.25',        default=0.25,      help="Anchored-tel filtering: Max minor p/q fraction in tels")
    parser.add_argument('-fc', type=float, required=False, metavar='0.15',        default=0.15,      help="Anchored-tel filtering: Minimum canonical frac in tels")
    #
    parser.add_argument('-sa', type=str,   required=False, metavar='largest',     default='largest', help="Subtel/tel anchoring strategy")
    parser.add_argument('-sm', type=str,   required=False, metavar='mapq',        default='mapq',    help="Repeated matches trimming strategy")
    #
    parser.add_argument('-cd', type=int,   required=False, metavar='2000',        default=2000,      help="Maximum subtel distance to cluster anchored tels together")
    parser.add_argument('-cr', type=int,   required=False, metavar='3',           default=3,         help="Minimum number of reads per anchored-tel cluster")
    parser.add_argument('-cm', type=str,   required=False, metavar='p90',         default='p90',     help="Method for computing chr TL: mean / median / max / p90")
    #
    parser.add_argument('-vtm', type=int,  required=False, metavar='20000',       default=20000,     help="Violin plot (tel length) max value")
    parser.add_argument('-vtt', type=int,  required=False, metavar='5000',        default=5000,      help="Violin plot (tel length) tick size")
    parser.add_argument('-vrm', type=int,  required=False, metavar='50000',       default=50000,     help="Violin plot (read length) max value")
    parser.add_argument('-vrt', type=int,  required=False, metavar='10000',       default=10000,     help="Violin plot (read length) tick size")
    #
    parser.add_argument('-m',   type=str,  required=False, metavar='muscle',      default='muscle',  help="/path/to/muscle executable")
    parser.add_argument('-ar',  type=int,  required=False, metavar='3',           default=3,         help="Minimum number of reads per phased tel allele")
    parser.add_argument('-am',  type=str,  required=False, metavar='max',         default='max',     help="Method for computing chr TL: mean / median / max / p90")
    parser.add_argument('-att', type=float,required=False, metavar='0.180',       default=0.180,     help="Treecut value: main TVR clustering")
    parser.add_argument('-atp', type=float,required=False, metavar='0.080',       default=0.080,     help="Treecut value: merging TVRs with similar prefixes")
    parser.add_argument('-atm', type=float,required=False, metavar='0.070',       default=0.070,     help="Treecut value: grouping multimapped alleles")
    parser.add_argument('-ato', type=float,required=False, metavar='0.080',       default=0.080,     help="Treecut value: clustering orphaned alleles")
    parser.add_argument('-ats', type=float,required=False, metavar='0.250',       default=0.250,     help="Treecut value: declustering alleles using subtels")
    parser.add_argument('-atc', type=str,  required=False, metavar='treecuts.tsv',default='',        help="Custom TVR treecut vals for specific anchors")
    #
    parser.add_argument('--plot',          required=False, action='store_true',   default=False,     help="Create read plots")
    parser.add_argument('--plot-filt',     required=False, action='store_true',   default=False,     help="Create read plots (filtered reads)")
    parser.add_argument('--plot-filt-tvr', required=False, action='store_true',   default=False,     help="Plot denoised TVR instead of raw signal")
    #
    parser.add_argument('--debug',               required=False, action='store_true', default=False, help="[DEBUG] Print extra info for each read as its processed")
    parser.add_argument('--debug-npy',           required=False, action='store_true', default=False, help="[DEBUG] Save .npy files and use existing .npy files")
    parser.add_argument('--debug-overwriteplot', required=False, action='store_true', default=False, help="[DEBUG] Do not regenerate plots that already exist")
    parser.add_argument('--debug-chr', type=str, required=False, metavar='chr_list',  default='',    help="[DEBUG] Only process: chr1p,chr1q,... (comma-delimited)")
    #
    parser.add_argument('--no-anchorfilt', required=False, action='store_true',   default=False,     help="Skip double-anchored read filtering")
    parser.add_argument('--no-orphans',    required=False, action='store_true',   default=False,     help="Skip orphan read clustering")
    parser.add_argument('--no-subtelclust',required=False, action='store_true',   default=False,     help="Skip subtelomere-based cluster splitting")
    parser.add_argument('--fast',          required=False, action='store_true',   default=False,     help="Use faster but less accurate pairwise alignment")
    #
    parser.add_argument('-p',   type=int,  required=False, metavar='4',           default=4,         help="Number of processes to use")

    args = parser.parse_args()

    #
    # initiate input / output files
    #
    INPUT_ALN  = args.i
    if exists_and_is_nonzero(INPUT_ALN) is False:
        print('Error reading -i input file')
        exit(1)
    if INPUT_ALN[-4:].lower() == '.sam':
        INPUT_TYPE = 'sam'
    elif INPUT_ALN[-4:].lower() == '.bam':
        INPUT_TYPE = 'bam'
    elif INPUT_ALN[-5:].lower() == '.cram':
        INPUT_TYPE = 'cram'
    elif INPUT_ALN[-2:].lower() == '.p':
        INPUT_TYPE = 'pickle'
    else:
        print('Error: unknown -i file type')
        exit(1)
    #
    OUT_DIR = args.o
    if OUT_DIR[-1] != '/':
        OUT_DIR += '/'
    makedir(OUT_DIR)
    OUT_PICKLE    = OUT_DIR + 'tel-data.p'
    OUT_PLOT_DIR  = OUT_DIR + 'read_plots/'
    OUT_CHR_TL    = OUT_DIR + 'tlens_by_chr.tsv'
    OUT_ALLELE_TL = OUT_DIR + 'tlens_by_allele.tsv'
    OUT_VIOLIN_TL = OUT_DIR + 'violin_tlens.png'
    OUT_VIOLIN_RL = OUT_DIR + 'violin_rlens.png'
    OUT_FILTSTATS = OUT_DIR + 'filtered_read_stats.tsv'
    #
    OUT_TVR_DIR  = OUT_DIR + 'tvr_clustering/'
    OUT_TVR_TEMP = OUT_TVR_DIR + 'temp/'
    makedir(OUT_TVR_DIR)
    makedir(OUT_TVR_TEMP)
    #
    READ_TYPE = args.r
    #
    PLOT_READS      = args.plot
    PLOT_FILT_READS = args.plot_filt
    PLOT_FILT_CVECS = args.plot_filt_tvr
    PRINT_DEBUG     = args.debug
    FAST_ALIGNMENT  = args.fast
    #
    DS_OR_MS = 'ds'
    RAND_SHUFFLE = 3
    if FAST_ALIGNMENT:
        DS_OR_MS = 'ms'
        RAND_SHUFFLE = 1
    #
    MUSCLE_EXE = args.m
    #
    if PLOT_READS:
        makedir(OUT_PLOT_DIR)
    #
    # for debugging purposes
    DONT_OVERWRITE_PLOTS = args.debug_overwriteplot   # (True = don't replot figures if they already exist)
    ALWAYS_REPROCESS = not(args.debug_npy)            # (True = don't write out .npy matrices, always recompute)
    # for debugging purposes (only cluster TVRs at specific arms)
    DEBUG_CHR_LIST = []
    if len(args.debug_chr):
        DEBUG_CHR_LIST = [n for n in args.debug_chr.split(',')]
    # merge all blank TVRs into "chrBq" and let them be clustered via subtel sequences instead
    MERGE_BLANK = True
    #
    orphans_chr = 'chrUq'
    orphans_pos = 0
    blank_chr   = 'chrBq'
    blank_pos   = 0
    fake_chr    = 'chrUq'
    fake_pos    = 0

    #
    # parse telomere kmer data
    #
    KMER_FILE = args.k
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

    #
    # parse custom treecut values, if provided
    #
    TREECUT_TSV = args.atc
    custom_treecut_vals = {}
    if len(TREECUT_TSV):
        f = open(TREECUT_TSV, 'r')
        for line in f:
            if len(line.strip()):
                splt = line.strip().split('\t')
                if ':' in splt[0]:
                    splt2  = splt[0].split(':')
                    tc_chr = splt2[0]
                    tc_pos = int(splt2[1])
                else:
                    tc_chr = splt[0].replace('_','')
                    tc_pos = None
                custom_treecut_vals[(tc_chr, tc_pos)] = float(splt[1])
        f.close()

    #
    # various parameters
    #
    MINIMUM_READ_LEN   = args.fl
    MINIMUM_TEL_BASES  = args.ft
    MAXIMUM_TEL_FRAC   = args.ff
    MAXIMUM_MINOR_PQ   = args.fp
    MIN_CANONICAL_FRAC = args.fc
    P_VS_Q_AMP_THRESH  = args.t
    #
    ANCHORING_STRATEGY  = args.sa
    MATCH_TRIM_STRATEGY = args.sm
    #
    ANCHOR_CLUSTER_DIST = args.cd
    MIN_READS_PER_CLUST = args.cr
    CHR_TL_METHOD       = args.cm
    #
    MIN_READS_PER_PHASE = args.ar
    ALLELE_TL_METHOD    = args.am
    TREECUT_DEFAULT     = args.att
    TREECUT_PREFIXMERGE = args.atp
    TREECUT_MULTIMAPPED = args.atm
    TREECUT_ORPHANS     = args.ato
    TREECUT_SUBTELS     = args.ats
    #
    (VIOLIN_TLEN_MAX, VIOLIN_TLEN_TICK) = (args.vtm, args.vtt)
    (VIOLIN_RLEN_MAX, VIOLIN_RLEN_TICK) = (args.vrm, args.vrt)
    #
    SKIP_FILT_DOUBLEANCHOR = args.no_anchorfilt
    SKIP_ORPHAN_CLUSTERING = args.no_orphans
    SKIP_SUBTEL_CLUSTERING = args.no_subtelclust
    #
    NUM_PROCESSES = args.p

    reads_skipped = {'trim_filter':0,
                     'min_readlen':0,
                     'unmapped':0,
                     'unknown_ref':0,
                     'no_chr_aln':0,
                     'min_telbases':0}

    ANCHORED_TEL_BY_CHR    = {}
    NONTEL_REFSPANS_BY_CHR = {}

    #
    # GET ALL READ DATA FROM INPUT ALIGNMENT
    #
    ALIGNMENTS_BY_RNAME = {}
    #
    if INPUT_TYPE == 'pickle':
        print('reading pickle input...')
        f = open(INPUT_ALN, 'rb')
        my_pickle = pickle.load(f)
        f.close()
        ANCHORED_TEL_BY_CHR = my_pickle['anchored-tels']
        NONTEL_REFSPANS_BY_CHR = my_pickle['non-tel-ref-spans']
        num_input_reads_for_filter_tsv = None
    else:
        if INPUT_TYPE == 'sam':
            samfile = pysam.AlignmentFile(INPUT_ALN, "r")
        elif INPUT_TYPE == 'bam':
            samfile = pysam.AlignmentFile(INPUT_ALN, "rb")
        elif INPUT_TYPE == 'cram':
            samfile = pysam.AlignmentFile(INPUT_ALN, "rc")
        #
        refseqs = samfile.references
        sys.stdout.write('getting all read data from input alignment...')
        sys.stdout.flush()
        tt = time.perf_counter()
        for aln in samfile.fetch(until_eof=True):
            sam_line    = str(aln).split('\t')
            #
            # pysam weirdness:
            #
            my_ref_ind  = sam_line[2].replace('#','')
            if my_ref_ind.isdigit():
                sam_line[2] = refseqs[int(my_ref_ind)]
            elif my_ref_ind == '-1':
                sam_line[2] = refseqs[-1]
            else:
                sam_line[2] = my_ref_ind
            #
            [rnm, ref_key, pos, read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat] = parse_read(sam_line)
            if rnm not in ALIGNMENTS_BY_RNAME:
                ALIGNMENTS_BY_RNAME[rnm] = []
            ALIGNMENTS_BY_RNAME[rnm].append([read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat])
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()
        samfile.close()
        #
        if len(ALIGNMENTS_BY_RNAME) == 0:
            print('Error: unable to get any reads from input ' + INPUT_TYPE)
            exit(1)

        #
        #
        # INITIAL FILTERING
        #
        #
        par_results = None
        if NUM_PROCESSES > 1:
            all_read_keys    = [[] for n in range(NUM_PROCESSES)]
            reads_to_process = [{} for n in range(NUM_PROCESSES)]
            k = 0
            for read_key in ALIGNMENTS_BY_RNAME.keys():
                all_read_keys[k % NUM_PROCESSES].append(read_key)
                reads_to_process[k % NUM_PROCESSES][read_key] = copy.deepcopy(ALIGNMENTS_BY_RNAME[read_key])
                k += 1
        #
        # execute parallel jobs
        #
        sys.stdout.write('initial read filtering...')
        sys.stdout.flush()
        tt = time.perf_counter()
        num_starting_reads = len(ALIGNMENTS_BY_RNAME.keys())
        num_input_reads_for_filter_tsv = num_starting_reads
        #
        par_params      = [CANONICAL_STRINGS, CANONICAL_STRINGS_REV, READ_TYPE, MATCH_TRIM_STRATEGY, INPUT_TYPE, PRINT_DEBUG]
        par_params_filt = [MINIMUM_READ_LEN, MINIMUM_TEL_BASES]
        if NUM_PROCESSES > 1:
            manager     = multiprocessing.Manager()
            par_results = manager.dict()
            processes   = []
            for i in range(NUM_PROCESSES):
                p = multiprocessing.Process(target=parallel_filtering_job, args=(reads_to_process[i], i, par_params, par_params_filt, par_results))
                processes.append(p)
            for i in range(NUM_PROCESSES):
                processes[i].start()
            for i in range(NUM_PROCESSES):
                processes[i].join()
        else:
            par_results = {}
            parallel_filtering_job(ALIGNMENTS_BY_RNAME, 0, par_params, par_params_filt, par_results)
        #
        # collect parallel job results
        #
        FILTERED_READS = []
        for i in range(NUM_PROCESSES):
            job_filtered_reads = par_results[(i,0)]
            job_filt_strings   = par_results[(i,1)]
            job_nontel_spans   = par_results[(i,2)]
            #
            FILTERED_READS.extend(job_filtered_reads)
            for filt_string,count in job_filt_strings.items():
                if filt_string not in reads_skipped:
                    reads_skipped[filt_string] = 0
                reads_skipped[filt_string] += count
            for k in job_nontel_spans.keys():
                if k not in NONTEL_REFSPANS_BY_CHR:
                    NONTEL_REFSPANS_BY_CHR[k] = []
                NONTEL_REFSPANS_BY_CHR[k].extend(job_nontel_spans[k])
        #
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        num_ending_reads = len(FILTERED_READS)
        sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads\n')
        sys.stdout.flush()

        #
        #
        # TELOGATOR 1.0 - IDENTIFY ANCHORED TELOMERES
        #
        #
        par_results = None
        if NUM_PROCESSES > 1:
            all_indices      = [[] for n in range(NUM_PROCESSES)]
            reads_to_process = [[] for n in range(NUM_PROCESSES)]
            k = 0
            for i in range(len(FILTERED_READS)):
                all_indices[k % NUM_PROCESSES].append(i)
                reads_to_process[k % NUM_PROCESSES].append(copy.deepcopy(FILTERED_READS[i]))
                k += 1
        #
        # execute parallel jobs
        #
        sys.stdout.write('identify anchorable telomeres...')
        sys.stdout.flush()
        num_starting_reads = len(FILTERED_READS)
        num_ending_reads   = 0
        tt = time.perf_counter()
        #
        par_params      = [KMER_LIST, KMER_LIST_REV, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH, ANCHORING_STRATEGY, PLOT_READS, INPUT_TYPE, OUT_PLOT_DIR, PRINT_DEBUG, PLOT_FILT_READS]
        par_params_filt = [MAXIMUM_TEL_FRAC, MAXIMUM_MINOR_PQ, MAXIMUM_UNEXPLAINED_FRAC, MAX_NONTEL_MEDIAN_KMER_DENSITY, READ_TYPE]
        if NUM_PROCESSES > 1:
            manager     = multiprocessing.Manager()
            par_results = manager.dict()
            processes   = []
            for i in range(NUM_PROCESSES):
                p = multiprocessing.Process(target=parallel_anchored_tel_job, args=(reads_to_process[i], all_indices[i], par_params, par_params_filt, par_results))
                processes.append(p)
            for i in range(NUM_PROCESSES):
                processes[i].start()
            for i in range(NUM_PROCESSES):
                processes[i].join()
        else:
            par_results = {}
            parallel_anchored_tel_job(FILTERED_READS, list(range(len(FILTERED_READS))), par_params, par_params_filt, par_results)
        #
        # collect parallel job results
        #
        for i in range(num_starting_reads):
            anchored_tel_dat = par_results[(i,0)]
            my_filt_string   = par_results[(i,1)]
            nontel_spans     = par_results[(i,2)]
            #
            if len(my_filt_string):
                if my_filt_string not in reads_skipped:
                    reads_skipped[my_filt_string] = 0
                reads_skipped[my_filt_string] += 1
            else:
                adjacent_chr = anchored_tel_dat[0]
                if adjacent_chr not in ANCHORED_TEL_BY_CHR:
                    ANCHORED_TEL_BY_CHR[adjacent_chr] = []
                ANCHORED_TEL_BY_CHR[adjacent_chr].append(anchored_tel_dat[1:])
                num_ending_reads += 1
            #
            for k in nontel_spans.keys():
                if k not in NONTEL_REFSPANS_BY_CHR:
                    NONTEL_REFSPANS_BY_CHR[k] = []
                NONTEL_REFSPANS_BY_CHR[k].extend(nontel_spans[k])
        #
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads\n')
        sys.stdout.flush()

        #
        # sort non-tel ref spans
        #
        for k in NONTEL_REFSPANS_BY_CHR.keys():
            NONTEL_REFSPANS_BY_CHR[k] = sorted(NONTEL_REFSPANS_BY_CHR[k])
        #
        # save output pickle
        #
        f = open(OUT_PICKLE, 'wb')
        pickle.dump({'anchored-tels':ANCHORED_TEL_BY_CHR, 'non-tel-ref-spans':NONTEL_REFSPANS_BY_CHR}, f)
        f.close()

    #
    # remove telomeres that we thought were single-anchored but other reads tell us it is actually double-anchored
    #
    if not SKIP_FILT_DOUBLEANCHOR:
        sys.stdout.write('removing anchored tels that double-anchored according to other reads...')
        sys.stdout.flush()
        num_starting_reads = sum([len(ANCHORED_TEL_BY_CHR[k]) for k in ANCHORED_TEL_BY_CHR.keys()])
        tt = time.perf_counter()
        #
        gdat_params = [MIN_DOUBLE_ANCHOR_LEN, MIN_DOUBLE_ANCHOR_READS, PRINT_DEBUG]
        del_keys    = get_double_anchored_tels(ANCHORED_TEL_BY_CHR, NONTEL_REFSPANS_BY_CHR, gdat_params)
        #
        num_ending_reads = num_starting_reads - len(del_keys)
        reads_skipped['double_anchored'] = len(del_keys)
        del_keys2 = []
        for (k,di) in del_keys:
            del ANCHORED_TEL_BY_CHR[k][di]
            if len(ANCHORED_TEL_BY_CHR[k]) == 0:
                del_keys2.append(k)
        for k in del_keys2:
            del ANCHORED_TEL_BY_CHR[k]
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads\n')
        sys.stdout.flush()

    #
    # remove reads where telomere region is insufficiently canonical
    #
    sys.stdout.write('removing reads with tel regions that are insufficiently canonical...')
    sys.stdout.flush()
    num_starting_reads = sum([len(ANCHORED_TEL_BY_CHR[k]) for k in ANCHORED_TEL_BY_CHR.keys()])
    tt = time.perf_counter()
    #
    gtbct_params = [MINIMUM_TEL_BASES, MIN_CANONICAL_FRAC, CANONICAL_STRINGS, CANONICAL_STRINGS_REV, READ_TYPE]
    del_keys     = get_tels_below_canonical_thresh(ANCHORED_TEL_BY_CHR, gtbct_params)
    #
    num_ending_reads = num_starting_reads - len(del_keys)
    reads_skipped['min_canonical'] = len(del_keys)
    del_keys2 = []
    for (k,di) in del_keys:
        del ANCHORED_TEL_BY_CHR[k][di]
        if len(ANCHORED_TEL_BY_CHR[k]) == 0:
            del_keys2.append(k)
    for k in del_keys2:
        del ANCHORED_TEL_BY_CHR[k]
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads\n')
    sys.stdout.flush()

    #
    # print stats on reads that were filtered
    #
    with open(OUT_FILTSTATS, 'w') as f:
        if num_input_reads_for_filter_tsv is not None:
            f.write(str(num_input_reads_for_filter_tsv) + '\t' + 'total_input_reads' + '\n')
        f.write(str(num_ending_reads) + '\t' + 'total_after_filtering' + '\n')
        if sum(reads_skipped.values()) > 0:
            print()
            print('reads filtered:')
            sk = sorted(reads_skipped.keys())
            max_str_len = 0
            for k in sk:
                max_str_len = max(len(str(reads_skipped[k])), max_str_len)
            for k in sk:
                if reads_skipped[k] > 0:
                    print(' - ' + ' '*(max_str_len - len(str(reads_skipped[k]))) + str(reads_skipped[k]) + ' ' + k)
                    f.write(str(reads_skipped[k]) + '\t' + k + '\n')
            print()

    #
    # cluster subtel/tel boundaries by position and compute TL statistics
    #
    sorted_ref_keys = []
    for k in sorted(ANCHORED_TEL_BY_CHR.keys()):
        if k[:3] == 'chr':
            my_chr = k.replace('_','')
            my_ind = 1
        elif k[:3] == 'alt':
            my_chr = ''.join(k.split('_')[:-1])
            my_chr = my_chr.replace('alt', 'chr')
            my_ind = 2
        else:
            print('skipping telomeres anchored to unknown ref:', k)
            continue
        sorted_ref_keys.append((my_ind, LEXICO_2_IND[my_chr[:-1]], my_chr, k))
    sorted_ref_keys = sorted(sorted_ref_keys)
    #
    print('clustering anchored tels by position & filtering by read count...')
    CHR_TEL_DAT = []
    TEL_LEN_OUT = {}    # indexed by chr
    READLEN_OUT = {}    # indexed by chr
    clust_num   = 0
    tel_composition_data = []
    for ki in range(len(sorted_ref_keys)):
        (my_chr, k) = (sorted_ref_keys[ki][2], sorted_ref_keys[ki][3])
        sort_list = sorted([(ANCHORED_TEL_BY_CHR[k][n][1], n) for n in range(len(ANCHORED_TEL_BY_CHR[k]))])
        clusters  = cluster_list(sort_list, ANCHOR_CLUSTER_DIST, which_val=0)
        num_reads = sum([len(n) for n in clusters])
        sys.stdout.write(' - ' + my_chr + ': ' + str(len(clusters)) + ' anchor / ' + str(num_reads) + ' read')
        num_clust = 0
        num_reads = 0
        for cl in clusters:
            if len(cl) < MIN_READS_PER_CLUST:
                continue
            #if my_chr not in DEBUG_CHR_LIST:
            #   continue
            num_clust += 1
            num_reads += len(cl)
            pos_list = [n[0] for n in cl]
            ind_list = [n[1] for n in cl]
            my_pos   = int(np.median(pos_list))
            clust_num += 1
            #
            sorted_by_tlen = sorted([(ANCHORED_TEL_BY_CHR[k][i][3], len(ANCHORED_TEL_BY_CHR[k][i][6]), ANCHORED_TEL_BY_CHR[k][i][0]) for i in ind_list])
            my_tlens       = [n_sbtl[0] for n_sbtl in sorted_by_tlen]
            my_rlens       = [n_sbtl[1] for n_sbtl in sorted_by_tlen]
            my_rnames      = [n_sbtl[2] for n_sbtl in sorted_by_tlen]
            #
            consensus_tl = choose_tl_from_observations(my_tlens, CHR_TL_METHOD)
            #
            CHR_TEL_DAT.append([my_chr,
                                str(my_pos),
                                str(int(consensus_tl)),
                                ','.join([str(n) for n in my_tlens]),
                                ','.join([str(n) for n in my_rlens])])
            #
            if my_chr not in TEL_LEN_OUT:
                TEL_LEN_OUT[my_chr] = []
                READLEN_OUT[my_chr] = []
            TEL_LEN_OUT[my_chr].extend(my_tlens)
            READLEN_OUT[my_chr].extend(my_rlens)
            #
            clust_tcdat = []
            gtc_params  = [my_chr, clust_num, KMER_LIST, KMER_LIST_REV, KMER_ISSUBSTRING, READ_TYPE]
            for i in ind_list:
                clust_tcdat.append(get_telomere_composition(ANCHORED_TEL_BY_CHR[k][i], gtc_params))
            clust_tcout = [my_chr, my_pos, clust_num, ind_list, my_rlens, my_rnames, clust_tcdat]
            tel_composition_data.append(copy.deepcopy(clust_tcout))
            #
        if num_clust != len(clusters):
            sys.stdout.write(' --> ' + str(num_clust) + ' anchor / ' + str(num_reads) + ' read\n')
        else:
            sys.stdout.write('\n')
        sys.stdout.flush()

    #
    # write out chromosome TLs
    #
    print('writing TL results to "' + OUT_CHR_TL.split('/')[-1] + '"...')
    f = open(OUT_CHR_TL, 'w')
    f.write('#chr' + '\t' +
            'position' + '\t' +
            'TL_' + CHR_TL_METHOD + '\t' +
            'read_TLs' + '\t' +
            'read_lengths' + '\n')
    for i in range(len(CHR_TEL_DAT)):
        f.write('\t'.join(CHR_TEL_DAT[i]) + '\n')
    f.close()

    #
    # violin plots
    #
    print('creating violin plots...')
    tlen_violin_params = {'p_ymax':VIOLIN_TLEN_MAX,
                          'q_ymax':VIOLIN_TLEN_MAX,
                          'y_step':VIOLIN_TLEN_TICK}
    rlen_violin_params = {'p_color':'gray',
                          'q_color':'gray',
                          'y_label':'<-- q     read length     p -->',
                          'p_ymax':VIOLIN_RLEN_MAX,
                          'q_ymax':VIOLIN_RLEN_MAX,
                          'y_step':VIOLIN_RLEN_TICK}
    tel_len_violin_plot(TEL_LEN_OUT, OUT_VIOLIN_TL, custom_plot_params=tlen_violin_params)
    tel_len_violin_plot(READLEN_OUT, OUT_VIOLIN_RL, custom_plot_params=rlen_violin_params)

    ####################################################
    #                                                  #
    # TELOGATOR 2.0: CLUSTERING ALLELES BY TVR PATTERN #
    #                                                  #
    ####################################################

    # first lets check to make sure the muscle executable works
    check_muscle_version(MUSCLE_EXE)

    #####
    ##### for HiFi reads, remove noisy telomeres
    #####
    ####if READ_TYPE in ['hifi']:
    ####    sys.stdout.write('[HiFi only] filtering reads with noisy telomeres...')
    ####    sys.stdout.flush()
    ####    tt = time.perf_counter()
    ####    num_starting_reads = 0
    ####    for tci in range(len(tel_composition_data)):
    ####        [my_chr, my_pos, clust_num, ind_list, my_rlens, my_rnames, kmer_hit_dat] = tel_composition_data[tci]
    ####        num_starting_reads += len(kmer_hit_dat)
    ####        denoise_clust_dat   = filter_by_denoise_frac(kmer_hit_dat, KMER_METADATA, my_chr)
    ####        del_keys = []
    ####        for dci in range(len(denoise_clust_dat)):
    ####            if denoise_clust_dat[dci][0] > MAX_NOISE_BASES or denoise_clust_dat[dci][1] > MAX_NOISE_FRAC:
    ####                del_keys.append(dci)
    ####        for di in sorted(del_keys,reverse=True):
    ####            del tel_composition_data[tci][3][di]
    ####            del tel_composition_data[tci][4][di]
    ####            del tel_composition_data[tci][5][di]
    ####            del tel_composition_data[tci][6][di]
    ####    for tci in range(len(tel_composition_data)-1,-1,-1):
    ####        if len(tel_composition_data[tci][3]) == 0:
    ####            del tel_composition_data[tci]
    ####    num_ending_reads = 0
    ####    for tci in range(len(tel_composition_data)):
    ####        num_ending_reads += len(tel_composition_data[tci][6])
    ####    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    ####    sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads\n')
    ####    sys.stdout.flush()

    #
    #
    #
    print('pairwise comparing TVRs at each cluster...')
    allele_tel_dat_temp = []
    ALLELE_TEL_DAT      = []
    orphaned_read_dat   = []
    temp_read_dat       = {}    # contains read fastas for downstream purposes. indexed by (chr, pos, allele_i)
    #
    gatd_params = [ALLELE_TL_METHOD, MIN_READS_PER_PHASE]
    mtp_params  = [KMER_METADATA, KMER_COLORS, MIN_READS_PER_PHASE, PLOT_FILT_CVECS, DUMMY_TEL_MAPQ, DONT_OVERWRITE_PLOTS]
    #
    for tc_dat in tel_composition_data:
        [my_chr, my_pos, clust_num, ind_list, my_rlens, my_rnames, kmer_hit_dat] = tc_dat
        if len(DEBUG_CHR_LIST) and my_chr not in DEBUG_CHR_LIST:
            print(' - skipping', my_chr, '(not in debug list)')
            continue
        sys.stdout.write(' - ' + my_chr + ':' + str(my_pos) + ' ' + str(len(ind_list)) + ' reads')
        sys.stdout.flush()
        tt = time.perf_counter()
        #
        my_tc = TREECUT_DEFAULT
        if (my_chr,None) in custom_treecut_vals:
            my_tc = custom_treecut_vals[(my_chr,None)]
        elif (my_chr,my_pos) in custom_treecut_vals:
            my_tc = custom_treecut_vals[(my_chr,my_pos)]
        #
        if k[:3] == 'alt':
            plotname_chr = 'alt-' + my_chr
        else:
            plotname_chr = my_chr
        zfcn = str(clust_num).zfill(3)
        telcompplot_fn = OUT_TVR_DIR  + 'tvr-reads-'     + zfcn + '_' + plotname_chr + '.png'
        telcompcons_fn = OUT_TVR_DIR  + 'tvr-consensus-' + zfcn + '_' + plotname_chr + '.png'
        dendrogram_fn  = OUT_TVR_TEMP + 'dendrogram-'    + zfcn + '_' + plotname_chr + '.png'
        dend_prefix_fn = OUT_TVR_TEMP + 'p_dendrogram-'  + zfcn + '_' + plotname_chr + '.png'
        dist_matrix_fn = OUT_TVR_TEMP + 'cluster-'       + zfcn + '_' + plotname_chr + '.npy'
        dist_prefix_fn = OUT_TVR_TEMP + 'p_cluster-'     + zfcn + '_' + plotname_chr + '.npy'
        consensus_fn   = OUT_TVR_TEMP + 'consensus-seq-' + zfcn + '_' + plotname_chr + '.fa'
        #
        if ALWAYS_REPROCESS:
            dist_matrix_fn = None
            dist_prefix_fn = None
            consensus_fn = None
        #
        read_clust_dat = cluster_tvrs(kmer_hit_dat, KMER_METADATA, my_chr, my_pos, my_tc, TREECUT_PREFIXMERGE,
                                      aln_mode=DS_OR_MS,
                                      alignment_processes=NUM_PROCESSES,
                                      rand_shuffle_count=RAND_SHUFFLE,
                                      dist_in=dist_matrix_fn,
                                      dist_in_prefix=dist_prefix_fn,
                                      fig_name=dendrogram_fn,
                                      fig_prefix_name=dend_prefix_fn,
                                      save_msa=consensus_fn,
                                      muscle_dir=OUT_TVR_TEMP,
                                      muscle_exe=MUSCLE_EXE,
                                      PRINT_DEBUG=PRINT_DEBUG)
        #
        # lets retain reads that didn't make it into a cluster
        #
        read_was_clustered = {n:-1 for n in range(len(kmer_hit_dat))}
        read_subtel_adj    = {n:0 for n in range(len(kmer_hit_dat))}
        for allele_i in range(len(read_clust_dat[0])):
            allele_readcount = len(read_clust_dat[0][allele_i])
            for rcd_i, read_i in enumerate(read_clust_dat[0][allele_i]):
                if allele_readcount >= MIN_READS_PER_PHASE:
                    read_was_clustered[read_i] = allele_i
                read_subtel_adj[read_i]    = kmer_hit_dat[read_i][1] - read_clust_dat[3][read_i]
        # save fasta sequences for final clustering of all alleles
        for khd_i, allele_i in read_was_clustered.items():
            # refine initial fasta sequences based on sub/tel lens determined during clustering
            my_tel_sequence_len = read_subtel_adj[khd_i]
            my_sub_sequence_len = len(kmer_hit_dat[khd_i][6][2][1]) - my_tel_sequence_len
            my_tel_type         = kmer_hit_dat[khd_i][6][0][0].split('_')[2][-1]
            #print(kmer_hit_dat[khd_i][4], my_tel_sequence_len, my_sub_sequence_len, my_tel_type)
            if my_chr[-1] == 'p':
                if my_tel_type == 'q':
                    kmer_hit_dat[khd_i][6][2] = (kmer_hit_dat[khd_i][6][2][0], RC(kmer_hit_dat[khd_i][6][2][1]))
                kmer_hit_dat[khd_i][6][0] = (kmer_hit_dat[khd_i][6][0][0], kmer_hit_dat[khd_i][6][2][1][:my_tel_sequence_len])
                kmer_hit_dat[khd_i][6][1] = (kmer_hit_dat[khd_i][6][1][0], kmer_hit_dat[khd_i][6][2][1][my_tel_sequence_len:])
            else:
                if my_tel_type == 'p':
                    kmer_hit_dat[khd_i][6][2] = (kmer_hit_dat[khd_i][6][2][0], RC(kmer_hit_dat[khd_i][6][2][1]))
                kmer_hit_dat[khd_i][6][0] = (kmer_hit_dat[khd_i][6][0][0], kmer_hit_dat[khd_i][6][2][1][my_sub_sequence_len:])
                kmer_hit_dat[khd_i][6][1] = (kmer_hit_dat[khd_i][6][1][0], kmer_hit_dat[khd_i][6][2][1][:my_sub_sequence_len])
            #
            if allele_i == -1:
                orphaned_read_dat.append((my_chr, my_pos, my_rlens[khd_i], copy.deepcopy(kmer_hit_dat[khd_i])))
            else:
                my_key = (my_chr, my_pos, allele_i)
                if my_key not in temp_read_dat:
                    temp_read_dat[my_key] = []
                temp_read_dat[my_key].append([(n[0], n[1]) for n in kmer_hit_dat[khd_i][6]])
        #
        # values for output tsv
        #
        allele_tel_dat_temp.extend(get_allele_tsv_dat(kmer_hit_dat, read_clust_dat, my_chr, my_pos, my_rlens, gatd_params))
        #
        # plotting!
        #
        make_tvr_plots(kmer_hit_dat, read_clust_dat, my_chr, my_pos, telcompplot_fn, telcompcons_fn, mtp_params)
        #
        # fasta data for subtel / tel, can be used for debugging
        #
        for allele_i in range(len(read_clust_dat[0])):
            allele_fasta_all = [kmer_hit_dat[n][6] for n in read_clust_dat[0][allele_i]]
            allele_fasta_tel = [n[0] for n in allele_fasta_all]
            #print(allele_fasta_tel)
            #exit(1)
        #
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()

    #
    # analyze orphaned reads
    #
    #   kmer_dat[i]          = [[kmer1_hits, kmer2_hits, ...], tlen, tel_anchor_dist, read_orientation, readname, anchor_mapq, fasta_dat]
    #   orphaned_read_dat[i] = [my_chr, my_pos, my_rlen, my_kmer_dat]
    #
    if SKIP_ORPHAN_CLUSTERING is False and len(orphaned_read_dat) >= MIN_READS_PER_PHASE:
        sys.stdout.write('clustering ' + str(len(orphaned_read_dat)) + ' orphan reads...')
        sys.stdout.flush()
        orphans_kmer_hit_dat = []
        orphans_rlens = [n[2] for n in orphaned_read_dat]
        for ori in range(len(orphaned_read_dat)):
            # reverse p-arm kmer hits so everything is in same orientation
            if orphaned_read_dat[ori][0][-1] == 'p':
                # reverse kmer hits
                or_tlen = orphaned_read_dat[ori][3][1]
                orphaned_read_dat[ori][3][0] = [[[or_tlen-n[1], or_tlen-n[0]] for n in m[::-1]] for m in orphaned_read_dat[ori][3][0]]
                # swap orientation label (not sure if this is ever used anywhere, honestly)
                orphaned_read_dat[ori][3][3] = ORR_SWAP[orphaned_read_dat[ori][3][3]]
                # swap read sequence fasta
                orphaned_read_dat[ori][3][6] = [(orphaned_read_dat[ori][3][6][0][0], RC(orphaned_read_dat[ori][3][6][0][1])),
                                                (orphaned_read_dat[ori][3][6][1][0], RC(orphaned_read_dat[ori][3][6][1][1][:SUBTEL_CLUSTER_SIZE[1]])),
                                                (orphaned_read_dat[ori][3][6][2][0], RC(orphaned_read_dat[ori][3][6][2][1]))]
            else:
                orphaned_read_dat[ori][3][6] = [(orphaned_read_dat[ori][3][6][0][0], orphaned_read_dat[ori][3][6][0][1]),
                                                (orphaned_read_dat[ori][3][6][1][0], orphaned_read_dat[ori][3][6][1][1][-SUBTEL_CLUSTER_SIZE[1]:]),
                                                (orphaned_read_dat[ori][3][6][2][0], orphaned_read_dat[ori][3][6][2][1])]
            orphans_kmer_hit_dat.append(copy.deepcopy(orphaned_read_dat[ori][3]))
            # embed original chromosome in the read name
            orphans_kmer_hit_dat[-1][4] = orphans_kmer_hit_dat[-1][4] + ' ' + orphaned_read_dat[ori][0] + ':' + str(orphaned_read_dat[ori][1])
        #
        orphans_telcompplot_fn = OUT_TVR_DIR  + 'tvr-reads-999_orphans.png'
        orphans_telcompcons_fn = OUT_TVR_DIR  + 'tvr-consensus-999_orphans.png'
        orphans_dendrogram_fn  = OUT_TVR_TEMP + 'dendrogram-999_orphans.png'
        orphans_dend_prefix_fn = OUT_TVR_TEMP + 'p_dendrogram-999_orphans.png'
        orphans_dist_matrix_fn = OUT_TVR_TEMP + 'cluster-999_orphans.npy'
        orphans_dist_prefix_fn = OUT_TVR_TEMP + 'p_cluster-999_orphans.npy'
        orphans_consensus_fn   = OUT_TVR_TEMP + 'consensus-seq-999_orphans.fa'
        #
        if ALWAYS_REPROCESS:
            orphans_dist_matrix_fn = None
            orphans_dist_prefix_fn = None
            orphans_consensus_fn = None
        #
        orphans_clust_dat = cluster_tvrs(orphans_kmer_hit_dat, KMER_METADATA, orphans_chr, orphans_pos, TREECUT_ORPHANS, TREECUT_PREFIXMERGE,
                                         aln_mode=DS_OR_MS,
                                         alignment_processes=NUM_PROCESSES,
                                         rand_shuffle_count=RAND_SHUFFLE,
                                         dist_in=orphans_dist_matrix_fn,
                                         dist_in_prefix=orphans_dist_prefix_fn,
                                         fig_name=orphans_dendrogram_fn,
                                         fig_prefix_name=orphans_dend_prefix_fn,
                                         save_msa=orphans_consensus_fn,
                                         muscle_dir=OUT_TVR_TEMP,
                                         muscle_exe=MUSCLE_EXE,
                                         PRINT_DEBUG=PRINT_DEBUG)
        #
        for allele_i in range(len(orphans_clust_dat[0])):
            allele_readcount = len(orphans_clust_dat[0][allele_i])
            if allele_readcount >= MIN_READS_PER_PHASE:
                for read_i in orphans_clust_dat[0][allele_i]:
                    my_key = (orphans_chr, orphans_pos, allele_i)
                    if my_key not in temp_read_dat:
                        temp_read_dat[my_key] = []
                    temp_read_dat[my_key].append(copy.deepcopy(orphaned_read_dat[read_i][3][6]))
        #
        allele_tel_dat_temp.extend(get_allele_tsv_dat(orphans_kmer_hit_dat, orphans_clust_dat, orphans_chr, orphans_pos, orphans_rlens, gatd_params))
        #
        make_tvr_plots(orphans_kmer_hit_dat, orphans_clust_dat, orphans_chr, orphans_pos, orphans_telcompplot_fn, orphans_telcompcons_fn, mtp_params)
        #
        sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
        sys.stdout.flush()

    #
    # collapse all blank TVRs into a single allele
    #
    del_inds_blank  = []
    blank_read_keys = []
    blank_entry = [blank_chr, str(blank_pos), '0', '0', '', '', '', '', '0', '', '', '']
    original_chr_mapping_of_blank_reads = {}
    if MERGE_BLANK:
        for i in range(len(allele_tel_dat_temp)):
            if len(allele_tel_dat_temp[i][9]) == 0:
                del_inds_blank.append(i)
                blank_entry[5]  += allele_tel_dat_temp[i][5] + ','
                blank_entry[6]  += allele_tel_dat_temp[i][6] + ','
                blank_entry[7]  += allele_tel_dat_temp[i][7] + ','
                blank_entry[10] += allele_tel_dat_temp[i][10] + ','
                blank_read_keys.append((allele_tel_dat_temp[i][0], int(allele_tel_dat_temp[i][1]), int(allele_tel_dat_temp[i][2])))
                for rname in allele_tel_dat_temp[i][10].split(','):
                    original_chr_mapping_of_blank_reads[rname] = (allele_tel_dat_temp[i][0], int(allele_tel_dat_temp[i][1]))
        for di in sorted(del_inds_blank, reverse=True):
            del allele_tel_dat_temp[di]
        if len(blank_entry[5]):
            if blank_entry[5][-1] == ',':
                blank_entry[5] = blank_entry[5][:-1]
            if blank_entry[6][-1] == ',':
                blank_entry[6] = blank_entry[6][:-1]
            if blank_entry[7][-1] == ',':
                blank_entry[7] = blank_entry[7][:-1]
            if blank_entry[10][-1] == ',':
                blank_entry[10] = blank_entry[10][:-1]

    #
    # organize read data and clustered alleles for the final stretch of processing
    #
    clustered_read_dat  = {}
    for atdat in allele_tel_dat_temp:
        ALLELE_TEL_DAT.append(copy.deepcopy(atdat))
        my_key = (atdat[0], int(atdat[1]), int(atdat[2]))
        clustered_read_dat[my_key] = copy.deepcopy(temp_read_dat[my_key])
    #
    # blank tvr/subtel sequences need to be manipulated so that ps and qs can be compared to each other later
    #
    if MERGE_BLANK and len(blank_read_keys):
        clustered_read_dat[(blank_chr, blank_pos, 0)] = []
        #for k in blank_read_keys:
        #   clustered_read_dat[(blank_chr, blank_pos, 0)] += [copy.deepcopy(n) for n in temp_read_dat[k]]
        for k in blank_read_keys:
            if k[0][-1] == 'p':
                clustered_read_dat[(blank_chr, blank_pos, 0)] += [[(n[0][0], RC(n[0][1])),
                                                                   (n[1][0], RC(n[1][1][:SUBTEL_CLUSTER_SIZE[1]])),
                                                                   (n[2][0], RC(n[2][1]))] for n in temp_read_dat[k]]
            else:
                clustered_read_dat[(blank_chr, blank_pos, 0)] += [[(n[0][0], n[0][1]),
                                                                   (n[1][0], n[1][1][-SUBTEL_CLUSTER_SIZE[1]:]),
                                                                   (n[2][0], n[2][1])] for n in temp_read_dat[k]]
        ALLELE_TEL_DAT.append(blank_entry)
    del allele_tel_dat_temp
    del temp_read_dat

    #
    #
    #
    if len(ALLELE_TEL_DAT) == 0:
        print('No telomere alleles found, something probably went wrong?...')
        exit(1)

    #
    # pairwise tvr comparison --> fill out allele_id field
    #
    sys.stdout.write('pairwise comparing all consensus TVR sequences...')
    sys.stdout.flush()
    tt = time.perf_counter()
    all_tvrs     = []
    tvr_labels   = []
    all_labels   = []
    blank_cluster = []
    for i in range(len(ALLELE_TEL_DAT)):
        if len(ALLELE_TEL_DAT[i][9]):
            all_tvrs.append(ALLELE_TEL_DAT[i][9])
            tvr_labels.append(ALLELE_TEL_DAT[i][0]+'_'+ALLELE_TEL_DAT[i][1]+'_'+ALLELE_TEL_DAT[i][2]+'_'+str(i))
        else:
            blank_cluster.append(i)
        all_labels.append(ALLELE_TEL_DAT[i][0]+'_'+ALLELE_TEL_DAT[i][1]+'_'+ALLELE_TEL_DAT[i][2]+'_'+str(i))
    if len(all_tvrs) == 1:      # this can happen when debugging individual arms
        alltvr_clustdat = [[0]]
    else:
        alltvr_dist_fn   = OUT_TVR_TEMP+'all-tvrs.npy'
        alltvr_plot_fn   = OUT_TVR_DIR+'all-tvrs.png'
        alltvr_dendro_fn = OUT_TVR_DIR+'all-tvrs_dendrogram.png'
        #
        if ALWAYS_REPROCESS:
            alltvr_dist_fn = None
        #
        alltvr_clustdat = cluster_consensus_tvrs(all_tvrs, KMER_METADATA, TREECUT_MULTIMAPPED,
                                                 alignment_processes=NUM_PROCESSES,
                                                 aln_mode='ds',
                                                 gap_bool=(True,False),                 # allow grouping prefixes (risky?)
                                                 rand_shuffle_count=RAND_SHUFFLE,
                                                 dist_in=alltvr_dist_fn,
                                                 fig_name=alltvr_plot_fn,
                                                 dendro_name=alltvr_dendro_fn,
                                                 samp_labels=tvr_labels,
                                                 linkage_method='complete',
                                                 job=(1,1),
                                                 dendrogram_height=12,
                                                 overwrite_figures=not(DONT_OVERWRITE_PLOTS))
    if MERGE_BLANK and len(blank_cluster):
        sort_alltvr_clustdat = sorted(alltvr_clustdat) + [blank_cluster]
    else:
        sort_alltvr_clustdat = sorted(alltvr_clustdat)
    sys.stdout.write(' (' + str(int(time.perf_counter() - tt)) + ' sec)\n')
    sys.stdout.flush()

    #
    # decluster alleles that have the same TVR if their subtelomeres are different enough
    #
    if SKIP_SUBTEL_CLUSTERING is False:
        print('comparing subtelomeres...')
    new_allele_tel_dat_entries = []
    num_allelles_added_because_subtel = 0
    for i in range(len(sort_alltvr_clustdat)):
        subtels_for_this_cluster = []
        subtel_labels = []
        chrs_hit = {}
        subtel_sizes = []
        for j in sort_alltvr_clustdat[i]:
            splt = all_labels[j].split('_')
            (my_c, my_p, my_a) = (splt[0], int(splt[1]), int(splt[2]))
            my_key = (my_c, my_p, my_a)
            chrs_hit[my_c] = True
            subtel_sizes.extend([len(n[1][1]) for n in clustered_read_dat[my_key]])
        subtel_size_for_this_cluster = max(min(min(subtel_sizes), SUBTEL_CLUSTER_SIZE[1]), SUBTEL_CLUSTER_SIZE[0])
        for j in sort_alltvr_clustdat[i]:
            splt = all_labels[j].split('_')
            (my_c, my_p, my_a) = (splt[0], int(splt[1]), int(splt[2]))
            my_key = (my_c, my_p, my_a)
            if my_c[-1] == 'p':
                subtels_for_this_cluster.extend([RC(n[1][1][:subtel_size_for_this_cluster]) for n in clustered_read_dat[my_key]])
            elif my_c[-1] == 'q':
                subtels_for_this_cluster.extend([n[1][1][-subtel_size_for_this_cluster:] for n in clustered_read_dat[my_key]])
            subtel_labels.extend(['_'.join(n[1][0].split('_')[3:]) for n in clustered_read_dat[my_key]])
        my_dendro_title = 'all-tvrs-cluster ' + str(i) + ' : ' + '/'.join(chrs_hit.keys())
        #
        if len(subtels_for_this_cluster) == 1:      # this can happen when debugging individual arms
            subtel_clustdat = [[0]]
        elif SKIP_SUBTEL_CLUSTERING:
            subtel_clustdat = [list(range(len(subtels_for_this_cluster)))]
        else:
            subtel_dist_fn   = OUT_TVR_TEMP+'subtels-'+str(i).zfill(3)+'.npy'
            subtel_dendro_fn = OUT_TVR_TEMP+'subtels_dendrogram-'+str(i).zfill(3)+'.png'
            #
            if ALWAYS_REPROCESS:
                subtel_dist_fn = None
            #
            subtel_clustdat = cluster_consensus_tvrs(subtels_for_this_cluster, KMER_METADATA, TREECUT_SUBTELS,
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
        usable_clusters = [n for n in subtel_clustdat if len(n) >= MIN_READS_PER_PHASE]
        if len(usable_clusters):
            subtel_clust_by_readname = {}
            for j in range(len(usable_clusters)):
                for k in usable_clusters[j]:
                    subtel_clust_by_readname[subtel_labels[k]] = j
            #
            # it's kind of silly how much I need to unpack and repack data in order to do this.
            #
            allele_tel_dat_by_subtel_clust = [[] for n in usable_clusters]
            for j in sort_alltvr_clustdat[i]:
                splt_tls   = ALLELE_TEL_DAT[j][5].split(',')
                splt_rlens = ALLELE_TEL_DAT[j][6].split(',')
                splt_mapq  = ALLELE_TEL_DAT[j][7].split(',')
                splt_rname = ALLELE_TEL_DAT[j][10].split(',')
                for k in range(len(splt_rname)):
                    if splt_rname[k] in subtel_clust_by_readname:
                        allele_tel_dat_by_subtel_clust[subtel_clust_by_readname[splt_rname[k]]].append((int(splt_tls[k]), splt_rlens[k], splt_mapq[k], splt_rname[k]))
            for j in range(len(allele_tel_dat_by_subtel_clust)):
                allele_tel_dat_by_subtel_clust[j] = sorted(allele_tel_dat_by_subtel_clust[j])
            #
            for j in sort_alltvr_clustdat[i]:
                splt = all_labels[j].split('_')
                if len(usable_clusters) == 1:
                    ALLELE_TEL_DAT[j][3] = str(i+1)
                else:
                    ALLELE_TEL_DAT[j][3] = str(i+1)+'a'
                ALLELE_TEL_DAT[j][4]  = str(int(choose_tl_from_observations([n[0] for n in allele_tel_dat_by_subtel_clust[0]], ALLELE_TL_METHOD)))
                ALLELE_TEL_DAT[j][5]  = ','.join([str(n[0]) for n in allele_tel_dat_by_subtel_clust[0]])
                ALLELE_TEL_DAT[j][6]  = ','.join([n[1] for n in allele_tel_dat_by_subtel_clust[0]])
                ALLELE_TEL_DAT[j][7]  = ','.join([n[2] for n in allele_tel_dat_by_subtel_clust[0]])
                ALLELE_TEL_DAT[j][10] = ','.join([n[3] for n in allele_tel_dat_by_subtel_clust[0]])
                if ALLELE_TEL_DAT[j][0] == blank_chr:
                    rname_list = [n[3] for n in allele_tel_dat_by_subtel_clust[0]]
                    ALLELE_TEL_DAT[j][11] = ','.join([original_chr_mapping_of_blank_reads[n][0] for n in rname_list])
                for k in range(1,len(allele_tel_dat_by_subtel_clust)):
                    new_allele_tel_dat_entries.append((j, [ALLELE_TEL_DAT[j][0],
                                                           ALLELE_TEL_DAT[j][1],
                                                           ALLELE_TEL_DAT[j][2],
                                                           str(i+1)+chr(97+k),
                                                           str(int(choose_tl_from_observations([n[0] for n in allele_tel_dat_by_subtel_clust[k]], ALLELE_TL_METHOD))),
                                                           ','.join([str(n[0]) for n in allele_tel_dat_by_subtel_clust[k]]),
                                                           ','.join([n[1] for n in allele_tel_dat_by_subtel_clust[k]]),
                                                           ','.join([n[2] for n in allele_tel_dat_by_subtel_clust[k]]),
                                                           ALLELE_TEL_DAT[j][8],
                                                           ALLELE_TEL_DAT[j][9],
                                                           ','.join([n[3] for n in allele_tel_dat_by_subtel_clust[k]]),
                                                           '']))
                    if new_allele_tel_dat_entries[-1][1][0] == blank_chr:
                        rname_list = [n[3] for n in allele_tel_dat_by_subtel_clust[k]]
                        new_allele_tel_dat_entries[-1][1][11] = ','.join([original_chr_mapping_of_blank_reads[n][0] for n in rname_list])
            num_allelles_added_because_subtel += len(usable_clusters) - 1
    #
    # insert subcluster splits
    if len(new_allele_tel_dat_entries):
        new_allele_tel_dat_entries = sorted(new_allele_tel_dat_entries, reverse=True)
        for n in new_allele_tel_dat_entries:
            ALLELE_TEL_DAT.insert(n[0]+1, n[1])
    #
    # at this point allele_id == 0 means cluster was declustered due to different subtels and no subcluster had enough reads
    if MERGE_BLANK:
        del_inds = []
        for i in range(len(ALLELE_TEL_DAT)):
            if ALLELE_TEL_DAT[i][3] == '0':
                del_inds.append(i)
        for di in sorted(del_inds, reverse=True):
            del ALLELE_TEL_DAT[di]

    #
    # print some TVR stats to console
    #
    tvr_nums = {}
    for n in ALLELE_TEL_DAT:
        if n[3].isdigit():
            tvr_nums[int(n[3])] = True
        else:
            tvr_nums[int(n[3][:-1])] = True
    num_unique_tvrs    = len(tvr_nums)
    num_unique_alleles = len(set([n[3] for n in ALLELE_TEL_DAT]))
    num_blank_alleles  = len([n[0] for n in ALLELE_TEL_DAT if n[0] == blank_chr])
    print(' -', num_unique_alleles, 'unique alleles')
    print(' -', num_unique_tvrs, 'unique TVRs')
    print(' -', num_blank_alleles, 'blank TVRs')
    print(' -', num_allelles_added_because_subtel, 'alleles split to different subtels')

    #
    # write out allele TLs
    #
    print('writing allele TL results to "' + OUT_ALLELE_TL.split('/')[-1] + '"...')
    f = open(OUT_ALLELE_TL, 'w')
    f.write('#chr' + '\t' +
            'position' + '\t' +
            'allele_num' + '\t' +
            'allele_id' + '\t' +
            'TL_' + ALLELE_TL_METHOD + '\t' +
            'read_TLs' + '\t' +
            'read_lengths' + '\t' +
            'read_mapq' + '\t' +
            'tvr_len' + '\t' +
            'tvr_consensus' + '\t' +
            'supporting_reads' + '\t' +
            'original_chr' + '\n')
    for i in range(len(ALLELE_TEL_DAT)):
        f.write('\t'.join(ALLELE_TEL_DAT[i]) + '\n')
    f.close()


if __name__ == '__main__':
    main()