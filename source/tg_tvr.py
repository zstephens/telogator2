import multiprocessing
import numpy as np
import matplotlib.pyplot as mpl

from Bio.Align import PairwiseAligner, substitution_matrices

from collections import Counter

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance  import squareform

from source.tg_muscle import get_muscle_msa
from source.tg_plot   import DEFAULT_DPI, MAX_PLOT_SIZE, plot_kmer_hits
from source.tg_reader import TG_Reader
from source.tg_util   import exists_and_is_nonzero, shuffle_seq

# we're going to pretend each kmer color is an amino acid, so alignment tools behave themselves
AMINO = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
AMINO_2_IND = {AMINO[i]:i for i in range(len(AMINO))}
UNKNOWN_LETTER = AMINO[0]

MAX_TVR_LEN       = 8000    # ignore variant repeats past this point when finding TVR boundary
MAX_TVR_LEN_SHORT = 3000    # when examining TVRs with very few variant repeats
TVR_BOUNDARY_BUFF = 30      # add this many bp to detected TVR boundary

# if tvr + tel region has at least this many unknown characters, use minimum-dens pos instead of first-below-dens ps
TOO_MUCH_UNKNOWN_IN_TVRTEL = 1000
# density parameters for identifing subtel / tvr boundaries (individual reads)
UNKNOWN_WIN_SIZE = 100
UNKNOWN_END_DENS = 0.120
# density parameters for recovering variant repeats from subtel consensus
UNKNOWN_WIN_SIZE_CONS = 50
UNKNOWN_END_DENS_CONS = 0.800
# density parameters for discerning canonical regions from sequencing artifacts
CANON_WIN_SIZE = 100
CANON_END_DENS = 0.700
# parameters for determining tvr / canonical boundaries: (denoise_region_size, cum_thresh, min_hits)
TVR_CANON_FILT_PARAMS_STRICT  = (10, 0.05, 100)
TVR_CANON_FILT_PARAMS_LENIENT = ( 5, 0.20,  50)

#
# scoring matrix parameters
#
MATCH_NORMAL  = 5
XMATCH_NORMAL = -4
GAP_OPEN = -4
GAP_EXT  = -4
#
MATCH_CANON  = 0
XMATCH_CANON = -4
MATCH_CANON_PREFIX = 1  # when doing prefix merging we want to consider canonical at least a little
#
MATCH_UNKNOWN  = 2
XMATCH_UNKNOWN = -4
MATCH_CANON_UNKNOWN = -2

# when comparing sequences of different lengths, choose min(len(seq1), len(seq2)), but only go this low
MIN_VIABLE_SEQ_LEN = 1000
# the log-based distance function can go to infinity so lets set an upper bound on it
MAX_SEQ_DIST = 10.0
# similarly, lets choose a small number as the minimum to prevent numerical weirdness from giving us negative values
MIN_SEQ_DIST = 0.0001
# to prevent pesky div-by-zeros in edge cases
MIN_MSD      = 3.0

#
#
#
def get_scoring_matrix(canonical_letter, canonical_score=MATCH_CANON):
    all_amino = tuple(AMINO)
    d = {}
    for c1 in all_amino:
        for c2 in all_amino:
            if c1 == c2:
                d[c1,c2] = float(MATCH_NORMAL)
            else:
                d[c1,c2] = float(XMATCH_NORMAL)
    for c1 in all_amino:
        d[c1,canonical_letter] = float(XMATCH_CANON)
        d[canonical_letter,c1] = float(XMATCH_CANON)
    for c1 in all_amino:
        d[c1,UNKNOWN_LETTER] = float(XMATCH_UNKNOWN)
        d[UNKNOWN_LETTER,c1] = float(XMATCH_UNKNOWN)
    d[canonical_letter,canonical_letter] = float(canonical_score)      # matching canonical
    d[UNKNOWN_LETTER,UNKNOWN_LETTER]     = float(MATCH_UNKNOWN)        # matching unknown
    d[canonical_letter,UNKNOWN_LETTER]   = float(MATCH_CANON_UNKNOWN)  # canon = unknown
    d[UNKNOWN_LETTER,canonical_letter]   = float(MATCH_CANON_UNKNOWN)  # canon = unknown
    return substitution_matrices.Array(alphabet=None, dims=None, data=d)

#
#
#
def parallel_alignment_job(sequences,
                           ij,
                           pq,
                           results_dict,
                           scoring_matrix=None,
                           gap_bool=(True,True),
                           adjust_lens=True,
                           min_viable=True,
                           randshuffle=3):
    #
    # create aligner object
    #
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = float(MATCH_NORMAL)
    aligner.mismatch_score = float(XMATCH_NORMAL)
    aligner.open_gap_score = float(GAP_OPEN)
    aligner.extend_gap_score = float(GAP_EXT)
    if scoring_matrix is not None:
        aligner.alphabet = ''.join(AMINO)
        aligner.substitution_matrix = scoring_matrix
    if gap_bool[0]:
        aligner.left_open_gap_score   = float(GAP_OPEN)
        aligner.left_extend_gap_score = float(GAP_EXT)
    else:
        aligner.left_open_gap_score   = 0.0
        aligner.left_extend_gap_score = 0.0
    if gap_bool[1]:
        aligner.right_open_gap_score   = float(GAP_OPEN)
        aligner.right_extend_gap_score = float(GAP_EXT)
    else:
        aligner.right_open_gap_score   = 0.0
        aligner.right_extend_gap_score = 0.0
    #
    # do the alignments
    #
    for (i,j) in ij:
        #
        if adjust_lens:
            min_len = min(len(sequences[i]), len(sequences[j]))
            if min_viable:
                min_len = max(min_len, MIN_VIABLE_SEQ_LEN)
            if pq == 'p':
                seq_i = sequences[i][-min_len:]
                seq_j = sequences[j][-min_len:]
            elif pq == 'q':
                seq_i = sequences[i][:min_len]
                seq_j = sequences[j][:min_len]
        else:
            seq_i = sequences[i]
            seq_j = sequences[j]
        #
        # aln score
        #
        aln_score = aligner.score(seq_i, seq_j)
        #
        # iden score
        #
        c1 = Counter(seq_i)
        c2 = Counter(seq_j)
        if scoring_matrix is None:
            iden_score = MATCH_NORMAL*((len(seq_i)+len(seq_j))/2.)
        else:
            is1 = sum([c1[n]*scoring_matrix[n,n] for n in c1.keys()])
            is2 = sum([c2[n]*scoring_matrix[n,n] for n in c2.keys()])
            iden_score = (is1+is2)/2.
        #
        # rand score
        #
        rand_scores = []
        for k in range(randshuffle):
            rand_scores.append(aligner.score(shuffle_seq(seq_i), shuffle_seq(seq_j)))
        rand_score_shuffle = np.mean(rand_scores)
        #
        # convert to distance, bounds checking
        #
        if rand_score_shuffle >= aln_score:
            my_dist_shuffle = MAX_SEQ_DIST
        else:
            my_dist_shuffle = min(-np.log((aln_score - rand_score_shuffle)/(iden_score - rand_score_shuffle)), MAX_SEQ_DIST)
        my_dist_shuffle = max(my_dist_shuffle, MIN_SEQ_DIST)
        #
        # output
        #
        #print((i,j), aln_score, iden_score, rand_scores, '{0:0.3f}'.format(my_dist_shuffle))
        results_dict[(i,j)] = my_dist_shuffle

#
#
#
def filter_by_denoise_frac(kmer_dat, repeats_metadata, my_chr):
    #
    n_reads = len(kmer_dat)
    pq      = my_chr[-1]
    [kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
    denoise_chars = []
    for i in range(len(kmer_list)):
        if 'denoise' in kmer_flags[i]:
            denoise_chars.append(kmer_letters[i])
    denoise_chars = list(set(denoise_chars))
    #
    # create color vector
    #
    my_out = []
    for i in range(n_reads):
        [my_kmer_hits, my_tlen, my_dbta, my_orr, my_rname, my_mapq, my_fastadat] = kmer_dat[i]
        my_letters = [UNKNOWN_LETTER for n in range(my_tlen)]
        for ki in range(len(my_kmer_hits)):
            if len(my_kmer_hits[ki]):
                for kmer_span in my_kmer_hits[ki]:
                    for j in range(kmer_span[0], kmer_span[1]):
                        my_letters[j] = kmer_letters[ki]
        my_cvec = ''.join(my_letters)
        #
        if pq == 'p':
            (seq_left, seq_right) = find_density_boundary(my_cvec[::-1], UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
        elif pq == 'q':
            (seq_left, seq_right) = find_density_boundary(my_cvec, UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
        #
        my_count_unknown = seq_right.count(UNKNOWN_LETTER)
        my_count_denoise = 0
        for dc in denoise_chars:
            my_count_denoise += seq_right.count(dc)
        my_num = my_count_unknown + my_count_denoise
        my_den = float(max([1, len(seq_right)]))
        my_out.append((my_num, my_num/my_den))
    return my_out

#
#   kmer_dat[i] = [[kmer1_hits, kmer2_hits, ...], tlen, tel_anchor_dist, read_orientation, readname, anchor_mapq, fasta_dat]
#
#   repeats_metadata = [kmer_list, kmer_colors, kmer_letters, kmer_flags]
#
def cluster_tvrs(kmer_dat,
                 repeats_metadata,
                 my_chr,
                 my_pos,
                 tree_cut,
                 tree_cut_prefix,
                 aln_mode='ds',
                 rand_shuffle_count=3,
                 dist_in=None,
                 dist_in_prefix=None,
                 fig_name=None,
                 fig_prefix_name=None,
                 muscle_dir='',
                 save_msa=None,
                 tvr_truncate=3000,
                 alignment_processes=4,
                 muscle_exe='muscle',
                 PRINT_DEBUG=False):
    #
    [kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
    #
    n_reads = len(kmer_dat)
    pq      = my_chr[-1]
    #
    canonical_letter = None
    denoise_chars    = []
    tvr_letters      = []
    nofilt_letters   = []
    for i in range(len(kmer_list)):
        if 'canonical' in kmer_flags[i]:
            canonical_letter = kmer_letters[i]
        if 'denoise' in kmer_flags[i]:
            denoise_chars.append(kmer_letters[i])
        if 'tvr' in kmer_flags[i]:
            tvr_letters.append(kmer_letters[i])
        if 'nofilt' in kmer_flags[i]:
            nofilt_letters.append(kmer_letters[i])
        if kmer_letters[i] == UNKNOWN_LETTER:
            print('Error: character A is reserved for unknown sequence')
            exit(1)
    if canonical_letter is None:
        print('Error: cluster_tvrs() received a kmer list that does not have any canonical')
        exit(1)

    #
    # when generating consensus sequence for cluster: in ties, prioritize canonical, deprioritize unknown
    #
    char_score_adj = {canonical_letter:1, UNKNOWN_LETTER:-1}
    noncanon_cheat = (canonical_letter, 3)
    #
    # create color vector
    #
    all_colorvecs = []
    for i in range(n_reads):
        [my_kmer_hits, my_tlen, my_dbta, my_orr, my_rname, my_mapq, my_fastadat] = kmer_dat[i]
        my_letters = [UNKNOWN_LETTER for n in range(my_tlen)]
        for ki in range(len(my_kmer_hits)):
            if len(my_kmer_hits[ki]):
                for kmer_span in my_kmer_hits[ki]:
                    for j in range(kmer_span[0], kmer_span[1]):
                        my_letters[j] = kmer_letters[ki]
        all_colorvecs.append(''.join(my_letters))

    #
    # identify subtel / tvr boundary
    #
    subtel_regions = []
    tvrtel_regions = []
    cleaned_colorvecs = []
    colorvecs_for_msa = []
    err_end_lens      = []
    for i in range(len(all_colorvecs)):
        #
        if pq == 'p':
            # identify subtel/tvr boundary based on density of unknown characters
            (seq_left, seq_right) = find_density_boundary(all_colorvecs[i][::-1], UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
            ####if seq_right.count(UNKNOWN_LETTER) >= TOO_MUCH_UNKNOWN_IN_TVRTEL:
            ####    (seq_left, seq_right) = find_density_boundary(all_colorvecs[i][::-1], UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below', use_lowest_dens=True)
            # denoise tvr+tel section
            seq_right_denoise = denoise_colorvec(seq_right, chars_to_delete=denoise_chars, char_to_merge=canonical_letter)
            # remove ends of reads that might be sequencing artifacts, based on density of canonical characters
            (err_left, err_right) = find_density_boundary(seq_right_denoise[::-1], canonical_letter, CANON_WIN_SIZE, CANON_END_DENS, thresh_dir='above')
            #
            cleaned_colorvecs.append(err_right + seq_left[::-1])    # the entire denoised read (for plotting)
            err_end_lens.append(len(err_left))                      # needed for adjusting offsets in plots
            #
            #colorvecs_for_msa.append(seq_right_denoise[::-1])      # I don't remember why I was using denoised tvr for msa...
            colorvecs_for_msa.append(err_right)
            #
            subtel_regions.append(seq_left[::-1])                   # subtel regions (currently not used for anything)
            tvrtel_regions.append(err_right)                        # tvr sequence used for clustering
            if len(tvrtel_regions[-1]) > tvr_truncate:
                tvrtel_regions[-1] = tvrtel_regions[-1][-tvr_truncate:]
            #cleaned_colorvecs[-1] = tvrtel_regions[-1]
        #
        elif pq == 'q':
            # identify subtel/tvr boundary based on density of unknown characters
            (seq_left, seq_right) = find_density_boundary(all_colorvecs[i], UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
            ####if seq_right.count(UNKNOWN_LETTER) >= TOO_MUCH_UNKNOWN_IN_TVRTEL:
            ####    (seq_left, seq_right) = find_density_boundary(all_colorvecs[i], UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below', use_lowest_dens=True)
            # denoise tvr+tel section
            seq_right_denoise = denoise_colorvec(seq_right, chars_to_delete=denoise_chars, char_to_merge=canonical_letter)
            # remove ends of reads that might be sequencing artifacts, based on density of canonical characters
            (err_left, err_right) = find_density_boundary(seq_right_denoise[::-1], canonical_letter, CANON_WIN_SIZE, CANON_END_DENS, thresh_dir='above')
            #
            cleaned_colorvecs.append(seq_left + err_right[::-1])    # the entire denoised read (for plotting)
            err_end_lens.append(len(err_left))                      # needed for adjusting offsets in plots
            #
            #colorvecs_for_msa.append(seq_right_denoise)            # I don't remember why I was using denoised tvr for msa...
            colorvecs_for_msa.append(err_right[::-1])
            #
            subtel_regions.append(seq_left)                         # subtel regions (currently not used for anything)
            tvrtel_regions.append(err_right[::-1])                  # tvr sequence used for clustering
            if len(tvrtel_regions[-1]) > tvr_truncate:
                tvrtel_regions[-1] = tvrtel_regions[-1][:tvr_truncate]
            #cleaned_colorvecs[-1] = tvrtel_regions[-1]

    #
    # PAIRWISE ALIGNMENT OF ALL SEQUENCES
    #
    if dist_in is None or exists_and_is_nonzero(dist_in) is False:
        all_indices = [[] for n in range(alignment_processes)]
        k = 0
        for i in range(n_reads):
            for j in range(i+1, n_reads):
                all_indices[k % alignment_processes].append((i,j))
                k += 1
        #
        # align tvr + tel regions
        #
        manager     = multiprocessing.Manager()
        tvrtel_dist = manager.dict()
        processes   = []
        for i in range(alignment_processes):
            if aln_mode == 'ms':
                p = multiprocessing.Process(target=parallel_alignment_job,
                                            args=(tvrtel_regions, all_indices[i], pq, tvrtel_dist, None, (True,True), True, True, rand_shuffle_count))
            elif aln_mode == 'ds':
                scoring_matrix = get_scoring_matrix(canonical_letter)
                p = multiprocessing.Process(target=parallel_alignment_job,
                                            args=(tvrtel_regions, all_indices[i], pq, tvrtel_dist, scoring_matrix, (True,True), True, True, rand_shuffle_count))
            processes.append(p)
        for i in range(alignment_processes):
            processes[i].start()
        for i in range(alignment_processes):
            processes[i].join()
        #
        # combine distances for final distance matrix
        #
        dist_matrix = np.zeros((n_reads,n_reads))
        for i in range(n_reads):
            for j in range(i+1,n_reads):
                ij_dist = tvrtel_dist[(i,j)]
                dist_matrix[i,j] = ij_dist
                dist_matrix[j,i] = ij_dist
        dist_norm    = max(np.max(dist_matrix), MIN_MSD)
        dist_matrix /= dist_norm
        #
        if dist_in is not None:
            np.save(dist_in, dist_matrix)
        #
    else:
        dist_matrix = np.load(dist_in, allow_pickle=True)

    #
    # hierarchal clustering + dendrogram plotting
    #
    dist_array = squareform(dist_matrix)
    Zread      = linkage(dist_array, method='ward')
    #
    if fig_name is not None:
        dendrogram_width = max(8, 0.12 * n_reads)
        dendrogram_width = min(dendrogram_width, MAX_PLOT_SIZE / DEFAULT_DPI)
        fig = mpl.figure(3, figsize=(dendrogram_width,6), dpi=DEFAULT_DPI)
        mpl.rcParams.update({'font.size': 20, 'lines.linewidth':2.0})
        dendrogram(Zread, color_threshold=tree_cut)
        mpl.axhline(y=[tree_cut], linestyle='dashed', color='black', alpha=0.7)
        mpl.xlabel('read #')
        mpl.ylabel('distance')
        mpl.title(my_chr + ' : ' + str(my_pos))
        mpl.tight_layout()
        mpl.savefig(fig_name)
        mpl.close(fig)
    #
    assignments = fcluster(Zread, tree_cut, 'distance').tolist()
    by_class = {}
    for i in range(len(assignments)):
        if assignments[i] not in by_class:
            by_class[assignments[i]] = []
        by_class[assignments[i]].append(i)
    out_clust = sorted([(len(by_class[k]), sorted(by_class[k])) for k in by_class.keys()], reverse=True)
    out_clust = [n[1] for n in out_clust]
    for i in range(len(out_clust)): # sort by length
        out_clust[i] = sorted([(kmer_dat[n][1],n) for n in out_clust[i]], reverse=True)
        out_clust[i] = [n[1] for n in out_clust[i]]

    #
    # do MSA of TVRs (and also subtels) to get a consensus sequences
    #
    out_consensus    = []
    subtel_consensus = []
    #
    if save_msa is not None and exists_and_is_nonzero(save_msa):
        my_reader = TG_Reader(save_msa, verbose=False)
        while True:
            read_dat = my_reader.get_next_read()
            if not read_dat[0]:
                break
            if read_dat[0][:3] == 'tvr':
                out_consensus.append(read_dat[1])
            elif read_dat[0][:6] == 'subtel':
                subtel_consensus.append(read_dat[1])
        my_reader.close()
        if len(out_consensus) != len(out_clust):    # msa we read from file has different number of clusters than we currently have, abort!
            if PRINT_DEBUG:
                print('mismatch #clusters from input consensus:', len(out_consensus), '!=', len(out_clust))
            out_consensus = []
    #
    buffered_tvrs = {}
    buffered_subs = {}
    #
    if save_msa is None or len(out_consensus) == 0:
        for i in range(len(out_clust)):
            max_tvr_len = max([len(colorvecs_for_msa[n]) for n in out_clust[i]])
            max_sub_len = max([len(subtel_regions[n]) for n in out_clust[i]])
            for ci in out_clust[i]:
                tvr_buff_seq = canonical_letter*(max_tvr_len - len(colorvecs_for_msa[ci]))
                sub_buff_seq = UNKNOWN_LETTER*(max_sub_len - len(subtel_regions[ci]))
                if pq == 'p':
                    buffered_tvrs[ci] = tvr_buff_seq + colorvecs_for_msa[ci]
                    buffered_subs[ci] = subtel_regions[ci] + sub_buff_seq
                elif pq == 'q':
                    buffered_tvrs[ci] = colorvecs_for_msa[ci] + tvr_buff_seq
                    buffered_subs[ci] = sub_buff_seq + subtel_regions[ci]
        #
        for i in range(len(out_clust)):
            if len(out_clust[i]) == 1:
                out_consensus.append(buffered_tvrs[out_clust[i][0]])
                subtel_consensus.append(buffered_subs[out_clust[i][0]])
            else:
                clust_seq = [buffered_tvrs[n] for n in out_clust[i]]
                [msa_seq, consensus_seq] = get_muscle_msa(clust_seq, muscle_exe, working_dir=muscle_dir, char_score_adj=char_score_adj, noncanon_cheat=noncanon_cheat)
                out_consensus.append(consensus_seq)
                clust_seq = [buffered_subs[n] for n in out_clust[i]]
                if clust_seq[0] == '':
                    subtel_consensus.append('')
                else:
                    [msa_seq, consensus_seq] = get_muscle_msa(clust_seq, muscle_exe, working_dir=muscle_dir)
                    subtel_consensus.append(consensus_seq)
        #
        if save_msa is not None:
            f = open(save_msa,'w')
            for i in range(len(out_consensus)):
                f.write('>tvr-' + str(i+1).zfill(2) + '\n')
                f.write(out_consensus[i] + '\n')
            for i in range(len(subtel_consensus)):
                f.write('>subtel-' + str(i+1).zfill(2) + '\n')
                f.write(subtel_consensus[i] + '\n')
            f.close()

    #
    # check subtel for variant repeats --> add to tvr
    #
    subtel_recovery = []
    for i in range(len(out_consensus)):
        if pq == 'p':
            (scons_left, scons_right) = find_density_boundary(subtel_consensus[i], UNKNOWN_LETTER, 50, 0.800, thresh_dir='above', debug_plot=False)
            out_consensus[i] = out_consensus[i] + scons_left
            subtel_recovery.append(scons_left)
        elif pq == 'q':
            (scons_left, scons_right) = find_density_boundary(subtel_consensus[i][::-1], UNKNOWN_LETTER, 50, 0.800, thresh_dir='above', debug_plot=False)
            out_consensus[i] = scons_left[::-1] + out_consensus[i]
            subtel_recovery.append(scons_left[::-1])

    #####
    ##### prune subtel from consensus
    #####
    ####for i in range(len(out_consensus)):
    ####    if pq == 'p':
    ####        (cons_left, cons_right) = find_density_boundary(out_consensus[i][::-1], UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
    ####        out_consensus[i] = cons_right[::-1]
    ####    elif pq == 'q':
    ####        (cons_left, cons_right) = find_density_boundary(out_consensus[i], UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
    ####        out_consensus[i] = cons_right

    #
    # identify tvr/tel boundary from consensus sequences [STRICT]
    #
    out_tvr_tel_boundaries = []
    for i in range(len(out_consensus)):
        if len(out_consensus[i]) > MAX_TVR_LEN:
            if pq == 'p':
                current_cons = canonical_letter*(len(out_consensus[i])-MAX_TVR_LEN) + out_consensus[i][-MAX_TVR_LEN:]
            elif pq == 'q':
                current_cons = out_consensus[i][:MAX_TVR_LEN] + canonical_letter*(len(out_consensus[i])-MAX_TVR_LEN)
        else:
            current_cons = out_consensus[i]
        denoised_consensus = denoise_colorvec(current_cons,
                                              chars_to_delete=[n for n in tvr_letters if n not in nofilt_letters],
                                              min_size=TVR_CANON_FILT_PARAMS_STRICT[0],
                                              char_to_merge=canonical_letter)
        if pq == 'q':
            denoised_consensus = denoised_consensus[::-1]
        tel_boundary = find_cumulative_boundary(denoised_consensus, tvr_letters,
                                                cum_thresh=TVR_CANON_FILT_PARAMS_STRICT[1],
                                                min_hits=TVR_CANON_FILT_PARAMS_STRICT[2])
        #
        # failed to find tel boundary, try again with [LENIENT] params
        #
        if tel_boundary == len(out_consensus[i])+1:
            if len(out_consensus[i]) > MAX_TVR_LEN_SHORT:
                if pq == 'p':
                    current_cons = canonical_letter*(len(out_consensus[i])-MAX_TVR_LEN_SHORT) + out_consensus[i][-MAX_TVR_LEN_SHORT:]
                elif pq == 'q':
                    current_cons = out_consensus[i][:MAX_TVR_LEN_SHORT] + canonical_letter*(len(out_consensus[i])-MAX_TVR_LEN_SHORT)
            else:
                current_cons = out_consensus[i]
            denoised_consensus = denoise_colorvec(current_cons,
                                                  chars_to_delete=[n for n in tvr_letters if n not in nofilt_letters],
                                                  min_size=TVR_CANON_FILT_PARAMS_LENIENT[0],
                                                  char_to_merge=canonical_letter)
            if pq == 'q':
                denoised_consensus = denoised_consensus[::-1]
            tel_boundary = find_cumulative_boundary(denoised_consensus, tvr_letters,
                                                    cum_thresh=TVR_CANON_FILT_PARAMS_LENIENT[1],
                                                    min_hits=TVR_CANON_FILT_PARAMS_LENIENT[2])
        #
        # if we did find a boundary, buffer slightly to include variant repeats at the edge
        # - then adjust to nearest canonical / variant repeat
        #
        if tel_boundary != len(out_consensus[i])+1:
            tel_boundary = max(0, tel_boundary - TVR_BOUNDARY_BUFF)
            if denoised_consensus[tel_boundary] == canonical_letter:
                while denoised_consensus[tel_boundary] == canonical_letter:
                    tel_boundary += 1
            else:
                while denoised_consensus[tel_boundary] not in [canonical_letter, UNKNOWN_LETTER]:
                    tel_boundary -= 1
                tel_boundary += 1
            #print(denoised_consensus[tel_boundary-100:tel_boundary])
            #print(denoised_consensus[tel_boundary:tel_boundary+100])
            my_telbound = len(out_consensus[i]) - tel_boundary
        else:
            my_telbound = len(out_consensus[i]) - tel_boundary + 1
        out_tvr_tel_boundaries.append(my_telbound)

    #
    # merge clusters of TVRs if they are a prefix of a larger TVR
    # --- update out_clust
    # --- update out_consensus
    # --- update out_tvr_tel_boundaries
    # --- update sub_recovery_adj
    #
    if len(out_clust) > 1:
        if dist_in_prefix is None or exists_and_is_nonzero(dist_in_prefix) is False:
            pref_indices = [[] for n in range(alignment_processes)]
            k = 0
            n_consensus = len(out_consensus)
            for i in range(n_consensus):
                for j in range(i+1, n_consensus):
                    pref_indices[k % alignment_processes].append((i,j))
                    k += 1
            #
            manager   = multiprocessing.Manager()
            pref_dist = manager.dict()
            processes = []
            for i in range(alignment_processes):
                if aln_mode == 'ms':
                    p = multiprocessing.Process(target=parallel_alignment_job,
                                                args=(out_consensus, pref_indices[i], pq, pref_dist, None, (True,True), True, False, rand_shuffle_count))
                elif aln_mode == 'ds':
                    scoring_matrix = get_scoring_matrix(canonical_letter, MATCH_CANON_PREFIX)
                    p = multiprocessing.Process(target=parallel_alignment_job,
                                                args=(out_consensus, pref_indices[i], pq, pref_dist, scoring_matrix, (True,True), True, False, rand_shuffle_count))
                processes.append(p)
            for i in range(alignment_processes):
                processes[i].start()
            for i in range(alignment_processes):
                processes[i].join()
            #
            dist_matrix_prefix = np.zeros((n_consensus, n_consensus))
            for i in range(n_consensus):
                for j in range(i+1, n_consensus):
                    ij_dist = pref_dist[(i,j)]
                    dist_matrix_prefix[i,j] = ij_dist
                    dist_matrix_prefix[j,i] = ij_dist
            dist_norm = max(np.max(dist_matrix_prefix), MIN_MSD)
            dist_matrix_prefix /= dist_norm
            if dist_in_prefix is not None:
                np.save(dist_in_prefix, dist_matrix_prefix)
        else:
            dist_matrix_prefix = np.load(dist_in_prefix, allow_pickle=True)
        #
        dist_array_prefix  = squareform(dist_matrix_prefix)
        Z_prefix           = linkage(dist_array_prefix, method='complete')
        assignments_prefix = fcluster(Z_prefix, tree_cut_prefix, 'distance').tolist()
        merge_clust = [[] for n in range(max(assignments_prefix))]
        for i in range(len(assignments_prefix)):
            merge_clust[assignments_prefix[i]-1].append(i)
        #
        if fig_prefix_name is not None:
            #pref_clust_labels = [','.join([str(n) for n in m]) for m in out_clust]
            pref_clust_labels = [str(len(n)) + ' reads' for n in out_clust]
            dendrogram_width = max(8, 0.12 * len(out_clust))
            dendrogram_width = min(dendrogram_width, MAX_PLOT_SIZE / DEFAULT_DPI)
            fig = mpl.figure(3, figsize=(dendrogram_width,6), dpi=DEFAULT_DPI)
            mpl.rcParams.update({'font.size': 16, 'font.weight':'bold', 'lines.linewidth':2.0})
            dendrogram(Z_prefix, color_threshold=tree_cut_prefix, labels=pref_clust_labels)
            mpl.axhline(y=[tree_cut_prefix], linestyle='dashed', color='black', alpha=0.7)
            mpl.xlabel('cluster #')
            mpl.ylabel('distance')
            mpl.title(my_chr + ' : ' + str(my_pos))
            mpl.tight_layout()
            mpl.savefig(fig_prefix_name)
            mpl.close(fig)
        #
        merge_clust = sorted([(sum([len(out_clust[n]) for n in m]), m) for m in merge_clust], reverse=True) # sort by readcount
        merge_clust = [n[1] for n in merge_clust]
        new_out_clust              = []
        new_out_consensus          = []
        new_out_tvr_tel_boundaries = []
        new_subtel_recovery        = []
        for mc in merge_clust:
            ci_with_longest_tvr = sorted([(len(out_consensus[ci]), ci) for ci in mc])[-1][1]
            new_out_clust.append([])
            for ci in mc:
                for ri in out_clust[ci]:
                    new_out_clust[-1].append((len(colorvecs_for_msa[ri]), ri))  # gonna sort by length
            new_out_clust[-1] = [n[1] for n in sorted(new_out_clust[-1], reverse=True)]
            new_out_consensus.append(out_consensus[ci_with_longest_tvr])
            new_out_tvr_tel_boundaries.append(out_tvr_tel_boundaries[ci_with_longest_tvr])
            new_subtel_recovery.append(subtel_recovery[ci_with_longest_tvr])
        out_clust              = new_out_clust
        out_consensus          = new_out_consensus
        out_tvr_tel_boundaries = new_out_tvr_tel_boundaries
        subtel_recovery        = new_subtel_recovery

    #
    # lets use tvr/subtel boundary as offset instead of msa offset
    #
    all_subtel_lens = [len(subtel_regions[n]) for n in range(len(subtel_regions))]
    longest_subtel  = max(all_subtel_lens)
    out_adj         = []
    for i in range(len(out_clust)):     # this is strangely complicated...
        my_subtel_lens    = [len(subtel_regions[n]) for n in out_clust[i]]
        longest_subtel_cl = max(my_subtel_lens)
        clust_subtel_adj  = longest_subtel - longest_subtel_cl
        if pq == 'p':
            sub_recovery_adj = -len(subtel_recovery[i])
        elif pq == 'q':
            sub_recovery_adj = len(subtel_recovery[i])
        out_adj.append([longest_subtel_cl - n + clust_subtel_adj + sub_recovery_adj for n in my_subtel_lens])

    #
    # output mapping quality as well (not sure what I'll do with this yet, maybe filtering?)
    #
    out_mapq = []
    for n in out_clust:
        out_mapq.append([kmer_dat[m][5] for m in n])
    #
    if PRINT_DEBUG:
        print()
        print('cluster_tvr() results:')
        print(out_clust)
        print(out_mapq)
        print(out_adj)
        print(all_subtel_lens)
        print('len(out_consensus):', len(out_consensus))
        print('len(cleaned_colorvecs):', len(cleaned_colorvecs))
        print(err_end_lens)
        print(out_tvr_tel_boundaries)
    #
    return [out_clust,
            out_mapq,
            out_adj,
            all_subtel_lens,
            out_consensus,
            cleaned_colorvecs,
            err_end_lens,
            out_tvr_tel_boundaries]

#
#
#
def find_density_boundary(sequence, which_letter, win_size, dens_thresh, thresh_dir='below', starting_coord=0, use_lowest_dens=False, debug_plot=False):
    my_unknown = np.zeros(len(sequence))
    for j in range(starting_coord, len(sequence)):
        if sequence[j] == which_letter:
            my_unknown[j] = 1
    my_unknown_cum  = np.cumsum(my_unknown)
    my_unknown_dens = []
    first_pos_below_thresh = None
    pos_with_lowest_dens   = None   # use pos of min dense in the event that we never go below thresh
    if thresh_dir == 'below':
        lowest_dens_thus_far = 1.0
    elif thresh_dir == 'above':
        lowest_dens_thus_far = 0.0
    else:
        print('Error: find_density_boundary thresh_dir must be above or below')
        exit(1)
    for j in range(starting_coord, len(sequence) - win_size):
        my_unknown_dens.append(float(my_unknown_cum[j+win_size] - my_unknown_cum[j]) / win_size)
        if thresh_dir == 'below':
            if first_pos_below_thresh is None and my_unknown_dens[-1] <= dens_thresh:
                first_pos_below_thresh = j
            if my_unknown_dens[-1] < lowest_dens_thus_far:
                lowest_dens_thus_far = my_unknown_dens[-1]
                pos_with_lowest_dens = j
        else:
            if first_pos_below_thresh is None and my_unknown_dens[-1] >= dens_thresh:
                first_pos_below_thresh = j
            if my_unknown_dens[-1] > lowest_dens_thus_far:
                lowest_dens_thus_far = my_unknown_dens[-1]
                pos_with_lowest_dens = j
    #
    if debug_plot:
        fig = mpl.figure(0)
        mpl.plot(list(range(len(my_unknown_dens))), my_unknown_dens)
        mpl.plot([first_pos_below_thresh, first_pos_below_thresh], [0,1], '--r')
        mpl.plot([pos_with_lowest_dens, pos_with_lowest_dens], [0,1], '--b')
        mpl.title(which_letter + ' ' + thresh_dir)
        mpl.show()
        mpl.close(fig)
    #
    if first_pos_below_thresh is None:
        first_pos_below_thresh = pos_with_lowest_dens
    #
    if use_lowest_dens:
        return (sequence[:pos_with_lowest_dens], sequence[pos_with_lowest_dens:])
    else:
        return (sequence[:first_pos_below_thresh], sequence[first_pos_below_thresh:])

#
#
#
def find_cumulative_boundary(sequence, which_letters, cum_thresh=0.05, min_hits=100, debug_plot=False):
    hits = [1*(n in which_letters) for n in sequence]
    hits_cum = np.cumsum(hits)
    if hits_cum[-1] < min_hits: # not enough hits to even bother trying
        return len(sequence)
    #
    hits_cum = hits_cum/hits_cum[-1]
    first_pos_below_thresh = len(sequence)
    for i in range(len(hits_cum)):
        if hits_cum[i] >= cum_thresh:
            first_pos_below_thresh = i
            break
    #
    if debug_plot:
        fig = mpl.figure(0)
        mpl.plot(list(range(len(hits_cum))), hits_cum)
        mpl.plot([first_pos_below_thresh, first_pos_below_thresh], [0,1], '--r')
        mpl.show()
        mpl.close(fig)
    #
    return first_pos_below_thresh

#
#
#
def denoise_colorvec(v, min_size=10, max_gap_fill=50, chars_to_delete=[], char_to_merge=''):
    blocks = []
    current_block = v[0]
    current_start = 0
    v += UNKNOWN_LETTER
    for i in range(1,len(v)):
        if v[i] != current_block:
            if current_block != UNKNOWN_LETTER:
                blocks.append((current_start, i, current_block))
            current_block = v[i]
            current_start = i
    #
    del_list = []
    for i in range(len(blocks)):
        if blocks[i][1] - blocks[i][0] < min_size and blocks[i][2] in chars_to_delete:
            del_list.append(i)
    del_list = sorted(del_list, reverse=True)
    for di in del_list:
        del blocks[di]
    #
    for i in range(len(blocks)-1,0,-1):
        if blocks[i][2] == blocks[i-1][2] and blocks[i][2] == char_to_merge:
            our_gap = blocks[i][0] - blocks[i-1][1]
            if our_gap <= max_gap_fill:
                blocks[i-1] = (blocks[i-1][0], blocks[i][1], blocks[i][2])
                del blocks[i]
    #
    v_out = [UNKNOWN_LETTER for n in v]
    for block in blocks:
        for i in range(block[0],block[1]):
            v_out[i] = block[2]
    return ''.join(v_out)

#
#
#
def convert_colorvec_to_kmerhits(colorvecs, repeats_metadata):
    #
    [kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
    #
    amino_2_kmer_ind = {}
    for i in range(len(kmer_colors)):
        amino_2_kmer_ind[kmer_letters[i]] = i
    #
    out_kmerhits = []
    for i in range(len(colorvecs)):
        current_block = colorvecs[i][0]
        current_start = 0
        out_kmerhits.append([[] for n in range(len(kmer_colors))])
        colorvecs[i] += UNKNOWN_LETTER
        for j in range(1,len(colorvecs[i])):
            if colorvecs[i][j] != current_block:
                if current_block != UNKNOWN_LETTER:
                    my_ind = amino_2_kmer_ind[current_block]
                    out_kmerhits[-1][my_ind].append((current_start, j))
                current_block = colorvecs[i][j]
                current_start = j
    return out_kmerhits

#
# helper function for quickly plotting a set of tvrs (assumes Q orientation)
#
def plot_some_tvrs(tvrs, labels, repeats_metadata, plot_fn, custom_plot_params={}):
    [KMER_LIST, KMER_COLORS, KMER_LETTER, KMER_FLAGS] = repeats_metadata
    redrawn_tvrs = convert_colorvec_to_kmerhits(tvrs, repeats_metadata)
    clust_khd = []
    for i in range(len(redrawn_tvrs)):
        clust_khd.append([redrawn_tvrs[i], len(tvrs[i]), 0, 'FWD', labels[i], 60, None])
    clust_dat = []
    clust_dat.append([list(range(len(clust_khd)))])
    clust_dat.append([[60]*len(clust_dat[0][0])])
    clust_dat.append([[0]*len(clust_dat[0][0])])
    clust_dat.append([0])
    plot_kmer_hits(clust_khd, KMER_COLORS, '', 0, plot_fn, clust_dat=clust_dat, plot_params=custom_plot_params)

#
#
#
def cluster_consensus_tvrs(sequences,
                           repeats_metadata,
                           tree_cut,
                           dist_in=None,
                           fig_name=None,
                           dendro_name=None,
                           samp_labels=None,
                           aln_mode='ms',
                           gap_bool=(True,True),
                           adjust_lens=False,
                           rand_shuffle_count=3,
                           linkage_method='complete',
                           normalize_dist_matrix=True,
                           alignment_processes=8,
                           job=(1,1),
                           dendrogram_title=None,
                           dendrogram_height=12,
                           overwrite_figures=True):
    #
    n_seq = len(sequences)
    #
    [kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
    #
    canonical_letter = None
    for i in range(len(kmer_list)):
        if 'canonical' in kmer_flags[i]:
            canonical_letter = kmer_letters[i]
        if kmer_letters[i] == UNKNOWN_LETTER:
            print('Error: character A is reserved for unknown sequence')
            exit(1)
    if canonical_letter is None:
        print('Error: cluster_consensus_tvr() received a kmer list that does not have any canonical')
        exit(1)
    #
    if dist_in is None or exists_and_is_nonzero(dist_in) is False:
        dist_matrix = np.zeros((n_seq,n_seq))
        all_indices = [[] for n in range(alignment_processes)]
        k = 0
        for i in range(n_seq):
            for j in range(i+1,n_seq):
                all_indices[k % alignment_processes].append((i,j))
                k += 1
        #
        #   even more parallelization! Any problem can be solved by throwing tons of CPU at it.
        #
        if job[1] > 1:
            my_job = job[0] - 1
            chunks = job[1]
            for i in range(alignment_processes):
                chunksize = int(len(all_indices[i])/chunks)
                chunks_by_job = []
                for j in range(chunks):
                    if j == chunks-1:
                        chunks_by_job.append(all_indices[i][j*chunksize:])
                    else:
                        chunks_by_job.append(all_indices[i][j*chunksize:(j+1)*chunksize])
                all_indices[i] = [n for n in chunks_by_job[my_job]]
        #
        manager     = multiprocessing.Manager()
        return_dict = manager.dict()
        processes   = []
        for i in range(alignment_processes):
            if aln_mode == 'ms':
                p = multiprocessing.Process(target=parallel_alignment_job,
                                            args=(sequences, all_indices[i], 'q', return_dict, None, gap_bool, adjust_lens, False, rand_shuffle_count))
            elif aln_mode == 'ds':
                scoring_matrix = get_scoring_matrix(canonical_letter)
                p = multiprocessing.Process(target=parallel_alignment_job,
                                            args=(sequences, all_indices[i], 'q', return_dict, scoring_matrix, gap_bool, adjust_lens, False, rand_shuffle_count))
            processes.append(p)
        for i in range(alignment_processes):
            processes[i].start()
        for i in range(alignment_processes):
            processes[i].join()
        #
        dist_matrix = np.zeros((n_seq,n_seq))
        for (i,j) in return_dict.keys():
            dist_matrix[i,j] = return_dict[(i,j)]
            dist_matrix[j,i] = return_dict[(i,j)]
        if job[1] > 1:
            partial_dist_fn = dist_in[:-4] + '_job' + str(job[0]).zfill(3) + '.npy'
            np.save(partial_dist_fn, dist_matrix)
        else:
            if normalize_dist_matrix:
                dist_norm    = max(np.max(dist_matrix), MIN_MSD)
                dist_matrix /= dist_norm
            if dist_in is not None:
                np.save(dist_in, dist_matrix)
    else:
        dist_matrix = np.load(dist_in, allow_pickle=True)
    #
    if samp_labels is None:
        samp_labels = [str(n) for n in range(len(sequences))]
    #
    if job[1] == 1 or job[0] == 0:
        d_arr = squareform(dist_matrix)
        Zread = linkage(d_arr, method=linkage_method)
        #
        mpl.rcParams.update({'font.size': 16, 'font.weight':'bold'})
        fig = mpl.figure(figsize=(8,dendrogram_height))
        dendro_dat = dendrogram(Zread, orientation='left', labels=samp_labels, color_threshold=tree_cut)
        labels_fromtop = dendro_dat['ivl'][::-1]
        # ugliness to keep dendrogram ordering consistent with other figures
        reorder_map = {n:samp_labels.index(m) for n,m in enumerate(labels_fromtop)}
        replot_labels = [f'({labels_fromtop.index(n)}) {n}' for n in samp_labels]
        mpl.close(fig)
        fig = mpl.figure(figsize=(8,dendrogram_height))
        dendrogram(Zread, orientation='left', labels=replot_labels, color_threshold=tree_cut)
        #
        mpl.axvline(x=[tree_cut], linestyle='dashed', color='black', alpha=0.7)
        mpl.xlabel('distance')
        if dendrogram_title is not None:
            mpl.title(dendrogram_title)
        mpl.tight_layout()
        if dendro_name is not None:
            mpl.savefig(dendro_name, dpi=200)  # default figure dpi = 100
        mpl.close(fig)
        #
        if fig_name is not None:
            if exists_and_is_nonzero(fig_name) is False or overwrite_figures is True:
                reordered_sequences = [sequences[reorder_map[n]] for n in range(len(sequences))]
                plot_some_tvrs(reordered_sequences, labels_fromtop, repeats_metadata, fig_name, custom_plot_params={'custom_title':''})
        #
        assignments = fcluster(Zread, tree_cut, 'distance').tolist()
        by_class = {}
        for i in range(len(assignments)):
            if assignments[i] not in by_class:
                by_class[assignments[i]] = []
            by_class[assignments[i]].append(i)
        out_clust = sorted([(len(by_class[k]), sorted(by_class[k])) for k in by_class.keys()], reverse=True)
        out_clust = [n[1] for n in out_clust]
        #
        return out_clust
    #
    return None

#
#
#
def make_tvr_plots(kmer_hit_dat, read_clust_dat, my_chr, my_pos, telcompplot_fn, telcompcons_fn, mtp_params):
    #
    [KMER_METADATA, KMER_COLORS, MIN_READS_PER_PHASE, PLOT_FILT_CVECS, DUMMY_TEL_MAPQ, DO_NOT_OVERWRITE] = mtp_params
    #
    # adjust kmer_hit_dat based on the filters and etc that were applied during clustering
    #
    my_consensus_vecs = read_clust_dat[4]
    my_color_vectors  = read_clust_dat[5]
    my_end_err_lens   = read_clust_dat[6]
    my_tvr_tel_bounds = read_clust_dat[7]
    redrawn_consensus = convert_colorvec_to_kmerhits(my_consensus_vecs, KMER_METADATA)
    # do we want to plot denoised tvrs of individual reads?
    if PLOT_FILT_CVECS:
        redrawn_kmerhits = convert_colorvec_to_kmerhits(my_color_vectors, KMER_METADATA)
        for rdki in range(len(redrawn_kmerhits)):
            kmer_hit_dat[rdki][0]  = redrawn_kmerhits[rdki] # replace kmer hit tuples for plotting
            kmer_hit_dat[rdki][1] -= my_end_err_lens[rdki]  # subtract size of artifacts at end of reads
    #
    consensus_kmer_hit_dat = []
    consensus_clust_dat    = [[],[],[],[0]] # fake data so that plot_kmer_hits doesn't complain
    consensus_tvr_tel_pos  = []
    for rdki in range(len(redrawn_consensus)):
        cons_readcount = len(read_clust_dat[0][rdki])
        cons_readname  = 'consensus-' + str(rdki) + ' [' + str(cons_readcount)
        cons_tvrlen    = my_tvr_tel_bounds[rdki]
        if cons_readcount == 1:
            cons_readname += ' read]'
        else:
            cons_readname += ' reads]'
        if cons_readcount >= MIN_READS_PER_PHASE:
            consensus_kmer_hit_dat.append([redrawn_consensus[rdki], len(my_consensus_vecs[rdki]), 0, 'FWD', cons_readname, DUMMY_TEL_MAPQ, None])
            consensus_clust_dat[0].append([rdki])
            consensus_clust_dat[1].append([DUMMY_TEL_MAPQ])
            consensus_clust_dat[2].append([0])
            consensus_tvr_tel_pos.append(cons_tvrlen)
    #
    # TVR plotting (clustered reads + consensus for each allele)
    #
    custom_plot_params = {'xlim':[-1000,15000]}
    if telcompplot_fn is not None:
        if DO_NOT_OVERWRITE is False or exists_and_is_nonzero(telcompplot_fn) is False:
            plot_kmer_hits(kmer_hit_dat, KMER_COLORS, my_chr, my_pos, telcompplot_fn,
                           clust_dat=read_clust_dat,
                           plot_params=custom_plot_params)
    if telcompcons_fn is not None and len(consensus_clust_dat[0]):
        if DO_NOT_OVERWRITE is False or exists_and_is_nonzero(telcompcons_fn) is False:
            plot_kmer_hits(consensus_kmer_hit_dat, KMER_COLORS, my_chr, my_pos, telcompcons_fn,
                           clust_dat=consensus_clust_dat,
                           draw_boundaries=consensus_tvr_tel_pos,
                           plot_params=custom_plot_params)
