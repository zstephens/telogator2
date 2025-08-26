import numpy as np
import matplotlib.pyplot as mpl

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance  import squareform

from source.tg_align  import get_aligner_object, get_dist_matrix_parallel, get_final_tvr_consensus, get_scoring_matrix, iterative_refinement, MAX_MSA_READCOUNT, MAX_SEQ_DIST, progressive_alignment, UNKNOWN
from source.tg_plot   import DEFAULT_DPI, MAX_PLOT_SIZE, plot_some_tvrs
from source.tg_reader import TG_Reader
from source.tg_util   import exists_and_is_nonzero, UNCLUST_CHR, UNCLUST_POS


MAX_TVR_LEN       = 8000    # ignore variant repeats past this point when finding TVR boundary
MAX_TVR_LEN_SHORT = 3000    # when examining TVRs with very few variant repeats
TVR_BOUNDARY_BUFF = 10      # add this many bp to detected TVR boundary
TVR_BOUNDARY_ADJ_STEP = 50  # how far are we willing to look for more variant repeats to extend TVR boundary?

# if tvr + tel region has at least this many unknown characters, use minimum-dens pos instead of first-below-dens ps
TOO_MUCH_UNKNOWN_IN_TVRTEL = 1000
# density parameters for identifing subtel / tvr boundaries (individual reads)
UNKNOWN_WIN_SIZE = 100
UNKNOWN_END_DENS = 0.125
UNKNOWN_WIN_SIZE_FINE = 50
UNKNOWN_END_DENS_FINE = 0.500
# density parameters for discerning canonical regions from sequencing artifacts
CANON_WIN_SIZE = 100
CANON_END_DENS = 0.700
# parameters for determining tvr / canonical boundaries: (denoise_region_size, cum_thresh, min_hits)
TVR_CANON_FILT_PARAMS_STRICT  = (10, 0.125, 100)
TVR_CANON_FILT_PARAMS_LENIENT = ( 5, 0.250,  50)


def cluster_tvrs(kmer_dat,
                 repeats_metadata,
                 tree_cut,
                 my_chr=UNCLUST_CHR,
                 my_pos=UNCLUST_POS,
                 aln_mode='ds',
                 rand_shuffle_count=3,
                 dist_in=None,
                 fig_name=None,
                 save_msa=None,
                 tvr_truncate=3000,
                 num_processes=4,
                 clustering_only=False,
                 print_matrix_progress=False):
    #
    #   kmer_dat[i] = [[kmer1_hits, kmer2_hits, ...], tlen, tel_anchor_dist, read_orientation, readname, anchor_mapq, fasta_dat]
    #
    #   repeats_metadata = [kmer_list, kmer_colors, kmer_letters, kmer_flags]
    #
    [kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
    #
    n_reads = len(kmer_dat)
    pq = my_chr[-1]
    if pq != 'q':
        print('Error: cluster_tvrs() now only supports telomeres in q configuration')
        exit(1)
    #
    canonical_letter = None
    denoise_letters  = []
    tvr_letters      = []
    nofilt_letters   = []
    dubious_letters  = [UNKNOWN]
    subfilt_letters  = []
    for i in range(len(kmer_list)):
        if 'canonical' in kmer_flags[i]:
            canonical_letter = kmer_letters[i]
        if 'denoise' in kmer_flags[i]:
            denoise_letters.append(kmer_letters[i])
        if 'tvr' in kmer_flags[i]:
            tvr_letters.append(kmer_letters[i])
        if 'nofilt' in kmer_flags[i]:
            nofilt_letters.append(kmer_letters[i])
        if 'dubious' in kmer_flags[i]:
            dubious_letters.append(kmer_letters[i])
        if 'subtel-filt' in kmer_flags[i]:
            subfilt_letters.append(kmer_letters[i])
        if kmer_letters[i] == UNKNOWN:
            print(f'Error: character {UNKNOWN} is reserved for unknown sequence')
            exit(1)
    if canonical_letter is None:
        print('Error: cluster_tvrs() received a kmer list that does not have any kmers marked as canonical')
        exit(1)
    denoise_letters = list(set(denoise_letters))
    tvr_letters     = list(set(tvr_letters))
    nofilt_letters  = list(set(nofilt_letters))
    dubious_letters = list(set(dubious_letters))
    subfilt_letters = list(set(subfilt_letters))
    letters_worth_chasing = [n for n in tvr_letters if n not in dubious_letters]

    #
    # create color vector
    #
    all_colorvecs = []
    for i in range(n_reads):
        [my_kmer_hits, my_tlen, my_dbta, my_orr, my_rname, my_mapq, my_fastadat] = kmer_dat[i]
        my_letters = [UNKNOWN for n in range(my_tlen)]
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
    tvrtels_for_clustering = []
    cleaned_colorvecs = []
    colorvecs_for_msa = []
    err_end_lens      = []
    for i in range(len(all_colorvecs)):
        # identify subtel/tvr boundary based on density of unknown characters
        (seq_left, seq_right) = find_density_boundary(all_colorvecs[i], [UNKNOWN]+subfilt_letters, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below', debug_plot=False, readname=kmer_dat[i][4])
        # lets walk backwards with a smaller window to see if there are some tvrs we missed
        (recovered_tvr, rest_is_subtel) = find_density_boundary(seq_left[::-1], [UNKNOWN]+subfilt_letters, UNKNOWN_WIN_SIZE_FINE, UNKNOWN_END_DENS_FINE, thresh_dir='above', debug_plot=False, readname=kmer_dat[i][4])
        seq_right = recovered_tvr[::-1] + seq_right
        seq_left = rest_is_subtel[::-1]
        # adjust subtel boundary so that tvr does not begin with an unknown character
        adj = 0
        while seq_right[adj] == UNKNOWN and adj < len(seq_right)-1:
            adj += 1
        seq_left = seq_left + seq_right[:adj]
        seq_right = seq_right[adj:]
        # denoise tvr+tel section
        seq_right_denoise = denoise_colorvec(seq_right, replace_char=canonical_letter, chars_to_delete=denoise_letters, chars_to_merge=[canonical_letter]+tvr_letters)
        # remove ends of reads that might be sequencing artifacts, based on density of canonical characters
        (err_left, err_right) = find_density_boundary(seq_right_denoise[::-1], [canonical_letter], CANON_WIN_SIZE, CANON_END_DENS, thresh_dir='above')
        #
        cleaned_colorvecs.append(seq_left + err_right[::-1])    # the entire denoised read (for plotting)
        err_end_lens.append(len(err_left))                      # needed for adjusting offsets in plots
        colorvecs_for_msa.append(err_right[::-1])
        #
        subtel_regions.append(seq_left)                         # subtel regions
        tvrtel_regions.append(seq_right)                        # tvrtel regions
        tvrtels_for_clustering.append(err_right[::-1])          # tvr sequence used for clustering
        if len(tvrtels_for_clustering[-1]) > tvr_truncate:
            tvrtels_for_clustering[-1] = tvrtels_for_clustering[-1][:tvr_truncate]

    #
    # PAIRWISE ALIGNMENT OF ALL SEQUENCES
    #
    if dist_in is None or exists_and_is_nonzero(dist_in) is False:
        if aln_mode == 'ms':
            scoring_matrix = None
        elif aln_mode == 'ds':
            scoring_matrix = get_scoring_matrix(canonical_letter, which_type='tvr')
        aligner = get_aligner_object(scoring_matrix=scoring_matrix, gap_bool=(True, True), which_type='tvr')
        dist_matrix = get_dist_matrix_parallel(tvrtels_for_clustering, aligner, True, True, rand_shuffle_count, max_workers=num_processes, print_progress=print_matrix_progress)
        dist_matrix /= MAX_SEQ_DIST
        if dist_in is not None:
            np.savez_compressed(dist_in, dist=dist_matrix)
    #
    else:
        dist_matrix = np.load(dist_in)['dist']

    #
    # hierarchal clustering + dendrogram plotting
    #
    out_clust = [list(range(len(tvrtels_for_clustering)))]
    if len(tvrtels_for_clustering) > 1:
        dist_array = squareform(dist_matrix)
        Zread = linkage(dist_array, method='ward')
        max_distance = np.max(Zread[:,2])
        if max_distance > 1:
            Zread[:,2] /= max_distance # normalize linkage distances so that dendrogram height is 1.
        #
        if fig_name is not None:
            dendrogram_width = max(8, 0.12 * n_reads)
            dendrogram_width = min(dendrogram_width, MAX_PLOT_SIZE / DEFAULT_DPI)
            mpl.rcParams.update({'font.size': 20, 'lines.linewidth':2.0})
            fig = mpl.figure(3, figsize=(dendrogram_width,6), dpi=DEFAULT_DPI, layout="constrained")
            dendrogram(Zread, color_threshold=tree_cut)
            mpl.axhline(y=[tree_cut], linestyle='dashed', color='black', alpha=0.7)
            mpl.xlabel('read #')
            mpl.ylabel('distance')
            mpl.title(my_chr + ' : ' + str(my_pos))
            mpl.savefig(fig_name)
            mpl.close(fig)
            mpl.rcdefaults()
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
    # if all we want is the cluster assignments, we can stop now
    #
    if clustering_only:
        return [out_clust, None, None, None, None, None, None, None]

    #
    # do MSA of TVRs (and also subtels) to get a consensus sequences
    #
    out_consensus = []
    if save_msa is not None and exists_and_is_nonzero(save_msa):
        my_reader = TG_Reader(save_msa, verbose=False)
        while True:
            read_dat = my_reader.get_next_read()
            if not read_dat[0]:
                break
            if read_dat[0][:3] == 'tvr':
                out_consensus.append(read_dat[1])
            elif read_dat[0][:6] == 'subtel':
                pass # only keeping this around for reprocessing previous runs that used to do subtel consensuses
        my_reader.close()
        if len(out_consensus) != len(out_clust):
            # msa we read from file has different number of clusters than we currently have, abort!
            out_consensus = []
    #
    if save_msa is None or len(out_consensus) == 0:
        msa_matrix = get_scoring_matrix(canonical_letter, which_type='msa')
        refine_matrix = get_scoring_matrix(canonical_letter, which_type='msa_refinement')
        msa_aligner = get_aligner_object(scoring_matrix=msa_matrix, gap_bool=(True,False), which_type='msa')
        refinement_aligner = get_aligner_object(scoring_matrix=refine_matrix, gap_bool=(True,False), which_type='msa')
        for i in range(len(out_clust)):
            clust_inds = [n[1] for n in sorted([(len(colorvecs_for_msa[ci]), ci) for ci in out_clust[i]], reverse=True)[:MAX_MSA_READCOUNT]]
            my_dists = dist_matrix[np.ix_(clust_inds, clust_inds)]
            my_seqs_tvr = [colorvecs_for_msa[ci] for ci in clust_inds]
            initial_msa = progressive_alignment(my_seqs_tvr, my_dists, msa_aligner, num_processes=num_processes)
            refined_msa = iterative_refinement(initial_msa, refinement_aligner)
            out_consensus.append(get_final_tvr_consensus(refined_msa,
                                                         default_char=canonical_letter,
                                                         untrustworthy_chars=[canonical_letter]+dubious_letters,
                                                         tiebreak_adj={canonical_letter:1, UNKNOWN:-1}))
        if save_msa is not None:
            f = open(save_msa,'w')
            for i in range(len(out_consensus)):
                f.write('>tvr-' + str(i+1).zfill(3) + '\n')
                f.write(out_consensus[i] + '\n')
            f.close()

    #
    # identify tvr/tel boundary from consensus sequences [STRICT]
    #
    out_tvr_tel_boundaries = []
    for i in range(len(out_consensus)):
        if len(out_consensus[i]) > MAX_TVR_LEN:
            current_cons = out_consensus[i][:MAX_TVR_LEN] + canonical_letter*(len(out_consensus[i])-MAX_TVR_LEN)
        else:
            current_cons = out_consensus[i]
        denoised_consensus = denoise_colorvec(current_cons,
                                              replace_char=canonical_letter,
                                              chars_to_delete=[n for n in tvr_letters if n not in nofilt_letters],
                                              min_size=TVR_CANON_FILT_PARAMS_STRICT[0],
                                              chars_to_merge=[canonical_letter]+tvr_letters)
        denoised_consensus = denoised_consensus[::-1]
        tel_boundary = find_cumulative_boundary(denoised_consensus,tvr_letters,
                                                cum_thresh=TVR_CANON_FILT_PARAMS_STRICT[1],
                                                min_hits=TVR_CANON_FILT_PARAMS_STRICT[2])
        #
        # failed to find tel boundary, try again with [LENIENT] params
        #
        if tel_boundary == len(out_consensus[i])+1:
            if len(out_consensus[i]) > MAX_TVR_LEN_SHORT:
                current_cons = out_consensus[i][:MAX_TVR_LEN_SHORT] + canonical_letter*(len(out_consensus[i])-MAX_TVR_LEN_SHORT)
            else:
                current_cons = out_consensus[i]
            denoised_consensus = denoise_colorvec(current_cons,
                                                  replace_char=canonical_letter,
                                                  chars_to_delete=[n for n in tvr_letters if n not in nofilt_letters],
                                                  min_size=TVR_CANON_FILT_PARAMS_LENIENT[0],
                                                  chars_to_merge=[canonical_letter]+tvr_letters)
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
            if tel_boundary > TVR_BOUNDARY_ADJ_STEP+1:  # do we actually have room to adjust it?
                found_letters_ahead = -1
                for j in range(0,TVR_BOUNDARY_ADJ_STEP+1):
                    if denoised_consensus[tel_boundary-j] in letters_worth_chasing:
                        found_letters_ahead = j
                if found_letters_ahead >= 0:
                    tel_boundary -= found_letters_ahead
                    while denoised_consensus[tel_boundary] in letters_worth_chasing:
                        if tel_boundary <= 0:
                            tel_boundary = 0
                            break
                        found_new_pos = False
                        for j in range(1,TVR_BOUNDARY_ADJ_STEP+1):
                            if denoised_consensus[tel_boundary-j] in letters_worth_chasing:
                                tel_boundary -= j
                                found_new_pos = True
                                break
                        if found_new_pos is False:
                            break
                else:
                    while denoised_consensus[tel_boundary] == canonical_letter:
                        tel_boundary += 1
            my_telbound = len(out_consensus[i]) - tel_boundary
        else:
            my_telbound = len(out_consensus[i]) - tel_boundary + 1
        out_tvr_tel_boundaries.append(my_telbound)

    #
    # lets use tvr/subtel boundary as offset instead of msa offset
    #
    all_subtel_lens = [len(subtel_regions[n]) for n in range(len(subtel_regions))]
    all_tvrtel_lens = [len(tvrtel_regions[n]) for n in range(len(tvrtel_regions))]
    longest_subtel  = max(all_subtel_lens)
    out_adj         = []
    for i in range(len(out_clust)):     # this is strangely complicated...
        my_subtel_lens    = [len(subtel_regions[n]) for n in out_clust[i]]
        longest_subtel_cl = max(my_subtel_lens)
        clust_subtel_adj  = longest_subtel - longest_subtel_cl
        out_adj.append([longest_subtel_cl - n + clust_subtel_adj for n in my_subtel_lens])

    return [out_clust,
            all_tvrtel_lens,
            out_adj,
            all_subtel_lens,
            out_consensus,
            cleaned_colorvecs,
            err_end_lens,
            out_tvr_tel_boundaries]


def cluster_consensus_tvrs(sequences,
                           repeats_metadata,
                           tree_cut,
                           dist_in=None,
                           fig_name=None,
                           fig_plot_params={},
                           dendro_name=None,
                           samp_labels=None,
                           aln_mode='ms',
                           gap_bool=(True,True),
                           adjust_lens=False,
                           tvr_truncate=-1,
                           rand_shuffle_count=3,
                           linkage_method='complete',
                           normalize_dist_matrix=True,
                           num_processes=8,
                           dendrogram_title=None,
                           dendrogram_height=12,
                           dendrogram_allblack=False,
                           dendrogram_xlim=None,
                           overwrite_figures=True,
                           print_matrix_progress=False):
    #
    n_seq = len(sequences)
    #
    [kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
    #
    canonical_letter = None
    for i in range(len(kmer_list)):
        if 'canonical' in kmer_flags[i]:
            canonical_letter = kmer_letters[i]
        if kmer_letters[i] == UNKNOWN:
            print('Error: character A is reserved for unknown sequence')
            exit(1)
    if canonical_letter is None:
        print('Error: cluster_consensus_tvr() received a kmer list that does not have any canonical')
        exit(1)
    #
    if dist_in is None or exists_and_is_nonzero(dist_in) is False:
        if aln_mode == 'ms':
            scoring_matrix = None
        elif aln_mode == 'ds':
            scoring_matrix = get_scoring_matrix(canonical_letter, which_type='consensus')
        aligner = get_aligner_object(scoring_matrix=scoring_matrix, gap_bool=gap_bool, which_type='tvr')
        sequences_for_dist = sequences
        if tvr_truncate > 0:
            sequences_for_dist = [n[:tvr_truncate] for n in sequences]
        dist_matrix = get_dist_matrix_parallel(sequences_for_dist, aligner, adjust_lens, False, rand_shuffle_count, max_workers=num_processes, print_progress=print_matrix_progress)
        if normalize_dist_matrix:
            dist_matrix /= MAX_SEQ_DIST
        if dist_in is not None:
            np.savez_compressed(dist_in, dist=dist_matrix)
    else:
        dist_matrix = np.load(dist_in)['dist']
    #
    if samp_labels is None:
        samp_labels = [str(n) for n in range(len(sequences))]
    #
    out_clust = [list(range(n_seq))]
    if n_seq > 1:
        d_arr = squareform(dist_matrix)
        Zread = linkage(d_arr, method=linkage_method)
        max_distance = np.max(Zread[:,2])
        if max_distance > 1:
            Zread[:,2] /= max_distance # normalize linkage distances so that dendrogram height is 1.
        #
        mpl.rcParams.update({'font.size': 16, 'font.weight':'bold'})
        #
        dendrogram_height = min(dendrogram_height, MAX_PLOT_SIZE / DEFAULT_DPI)
        fig = mpl.figure(figsize=(8,dendrogram_height), dpi=DEFAULT_DPI, layout="constrained")
        dendro_dat = dendrogram(Zread, orientation='left', labels=samp_labels, color_threshold=tree_cut)
        labels_fromtop = dendro_dat['ivl'][::-1]
        # ugliness to keep dendrogram ordering consistent with other figures
        reorder_map = {n:samp_labels.index(m) for n,m in enumerate(labels_fromtop)}
        replot_labels = [f'({labels_fromtop.index(n)}) {n}' for n in samp_labels]
        mpl.close(fig)
        #
        fig = mpl.figure(figsize=(8,dendrogram_height))
        if dendrogram_allblack:
            dendrogram(Zread, orientation='left', labels=replot_labels, color_threshold=tree_cut, above_threshold_color='black', link_color_func=lambda k:'black')
        else:
            dendrogram(Zread, orientation='left', labels=replot_labels, color_threshold=tree_cut)
        if dendrogram_xlim is not None:
            mpl.xlim(dendrogram_xlim[0], dendrogram_xlim[1])
        #
        mpl.axvline(x=[tree_cut], linestyle='dashed', color='black', alpha=0.7)
        mpl.xlabel('distance')
        if dendrogram_title is not None:
            mpl.title(dendrogram_title)
        #
        if dendro_name is not None and (exists_and_is_nonzero(dendro_name) is False or overwrite_figures is True):
            mpl.savefig(dendro_name, dpi=200)  # default figure dpi = 100
        mpl.close(fig)
        mpl.rcdefaults()
        #
        if fig_name is not None:
            if exists_and_is_nonzero(fig_name) is False or overwrite_figures is True:
                reordered_sequences = [sequences[reorder_map[n]] for n in range(len(sequences))]
                plot_some_tvrs(reordered_sequences, labels_fromtop, repeats_metadata, fig_name, custom_plot_params=fig_plot_params)
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


def find_density_boundary(sequence, which_letter, win_size, dens_thresh, thresh_dir='below', starting_coord=0, use_lowest_dens=False, debug_plot=False, readname=''):
    my_unknown = np.zeros(len(sequence))
    for j in range(starting_coord, len(sequence)):
        if sequence[j] in which_letter:
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
        print('Error: find_density_boundary() thresh_dir must be above or below')
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
        fig = mpl.figure(0,figsize=(12,3))
        mpl.plot(list(range(len(my_unknown_dens))), my_unknown_dens)
        mpl.plot([first_pos_below_thresh, first_pos_below_thresh], [0,1], '--r')
        mpl.plot([pos_with_lowest_dens, pos_with_lowest_dens], [0,1], '--b')
        mpl.yticks([0,dens_thresh,1],['0',f'{dens_thresh:.3f}','1'])
        mpl.grid(which='both', linestyle='--', alpha=0.6)
        mpl.axis([0,2000,0,1])
        mpl.title(readname + ' ' + which_letter + ' ' + thresh_dir)
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


def denoise_colorvec(v, replace_char=UNKNOWN, min_size=10, max_gap_fill=50, chars_to_delete=[], chars_to_merge=[]):
    if len(v) == 0:
        return ''
    blocks = []
    current_block = v[0]
    current_start = 0
    v += replace_char
    for i in range(1,len(v)):
        if v[i] != current_block:
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
        if blocks[i][2] == blocks[i-1][2] and blocks[i][2] in chars_to_merge:
            our_gap = blocks[i][0] - blocks[i-1][1]
            if our_gap <= max_gap_fill:
                blocks[i-1] = (blocks[i-1][0], blocks[i][1], blocks[i][2])
                del blocks[i]
    #
    v_out = [replace_char for n in v]
    for block in blocks:
        for i in range(block[0],block[1]):
            v_out[i] = block[2]
    return ''.join(v_out)


def quick_get_tvrtel_lens(kmer_dat, repeats_metadata):
    [kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
    subfilt_letters  = []
    for i in range(len(kmer_list)):
        if 'subtel-filt' in kmer_flags[i]:
            subfilt_letters.append(kmer_letters[i])
    subfilt_letters = list(set(subfilt_letters))
    #
    all_colorvecs = []
    for i in range(len(kmer_dat)):
        [my_kmer_hits, my_tlen, my_dbta, my_orr, my_rname, my_mapq, my_fastadat] = kmer_dat[i]
        my_letters = [UNKNOWN for n in range(my_tlen)]
        for ki in range(len(my_kmer_hits)):
            if len(my_kmer_hits[ki]):
                for kmer_span in my_kmer_hits[ki]:
                    for j in range(kmer_span[0], kmer_span[1]):
                        my_letters[j] = kmer_letters[ki]
        all_colorvecs.append(''.join(my_letters))
    #
    all_tvrtel_lens = []
    for i in range(len(all_colorvecs)):
        # identify subtel/tvr boundary based on density of unknown characters
        (seq_left, seq_right) = find_density_boundary(all_colorvecs[i], [UNKNOWN]+subfilt_letters, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below', debug_plot=False, readname=kmer_dat[i][4])
        # lets walk backwards with a smaller window to see if there are some tvrs we missed
        (recovered_tvr, rest_is_subtel) = find_density_boundary(seq_left[::-1], [UNKNOWN]+subfilt_letters, UNKNOWN_WIN_SIZE_FINE, UNKNOWN_END_DENS_FINE, thresh_dir='above', debug_plot=False, readname=kmer_dat[i][4])
        seq_right = recovered_tvr[::-1] + seq_right
        seq_left = rest_is_subtel[::-1]
        # adjust subtel boundary so that tvr does not begin with an unknown character
        adj = 0
        while seq_right[adj] == UNKNOWN and adj < len(seq_right)-1:
            adj += 1
        seq_left = seq_left + seq_right[:adj]
        seq_right = seq_right[adj:]
        all_tvrtel_lens.append(len(seq_right))
    #
    return all_tvrtel_lens
