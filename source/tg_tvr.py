import multiprocessing
import numpy as np
import matplotlib.pyplot as mpl

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance  import squareform

from source.tg_align  import UNKNOWN_LETTER, get_aligner_object, get_scoring_matrix, MATCH_CANON_PREFIX, tvr_distance
from source.tg_muscle import get_muscle_msa
from source.tg_plot   import DEFAULT_DPI, MAX_PLOT_SIZE, plot_some_tvrs
from source.tg_reader import TG_Reader
from source.tg_util   import exists_and_is_nonzero


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
# to prevent pesky div-by-zeros in edge cases
MIN_MSD = 3.0


def parallel_alignment_job(sequences, ij, pq, results_dict,
                           scoring_matrix=None,
                           gap_bool=(True,True),
                           adjust_lens=True,
                           min_viable=True,
                           randshuffle=3):
    aligner = get_aligner_object(scoring_matrix=scoring_matrix, gap_bool=gap_bool)
    for (i,j) in ij:
        results_dict[(i,j)] = tvr_distance(sequences[i], sequences[j], aligner, adjust_lens=adjust_lens, min_viable=min_viable, randshuffle=randshuffle, pq=pq)


def parallel_msa_job(my_inds, out_clust, buffered_tvrs, buffered_subs, muscle_params, results_dict):
    [muscle_exe, muscle_prefix, char_score_adj, noncanon_cheat] = muscle_params
    for i in my_inds:
        if len(out_clust[i]) == 1:
            results_dict[('tvr',i)] = buffered_tvrs[out_clust[i][0]]
            results_dict[('sub',i)] = buffered_subs[out_clust[i][0]]
        else:
            clust_seq = [buffered_tvrs[n] for n in out_clust[i]]
            my_prefix = muscle_prefix + '_'*(len(muscle_prefix) > 0 and muscle_prefix[-1] != '/') + 'tvr' + str(i).zfill(5)
            [msa_seq, consensus_seq] = get_muscle_msa(clust_seq, muscle_exe, tempfile_prefix=my_prefix, char_score_adj=char_score_adj, noncanon_cheat=noncanon_cheat)
            results_dict[('tvr',i)] = consensus_seq
            clust_seq = [buffered_subs[n] for n in out_clust[i]]
            if clust_seq[0] == '':
                results_dict[('sub',i)] = ''
            else:
                my_prefix = muscle_prefix + '_'*(len(muscle_prefix) > 0 and muscle_prefix[-1] != '/') + 'sub' + str(i).zfill(5)
                [msa_seq, consensus_seq] = get_muscle_msa(clust_seq, muscle_exe, tempfile_prefix=my_prefix)
                results_dict[('sub',i)] = consensus_seq


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
                 save_msa=None,
                 tvr_truncate=3000,
                 alignment_processes=4,
                 muscle_exe='muscle',
                 muscle_dir='',
                 PRINT_DEBUG=False):
    #
    #   kmer_dat[i] = [[kmer1_hits, kmer2_hits, ...], tlen, tel_anchor_dist, read_orientation, readname, anchor_mapq, fasta_dat]
    #
    #   repeats_metadata = [kmer_list, kmer_colors, kmer_letters, kmer_flags]
    #
    [kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
    #
    muscle_prefix = muscle_dir
    if save_msa is not None and '.' in save_msa:
        muscle_prefix = '.'.join(save_msa.split('.')[:-1])
    #
    n_reads = len(kmer_dat)
    pq      = my_chr[-1]
    #
    canonical_letter = None
    denoise_letters  = []
    tvr_letters      = []
    nofilt_letters   = []
    dubious_letters  = [UNKNOWN_LETTER]
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
        if kmer_letters[i] == UNKNOWN_LETTER:
            print(f'Error: character {UNKNOWN_LETTER} is reserved for unknown sequence')
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
    # when generating consensus sequence for cluster: in ties, prioritize canonical, deprioritize unknown
    #
    char_score_adj = {canonical_letter:1, UNKNOWN_LETTER:-1}
    noncanon_cheat = ([canonical_letter]+dubious_letters, 3)
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
            (seq_left, seq_right) = find_density_boundary(all_colorvecs[i][::-1], [UNKNOWN_LETTER], UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
            # denoise tvr+tel section
            seq_right_denoise = denoise_colorvec(seq_right, replace_char=canonical_letter, chars_to_delete=denoise_letters, chars_to_merge=[canonical_letter]+tvr_letters)
            # remove ends of reads that might be sequencing artifacts, based on density of canonical characters
            (err_left, err_right) = find_density_boundary(seq_right_denoise[::-1], [canonical_letter], CANON_WIN_SIZE, CANON_END_DENS, thresh_dir='above')
            #
            cleaned_colorvecs.append(err_right + seq_left[::-1])    # the entire denoised read (for plotting)
            err_end_lens.append(len(err_left))                      # needed for adjusting offsets in plots
            colorvecs_for_msa.append(err_right)
            #
            subtel_regions.append(seq_left[::-1])                   # subtel regions (currently not used for anything)
            tvrtel_regions.append(err_right)                        # tvr sequence used for clustering
            if len(tvrtel_regions[-1]) > tvr_truncate:
                tvrtel_regions[-1] = tvrtel_regions[-1][-tvr_truncate:]
        #
        elif pq == 'q':
            # identify subtel/tvr boundary based on density of unknown characters
            (seq_left, seq_right) = find_density_boundary(all_colorvecs[i], [UNKNOWN_LETTER]+subfilt_letters, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below', debug_plot=False, readname=kmer_dat[i][4])
            # lets walk backwards with a smaller window to see if there are some tvrs we missed
            (recovered_tvr, rest_is_subtel) = find_density_boundary(seq_left[::-1], [UNKNOWN_LETTER]+subfilt_letters, UNKNOWN_WIN_SIZE_FINE, UNKNOWN_END_DENS_FINE, thresh_dir='above', debug_plot=False, readname=kmer_dat[i][4])
            seq_right = recovered_tvr[::-1] + seq_right
            seq_left = rest_is_subtel[::-1]
            # adjust subtel boundary so that tvr does not begin with an unknown character
            adj = 0
            while seq_right[adj] == UNKNOWN_LETTER and adj < len(seq_right)-1:
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
            subtel_regions.append(seq_left)                         # subtel regions (currently not used for anything)
            tvrtel_regions.append(err_right[::-1])                  # tvr sequence used for clustering
            if len(tvrtel_regions[-1]) > tvr_truncate:
                tvrtel_regions[-1] = tvrtel_regions[-1][:tvr_truncate]

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
    out_clust = [list(range(len(tvrtel_regions)))]
    if len(tvrtel_regions) > 1:
        dist_array = squareform(dist_matrix)
        Zread = linkage(dist_array, method='ward')
        #
        if fig_name is not None:
            dendrogram_width = max(8, 0.12 * n_reads)
            dendrogram_width = min(dendrogram_width, MAX_PLOT_SIZE / DEFAULT_DPI)
            mpl.rcParams.update({'font.size': 20, 'lines.linewidth':2.0})
            fig = mpl.figure(3, figsize=(dendrogram_width,6), dpi=DEFAULT_DPI)
            dendrogram(Zread, color_threshold=tree_cut)
            mpl.axhline(y=[tree_cut], linestyle='dashed', color='black', alpha=0.7)
            mpl.xlabel('read #')
            mpl.ylabel('distance')
            mpl.title(my_chr + ' : ' + str(my_pos))
            #mpl.tight_layout() # for some reason this was using 200 dpi instead of 100 and was returning ValueErrors for figures being too large
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
        manager     = multiprocessing.Manager()
        msa_results = manager.dict()
        processes   = []
        muscle_params = [muscle_exe, muscle_prefix, char_score_adj, noncanon_cheat]
        for i in range(alignment_processes):
            my_inds = list(range(i,len(out_clust),alignment_processes))
            p = multiprocessing.Process(target=parallel_msa_job,
                                        args=(my_inds, out_clust, buffered_tvrs, buffered_subs, muscle_params, msa_results))
            processes.append(p)
        for i in range(alignment_processes):
            processes[i].start()
        for i in range(alignment_processes):
            processes[i].join()
        #
        out_consensus    = ['' for n in out_clust]
        subtel_consensus = ['' for n in out_clust]
        msa_keys = sorted(msa_results.keys())
        for k in msa_keys:
            if k[0] == 'tvr':
                out_consensus[k[1]] = msa_results[k]
            elif k[0] == 'sub':
                subtel_consensus[k[1]] = msa_results[k]
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
                                              replace_char=canonical_letter,
                                              chars_to_delete=[n for n in tvr_letters if n not in nofilt_letters],
                                              min_size=TVR_CANON_FILT_PARAMS_STRICT[0],
                                              chars_to_merge=[canonical_letter]+tvr_letters)
        if pq == 'q':
            denoised_consensus = denoised_consensus[::-1]
        tel_boundary = find_cumulative_boundary(denoised_consensus,tvr_letters,
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
                                                  replace_char=canonical_letter,
                                                  chars_to_delete=[n for n in tvr_letters if n not in nofilt_letters],
                                                  min_size=TVR_CANON_FILT_PARAMS_LENIENT[0],
                                                  chars_to_merge=[canonical_letter]+tvr_letters)
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
        for mc in merge_clust:
            ci_with_longest_tvr = sorted([(len(out_consensus[ci]), ci) for ci in mc])[-1][1]
            new_out_clust.append([])
            for ci in mc:
                for ri in out_clust[ci]:
                    new_out_clust[-1].append((len(colorvecs_for_msa[ri]), ri))  # gonna sort by length
            new_out_clust[-1] = [n[1] for n in sorted(new_out_clust[-1], reverse=True)]
            new_out_consensus.append(out_consensus[ci_with_longest_tvr])
            new_out_tvr_tel_boundaries.append(out_tvr_tel_boundaries[ci_with_longest_tvr])
        out_clust              = new_out_clust
        out_consensus          = new_out_consensus
        out_tvr_tel_boundaries = new_out_tvr_tel_boundaries

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
        out_adj.append([longest_subtel_cl - n + clust_subtel_adj for n in my_subtel_lens])

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
                           rand_shuffle_count=3,
                           linkage_method='complete',
                           normalize_dist_matrix=True,
                           alignment_processes=8,
                           job=(1,1),
                           dendrogram_title=None,
                           dendrogram_height=12,
                           dendrogram_allblack=False,
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
        #
        out_clust = [list(range(n_seq))]
        if n_seq > 1:
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
            if dendrogram_allblack:
                dendrogram(Zread, orientation='left', labels=replot_labels, color_threshold=tree_cut, above_threshold_color='black', link_color_func=lambda k:'black')
            else:
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
    #
    return None


def filter_by_denoise_frac(kmer_dat, repeats_metadata, my_chr):
    #
    n_reads = len(kmer_dat)
    pq      = my_chr[-1]
    [kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
    denoise_letters = []
    for i in range(len(kmer_list)):
        if 'denoise' in kmer_flags[i]:
            denoise_letters.append(kmer_letters[i])
    denoise_letters = list(set(denoise_letters))
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
            (seq_left, seq_right) = find_density_boundary(my_cvec[::-1], [UNKNOWN_LETTER], UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
        elif pq == 'q':
            (seq_left, seq_right) = find_density_boundary(my_cvec, [UNKNOWN_LETTER], UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
        #
        my_count_unknown = seq_right.count(UNKNOWN_LETTER)
        my_count_denoise = 0
        for dc in denoise_letters:
            my_count_denoise += seq_right.count(dc)
        my_num = my_count_unknown + my_count_denoise
        my_den = float(max([1, len(seq_right)]))
        my_out.append((my_num, my_num/my_den))
    return my_out


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


def denoise_colorvec(v, replace_char=UNKNOWN_LETTER, min_size=10, max_gap_fill=50, chars_to_delete=[], chars_to_merge=[]):
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
