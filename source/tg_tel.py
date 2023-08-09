import bisect
import copy

import numpy as np

from source.tg_kmer import get_telomere_base_count, get_telomere_kmer_density, get_telomere_regions
from source.tg_plot import plot_all_read_data
from source.tg_util import makedir, posmax, RC, repeated_matches_trimming

MIN_TEL_SCORE   = 100
BTWN_TEL_SCORE  = 0.8
MAX_EDGE_NONTEL = 1000

#
#
#
def parallel_filtering_job(read_subset, my_job_i, params, params_filt, out_dict):
    [CANONICAL_STRINGS, CANONICAL_STRINGS_REV, READ_TYPE, MATCH_TRIM_STRATEGY, INPUT_TYPE, PRINT_DEBUG] = params
    [MINIMUM_READ_LEN, MINIMUM_TEL_BASES] = params_filt
    filtered_reads_out = []
    my_nontel_spans    = {}
    filt_counts_out    = {'trim_filter':0,
                          'min_readlen':0,
                          'unmapped':0,
                          'unknown_ref':0,
                          'no_chr_aln':0,
                          'min_telbases':0}
    for readname in read_subset.keys():
        abns_k = repeated_matches_trimming(sorted(read_subset[readname]), strategy=MATCH_TRIM_STRATEGY, print_debug=PRINT_DEBUG)
        # did we lose all of our alignments during trimming?
        if len(abns_k) == 0:
            filt_counts_out['trim_filter'] += 1
            continue
        # make sure string used for kmer matching is same orientation as the alignments
        # assuming softclipping was used. i.e. all alignments should have same sequence... (don't pick tel though)
        which_i = 0
        for i in range(len(abns_k)):
            if abns_k[i][2][:3] != 'tel':
                which_i = i
                break
        if abns_k[which_i][5] == 'FWD':
            rdat = abns_k[which_i][7]
        elif abns_k[which_i][5] == 'REV':
            rdat = RC(abns_k[which_i][7])
        # read len filter
        if len(rdat) < MINIMUM_READ_LEN:
            filt_counts_out['min_readlen'] += 1
            continue
        # check if we're unmapped
        refs_we_aln_to = [aln[2] for aln in abns_k]
        refs_we_aln_to = sorted(list(set(refs_we_aln_to)))
        if refs_we_aln_to == ['*']:
            filt_counts_out['unmapped'] += 1
            continue
        # check for alignments to unexpected reference contigs
        any_chr = any([n[:3] == 'chr' for n in refs_we_aln_to])
        any_tel = any([n[:3] == 'tel' for n in refs_we_aln_to])
        if any_chr is False and any_tel is False:
            filt_counts_out['unknown_ref'] += 1
            continue
        # we need at least 1 chr alignment to be anchorable anywhere
        if any_chr is False:
            filt_counts_out['no_chr_aln'] += 1
            continue
        # minimum tel content
        tel_bc = get_telomere_base_count(rdat, CANONICAL_STRINGS + CANONICAL_STRINGS_REV, mode=READ_TYPE)
        if tel_bc < MINIMUM_TEL_BASES:
            filt_counts_out['min_telbases'] += 1
            if INPUT_TYPE != 'pickle':  # use non-tel reads for downstream filtering of double-anchored tels
                for aln in abns_k:
                    if aln[2][:3] != 'tel':
                        if aln[2] not in my_nontel_spans:
                            my_nontel_spans[aln[2]] = []
                        my_nontel_spans[aln[2]].append(tuple(sorted(aln[3:5])))
            continue
        #
        filtered_reads_out.append([readname, rdat, copy.deepcopy(abns_k)])
    #
    # output
    #
    out_dict[(my_job_i,0)] = filtered_reads_out
    out_dict[(my_job_i,1)] = copy.deepcopy(filt_counts_out)
    out_dict[(my_job_i,2)] = copy.deepcopy(my_nontel_spans)

#
#
#
def get_anchored_tel(p_vs_q_power, tel_regions, abns_k, rdat, TEL_WINDOW_SIZE, FILT_PARAMS, ANCHORING_STRATEGY):
    #
    [MAXIMUM_TEL_FRAC, MAXIMUM_MINOR_PQ, MAXIMUM_UNEXPLAINED_FRAC, MAX_NONTEL_MEDIAN_KMER_DENSITY, READ_TYPE] = FILT_PARAMS
    #
    score_scalars = np.ones(len(tel_regions))
    for i in range(1, len(tel_regions)-1):
        if tel_regions[i][2] is None and tel_regions[i-1][2] == tel_regions[i+1][2]:
            score_scalars[i] = BTWN_TEL_SCORE
    #
    my_tel_len_p  = 0
    my_tel_len_q  = 0
    my_plot_title = ''
    my_tel_aln    = []
    my_tel_len    = None
    my_tel_conc   = [{'p':0, 'q':0, None:0}, {'p':0, 'q':0, None:0}]
    #
    #   get tel lengths, from the left (p-arm)
    #
    if tel_regions[0][2] is None and tel_regions[0][1] - tel_regions[0][0] > MAX_EDGE_NONTEL:
        pass
    else:
        my_score = []
        for i in range(len(tel_regions)):
            my_s = (tel_regions[i][1] - tel_regions[i][0]) * score_scalars[i]
            if tel_regions[i][2] == 'p':
                my_score.append(my_s)
            elif tel_regions[i][2] == 'q':  # not too sure how to handle this
                my_score.append(my_s)       # should this be -my_s ??
            elif tel_regions[i][2] is None:
                my_score.append(-my_s)
        cum_score = np.cumsum(my_score)
        max_i     = posmax(cum_score)
        if cum_score[max_i] >= MIN_TEL_SCORE:
            my_tel_len_p   = int(tel_regions[max_i][1] + TEL_WINDOW_SIZE/2)
            #print('P-ARM TEL:', my_tel_len_p, [int(n) for n in cum_score.tolist()], max_i, '\n')
            for i in range(0, max_i+1):
                #print('tel_regions['+str(i)+'] =', tel_regions[i])
                my_tel_conc[0][tel_regions[i][2]] += abs(tel_regions[i][1] - tel_regions[i][0])
    #
    #   get tel lengths, from the right (q-arm)
    #
    if tel_regions[-1][2] is None and tel_regions[-1][1] - tel_regions[-1][0] > MAX_EDGE_NONTEL:
        pass
    else:
        my_score = []
        for i in range(len(tel_regions)):
            my_s = (tel_regions[i][1] - tel_regions[i][0]) * score_scalars[i]
            if tel_regions[i][2] == 'q':
                my_score.append(my_s)
            elif tel_regions[i][2] == 'p':  # not too sure how to handle this
                my_score.append(my_s)       # should this be -my_s ??
            elif tel_regions[i][2] is None:
                my_score.append(-my_s)
        cum_score = np.cumsum(my_score[::-1])[::-1]
        max_i     = posmax(cum_score)
        if cum_score[max_i] >= MIN_TEL_SCORE:
            my_tel_len_q   = int(tel_regions[-1][1] - tel_regions[max_i][0] + TEL_WINDOW_SIZE/2)
            #print('Q-ARM TEL:', my_tel_len_q, [int(n) for n in cum_score.tolist()], max_i, '\n')
            for i in range(max_i, len(tel_regions)):
                #print('tel_regions['+str(i)+'] =', tel_regions[i])
                my_tel_conc[1][tel_regions[i][2]] += abs(tel_regions[i][1] - tel_regions[i][0])
    #
    # sort out which tel to pick
    #
    which_to_choose = None
    if my_tel_len_p > 0 and my_tel_len_q > 0:
        tel_frac_p = my_tel_len_p / float(len(rdat))
        tel_frac_q = my_tel_len_q / float(len(rdat))
        # try to recover anchors in cases where we cruised through short subtel sequence
        if tel_frac_p > MAXIMUM_TEL_FRAC and tel_frac_q <= MAXIMUM_TEL_FRAC:
            which_to_choose = 'q'
        elif tel_frac_q > MAXIMUM_TEL_FRAC and tel_frac_p <= MAXIMUM_TEL_FRAC:
            which_to_choose = 'p'
        # otherwise pick whichever one is longest
        elif my_tel_len_p >= my_tel_len_q:
            which_to_choose = 'p'
        else:
            which_to_choose = 'q'
    else:
        if my_tel_len_p >= my_tel_len_q:
            which_to_choose = 'p'
        else:
            which_to_choose = 'q'
    #
    if which_to_choose == 'p':
        my_tel_len_q  = 0
        my_tel_len    = my_tel_len_p
        #my_plot_title = ' : p-arm tel = ' + str(my_tel_len_p)
        my_tel_aln    = [[0, my_tel_len_p, 'tel-p', None, None, 'FWD']]
        my_tel_conc   = my_tel_conc[0]
    elif which_to_choose == 'q':
        my_tel_len_p  = 0
        my_tel_len    = my_tel_len_q
        #my_plot_title = ' : q-arm tel = ' + str(my_tel_len_q)
        my_tel_aln    = [[len(rdat)-my_tel_len_q, len(rdat), 'tel-q', None, None, 'REV']]
        my_tel_conc   = my_tel_conc[1]
    tlen_vals = [my_tel_len_p, my_tel_len_q]
    #
    # grab alignment data for anchoring
    #
    alns_with_tel = []
    tel_ref_span  = None    # estimated ref coords spanned by tel
    if my_tel_len_p > 0:
        my_tel_pos = my_tel_len_p
        for i in range(len(abns_k)):
            if abns_k[i][2][:3] == 'tel' or abns_k[i][1] <= my_tel_pos: # tel supersedes
                continue
            if abns_k[i][0] <= my_tel_pos and abns_k[i][1] > my_tel_pos:
                alns_with_tel.append(copy.deepcopy(abns_k[i]))
                adj = my_tel_pos - abns_k[i][0]
                alns_with_tel[-1][0] += adj
                if alns_with_tel[-1][5] == 'FWD':
                    alns_with_tel[-1][3] += adj
                    tel_ref_span = (alns_with_tel[-1][3] - my_tel_len_p, alns_with_tel[-1][3])
                elif alns_with_tel[-1][5] == 'REV':
                    alns_with_tel[-1][3] -= adj
                    tel_ref_span = (alns_with_tel[-1][3], alns_with_tel[-1][3] + my_tel_len_p)
            else:
                alns_with_tel.append(copy.deepcopy(abns_k[i]))
    #
    elif my_tel_len_q > 0:
        my_tel_pos = len(rdat)-my_tel_len_q
        for i in range(len(abns_k)):
            if abns_k[i][2][:3] == 'tel' or abns_k[i][0] >= my_tel_pos: # tel supersedes
                continue
            if abns_k[i][0] <= my_tel_pos and abns_k[i][1] > my_tel_pos:
                alns_with_tel.append(copy.deepcopy(abns_k[i]))
                adj = abns_k[i][1] - my_tel_pos
                alns_with_tel[-1][1] -= adj
                if alns_with_tel[-1][5] == 'FWD':
                    alns_with_tel[-1][4] -= adj
                    tel_ref_span = (alns_with_tel[-1][4], alns_with_tel[-1][4] + my_tel_len_q)
                elif alns_with_tel[-1][5] == 'REV':
                    alns_with_tel[-1][4] += adj
                    tel_ref_span = (alns_with_tel[-1][4] - my_tel_len_q, alns_with_tel[-1][4])
            else:
                alns_with_tel.append(copy.deepcopy(abns_k[i]))
    #
    # FILTERING
    #
    my_filt_string = ''
    #
    # too much tel
    #
    my_tel_frac = (my_tel_len_p + my_tel_len_q) / float(len(rdat))
    if my_tel_frac > MAXIMUM_TEL_FRAC:
        my_filt_string = 'tel_frac'
    #
    # no tel at all (read did not terminate in telomere)
    #
    if len(my_filt_string) == 0:
        if my_tel_len_p == 0 and my_tel_len_q == 0:
            my_filt_string = 'no_tel_end'
    #
    # too much mixture of p and q (skip for ont)
    #
    if len(my_filt_string) == 0 and READ_TYPE != 'ont':
        tel_types     = sorted([my_tel_conc['p'], my_tel_conc['q']])
        my_minor_frac = float(tel_types[0])/sum(tel_types)
        if my_minor_frac > MAXIMUM_MINOR_PQ:
            my_filt_string = 'pq_mix'
    #
    # remove reads with too much unexplained seq (unaligned regions)
    #
    if len(my_filt_string) == 0:
        aln_cov = np.zeros(len(rdat))
        for aln in abns_k:
            if aln[2][:3] != 'tel':
                aln_cov[aln[0]:aln[1]] = 1
        (p1, p2) = (my_tel_len_p, len(rdat) - my_tel_len_q)
        my_unexplained_frac = 1.0 - (np.sum(aln_cov[p1:p2]) / float(p2 - p1))
        if my_unexplained_frac > MAXIMUM_UNEXPLAINED_FRAC:
            my_filt_string = 'unexplained_seq'
    #
    # remove reads where the non-telomere regions look too much like telomeres
    #
    if len(my_filt_string) == 0:
        my_nontel_kmer_dens = abs(np.median(p_vs_q_power[my_tel_len_p:len(p_vs_q_power)-my_tel_len_q]))
        if my_nontel_kmer_dens > MAX_NONTEL_MEDIAN_KMER_DENSITY:
            my_filt_string = 'nontel_kmer_dens'
    #
    # remove reads where we have no viable subtel alignment to anchor tel to
    #
    if len(my_filt_string) == 0:
        if len(alns_with_tel) == 0:
            my_filt_string = 'no_viable_subtel'
    #
    # ANCHOR TELOMERE TO SUBTELOMERE ALIGNMENT
    #
    dist_to_nearest_aln = None
    adjacent_chr        = None
    adjacent_pos        = None
    adjacent_span       = None
    my_tel_type         = None
    #
    if len(my_filt_string) == 0:
        if my_tel_len_p > 0:
            alns_with_tel = my_tel_aln + alns_with_tel
            if ANCHORING_STRATEGY == 'closest':
                ati = 1
            elif ANCHORING_STRATEGY == 'largest':
                ati = sorted([(alns_with_tel[n][1] - alns_with_tel[n][0], n) for n in range(len(alns_with_tel)) if alns_with_tel[n][2][:3] != 'tel'])[-1][1]
            #
            dist_to_nearest_aln = alns_with_tel[ati][0] - alns_with_tel[0][1]
            adjacent_chr        = alns_with_tel[ati][2]
            adjacent_pos        = alns_with_tel[ati][3]
            adjacent_span       = sorted([alns_with_tel[ati][4], alns_with_tel[ati][3]])
            my_tel_type         = 'p'
            if tel_ref_span is None:
                if alns_with_tel[ati][5] == 'FWD':
                    tel_ref_span = (alns_with_tel[ati][3] - my_tel_len_p - dist_to_nearest_aln, alns_with_tel[ati][3] - dist_to_nearest_aln)
                elif alns_with_tel[ati][5] == 'REV':
                    tel_ref_span = (alns_with_tel[ati][3] + dist_to_nearest_aln, alns_with_tel[ati][3] + my_tel_len_p + dist_to_nearest_aln)
        #
        elif my_tel_len_q > 0:
            alns_with_tel = alns_with_tel + my_tel_aln
            if ANCHORING_STRATEGY == 'closest':
                ati = -2
            elif ANCHORING_STRATEGY == 'largest':
                ati = sorted([(alns_with_tel[n][1] - alns_with_tel[n][0], n) for n in range(len(alns_with_tel)) if alns_with_tel[n][2][:3] != 'tel'])[-1][1]
            #
            dist_to_nearest_aln = alns_with_tel[-1][0] - alns_with_tel[ati][1]
            adjacent_chr        = alns_with_tel[ati][2]
            adjacent_pos        = alns_with_tel[ati][4]
            adjacent_span       = sorted([alns_with_tel[ati][4], alns_with_tel[ati][3]])
            my_tel_type         = 'q'
            if tel_ref_span is None:
                if alns_with_tel[ati][5] == 'FWD':
                    tel_ref_span = (alns_with_tel[ati][4] + dist_to_nearest_aln, alns_with_tel[ati][4] + my_tel_len_q + dist_to_nearest_aln)
                elif alns_with_tel[ati][5] == 'REV':
                    tel_ref_span = (alns_with_tel[ati][4] - my_tel_len_q - dist_to_nearest_aln, alns_with_tel[ati][4] - dist_to_nearest_aln)
    #
    return (alns_with_tel,
            tel_ref_span,
            dist_to_nearest_aln,
            adjacent_chr,
            adjacent_pos,
            adjacent_span,
            my_tel_len,
            my_tel_type,
            tlen_vals,
            my_filt_string)


#
# parallel jobs:
# --- input:  FILTERED_READS + list of indices to process
# --- output: anchored_tel_dat
#             filt_string
#             tuples to add to NONTEL_REF_SPANS
#
def parallel_anchored_tel_job(read_subset, my_indices, params, params_filt, out_dict):
    [KMER_LIST, KMER_LIST_REV, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH, ANCHORING_STRATEGY, PLOT_READS, INPUT_TYPE, OUT_PLOT_DIR, PRINT_DEBUG, PLOT_FILT_READS] = params
    [MAXIMUM_TEL_FRAC, MAXIMUM_MINOR_PQ, MAXIMUM_UNEXPLAINED_FRAC, MAX_NONTEL_MEDIAN_KMER_DENSITY, READ_TYPE]                                              = params_filt    # just so I know the expected order
    for ri in range(len(read_subset)):
        [readname, rdat, abns_k] = read_subset[ri]
        my_index = my_indices[ri]
        #
        (td_p_e0, td_p_e1) = get_telomere_kmer_density(rdat, KMER_LIST,     TEL_WINDOW_SIZE)
        (td_q_e0, td_q_e1) = get_telomere_kmer_density(rdat, KMER_LIST_REV, TEL_WINDOW_SIZE)
        (p_vs_q_power, tel_regions) = get_telomere_regions(td_p_e0, td_p_e1, td_q_e0, td_q_e1, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH)
        #
        anchored_teldat = get_anchored_tel(p_vs_q_power, tel_regions, abns_k, rdat, TEL_WINDOW_SIZE, params_filt, ANCHORING_STRATEGY)
        #
        alns_with_tel       = anchored_teldat[0]
        tel_ref_span        = anchored_teldat[1]
        dist_to_nearest_aln = anchored_teldat[2]
        adjacent_chr        = anchored_teldat[3]
        adjacent_pos        = anchored_teldat[4]
        adjacent_span       = anchored_teldat[5]
        my_tel_len          = anchored_teldat[6]
        my_tel_type         = anchored_teldat[7]
        tlen_vals           = anchored_teldat[8]
        my_filt_string      = anchored_teldat[9]
        #
        # plot read data if desired
        #
        if PLOT_READS:
            if len(my_filt_string) == 0:
                my_plot_dir = OUT_PLOT_DIR + adjacent_chr + '/'
                makedir(my_plot_dir)
                plot_title = readname + ' ' + my_tel_type + '-arm tel= ' + str(my_tel_len)
                fig_name   = my_plot_dir + 'read_' + str(my_index) + '.png'
                dens_data  = [td_p_e0, td_p_e1, td_q_e0, td_q_e1, p_vs_q_power]
                plot_all_read_data(dens_data, tlen_vals, abns_k, TEL_WINDOW_SIZE, plot_title, fig_name)
            elif PLOT_FILT_READS:
                my_plot_dir = OUT_PLOT_DIR + 'filt_' + my_filt_string + '/'
                makedir(my_plot_dir)
                plot_title = readname
                fig_name   = my_plot_dir + 'read_' + str(my_index) + '.png'
                dens_data  = [td_p_e0, td_p_e1, td_q_e0, td_q_e1, p_vs_q_power]
                plot_all_read_data(dens_data, [0,0], abns_k, TEL_WINDOW_SIZE, plot_title, fig_name)

        if PRINT_DEBUG:
            print(my_indices[ri], readname)
            print(adjacent_chr, adjacent_pos, my_tel_len, tel_ref_span, dist_to_nearest_aln)
        #
        my_nontel_spans = {}
        if INPUT_TYPE != 'pickle':
            if len(my_filt_string) == 0:
                for aln in alns_with_tel:
                    if aln[2][:3] != 'tel':
                        if aln[2] not in my_nontel_spans:
                            my_nontel_spans[aln[2]] = []
                        my_nontel_spans[aln[2]].append(tuple(sorted(aln[3:5])))
                    if PRINT_DEBUG:
                        print(aln[:7])
            elif my_filt_string == 'no_anchored_tel':
                for aln in abns_k:
                    if aln[2][:3] != 'tel':
                        if aln[2] not in my_nontel_spans:
                            my_nontel_spans[aln[2]] = []
                        my_nontel_spans[aln[2]].append(tuple(sorted(aln[3:5])))
        #
        # output
        #
        out_dict[(my_index,0)] = [adjacent_chr,
                                  readname,
                                  adjacent_pos,
                                  adjacent_span,
                                  my_tel_len,
                                  my_tel_type,
                                  tel_ref_span,
                                  copy.deepcopy(rdat),
                                  [n[:7] for n in alns_with_tel]]
        out_dict[(my_index,1)] = my_filt_string
        out_dict[(my_index,2)] = copy.deepcopy(my_nontel_spans)

#
#
#
def get_double_anchored_tels(ANCHORED_TEL_BY_CHR, NONTEL_REFSPANS_BY_CHR, gdat_params):
    [MIN_DOUBLE_ANCHOR_LEN, MIN_DOUBLE_ANCHOR_READS, PRINT_DEBUG] = gdat_params
    del_keys = []
    for k in sorted(ANCHORED_TEL_BY_CHR.keys()):
        del_list = []
        if k in NONTEL_REFSPANS_BY_CHR:
            for i in range(len(ANCHORED_TEL_BY_CHR[k])):
                tel_ref_span   = ANCHORED_TEL_BY_CHR[k][i][5]
                relevant_spans = []
                if PRINT_DEBUG:
                    print('double-anchor filt:', ANCHORED_TEL_BY_CHR[k][i][:6], tel_ref_span)
                bi  = bisect.bisect(NONTEL_REFSPANS_BY_CHR[k], tel_ref_span)
                bi2 = min([bi, len(NONTEL_REFSPANS_BY_CHR[k]) - 1])
                while True:
                    if bi2 < 0 or tel_ref_span[0] - NONTEL_REFSPANS_BY_CHR[k][bi2][1] > 10*MIN_DOUBLE_ANCHOR_LEN:
                        break
                    relevant_spans.append(NONTEL_REFSPANS_BY_CHR[k][bi2])
                    bi2 -= 1
                bi2 = min([bi, len(NONTEL_REFSPANS_BY_CHR[k]) - 1])
                while True:
                    if bi2 >= len(NONTEL_REFSPANS_BY_CHR[k]) or NONTEL_REFSPANS_BY_CHR[k][bi2][0] - tel_ref_span[1] > 10*MIN_DOUBLE_ANCHOR_LEN:
                        break
                    relevant_spans.append(NONTEL_REFSPANS_BY_CHR[k][bi2])
                    bi2 += 1
                relevant_spans = sorted(relevant_spans)
                n_spanning_reads = 0
                for span in relevant_spans:
                    if tel_ref_span[0] - span[0] >= MIN_DOUBLE_ANCHOR_LEN and span[1] - tel_ref_span[1] >= MIN_DOUBLE_ANCHOR_LEN:
                        n_spanning_reads += 1
                if n_spanning_reads >= MIN_DOUBLE_ANCHOR_READS:
                    del_list.append(i)
        for di in sorted(del_list, reverse=True):
            del_keys.append((k,di))
    return del_keys

#
#
#
def get_tels_below_canonical_thresh(ANCHORED_TEL_BY_CHR, gtbct_params):
    [MINIMUM_TEL_BASES, MIN_CANONICAL_FRAC, KMER_LIST, KMER_LIST_REV, READ_TYPE] = gtbct_params
    del_keys = []
    for k in sorted(ANCHORED_TEL_BY_CHR.keys()):
        del_list = []
        for i in range(len(ANCHORED_TEL_BY_CHR[k])):
            my_tlen = ANCHORED_TEL_BY_CHR[k][i][3]
            my_type = ANCHORED_TEL_BY_CHR[k][i][4]
            my_rdat = ANCHORED_TEL_BY_CHR[k][i][6]
            if k[-1] == 'p':
                kmers_to_use = KMER_LIST
                if my_type == 'p':
                    my_telseq = my_rdat[:my_tlen]
                elif my_type == 'q':
                    my_telseq = RC(my_rdat[-my_tlen:])
            elif k[-1] == 'q':
                kmers_to_use = KMER_LIST_REV
                if my_type == 'p':
                    my_telseq = RC(my_rdat[:my_tlen])
                elif my_type == 'q':
                    my_telseq = my_rdat[-my_tlen:]
            my_canonical_bases = get_telomere_base_count(my_telseq, kmers_to_use, mode=READ_TYPE)
            my_canonical_frac  = my_canonical_bases / float(len(my_telseq))
            if my_canonical_bases < MINIMUM_TEL_BASES or my_canonical_frac < MIN_CANONICAL_FRAC:
                del_list.append(i)
        for di in sorted(del_list, reverse=True):
            del_keys.append((k,di))
    return del_keys

#
#
#
def get_allele_tsv_dat(kmer_hit_dat, read_clust_dat, my_chr, my_pos, my_rlens, gatd_params):
    [ALLELE_TL_METHOD, MIN_READS_PER_PHASE] = gatd_params
    out_dat = []
    for allele_i in range(len(read_clust_dat[0])):
        allele_tvrlen   = read_clust_dat[7][allele_i]
        allele_cons_out = ''
        if allele_tvrlen > 0:
            if my_chr[-1] == 'p':   # p will be reversed so its in subtel --> tvr --> tel orientation
                allele_cons_out = read_clust_dat[4][allele_i][-allele_tvrlen:][::-1]
            elif my_chr[-1] == 'q':
                allele_cons_out = read_clust_dat[4][allele_i][:allele_tvrlen]
        #
        # kmer_hit_dat[n][1]   = tlen + all the extra subtel buffers
        # kmer_hit_dat[n][4]   = readname (needed for downstream analyses)
        # read_clust_dat[3][n] = the length of the subtel region present before tvr/tel region
        # read_clust_dat[6][n] = the length of filtered regions past canonical telomere region
        #
        # the difference of these two will be the actual size of the (tvr + tel) region in the read
        # -- also subtract length of filtered error regions
        #
        allele_readcount = len(read_clust_dat[0][allele_i])
        allele_tlen_mapq = sorted([(kmer_hit_dat[n][1] - read_clust_dat[3][n] - read_clust_dat[6][n], my_rlens[n], kmer_hit_dat[n][4], kmer_hit_dat[n][5]) for n in read_clust_dat[0][allele_i]])
        # subtracting tvr so that "actual" TL is output. values can be negative
        allele_tlens     = [n[0]-len(allele_cons_out) for n in allele_tlen_mapq]
        allele_tlen_str  = ','.join([str(n) for n in allele_tlens])
        rlen_str         = ','.join([str(n[1]) for n in allele_tlen_mapq])
        rname_str        = ','.join([str(n[2]) for n in allele_tlen_mapq])
        mapq_str         = ','.join([str(n[3]) for n in allele_tlen_mapq])
        #
        if my_chr in ['chrUq', 'chrBq']:
            rname_str = ','.join([' '.join(n[2].split(' ')[:-1]) for n in allele_tlen_mapq])
        #
        consensus_tl_allele = choose_tl_from_observations(allele_tlens, ALLELE_TL_METHOD)
        #
        if allele_readcount >= MIN_READS_PER_PHASE:
            out_dat.append([my_chr,
                            str(my_pos),
                            str(allele_i),
                            '0',                            # allele id (will be filled out later)
                            str(int(consensus_tl_allele)),
                            allele_tlen_str,
                            rlen_str,
                            mapq_str,
                            str(len(allele_cons_out)),
                            allele_cons_out,
                            rname_str,
                            ''])                            # original chr mapping of blank tvrs (will be filled out later)
    return out_dat

#
#
#
def choose_tl_from_observations(allele_tlens, ALLELE_TL_METHOD):
    consensus_tl_allele = None
    if ALLELE_TL_METHOD == 'mean':
        consensus_tl_allele = np.mean(allele_tlens)
    elif ALLELE_TL_METHOD == 'median':
        consensus_tl_allele = np.median(allele_tlens)
    elif ALLELE_TL_METHOD == 'max':
        consensus_tl_allele = np.max(allele_tlens)
    elif ALLELE_TL_METHOD[0] == 'p':
        my_percentile = int(ALLELE_TL_METHOD[1:])
        consensus_tl_allele = np.percentile(allele_tlens, my_percentile)
    elif ALLELE_TL_METHOD[-7:] == 'fromtop':
        my_fromtop = min(len(allele_tlens), int(ALLELE_TL_METHOD[:-7]))
        consensus_tl_allele = allele_tlens[-my_fromtop]
    return consensus_tl_allele
