import numpy as np

from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

from source.tg_kmer import get_nonoverlapping_kmer_hits, get_telomere_base_count, get_telomere_kmer_density, get_telomere_regions
from source.tg_plot import plot_tel_signal
from source.tg_util import posmax, RC


def choose_tl_from_observations(allele_tlens, ALLELE_TL_METHOD, skip_negative_vals=False):
    if skip_negative_vals:
        atls = [n for n in allele_tlens if n > 0]
    else:
        atls = allele_tlens
    if len(atls) == 0:
        return None
    consensus_tl_allele = None
    if ALLELE_TL_METHOD == 'mean':
        consensus_tl_allele = np.mean(atls)
    elif ALLELE_TL_METHOD == 'median':
        consensus_tl_allele = np.median(atls)
    elif ALLELE_TL_METHOD == 'max':
        consensus_tl_allele = np.max(atls)
    elif ALLELE_TL_METHOD[0] == 'p':
        my_percentile = int(ALLELE_TL_METHOD[1:])
        consensus_tl_allele = np.percentile(atls, my_percentile)
    elif ALLELE_TL_METHOD[-7:] == 'fromtop':
        my_fromtop = min(len(atls), int(ALLELE_TL_METHOD[:-7]))
        consensus_tl_allele = atls[-my_fromtop]
    return int(consensus_tl_allele)


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
        consensus_atl = choose_tl_from_observations(allele_tlens, ALLELE_TL_METHOD)
        #
        if allele_readcount >= MIN_READS_PER_PHASE:
            out_dat.append([my_chr,
                            str(my_pos),
                            '-',                            # ref builds (will be filled out later)
                            '0',                            # allele id (will be filled out later)
                            str(int(consensus_atl)),
                            allele_tlen_str,
                            rlen_str,
                            mapq_str,
                            str(len(allele_cons_out)),
                            allele_cons_out,
                            rname_str])
    return out_dat


def merge_allele_tsv_dat(dat1, dat2, ALLELE_TL_METHOD):
    # (read_tls, read_lens, read_mapq, read_names)
    teldat1 = [[int(n) for n in dat1[5].split(',')], dat1[6].split(','), dat1[7].split(','), dat1[10].split(',')]
    teldat2 = [[int(n) for n in dat2[5].split(',')], dat2[6].split(','), dat2[7].split(','), dat2[10].split(',')]
    my_id = dat1[3] + ',' + dat2[3]
    my_chr, my_pos, my_ref = dat1[:3]
    my_tvr_len, my_tvr = dat1[8], dat1[9]
    which_ind = 0
    if len(teldat2[0]) > len(teldat1[0]):
        my_id = dat2[3] + ',' + dat1[3]
        my_chr, my_pos, my_ref = dat2[:3]
        my_tvr_len, my_tvr = dat2[8], dat2[9]
        which_ind = 1
    merged_dat = [teldat1[0] + teldat2[0], teldat1[1] + teldat2[1], teldat1[2] + teldat2[2], teldat1[3] + teldat2[3]]
    sorted_dat = sorted([col for col in zip(*merged_dat)])
    consensus_atl = choose_tl_from_observations([n[0] for n in sorted_dat], ALLELE_TL_METHOD)
    allele_tlen_str = ','.join([str(n[0]) for n in sorted_dat])
    rlen_str = ','.join([n[1] for n in sorted_dat])
    mapq_str = ','.join([n[2] for n in sorted_dat])
    rname_str = ','.join([n[3] for n in sorted_dat])
    return which_ind, [my_chr, my_pos, my_ref, my_id, str(int(consensus_atl)), allele_tlen_str, rlen_str, mapq_str, my_tvr_len, my_tvr, rname_str]


def split_allele_tsv_dat_by_readnames(dat, readname_list):
    # (read_tls, read_lens, read_mapq, read_names)
    teldat = [[int(n) for n in dat[5].split(',')], dat[6].split(','), dat[7].split(','), dat[10].split(',')]
    ind_is_in_readlist = [n in readname_list for n in teldat[3]]
    print(dat[0:3])
    print(ind_is_in_readlist.count(True), '/', len(ind_is_in_readlist))


def parse_tsv(fn, min_reads=3, min_tvr=100, min_atl=-2000, max_atl=20000, min_maxatl=100, print_warnings=False):
    out_dat = []
    fail_dict = {'interstitial':0,
                 'min_tvr':0,
                 'min_maxatl':0,
                 'no_atl_in_range':0,
                 'min_reads':0}
    with open(fn, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            splt = line.strip().split('\t')
            my_chr = splt[0].split(',')[0]
            #my_pos = int(splt[1].split(',')[0])
            my_refbuild = splt[2].split(',')[0]
            my_aid = splt[3].split(',')
            if any([n[-1] == 'i' for n in my_aid]):
                fail_dict['interstitial'] += 1
                continue
            else:
                my_aid = int(my_aid[0])
            tvr_len = int(splt[8])
            if tvr_len < min_tvr:
                fail_dict['min_tvr'] += 1
                continue
            my_tvr = splt[9]
            allele_tls = [int(n) for n in splt[5].split(',')]
            num_alleles_prefilt = len(allele_tls)
            allele_tls = [n for n in allele_tls if n < max_atl and n > min_atl]
            if not allele_tls:
                fail_dict['no_atl_in_range'] += 1
                continue
            if max(allele_tls) < min_maxatl:
                fail_dict['min_maxatl'] += 1
                continue
            if print_warnings and len(allele_tls) < num_alleles_prefilt:
                print(f'warning: skipping {num_alleles_prefilt-len(allele_tls)} reads with ATL >{max_atl} or <{min_atl} in allele {my_aid}')
            if len(allele_tls) < min_reads:
                fail_dict['min_reads'] += 1
                continue
            consensus_tl = int(splt[4])
            out_dat.append((my_chr, my_refbuild, my_aid, allele_tls, consensus_tl, my_tvr))
    return out_dat, fail_dict


def get_terminating_tl(rdat, pq, gtt_params, telplot_dat=None):
    [SIGNAL_KMERS, SIGNAL_KMERS_REV, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH, MIN_TEL_SCORE] = gtt_params
    window_adj = TEL_WINDOW_SIZE // 2
    #
    (td_p_e0, td_p_e1) = get_telomere_kmer_density(rdat, SIGNAL_KMERS,     TEL_WINDOW_SIZE)
    (td_q_e0, td_q_e1) = get_telomere_kmer_density(rdat, SIGNAL_KMERS_REV, TEL_WINDOW_SIZE)
    (p_vs_q_power, tel_regions) = get_telomere_regions(td_p_e0, td_p_e1, td_q_e0, td_q_e1, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH)
    #
    my_score = []
    my_tel_len = None
    edge_nontel = 0
    interstitial_regions = [n for n in tel_regions]
    if pq == 'p':
        if tel_regions[0][2] is None:
            edge_nontel = tel_regions[0][1] - tel_regions[0][0]
        for i in range(len(tel_regions)):
            my_s = (tel_regions[i][1] - tel_regions[i][0])
            if tel_regions[i][2] == 'p':
                my_score.append(my_s)
            elif tel_regions[i][2] == 'q':  # not too sure how to handle this
                my_score.append(my_s)       # should this be -my_s ??
            elif tel_regions[i][2] is None:
                my_score.append(-my_s)
        cum_score = np.cumsum(my_score)
        max_i = posmax(cum_score)
        if cum_score[max_i] >= MIN_TEL_SCORE:
            my_tel_len = int(tel_regions[max_i][1] + window_adj)
            interstitial_regions = tel_regions[max_i:]
    elif pq == 'q':
        if tel_regions[-1][2] is None:
            edge_nontel = tel_regions[-1][1] - tel_regions[-1][0]
        for i in range(len(tel_regions)):
            my_s = (tel_regions[i][1] - tel_regions[i][0])
            if tel_regions[i][2] == 'q':
                my_score.append(my_s)
            elif tel_regions[i][2] == 'p':  # not too sure how to handle this
                my_score.append(my_s)       # should this be -my_s ??
            elif tel_regions[i][2] is None:
                my_score.append(-my_s)
        cum_score = np.cumsum(my_score[::-1])[::-1]
        max_i = posmax(cum_score)
        if cum_score[max_i] >= MIN_TEL_SCORE:
            my_tel_len = int(tel_regions[-1][1] - tel_regions[max_i][0] + window_adj)
            interstitial_regions = tel_regions[:max_i]
    #
    # get largest interstitial telomere region
    #
    largest_interstitial_tel = None
    interstitial_tels = [n for n in interstitial_regions if n[2] is not None]
    if interstitial_tels:
        clustered_regions = []
        current_clust = [0]
        for i in range(len(interstitial_tels) - 1):
            my_size   = interstitial_tels[i][1]   - interstitial_tels[i][0]
            next_size = interstitial_tels[i+1][1] - interstitial_tels[i+1][0]
            our_dist  = interstitial_tels[i+1][0] - interstitial_tels[i][1]
            if max(my_size, next_size) > our_dist:
                current_clust.append(i+1)
            else:
                clustered_regions.append([n for n in current_clust])
                current_clust = [i+1]
        if current_clust:
            clustered_regions.append([n for n in current_clust])
        largest_interstitial_tel = [(interstitial_tels[min(n)][0] + window_adj, interstitial_tels[max(n)][1] + window_adj) for n in clustered_regions]
        largest_interstitial_tel = sorted([(n[1] - n[0], n) for n in largest_interstitial_tel], reverse=True)[0][1]
    #
    # plot telomere signals, if desired
    #
    if telplot_dat is not None:
        (telplot_title, telplot_fn) = telplot_dat
        density_data = [td_p_e0, td_p_e1, td_q_e0, td_q_e1, p_vs_q_power, len(rdat), TEL_WINDOW_SIZE]
        tl_vals_to_plot = None
        if my_tel_len is not None and my_tel_len > 0:
            tl_vals_to_plot = [0, my_tel_len]
        print(telplot_title, telplot_fn)
        plot_tel_signal(density_data, telplot_title, telplot_fn, tl_vals=tl_vals_to_plot, interstitial_tels=[largest_interstitial_tel])
    #
    if my_tel_len is None or my_tel_len <= 0:
        return (0, 0, largest_interstitial_tel)
    return (my_tel_len, edge_nontel, largest_interstitial_tel)


def gtrc_parallel_job(my_read_dat, tel_signal_plot_num, gtrc_params):
    (my_rnm, my_rdat, my_qdat) = my_read_dat
    [READ_TYPE, DUMMY_TEL_MAPQ,
     KMER_LIST, KMER_LIST_REV, KMER_ISSUBSTRING, CANONICAL_STRINGS, CANONICAL_STRINGS_REV,
     TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH,
     MIN_TEL_SCORE, FILT_TERM_TEL, FILT_TERM_NONTEL, FILT_TERM_SUBTEL,
     PLOT_TEL_SIGNALS, TEL_SIGNAL_DIR,
     MIN_INTERSTITIAL_TL, MIN_FUSION_ANCHOR] = gtrc_params
    #
    gtt_min_tel_score = MIN_TEL_SCORE
    if FILT_TERM_TEL > 0:
        gtt_min_tel_score = min(FILT_TERM_TEL, MIN_TEL_SCORE)
    gtt_params = [KMER_LIST, KMER_LIST_REV, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH, gtt_min_tel_score]
    #
    result_khd = None
    result_termtel = None
    result_nontelend = None
    result_isubtels = None
    result_ikhd = None
    result_ikhdrev = None
    skipped_termtel = False
    skipped_termunk = False
    skipped_termsub = False
    #
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
    (my_terminating_tel, my_nontel_end, interstitial_tel) = get_terminating_tl(my_rdat, 'q', gtt_params, telplot_dat=tel_signal_plot_dat)
    # some paranoid bounds checking
    my_terminating_tel = min(my_terminating_tel, len(my_rdat))
    my_nontel_end = min(my_nontel_end, len(my_rdat))
    #
    if my_terminating_tel == 0 and interstitial_tel is not None and interstitial_tel[1] - interstitial_tel[0] >= MIN_INTERSTITIAL_TL:
        isubtel_left  = my_rdat[:interstitial_tel[0]]
        isubtel_right = my_rdat[interstitial_tel[1]:]
        if len(isubtel_left) >= MIN_FUSION_ANCHOR and len(isubtel_right) >= MIN_FUSION_ANCHOR:
            result_isubtels = [(f'interstitial-left_{len(my_rdat)}_{my_rnm}', isubtel_left),
                               (f'interstitial-right_{len(my_rdat)}_{my_rnm}', isubtel_right)]
            result_ikhd = [get_nonoverlapping_kmer_hits(my_rdat, KMER_LIST, KMER_ISSUBSTRING),
                           len(my_rdat), 0, 'q', my_rnm.split(' ')[0], DUMMY_TEL_MAPQ, my_rdat]
            result_ikhdrev = [get_nonoverlapping_kmer_hits(my_rdat, KMER_LIST_REV, KMER_ISSUBSTRING),
                              len(my_rdat), 0, 'q', my_rnm.split(' ')[0], DUMMY_TEL_MAPQ, my_rdat]
    #
    if FILT_TERM_TEL > 0 and my_terminating_tel < FILT_TERM_TEL:
        skipped_termtel = True
    if FILT_TERM_NONTEL > 0 and my_nontel_end > FILT_TERM_NONTEL:
        skipped_termunk = True
    # too little subtel sequence?
    if FILT_TERM_SUBTEL > 0 and len(my_rdat) < my_terminating_tel + FILT_TERM_SUBTEL:
        skipped_termsub = True
    #
    if all([n is False for n in [skipped_termtel, skipped_termunk, skipped_termsub]]):
        my_subtel_end = max(len(my_rdat)-my_terminating_tel-FILT_TERM_SUBTEL, 0)
        my_teltvr_seq = my_rdat[my_subtel_end:]
        # if there's no terminating tel at all, then lets pretend entire read is tvr+tel
        if len(my_teltvr_seq) == 0:
            my_teltvr_seq = my_rdat
        #
        # a lot of these fields in kmer_hit_dat are holdovers from telogator1
        # and no longer used, but still need to be populated with values.
        #
        result_khd = [get_nonoverlapping_kmer_hits(my_teltvr_seq, KMER_LIST_REV, KMER_ISSUBSTRING),
                      len(my_teltvr_seq),   # atb, lets pretend entire read is tel
                      0,                    # my_dbta
                      'q',                  # my_type
                      my_rnm.split(' ')[0], # my_rnm
                      DUMMY_TEL_MAPQ,       # my_mapq
                      my_rdat]              # read sequence
        result_termtel = my_terminating_tel
        result_nontelend = my_nontel_end
    #
    return [result_khd,
            result_termtel,
            result_nontelend,
            result_isubtels,
            result_ikhd,
            result_ikhdrev,
            skipped_termtel,
            skipped_termunk,
            skipped_termsub]


def get_tel_repeat_comp_parallel(all_read_dat, gtrc_params, max_workers=4, max_pending=100):
    kmer_hit_dat = []
    all_terminating_tl = []
    all_nontel_end = []
    reads_removed_term_tel = 0
    reads_removed_term_unk = 0
    reads_removed_term_sub = 0
    interstitial_subtels = []
    interstitial_khd_by_readname = {}
    interstitial_khd_rev_by_readname = {}
    #
    tasks = iter(range(len(all_read_dat)))
    sorted_gtrc_results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        pending_futures = {}
        while True:
            while len(pending_futures) < max_pending:
                try:
                    i = next(tasks)
                    future = executor.submit(gtrc_parallel_job, all_read_dat[i], i, gtrc_params)
                    pending_futures[future] = i
                except StopIteration:
                    break
            if not pending_futures:
                sorted_gtrc_results = sorted(sorted_gtrc_results)
                for (i, [result_khd,
                         result_termtel,
                         result_nontelend,
                         result_isubtels,
                         result_ikhd,
                         result_ikhdrev,
                         skipped_termtel,
                         skipped_termunk,
                         skipped_termsub]) in sorted_gtrc_results:
                    if result_khd is not None:
                        kmer_hit_dat.append(result_khd)
                        all_terminating_tl.append(result_termtel)
                        all_nontel_end.append(result_nontelend)
                    if result_isubtels is not None:
                        my_rnm = all_read_dat[i][0]
                        interstitial_subtels.append(result_isubtels)
                        interstitial_khd_by_readname[my_rnm] = result_ikhd
                        interstitial_khd_rev_by_readname[my_rnm] = result_ikhdrev
                    reads_removed_term_tel += skipped_termtel * 1
                    reads_removed_term_unk += skipped_termunk * 1
                    reads_removed_term_sub += skipped_termsub * 1
                return [kmer_hit_dat,
                        all_terminating_tl,
                        all_nontel_end,
                        reads_removed_term_tel,
                        reads_removed_term_unk,
                        reads_removed_term_sub,
                        interstitial_subtels,
                        interstitial_khd_by_readname,
                        interstitial_khd_rev_by_readname]
            #
            done, _ = wait(pending_futures, return_when=FIRST_COMPLETED)
            for future in done:
                sorted_gtrc_results.append((pending_futures.pop(future), future.result()))
