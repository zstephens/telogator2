import numpy as np

from source.tg_kmer import get_telomere_kmer_density, get_telomere_regions
from source.tg_plot import plot_tel_signal
from source.tg_util import posmax

MIN_TEL_SCORE = 100


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
    my_chr, my_pos, my_ref, my_id = dat1[:4]
    my_tvr_len, my_tvr = dat1[8], dat1[9]
    which_ind = 0
    if len(teldat2[0]) > len(teldat1[0]):
        my_chr, my_pos, my_ref, my_id = dat2[:4]
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


def parse_tsv(fn, min_reads=3, min_tvr=100, min_atl=-2000, max_atl=20000, min_maxatl=100, print_warnings=False):
    out_dat = []
    fail_dict = {'interstitial':0,
                 'min_tvr':0,
                 'min_maxatl':0,
                 'min_reads':0}
    with open(fn, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            splt = line.strip().split('\t')
            my_chr = splt[0].split(',')[0]
            #my_pos = int(splt[1].split(',')[0])
            #my_t2t = splt[2].split(',')[0]
            my_aid = splt[3]
            if my_aid[-1] == 'i':
                fail_dict['interstitial'] += 1
                continue
            else:
                my_aid = int(my_aid)
            tvr_len = int(splt[8])
            if tvr_len < min_tvr:
                fail_dict['min_tvr'] += 1
                continue
            my_tvr = splt[9]
            allele_tls = [int(n) for n in splt[5].split(',')]
            num_alleles_prefilt = len(allele_tls)
            allele_tls = [n for n in allele_tls if n < max_atl and n > min_atl]
            if max(allele_tls) < min_maxatl:
                fail_dict['min_maxatl'] += 1
                continue
            if print_warnings and len(allele_tls) < num_alleles_prefilt:
                print(f'warning: skipping {num_alleles_prefilt-len(allele_tls)} reads with ATL >{max_atl} or <{min_atl} in allele {my_aid}')
            if len(allele_tls) < min_reads:
                fail_dict['min_reads'] += 1
                continue
            consensus_tl = int(splt[4])
            out_dat.append((my_chr, my_aid, allele_tls, consensus_tl, my_tvr))
    return out_dat, fail_dict


def get_terminating_tl(rdat, pq, gtt_params, telplot_dat=None):
    [SIGNAL_KMERS, SIGNAL_KMERS_REV, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH] = gtt_params
    #
    (td_p_e0, td_p_e1) = get_telomere_kmer_density(rdat, SIGNAL_KMERS,     TEL_WINDOW_SIZE)
    (td_q_e0, td_q_e1) = get_telomere_kmer_density(rdat, SIGNAL_KMERS_REV, TEL_WINDOW_SIZE)
    (p_vs_q_power, tel_regions) = get_telomere_regions(td_p_e0, td_p_e1, td_q_e0, td_q_e1, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH)
    #
    my_score = []
    my_tel_len = None
    edge_nontel = 0
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
        max_i     = posmax(cum_score)
        if cum_score[max_i] >= MIN_TEL_SCORE:
            my_tel_len = int(tel_regions[max_i][1] + TEL_WINDOW_SIZE/2)
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
        max_i     = posmax(cum_score)
        if cum_score[max_i] >= MIN_TEL_SCORE:
            my_tel_len = int(tel_regions[-1][1] - tel_regions[max_i][0] + TEL_WINDOW_SIZE/2)
    #
    if telplot_dat is not None:
        (telplot_title, telplot_fn) = telplot_dat
        density_data = [td_p_e0, td_p_e1, td_q_e0, td_q_e1, p_vs_q_power, len(rdat), TEL_WINDOW_SIZE]
        tl_vals_to_plot = None
        if my_tel_len is not None and my_tel_len > 0:
            tl_vals_to_plot = [0, my_tel_len]
        plot_tel_signal(density_data, telplot_title, telplot_fn, tl_vals=tl_vals_to_plot)
    #
    if my_tel_len is None or my_tel_len <= 0:
        return (0, 0)
    return (my_tel_len, edge_nontel)
