import pywt
import re

import numpy as np

from source.tg_util import exists_and_is_nonzero

INCLUDE_SUBTEL_BUFF = 500


def read_kmer_tsv(fn, READ_TYPE):
    if exists_and_is_nonzero(fn) is False:
        print('Error: ' + fn + ' not found.')
        exit(1)
    #
    KMER_LIST   = []
    KMER_COLORS = []
    KMER_LETTER = []
    KMER_FLAGS  = []
    CANONICAL_STRINGS = []
    #
    with open(fn,'r') as f:
        for line in f:
            if line[0] != '#' and len(line.strip()):
                splt = line.strip().split('\t')
                splt2 = splt[3].split(',')
                if 'pacbio_only' in splt2 and READ_TYPE not in ['hifi']:
                    continue
                if 'nanopore_only' in splt2 and READ_TYPE not in ['ont']:
                    continue
                KMER_LIST.append(splt[0])
                KMER_COLORS.append(splt[1])
                KMER_LETTER.append(splt[2])
                KMER_FLAGS.append(tuple(splt2))
                if 'canonical' in KMER_FLAGS[-1]:
                    CANONICAL_STRINGS.append(KMER_LIST[-1])
    #
    sorted_kmer_dat  = sorted(list(set([(len(KMER_LIST[n]), KMER_LIST[n], KMER_COLORS[n], KMER_LETTER[n], KMER_FLAGS[n]) for n in range(len(KMER_LIST))])), reverse=True)   # sort by length
    KMER_LIST        = [n[1] for n in sorted_kmer_dat]
    KMER_COLORS      = [n[2] for n in sorted_kmer_dat]
    KMER_LETTER      = [n[3] for n in sorted_kmer_dat]
    KMER_FLAGS       = [n[4] for n in sorted_kmer_dat]
    KMER_METADATA    = [KMER_LIST, KMER_COLORS, KMER_LETTER, KMER_FLAGS]
    KMER_ISSUBSTRING = []
    for i in range(len(KMER_LIST)):
        KMER_ISSUBSTRING.append([j for j in range(len(KMER_LIST)) if (j != i and KMER_LIST[i] in KMER_LIST[j])])
    #
    return (KMER_METADATA, KMER_ISSUBSTRING, CANONICAL_STRINGS)


def get_canonical_letter(fn, READ_TYPE):
    with open(fn,'r') as f:
        for line in f:
            if line[0] != '#' and len(line.strip()):
                splt = line.strip().split('\t')
                splt2 = splt[3].split(',')
                if 'pacbio_only' in splt2 and READ_TYPE not in ['hifi']:
                    continue
                if 'nanopore_only' in splt2 and READ_TYPE not in ['ont']:
                    continue
                if 'canonical' in splt2:
                    return splt[2]
    return None


def get_kmer_coverage(read_dat, kmer_list):
    # positions covered by at least one regex hit of any kmer in the list
    starts = []
    ends = []
    for kmer in kmer_list:
        for n in re.finditer(kmer, read_dat):
            starts.append(n.start(0))
            ends.append(n.end(0))
    delta = np.zeros(len(read_dat) + 1, dtype='<i8')
    np.add.at(delta, np.array(starts, dtype='<i8'), 1)
    np.add.at(delta, np.array(ends, dtype='<i8'), -1)
    return np.cumsum(delta[:-1]) > 0


def get_telomere_kmer_density(read_dat, kmer_list, tel_window, mode='hifi', smoothing=False):
    n_win = len(read_dat) - tel_window
    if n_win <= 0:
        return ([], [])
    tel_hit_cum_e0 = np.cumsum(get_kmer_coverage(read_dat, kmer_list), dtype='f8')
    tel_density_e0 = (tel_hit_cum_e0[tel_window:tel_window+n_win] - tel_hit_cum_e0[:n_win]) / tel_window
    tel_density_e1 = np.zeros(n_win)
    #
    if smoothing:
        tel_density_e0 = wavelet_smooth(tel_density_e0)
        tel_density_e1 = wavelet_smooth(tel_density_e1)
    #
    return (tel_density_e0, tel_density_e1)


def get_telomere_base_count(read_dat, kmer_list, mode='hifi'):
    return int(np.count_nonzero(get_kmer_coverage(read_dat, kmer_list)))


def get_telomere_regions(td_p_e0, td_p_e1, td_q_e0, td_q_e1, tel_window, pthresh, mode='hifi', smoothing=True):
    min_win_pos  = min([len(td_p_e0), len(td_p_e1), len(td_q_e0), len(td_q_e1)])
    max_win_pos  = max([len(td_p_e0), len(td_p_e1), len(td_q_e0), len(td_q_e1)])
    p_vs_q_power = np.zeros(max_win_pos)
    tel_regions  = [[0,None,None]]
    # somehow we received a sequence that's too short to work with (what went wrong?)
    if max_win_pos < tel_window:
        return ([], [[0,0,None]])
    #
    p_vs_q_power[:min_win_pos] = np.asarray(td_p_e0[:min_win_pos], dtype='f8') - np.asarray(td_q_e0[:min_win_pos], dtype='f8')
    if smoothing:
        p_vs_q_power = wavelet_smooth(p_vs_q_power)
    p_vs_q_power = p_vs_q_power[:-tel_window]
    #
    # 1 = p, 2 = q, 0 = neither; regions are runs of the same state
    state = np.where(p_vs_q_power >= pthresh, 1, np.where(p_vs_q_power <= -pthresh, 2, 0))
    run_starts = np.concatenate([[0], np.flatnonzero(state[1:] != state[:-1]) + 1]).tolist()
    run_ends = run_starts[1:] + [len(state)]
    state_label = {0:None, 1:'p', 2:'q'}
    tel_regions = [[run_starts[i], run_ends[i], state_label[int(state[run_starts[i]])]] for i in range(len(run_starts))]
    if run_ends[-1] - run_starts[-1] == 1:  # last region only 1bp wide
        tel_regions[-1][1] = len(state) - 1
    #
    return (p_vs_q_power, tel_regions)


def mad(x, c=1.4826):
    #
    # median absolute deviation, with some normalization
    #
    return c * np.median(np.abs(x - np.median(x)))


def wavelet_smooth(x, wavelet="db4", smoothing_level=5, fail_sigma=1e-3):
    #
    # adapted from: http://connor-johnson.com/2016/01/24/using-pywavelets-to-remove-high-frequency-noise/
    #
    coeff = pywt.wavedec(x, wavelet, mode="per")
    # not enough coeffs
    if len(coeff) < smoothing_level:
        return x
    # calculate a threshold
    sigma = mad(coeff[-smoothing_level])
    # if sigma below some value (i.e. signal is mostly 0s) return unsmoothed signal
    if sigma < fail_sigma:
        return x
    # changing this threshold also changes the behavior, but I have not played with this very much
    uthresh = sigma * np.sqrt(2 * np.log(len(x)))
    coeff[1:] = (pywt.threshold(i, value=uthresh, mode="soft") for i in coeff[1:])
    # reconstruct the signal using the thresholded coefficients
    y = pywt.waverec(coeff, wavelet, mode="per")
    return y


def get_nonoverlapping_kmer_hits(my_telseq, KMER_LIST, KMER_ISSUBSTRING, mode='hifi'):
    seq_len = len(my_telseq)
    raw_starts = []
    raw_ends = []
    raw_cumcov = []
    for kmer in KMER_LIST:
        hits = [(n.start(0), n.end(0)) for n in re.finditer(kmer, my_telseq)]
        my_s = np.array([n[0] for n in hits], dtype='<i8')
        my_e = np.array([n[1] for n in hits], dtype='<i8')
        delta = np.zeros(seq_len + 1, dtype='<i8')
        np.add.at(delta, my_s, 1)
        np.add.at(delta, my_e, -1)
        raw_starts.append(my_s)
        raw_ends.append(my_e)
        # cumulative count of covered positions, so coverage of [s,e) = cumcov[e] - cumcov[s]
        raw_cumcov.append(np.concatenate([[0], np.cumsum(np.cumsum(delta[:-1]) > 0)]))
    #
    # remove hits of kmers that overlap with a hit from any of their superstrings,
    #
    col_starts = []
    col_ends = []
    for kmer_list_i in range(len(KMER_LIST)):
        my_s = raw_starts[kmer_list_i]
        my_e = raw_ends[kmer_list_i]
        keep = np.ones(len(my_s), dtype=bool)
        for sub_i in KMER_ISSUBSTRING[kmer_list_i]:
            keep &= (raw_cumcov[sub_i][my_e] - raw_cumcov[sub_i][my_s]) == 0
        my_s = my_s[keep]
        my_e = my_e[keep]
        if len(my_s):
            new_block = np.concatenate([[True], my_s[1:] != my_e[:-1]])
            block_id = np.cumsum(new_block) - 1
            block_end = np.zeros(int(block_id[-1]) + 1, dtype='<i8')
            block_end[block_id] = my_e  # last write per block wins = end of the block's last hit
            my_s = my_s[new_block]
            my_e = block_end
        col_starts.append(my_s)
        col_ends.append(my_e)
    #
    # truncate hits overlapping the prefix of another hit:
    #
    all_starts = np.concatenate(col_starts) if col_starts else np.zeros(0, dtype='<i8')
    all_kmer = np.concatenate([np.full(len(n), i, dtype='<i8') for (i, n) in enumerate(col_starts)]) if col_starts else np.zeros(0, dtype='<i8')
    out_dat = []
    for kmer_list_i in range(len(KMER_LIST)):
        my_s = col_starts[kmer_list_i]
        my_e = col_ends[kmer_list_i].copy()
        other_starts = np.sort(all_starts[all_kmer != kmer_list_i])
        if len(my_s) and len(other_starts):
            next_i = np.searchsorted(other_starts, my_s, side='right')
            has_next = next_i < len(other_starts)
            next_start = np.where(has_next, other_starts[np.minimum(next_i, len(other_starts) - 1)], my_e)
            my_e = np.where(has_next & (next_start < my_e), next_start, my_e)
        out_dat.append(np.column_stack([my_s, my_e]).tolist() if len(my_s) else [])
    #
    return out_dat
