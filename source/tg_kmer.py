import pywt
import re

import numpy as np

from source.tg_util import exists_and_is_nonzero, RC

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


def get_telomere_kmer_density(read_dat, kmer_list, tel_window, mode='hifi', smoothing=False):
    re_hits = [[], []]
    for i in range(len(kmer_list)):
        re_hits[0].extend([(n.start(0), n.end(0), i, 0) for n in re.finditer(kmer_list[i], read_dat)])
    #
    tel_hit_p0 = np.zeros(len(read_dat))
    tel_hit_p1 = np.zeros(len(read_dat))
    for re_hit in re_hits[0]:
        tel_hit_p0[re_hit[0]:re_hit[1]] = 1
    for re_hit in re_hits[1]:
        tel_hit_p1[re_hit[0]:re_hit[1]] = 1
    tel_hit_cum_e0 = np.cumsum(tel_hit_p0)
    tel_hit_cum_e1 = np.cumsum(tel_hit_p1)
    #
    tel_density_e0 = []
    tel_density_e1 = []
    for i in range(len(read_dat) - tel_window):
        tel_density_e0.append(float(tel_hit_cum_e0[i+tel_window] - tel_hit_cum_e0[i]) / tel_window)
        tel_density_e1.append(float(tel_hit_cum_e1[i+tel_window] - tel_hit_cum_e1[i]) / tel_window)
    #
    if smoothing:
        tel_density_e0 = wavelet_smooth(tel_density_e0)
        tel_density_e1 = wavelet_smooth(tel_density_e1)
    #
    return (tel_density_e0, tel_density_e1)


def get_telomere_base_count(read_dat, kmer_list, mode='hifi'):
    re_hits = []
    for i in range(len(kmer_list)):
        re_hits.extend([(n.start(0), n.end(0), i, 0) for n in re.finditer(kmer_list[i], read_dat)])
    #
    tel_hit = np.zeros(len(read_dat),dtype='<i4')
    for re_hit in re_hits:
        tel_hit[re_hit[0]:re_hit[1]] = 1
    tel_hit = tel_hit.tolist()
    return tel_hit.count(1)


def get_telomere_regions(td_p_e0, td_p_e1, td_q_e0, td_q_e1, tel_window, pthresh, mode='hifi', smoothing=True):
    min_win_pos  = min([len(td_p_e0), len(td_p_e1), len(td_q_e0), len(td_q_e1)])
    max_win_pos  = max([len(td_p_e0), len(td_p_e1), len(td_q_e0), len(td_q_e1)])
    p_vs_q_power = np.zeros(max_win_pos)
    tel_regions  = [[0,None,None]]
    # somehow we received a sequence that's too short to work with (what went wrong?)
    if max_win_pos < tel_window:
        return ([], [[0,0,None]])
    #
    for i in range(min_win_pos):
        p_vs_q_power[i] = td_p_e0[i] - td_q_e0[i]
    if smoothing:
        p_vs_q_power = wavelet_smooth(p_vs_q_power)
    p_vs_q_power = p_vs_q_power[:-tel_window]
    #
    if p_vs_q_power[0] >= pthresh:
        tel_regions[-1][2] = 'p'
    elif p_vs_q_power[0] <= -pthresh:
        tel_regions[-1][2] = 'q'
    for i in range(1,len(p_vs_q_power)):
        if tel_regions[-1][2] == 'p' and p_vs_q_power[i] >= pthresh:
            tel_regions[-1][1] = i+1
        elif tel_regions[-1][2] == 'q' and p_vs_q_power[i] <= -pthresh:
            tel_regions[-1][1] = i+1
        elif tel_regions[-1][2] is None and abs(p_vs_q_power[i]) < pthresh:
            tel_regions[-1][1] = i+1
        else:
            if tel_regions[-1][1] is None:  # previous region was only 1bp wide, ugh
                tel_regions[-1][1] = tel_regions[-1][0]+1
            tel_regions.append([i,None,None])
            if p_vs_q_power[i] >= pthresh:
                tel_regions[-1][2] = 'p'
            elif p_vs_q_power[i] <= -pthresh:
                tel_regions[-1][2] = 'q'
    if tel_regions[-1][1] is None:
        tel_regions[-1][1] = len(p_vs_q_power)-1
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
    out_dat = []
    coord_hit_dict = []
    for kmer_list_i in range(len(KMER_LIST)):
        # get all hits
        raw_hits = [(n.start(0), n.end(0)) for n in re.finditer(KMER_LIST[kmer_list_i], my_telseq)]
        coord_hit_dict.append({})
        for kmer_span in raw_hits:
            for j in range(kmer_span[0],kmer_span[1]):
                coord_hit_dict[-1][j] = True
        out_dat.append([n for n in raw_hits])
    #
    # remove hits of kmers that overlap with a hit from any of their superstrings
    #
    for kmer_list_i in range(len(KMER_LIST)):
        del_list = []
        for ksi in range(len(out_dat[kmer_list_i])):
            kmer_span  = out_dat[kmer_list_i][ksi]
            are_we_hit = False
            for sub_i in KMER_ISSUBSTRING[kmer_list_i]:
                for j in range(kmer_span[0], kmer_span[1]):
                    if j in coord_hit_dict[sub_i]:
                        are_we_hit = True
                        break
                if are_we_hit:
                    break
            if are_we_hit:
                del_list.append(ksi)
        for di in sorted(del_list, reverse=True):
            del out_dat[kmer_list_i][di]
    #
    # collapse adjacent hits into larger blocks (so we have less polygons to plot)
    #
    for kmer_list_i in range(len(KMER_LIST)):
        collapsed_kmer_spans = [[n[0],n[1]] for n in out_dat[kmer_list_i]]
        for j in range(len(collapsed_kmer_spans)-1,0,-1):
            if collapsed_kmer_spans[j-1][1] == collapsed_kmer_spans[j][0]:
                collapsed_kmer_spans[j-1][1] = collapsed_kmer_spans[j][1]
                del collapsed_kmer_spans[j]
        out_dat[kmer_list_i] = [n for n in collapsed_kmer_spans]
    #
    # truncate hits overlapping the prefix of another hit
    #
    del_list = []
    for kmer_list_i in range(len(KMER_LIST)):
        for kmer_list_j in range(len(KMER_LIST)):
            if kmer_list_j == kmer_list_i:
                continue
            for ksi in range(len(out_dat[kmer_list_i])):
                for ksj in range(len(out_dat[kmer_list_j])):
                    kmer_span_i = out_dat[kmer_list_i][ksi]
                    kmer_span_j = out_dat[kmer_list_j][ksj]
                    if kmer_span_j[0] > kmer_span_i[1]:
                        break
                    if kmer_span_i[0] < kmer_span_j[0] and kmer_span_i[1] > kmer_span_j[0]:
                        out_dat[kmer_list_i][ksi] = [kmer_span_i[0], kmer_span_j[0]]
    #
    return out_dat
