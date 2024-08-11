from Bio.Align import PairwiseAligner, substitution_matrices
import numpy as np

from collections import Counter

from source.tg_util import shuffle_seq


# we're going to pretend each kmer color is an amino acid, so alignment tools behave themselves
AMINO = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
AMINO_2_IND = {AMINO[i]:i for i in range(len(AMINO))}
UNKNOWN_LETTER = 'A'
GAP_CHR = '-'

#
# scoring matrix parameters
#
MATCH_NORMAL  = 5
XMATCH_NORMAL = -4
GAP_OPEN = -4
GAP_EXT  = -4
GAP_MATCH = 5
GAP_XMATCH = -4
MATCH_GAP_UNKNOWN = -2
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


def get_scoring_matrix(canonical_letter, canonical_score=MATCH_CANON):
    m = substitution_matrices.Array(''.join(AMINO) + GAP_CHR, dims=2)
    for c1 in AMINO:
        for c2 in AMINO:
            if c1 == c2:
                m[c1,c2] = float(MATCH_NORMAL)
            else:
                m[c1,c2] = float(XMATCH_NORMAL)
    for c1 in AMINO:
        m[c1,canonical_letter] = float(XMATCH_CANON)
        m[canonical_letter,c1] = float(XMATCH_CANON)
        m[c1,UNKNOWN_LETTER] = float(XMATCH_UNKNOWN)
        m[UNKNOWN_LETTER,c1] = float(XMATCH_UNKNOWN)
        m[c1,GAP_CHR] = float(GAP_XMATCH)
        m[GAP_CHR,c1] = float(GAP_XMATCH)
    m[canonical_letter,canonical_letter] = float(canonical_score)      # matching canonical
    m[UNKNOWN_LETTER,UNKNOWN_LETTER]     = float(MATCH_UNKNOWN)        # matching unknown
    m[GAP_CHR,GAP_CHR]                   = float(GAP_MATCH)            # matching gap
    m[canonical_letter,UNKNOWN_LETTER]   = float(MATCH_CANON_UNKNOWN)  # canon = unknown
    m[UNKNOWN_LETTER,canonical_letter]   = float(MATCH_CANON_UNKNOWN)  # canon = unknown
    m[UNKNOWN_LETTER,GAP_CHR]            = float(MATCH_GAP_UNKNOWN)    # gap = unknown
    m[GAP_CHR,UNKNOWN_LETTER]            = float(MATCH_GAP_UNKNOWN)    # gap = unknown
    return m


def get_aligner_object(scoring_matrix=None, gap_bool=(True,True)):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.alphabet = ''.join(AMINO) + GAP_CHR
    aligner.match_score = float(MATCH_NORMAL)
    aligner.mismatch_score = float(XMATCH_NORMAL)
    aligner.open_gap_score = float(GAP_OPEN)
    aligner.extend_gap_score = float(GAP_EXT)
    #
    if scoring_matrix is not None:
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
    return aligner


def tvr_distance(tvr_i, tvr_j, aligner, adjust_lens=True, min_viable=True, randshuffle=3):
    if adjust_lens:
        min_len = min(len(tvr_i), len(tvr_j))
        if min_viable:
            min_len = max(min_len, MIN_VIABLE_SEQ_LEN)
        seq_i = tvr_i[:min_len]
        seq_j = tvr_j[:min_len]
    else:
        seq_i = tvr_i
        seq_j = tvr_j
    # aln score
    aln_score = aligner.score(seq_i, seq_j)
    # iden score
    c1 = Counter(seq_i)
    c2 = Counter(seq_j)
    if aligner.substitution_matrix is None:
        iden_score = aligner.match_score * ((len(seq_i)+len(seq_j)) / 2.)
    else:
        is1 = sum([c1[n]*aligner.substitution_matrix[n,n] for n in c1.keys()])
        is2 = sum([c2[n]*aligner.substitution_matrix[n,n] for n in c2.keys()])
        iden_score = (is1+is2) / 2.
    # rand score
    rand_scores = []
    for k in range(randshuffle):
        rand_scores.append(aligner.score(shuffle_seq(seq_i), shuffle_seq(seq_j)))
    rand_score_shuffle = np.mean(rand_scores)
    # convert to distance, bounds checking
    if rand_score_shuffle >= aln_score:
        my_dist_out = MAX_SEQ_DIST
    else:
        my_dist_out = min(-np.log((aln_score - rand_score_shuffle) / (iden_score - rand_score_shuffle)), MAX_SEQ_DIST)
    my_dist_out = max(my_dist_out, MIN_SEQ_DIST)
    #print(aln_score, iden_score, rand_score_shuffle, my_dist_out)
    return my_dist_out
